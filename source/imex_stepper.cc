//-----------------------------------------------------------
//
//    Copyright (C) 2015 - 2016 by the deal2lkit authors
//
//    This file is part of the deal2lkit library.
//
//    The deal2lkit library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal2lkit distribution.
//
//-----------------------------------------------------------


#include <deal2lkit/imex_stepper.h>
#ifdef D2K_WITH_SUNDIALS

#include <deal.II/base/utilities.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/lac/block_vector.h>
#ifdef DEAL_II_WITH_TRILINOS
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#endif
#ifdef DEAL_II_WITH_PETSC
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#endif
#include <deal.II/base/utilities.h>

#include <iostream>
#include <iomanip>

#include <math.h>

#ifdef DEAL_II_WITH_MPI
#include <nvector/nvector_parallel.h>
#endif

using namespace dealii;


D2K_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_MPI
template <typename VEC>
IMEXStepper<VEC>::IMEXStepper(std::string name,
                              MPI_Comm comm) :
  ParameterAcceptor(name),
  communicator(Utilities::MPI::duplicate_communicator(comm)),
  kinsol("KINSOL for IMEX",comm),
  pcout(std::cout,
        Utilities::MPI::this_mpi_process(communicator)==0)
{
  set_functions_to_trigger_an_assert();
}

#else
template <typename VEC>
IMEXStepper<VEC>::IMEXStepper(std::string name) :
  ParameterAcceptor(name),
  kinsol("KINSOL for IMEX"),
  pcout(std::cout)
{
  set_functions_to_trigger_an_assert();
}
#endif

template <typename VEC>
IMEXStepper<VEC>::~IMEXStepper()
{
#ifdef DEAL_II_WITH_MPI
  MPI_Comm_free(&communicator);
#endif
}

template <typename VEC>
double IMEXStepper<VEC>::get_alpha() const
{
  double alpha;
  if (initial_time == final_time || step_size == 0.0)
    alpha = 0.0;
  else
    alpha = 1./step_size;

  return alpha;
}

template <typename VEC>
void IMEXStepper<VEC>::set_initial_time(const double &t)
{
  initial_time = t;
}

template <typename VEC>
void IMEXStepper<VEC>::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &_step_size,
                "Step size", "1e-2", Patterns::Anything(),
                "Conditional statement can be used as well:\n"
                "(t<0.5?1e-2:1e-3)");

  add_parameter(prm, &abs_tol,
                "Absolute error tolerance", "1e-6",
                Patterns::Double());

  add_parameter(prm, &rel_tol,
                "Relative error tolerance", "1e-5",
                Patterns::Double());

  add_parameter(prm, &initial_time,
                "Initial time", "0.0",
                Patterns::Double());

  add_parameter(prm, &final_time,
                "Final time", "1.0",
                Patterns::Double());

  add_parameter(prm, &output_period,
                "Intervals between outputs", "1",
                Patterns::Integer());

  add_parameter(prm, &max_outer_non_linear_iterations,
                "Maximum number of outer nonlinear iterations", "5",
                Patterns::Integer(),
                "At each outer iteration the Jacobian is updated if it is set that the \n"
                "Jacobian is continuously updated and a cycle of inner iterations is \n"
                "perfomed.");

  add_parameter(prm, &max_inner_non_linear_iterations,
                "Maximum number of inner nonlinear iterations", "3",
                Patterns::Integer(),
                "At each inner iteration the Jacobian is NOT updated.");

  add_parameter(prm, &newton_alpha,
                "Newton relaxation parameter", "1",
                Patterns::Double());

  add_parameter(prm, &update_jacobian_continuously,
                "Update continuously Jacobian", "true",
                Patterns::Bool());

  add_parameter(prm, &n_max_backtracking,
                "Number of elements in backtracking sequence", "5",
                Patterns::Integer(1),
                "In the line seach method with backtracking the following alphas are\n"
                "tested: 1, 1/2, 1/4,..., 2^-i. This parameter sets the maximum i.");
  add_parameter(prm, &method,
                "Method used", "fixed_alpha",
                Patterns::Selection("fixed_alpha|LS_backtracking"),
                "Fixed alpha means that the parsed alpha is used in each Newton iteration\n"
                "LS_backtracking is the line search with backtracking method.");

  add_parameter(prm, &verbose,
                "Print useful informations", "false",
                Patterns::Bool());

  add_parameter(prm, &use_kinsol,
                "Use the KINSOL solver", "true",
                Patterns::Bool());

}

template <typename VEC>
void IMEXStepper<VEC>::compute_y_dot(const VEC &y, const VEC &prev, const double alpha, VEC &y_dot)
{
  y_dot = y;
  y_dot -= prev;
  y_dot *= alpha;
}

template <typename VEC>
unsigned int IMEXStepper<VEC>::solve_dae(VEC &solution, VEC &solution_dot)
{
  unsigned int step_number = 0;

  auto previous_solution = create_new_vector();
  auto residual_ = create_new_vector();
  auto rhs = create_new_vector();

  double t = initial_time;
  double alpha;
  step_size = evaluate_step_size(t);
  alpha = get_alpha();
  // check if it is a stationary problem
  if (initial_time == final_time)
    alpha = 0.0;
  else
    alpha = 1./step_size;

  compute_previous_solution(solution,solution_dot,alpha, *previous_solution);



  std::function<int(const VEC &, VEC &)> my_residual = [&] (const VEC &y, VEC &res)
  {
    double a = this->get_alpha();
    compute_y_dot(y,*previous_solution,a,solution_dot);
    int ret = this->residual(t,y,solution_dot,res);
    *residual_ = res;
    AssertThrow(!std::isnan(residual_->l2_norm()),ExcMessage("Residual contains one or more NaNs."));
    return ret;
  };

  std::function<int(const VEC &)> my_jac = [&] (const VEC &y)
  {
    double a = this->get_alpha();
    compute_y_dot(y,*previous_solution,a,solution_dot);
    return this->setup_jacobian(t,y,solution_dot,a);
  };

  std::function<int(const VEC &, VEC &)> my_solve = [&] (const VEC &, VEC &dst)
  {
    *rhs = *residual_;
    *rhs *= -1.0;
    return this->solve_jacobian_system(*rhs,dst);
  };


  kinsol.create_new_vector = create_new_vector;
  kinsol.residual = my_residual;
  kinsol.setup_jacobian = my_jac;
  kinsol.solve_linear_system = my_solve;
  kinsol.jacobian_vmult = jacobian_vmult;

  // Initialization of the state of the boolean variable
  // responsible to keep track of the requirement that the
  // system's Jacobian be updated.
  bool update_Jacobian = true;

  // silence a warning when using kinsol
  (void) update_Jacobian;


//   store initial conditions
  output_step(t, solution, solution_dot,  step_number);

  bool restart=false;

  auto L2 = create_new_vector();
  if (use_kinsol)
    {
      // call kinsol initialization. this is mandatory if I am doing multiple cycle in pi-DoMUS
      kinsol.initialize_solver(solution);
      *L2 = get_lumped_mass_matrix();
      kinsol.set_scaling_vectors(*L2, *L2);
    }

  compute_consistent_initial_conditions(initial_time,
                                        solution,
                                        solution_dot);
  // The overall cycle over time begins here.
  while (t<=final_time+1e-15)
    {
      pcout << "Solving for t = " << t
            << " (step size = "<< step_size<<")"
            << std::endl;

      if (use_kinsol)
        {
          kinsol.solve(solution);
          compute_y_dot(solution,*previous_solution,alpha,solution_dot);
        }
      else
        do_newton(t,alpha,update_Jacobian,*previous_solution,solution,solution_dot);

      restart = solver_should_restart(t,solution,solution_dot);

      while (restart)
        {

          previous_solution = create_new_vector();
          residual_ = create_new_vector();
          rhs = create_new_vector();
          L2 = create_new_vector();


          if (use_kinsol)
            {


              kinsol.initialize_solver(solution);
              *L2 = get_lumped_mass_matrix();
              kinsol.set_scaling_vectors(*L2, *L2);
            }
          compute_consistent_initial_conditions(t,
                                                solution,
                                                solution_dot);

          compute_previous_solution(solution,solution_dot,alpha,*previous_solution);

          restart = solver_should_restart(t,solution,solution_dot);

        }


      step_number += 1;

      if ((step_number % output_period) == 0)
        output_step(t, solution, solution_dot,  step_number);
      step_size = evaluate_step_size(t);
      t += step_size;

      if (initial_time == final_time)
        alpha = 0.0;
      else
        alpha = 1./step_size;

      *previous_solution = solution;
      update_Jacobian = update_jacobian_continuously;

    } // End of the cycle over time.
  return 0;
}



template <typename VEC>
double IMEXStepper<VEC>::
line_search_with_backtracking(const VEC &update,
                              const VEC &previous_solution,
                              const double &alpha,
                              const double &t,
                              VEC &solution,
                              VEC &solution_dot,
                              VEC &residual)
{
  auto first_trial = create_new_vector();
  auto first_residual = create_new_vector();

  *first_trial = solution;
  double n_alpha = 1.0;

  first_trial->sadd(1.0, n_alpha, update);

  solution_dot = *first_trial;
  solution_dot -= previous_solution;
  solution_dot *= alpha;

  this->residual(t, *first_trial, solution_dot, *first_residual);

  double first_res_norm = vector_norm(*first_residual);

  auto second_trial = create_new_vector();
  auto second_residual = create_new_vector();

  for (unsigned int i=1; i<=n_max_backtracking; ++i)
    {
      *second_trial = solution;

      n_alpha = 1.0/(std::pow(2.0,i));

      second_trial->sadd(1.0, n_alpha, update);

      solution_dot = *second_trial;
      solution_dot -= previous_solution;
      solution_dot *= alpha;

      this->residual(t, *second_trial, solution_dot, *second_residual);
      double second_res_norm = vector_norm(*second_residual);
      if (first_res_norm < second_res_norm)
        {
          solution = *first_trial;
          residual = *first_residual;

          solution_dot = *first_trial;
          solution_dot -= previous_solution;
          solution_dot *= alpha;
          n_alpha = 1.0/(std::pow(2.0,i-1));
          return n_alpha;

        }
      else if (i<(n_max_backtracking))
        {
          *first_trial = *second_trial;
          *first_residual = *second_residual;
          first_res_norm = second_res_norm;
        }
      else
        {
          solution = *second_trial;
          residual = *second_residual;
          /* solution dot is ok */
          return n_alpha;
        }
    }
  return n_alpha;
}


template <typename VEC>
void IMEXStepper<VEC>::
do_newton (const double t,
           const double alpha,
           const bool update_Jacobian,
           const VEC &previous_solution,
           VEC &solution,
           VEC &solution_dot)
{
  auto solution_update = create_new_vector();
  auto res = create_new_vector();
  auto rhs = create_new_vector();



  // Initialization of two counters for the monitoring of
  // progress of the nonlinear solver.
  unsigned int inner_iter = 0;
  unsigned int outer_iter = 0;
  unsigned int nonlin_iter = 0;
  this->residual(t, solution, solution_dot, *res);
  double res_norm = 0.0;
  double solution_norm = 0.0;

  if (abs_tol>0.0||rel_tol>0.0)
    res_norm = this->vector_norm(*res);
  // if (rel_tol>0.0)
  //   solution_norm = interface.vector_norm(solution);

  // The nonlinear solver iteration cycle begins here.
  // using a do while approach, we ensure that the system
  // is solved at least once
  do
    {
      outer_iter += 1;
      if (update_Jacobian == true)
        {
          setup_jacobian(t,
                         solution,
                         solution_dot,
                         alpha);
        }

      inner_iter = 0;
      while (inner_iter < max_inner_non_linear_iterations &&
             res_norm > abs_tol &&
             res_norm > rel_tol*solution_norm)
        {


          inner_iter += 1;

          *rhs = *res;
          *rhs *= -1.0;

          solve_jacobian_system(*rhs, *solution_update);


          if (method == "LS_backtracking")
            {
              newton_alpha = line_search_with_backtracking(*solution_update,
                                                           previous_solution,
                                                           alpha,
                                                           t,
                                                           solution,
                                                           solution_dot,
                                                           *res);
            }
          else if (method == "fixed_alpha")
            {
              solution.sadd(1.0,
                            newton_alpha, *solution_update);

            }

          compute_y_dot(solution,previous_solution,alpha, solution_dot);

          res_norm = vector_norm(*solution_update);

          AssertThrow(!std::isnan(res->l2_norm()),ExcMessage("Residual contains one or more NaNs."));

          if (rel_tol>0.0)
            {
              solution_norm = vector_norm(solution);

              if (verbose)
                {
                  pcout << std::endl
                        << "   "
                        << " iteration "
                        << nonlin_iter + inner_iter
                        << ":\n"
                        << std::setw(19) << std::scientific << res_norm
                        << "   update norm\n"
                        << std::setw(19) << std::scientific << solution_norm
                        << "   solution norm\n"
                        << std::setw(19) << newton_alpha
                        << "   newton alpha\n\n"
                        << std::endl;
                }
            }
          else if (verbose)
            {
              pcout << std::endl
                    << "   "
                    << " iteration "
                    << nonlin_iter + inner_iter
                    << ":\n"
                    << std::setw(19) << std::scientific << res_norm
                    << "   update norm\n"
                    << std::setw(19) << newton_alpha
                    << "   newton alpha\n\n"
                    << std::endl;
            }

          this->residual(t,solution,solution_dot,*res);
        }

      nonlin_iter += inner_iter;

      if (std::fabs(res_norm) < abs_tol ||
          std::fabs(res_norm) < rel_tol*solution_norm)
        {
          pcout << std::endl
                << "   "
                << std::setw(19) << std::scientific << res_norm
                << " (converged in "
                << nonlin_iter
                << " iterations)\n\n"
                << std::endl;
          break; // Break of the while cycle
        }
      else if (outer_iter == max_outer_non_linear_iterations)
        {
          pcout << std::endl
                << "   "
                << std::setw(19) << std::scientific << res_norm
                << " (not converged in "
                << std::setw(3) << nonlin_iter
                << " iterations)\n\n"
                << std::endl;
          AssertThrow(false,
                      ExcMessage ("No convergence in nonlinear solver"));
        }

    } // The nonlinear solver iteration cycle ends here.
  while (outer_iter < max_outer_non_linear_iterations &&
         res_norm > abs_tol &&
         res_norm > rel_tol*solution_norm);
}

template <typename VEC>
void IMEXStepper<VEC>::
compute_previous_solution(const VEC &sol,
                          const VEC &sol_dot,
                          const double &alpha,
                          VEC &prev)
{
  if (alpha > 0.0)
    {
      prev = sol_dot;
      prev /= (-1.0*alpha);
      prev += sol;
    }
  else
    {
      prev = sol;
      (void)sol_dot;
      (void)alpha;
    }
}

template <typename VEC>
double IMEXStepper<VEC>::evaluate_step_size(const double &t)
{
  std::string variables = "t";
  std::map<std::string,double> constants;
  // FunctionParser with 1 variables and 1 component:
  FunctionParser<1> fp(1);
  fp.initialize(variables,
                _step_size,
                constants);
  // Point at which we want to evaluate the function
  Point<1> time(t);
  // evaluate the expression at 'time':
  double result = fp.value(time);
  return result;
}

template <typename VEC>
void IMEXStepper<VEC>::compute_consistent_initial_conditions(const double &t,
    VEC &y,
    VEC &y_dot)
{
  auto previous_solution = create_new_vector();
  auto first_guess = create_new_vector();
  auto first_guess_dot = create_new_vector();
  auto update = create_new_vector();

  *first_guess     = y;
  *first_guess_dot = y_dot;



  step_size = 0;
  if (use_kinsol)
    {
      kinsol.solve(y);
    }
  else
    do_newton(t,0.0,true,*previous_solution,y,y_dot);

  *update = y;
  *update -= *first_guess;
  step_size = evaluate_step_size(t);
  *update /= step_size;
  *first_guess_dot += *update;
  y_dot = *first_guess_dot;


}


template<typename VEC>
void IMEXStepper<VEC>::set_functions_to_trigger_an_assert()
{

  create_new_vector = []() ->shared_ptr<VEC>
  {
    shared_ptr<VEC> p;
    AssertThrow(false, ExcPureFunctionCalled("Please implement create_new_vector function."));
    return p;
  };

  residual = [](const double,
                const VEC &,
                const VEC &,
                VEC &) ->int
  {
    int ret=0;
    AssertThrow(false, ExcPureFunctionCalled("Please implement residual function."));
    return ret;
  };

  setup_jacobian = [](const double,
                      const VEC &,
                      const VEC &,
                      const double) ->int
  {
    int ret=0;
    AssertThrow(false, ExcPureFunctionCalled("Please implement setup_jacobian function."));
    return ret;
  };

  solve_jacobian_system = [](const VEC &,
                             VEC &) ->int
  {
    int ret=0;
    AssertThrow(false, ExcPureFunctionCalled("Please implement solve_jacobian_system function."));
    return ret;
  };

  output_step = [](const double,
                   const VEC &,
                   const VEC &,
                   const unsigned int)
  {
    AssertThrow(false, ExcPureFunctionCalled("Please implement output_step function."));
  };

  solver_should_restart = [](const double,
                             VEC &,
                             VEC &) ->bool
  {
    bool ret=false;
    AssertThrow(false, ExcPureFunctionCalled("Please implement solver_should_restart function."));
    return ret;
  };

  get_lumped_mass_matrix = []() ->VEC &
  {
    shared_ptr<VEC> y;
    AssertThrow(false, ExcPureFunctionCalled("Please implement get_lumped_mass_matrix function."));
    return *y;
  };

  jacobian_vmult = [](const VEC &,
                      VEC &) ->int
  {
    int ret=0;
    AssertThrow(false, ExcPureFunctionCalled("Please implement jacobian_vmult function."));
    return ret;
  };

  vector_norm = [](const VEC &vector) ->double
  {
    return vector.l2_norm();
  };
}

D2K_NAMESPACE_CLOSE

template class deal2lkit::IMEXStepper<BlockVector<double> >;

#ifdef DEAL_II_WITH_MPI

#ifdef DEAL_II_WITH_TRILINOS
template class deal2lkit::IMEXStepper<TrilinosWrappers::MPI::Vector>;
template class deal2lkit::IMEXStepper<TrilinosWrappers::MPI::BlockVector>;
#endif

#ifdef DEAL_II_WITH_PETSC
template class deal2lkit::IMEXStepper<PETScWrappers::MPI::Vector>;
template class deal2lkit::IMEXStepper<PETScWrappers::MPI::BlockVector>;
#endif

#endif

#endif
