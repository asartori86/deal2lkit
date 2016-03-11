//-----------------------------------------------------------
//
//    Copyright (C) 2015 by the deal2lkit authors
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


#include <deal2lkit/sundials_interface.h>
#include <deal2lkit/imex_stepper.h>
#ifdef D2K_WITH_SUNDIALS

#include <deal.II/base/utilities.h>
#include <deal.II/lac/block_vector.h>
#ifdef DEAL_II_WITH_TRILINOS
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
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


template <typename VEC>
IMEXStepper<VEC>::IMEXStepper(SundialsInterface<VEC> &interface,
                              const double &step_size,
                              const double &initial_time,
                              const double &final_time) :
  ParameterAcceptor("IMEX Parameters"),
  interface(interface),
  step_size(step_size),
  initial_time(initial_time),
  final_time(final_time),
  pout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
{
  abs_tol = 1e-6;
  rel_tol = 1e-8;
  output_period = 1;
  newton_alpha = 1.0;
  max_outer_non_linear_iterations = 5;
  max_inner_non_linear_iterations = 3;
  verbose = false;
  n_max_backtracking = 5;
  method = "fixed alpha";
}

template <typename VEC>
void IMEXStepper<VEC>::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &step_size,
                "Step size", std::to_string(step_size), Patterns::Double());

  add_parameter(prm, &abs_tol,
                "Absolute error tolerance", std::to_string(abs_tol),
                Patterns::Double());

  add_parameter(prm, &rel_tol,
                "Relative error tolerance", std::to_string(rel_tol),
                Patterns::Double());

  add_parameter(prm, &initial_time,
                "Initial time", std::to_string(initial_time),
                Patterns::Double());

  add_parameter(prm, &final_time,
                "Final time", std::to_string(final_time),
                Patterns::Double());

  add_parameter(prm, &output_period,
                "Intervals between outputs", std::to_string(output_period),
                Patterns::Integer());

  add_parameter(prm, &max_outer_non_linear_iterations,
                "Maximum number of outer nonlinear iterations", std::to_string(max_outer_non_linear_iterations),
                Patterns::Integer(),
                "At each outer iteration the Jacobian is updated if it is set that the \n"
                "Jacobian is continuously updated and a cycle of inner iterations is \n"
                "perfomed.");

  add_parameter(prm, &max_inner_non_linear_iterations,
                "Maximum number of inner nonlinear iterations", std::to_string(max_inner_non_linear_iterations),
                Patterns::Integer(),
                "At each inner iteration the Jacobian is NOT updated.");

  add_parameter(prm, &newton_alpha,
                "Newton relaxation parameter", std::to_string(newton_alpha),
                Patterns::Double());

  add_parameter(prm, &update_jacobian_continuously,
                "Update continuously Jacobian", "true",
                Patterns::Bool());

  add_parameter(prm, &n_max_backtracking,
                "Number of elements in backtracking sequence", std::to_string(n_max_backtracking),
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

}


template <typename VEC>
unsigned int IMEXStepper<VEC>::start_ode(VEC &solution, VEC &solution_dot)
{
  AssertDimension(solution.size(), interface.n_dofs());

  unsigned int step_number = 0;


  auto previous_solution = interface.create_new_vector();
  auto solution_update = interface.create_new_vector();
  auto residual = interface.create_new_vector();
  auto rhs = interface.create_new_vector();

  *previous_solution = solution;

  double t = initial_time;
  double alpha;

  // check if it is a stationary problem
  if (initial_time == final_time)
    alpha = 0.0;
  else
    alpha = 1./step_size;

  interface.output_step( 0, solution, solution_dot, 0, step_size);

  // Initialization of the state of the boolean variable
  // responsible to keep track of the requirement that the
  // system's Jacobian be updated.
  bool update_Jacobian = true;

  bool restart=false;
  // The overall cycle over time begins here.
  for (; t<=final_time+1e-15; t+= step_size, ++step_number)
    {
      pout << "Time = " << t << std::endl;
      // Implicit Euler scheme.
	solution_dot = solution;
	solution_dot -= *previous_solution;
	solution_dot *= alpha;

      // Initialization of two counters for the monitoring of
      // progress of the nonlinear solver.
      unsigned int inner_iter = 0;
      unsigned int outer_iter = 0;
      unsigned int nonlin_iter = 0;
      interface.residual(t, solution, solution_dot, *residual);
      double res_norm = 0.0;
      double solution_norm = 0.0;

      if (abs_tol>0.0||rel_tol>0.0)
        res_norm = interface.vector_norm(*residual);
      // if (rel_tol>0.0)
      //   solution_norm = interface.vector_norm(solution);

      // The nonlinear solver iteration cycle begins here.
      while (outer_iter < max_outer_non_linear_iterations &&
             res_norm > abs_tol &&
             res_norm > rel_tol*solution_norm)
        {
          outer_iter += 1;
          if (update_Jacobian == true)
            {
              interface.setup_jacobian(t, solution, solution_dot,
                                       *residual, alpha);
            }

          inner_iter = 0;
          while (inner_iter < max_inner_non_linear_iterations &&
                 res_norm > abs_tol &&
                 res_norm > rel_tol*solution_norm)
            {


              inner_iter += 1;

              *rhs = *residual;
              *rhs *= -1.0;

              interface.solve_jacobian_system(t, solution, solution_dot,
                                              *residual, alpha,
                                              *rhs, *solution_update);


              if (method == "LS_backtracking")
                {
                  newton_alpha = line_search_with_backtracking(*solution_update,
                                                               *previous_solution,
                                                               alpha,
                                                               t,
                                                               solution,
                                                               solution_dot,
                                                               *residual);
                }
              else if (method == "fixed_alpha")
                {
                  solution.sadd(1.0,
                                newton_alpha, *solution_update);

                  // Implicit Euler scheme.
                  solution_dot = solution;
                  solution_dot -= *previous_solution;
                  solution_dot *= alpha;
                }

              res_norm = interface.vector_norm(*solution_update);

              if (rel_tol>0.0)
                {
                  solution_norm = interface.vector_norm(solution);

                  if (verbose)
                    {
                      pout << std::endl
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
                  pout << std::endl
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

              interface.residual(t,solution,solution_dot,*residual);
            }

          nonlin_iter += inner_iter;

          if (std::fabs(res_norm) < abs_tol ||
              std::fabs(res_norm) < rel_tol*solution_norm)
            {
              pout << std::endl
                   << "   "
                   << std::setw(19) << std::scientific << res_norm
                   << " (converged in "
                   << nonlin_iter
                   << " iterations)\n\n"
                   << std::endl;
              break; // Break of the while cycle ... after this a time advancement happens.
            }
          else if (outer_iter == max_outer_non_linear_iterations)
            {
              pout << std::endl
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

      restart = interface.solver_should_restart(t,step_number,step_size,solution,solution_dot);

      if (restart)
        {
          previous_solution = interface.create_new_vector();
          solution_update = interface.create_new_vector();
          residual = interface.create_new_vector();
          rhs = interface.create_new_vector();
	  //	  *previous_solution = solution;

          t -= step_size;
          --step_number;
	  /// fix solution after restart

	  //	  solution_dot *= 0; // alpha=0

      interface.residual(t, solution, solution_dot, *residual);

      	  double res_norm_restart = 0.0;
      	  double solution_norm_restart = 0.0;

      	  unsigned int outer_iter_restart = 0;
      	  unsigned int inner_iter_restart = 0;
      	  res_norm_restart = interface.vector_norm(*residual);
      	  while (outer_iter_restart < max_outer_non_linear_iterations &&
      		 res_norm_restart > abs_tol &&
      		 res_norm_restart > rel_tol*solution_norm)
      	    {
      	      outer_iter_restart += 1;
      	      if (update_Jacobian == true)
      		{
      		  interface.setup_jacobian(t, solution, solution_dot,
      					   *residual, 0);
      		}

      	      inner_iter_restart = 0;
      	      while (inner_iter_restart < max_inner_non_linear_iterations &&
      		     res_norm_restart > abs_tol &&
      		     res_norm_restart > rel_tol*solution_norm)
      		{


      		  inner_iter_restart += 1;

      		  *rhs = *residual;
      		  *rhs *= -1.0;

      		  interface.solve_jacobian_system(t, solution, solution_dot,
      						  *residual, 0,
      						  *rhs, *solution_update);


      		  newton_alpha = line_search_with_backtracking(*solution_update,
      							       *previous_solution,
      							       0,
      							       t,
      							       solution,
      							       solution_dot,
      							       *residual);


      		  res_norm_restart = interface.vector_norm(*solution_update);

      		  if (rel_tol>0.0)
      		    {
      		      solution_norm_restart = interface.vector_norm(solution);

                      pout << std::endl
                           << "   "
      			//                           << " iteration "
                        //   << nonlin_iter + inner_iter
                           << ":\n"
                           << std::setw(19) << std::scientific << res_norm_restart
                           << "   update norm\n"
                           << std::setw(19) << std::scientific << solution_norm_restart
                           << "   solution norm\n"
                           << std::setw(19) << newton_alpha
                           << "   newton alpha\n\n"
                           << std::endl;

      		    }
      		  else 
      		    {
      		      pout << std::endl
      			   << "   "
      			   << " iteration "
      			//	   << nonlin_iter + inner_iter
      			   << ":\n"
      			   << std::setw(19) << std::scientific << res_norm_restart
      			   << "   update norm\n"
      			   << std::setw(19) << newton_alpha
      			   << "   newton alpha\n\n"
      			   << std::endl;
      		    }

      		  interface.residual(t,solution,solution_dot,*residual);
      		}

      	      if (std::fabs(res_norm) < abs_tol ||
      		  std::fabs(res_norm) < rel_tol*solution_norm)
      		{
      		  pout << std::endl
      		       << "   "
      		       << std::setw(19) << std::scientific << res_norm_restart
      		       << " (converged in "
      		    //  << nonlin_iter
      		       << " iterations)\n\n"
      		       << std::endl;
      		  *previous_solution = solution;
      		  break; // Break of the while cycle ... after this a time advancement happens.
      		}
      	      else if (outer_iter == max_outer_non_linear_iterations)
      		{
      		  pout << std::endl
      		       << "   "
      		       << std::setw(19) << std::scientific << res_norm_restart
      		       << " (not converged in "
      		    //		       << std::setw(3) << nonlin_iter
      		       << " iterations)\n\n"
      		       << std::endl;
      		  AssertThrow(false,
      			      ExcMessage ("No convergence in nonlinear solver"));
      		}
      	    }





	  // solution *=0;
	  // solution_dot *=0;






























	  

















	  
        }
      else
        {

          if ((step_number % output_period) == 0)
            interface.output_step(t, solution, solution_dot,  step_number, step_size);
        }
      if(!restart)
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
  auto first_trial = interface.create_new_vector();
  auto first_residual = interface.create_new_vector();

  *first_trial = solution;
  double n_alpha = 1.0;

  first_trial->sadd(1.0, n_alpha, update);

  solution_dot = *first_trial;
  solution_dot -= previous_solution;
  solution_dot *= alpha;

  interface.residual(t, *first_trial, solution_dot, *first_residual);

  double first_res_norm = interface.vector_norm(*first_residual);

  auto second_trial = interface.create_new_vector();
  auto second_residual = interface.create_new_vector();

  for (unsigned int i=1; i<=n_max_backtracking; ++i)
    {
      *second_trial = solution;

      n_alpha = 1.0/(std::pow(2.0,i));

      second_trial->sadd(1.0, n_alpha, update);

      solution_dot = *second_trial;
      solution_dot -= previous_solution;
      solution_dot *= alpha;

      interface.residual(t, *second_trial, solution_dot, *second_residual);
      double second_res_norm = interface.vector_norm(*second_residual);
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


D2K_NAMESPACE_CLOSE

template class deal2lkit::IMEXStepper<BlockVector<double> >;

#ifdef DEAL_II_WITH_MPI

#ifdef DEAL_II_WITH_TRILINOS
template class deal2lkit::IMEXStepper<TrilinosWrappers::MPI::Vector>;
template class deal2lkit::IMEXStepper<TrilinosWrappers::MPI::BlockVector>;
#endif

#endif

#endif
