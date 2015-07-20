#!/bin/sh

PRG=$PWD/programs
CASA=$PWD

if [ ! -d programs ] 
then
	echo "create folder `$PRG`"
	mkdir $PRG
fi

# installing cmake
if [ ! -d $PRG/cmake ]
then
	echo "installing cmake"
	wget http://www.cmake.org/files/v3.3/cmake-3.3.0-rc3-Linux-x86_64.tar.gz > /dev/null
	tar xzf cmake-3.3.0-rc3-Linux-x86_64.tar.gz -C $PRG
	rm cmake-3.3.0-rc3-Linux-x86_64.tar.gz
	mv $PRG/cmake-* $PRG/cmake
fi

# installing numdiff
if [ ! -d $PRG/numdiff ]
then
	echo "installing numdiff"
	mkdir $PRG/numdiff-tmp
	cd $PRG/numdiff-tmp
	wget http://mirror.lihnidos.org/GNU/savannah//numdiff/numdiff-5.8.1.tar.gz
	tar xzf numdiff-5.8.1.tar.gz
	rm numdiff-5.8.1.tar.gz
	cd numdiff-5.8.1
	DST_INST=$PRG/numdiff
	./configure --prefix=$DST_INST > /dev/null
	make -j4 install > /dev/null
	cd $CASA
	rm -rf $PRG/numdiff-tmp
fi

# installing ninja
if [ ! -d $PRG/ninja ]
then
	echo "installing ninja"
	cd $PRG
	git clone git://github.com/martine/ninja.git
	cd ninja
	git checkout release
	./configure.py --bootstrap > /dev/null
	cd $CASA
fi

# astyle
if [ ! -d $PRG/astyle ]
then
	echo "Downloading and installing astyle."
	mkdir $PRG/astyle
	wget http://downloads.sourceforge.net/project/astyle/astyle/astyle%202.04/astyle_2.04_linux.tar.gz  > /dev/null
	tar xfz astyle_2.04_linux.tar.gz -C $PRG > /dev/null
	cd $PRG/astyle/build/gcc
	make -j4 > /dev/null
	cd $CASA
fi

# petsc
if [ ! -d $PRG/petsc ]
then
	unset PETSC_DIR
	unset PETSC_ARCH
	echo "installing petsc"
	cd $PRG
	git clone -b maint https://bitbucket.org/petsc/petsc petsc-tmp
	export PETSC_ARCH=arch
	DST=$PRG/petsc
	cd petsc-tmp
	./configure \
		--prefix=$DST \
		--with-debugging=0 \
		--with-shared-libraries \
		--download-mpich \
	  --download-fblas-lapack 1> /dev/null 2>&1
	#	--download-parmetis \
	#	--download-metis \
	#	--download-hypre \
	#	--download-mumps \
	#	--download-scalapack \
	#	--download-sundials \
	#	--download-superlu \
	#	--download-superlu_dist \
	#	--download-hdf5
	make PETSC_DIR=$PWD PETSC_ARCH=arch all
	make PETSC_DIR=$PWD PETSC_ARCH=arch install
	export PETSC_DIR=$PRG/petsc
	export PETSC_ARCH=
	cd $CASA 
	rm -rf $PRG/petsc-tmp
fi


# trilinos
if [ ! -d $PRG/trilinos ]
then
	echo "installing trilinos"
	DST_INST=$PRG/trilinos
  export PATH=$PRG/cmake/bin:$PATH
	cd $PRG
	git clone https://github.com/trilinos/trilinos.git trilinos-tmp
	cd trilinos-tmp
	mkdir build
	cd build
	cmake \
		-G Ninja \
		-D CMAKE_BUILD_TYPE:STRING=RELEASE \
		-D TPL_ENABLE_Boost:BOOL=OFF \
		-D TPL_ENABLE_BoostLib:BOOL=OFF \
		-D TPL_ENABLE_MPI:BOOL=ON \
		-D TPL_ENABLE_Netcdf:BOOL=OFF \
		-D CMAKE_INSTALL_PREFIX:PATH=$DST_INST \
		-D Trilinos_ENABLE_OpenMP:BOOL=OFF \
		-D BUILD_SHARED_LIBS:BOOL=ON \
		-D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
		-D CMAKE_VERBOSE_MAKEFILE:BOOL=FALSE \
		-D Trilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF \
		-D Trilinos_ENABLE_TESTS:BOOL=OFF \
		-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
		-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
		-D Trilinos_ENABLE_Epetra:BOOL=ON \
		-D Trilinos_ENABLE_NOX:BOOL=ON \
		-D Trilinos_ENABLE_EpetraExt:BOOL=ON \
		-D Trilinos_ENABLE_Tpetra:BOOL=ON \
		-D Trilinos_ENABLE_Kokkos:BOOL=ON \
		-D Trilinos_ENABLE_Sacado:BOOL=ON \
		-D Trilinos_ENABLE_Amesos:BOOL=ON \
		-D Trilinos_ENABLE_AztecOO:BOOL=ON \
		-D Trilinos_ENABLE_Ifpack:BOOL=ON \
		-D Trilinos_ENABLE_Rythmos:BOOL=ON \
		-D Trilinos_ENABLE_Piro:BOOL=ON \
		-D Trilinos_ENABLE_MOOCHO:BOOL=ON \
		-D Trilinos_ENABLE_ML:BOOL=ON \
		-D Trilinos_ENABLE_MueLu:BOOL=ON \
		-D Trilinos_ENABLE_Komplex:BOOL=ON \
		-D Trilinos_ENABLE_Thyra:BOOL=ON \
		-D Trilinos_ENABLE_TrilinosCouplings:BOOL=ON \
		-D Trilinos_ENABLE_Fortran:BOOL=OFF \
		-D MPI_BASE_DIR:PATH=$PETSC_DIR \
		-D CMAKE_CXX_COMPILER:PATH=/usr/bin/clang++-3.6 \
		-D CMAKE_CXX_FLAGS:STRING=-w \
		-D CMAKE_C_COMPILER:PATH=/usr/bin/clang-3.6 \
		-D CMAKE_C_FLAGS:STRING=-w \
		.. #> $CASA/trilinos_cmake.log 2>&1

	ninja -j4 
	ninja -j4 install > /dev/null
	cd $PRG
	rm -rf trilinos-tmp
	cd $CASA
	tar cfz $PRG/trilinos-serial-CI-build.tgz $PRG/trilinos
fi

# dealii
if [ ! -d $PRG/dealii ]
then
	echo "installing dealii"
	DST_INST=$PRG/dealii 
	cd $PRG
	git clone https://github.com/dealii/dealii.git dealii-tmp
	cd dealii-tmp
	mkdir build
	cd build
  export PATH=$PRG/cmake/bin:$PATH
  export PATH=$PWD/programs/ninja:$PATH
	export PETSC_DIR=$PRG/petsc
	export PETSC_ARCH=
	cmake \
		-G Ninja \
		-D CMAKE_INSTALL_PREFIX:PATH=$DST_INST \
		-D CMAKE_CXX_FLAGS:STRING=-w \
		-D DEAL_II_WITH_MPI:BOOL=ON	\
		-D DEAL_II_WITH_THREADS:BOOL=ON \
		.. #> $CASA/dealii_cmake.log 2>&1
	ninja -j3 
	ninja -j4 install > /dev/null
	cd $PRG
	rm -rf dealii-tmp
	cd $CASA
	tar cfz $PRG/dealii-trilinos-serial-CI-build.tgz $PRG/dealii
fi

