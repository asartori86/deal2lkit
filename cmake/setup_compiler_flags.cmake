## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Setup default compiler flags: This file sets up sensible default compiler
# flags for the various platforms, compilers and build targets supported by
# the deal.II library.
#
#
# ####################
# #     FAT NOTE:    #
# ####################
#
# All configuration in setup_compiler_flags.cmake and
# setup_compiler_flags_<compiler>.cmake shall ONLY modify:
#
#   D2K_CXX_FLAGS
#   D2K_CXX_FLAGS_DEBUG
#   D2K_CXX_FLAGS_RELEASE
#   D2K_LINKER_FLAGS
#   D2K_LINKER_FLAGS_DEBUG
#   D2K_LINKER_FLAGS_RELEASE
#   D2K_DEFINITIONS
#   D2K_DEFINITIONS_DEBUG
#   D2K_DEFINITIONS_RELEASE
#   D2K_USER_DEFINITIONS
#   D2K_USER_DEFINITIONS_DEBUG
#   D2K_USER_DEFINITIONS_RELEASE
#
# All modifications shall be guarded with the ENABLE_IF_SUPPORTED
# or ENABLE_IF_LINKS macro, e.g.
#
#   ENABLE_IF_SUPPORTED(D2K_CXX_FLAGS "-fpic")
#   ENABLE_IF_LINKS(D2K_LINKER_FLAGS "-Wl,--as-needed")
#
# Checks for compiler features (such as C++11 support) and compiler
# specific bugs that
#   - usually set up further configuration (such as preprocessor
#     definitions)
#   - disable a specific flag for a specific compiler version.
#
# belong the corresponding file:
#
#   ./cmake/checks/check_01_compiler_features.cmake
#   ./cmake/checks/check_01_cpu_features.cmake
#   ./cmake/checks/check_01_cxx_features.cmake
#   ./cmake/checks/check_01_system_features.cmake
#   ./cmake/checks/check_02_compiler_bugs.cmake
#


########################################################################
#                                                                      #
#                            Sanity checks:                            #
#                                                                      #
########################################################################

#
# Check the user provided CXX flags:
#

IF(NOT "${D2K_CXX_FLAGS_SAVED}" STREQUAL "${CACHED_D2K_CXX_FLAGS_SAVED}"
   OR NOT "${D2K_LINKER_FLAGS_SAVED}" STREQUAL "${CACHED_D2K_LINKER_FLAGS_SAVED}")
  # Rerun this test if cxx flags changed:
  UNSET(D2K_HAVE_USABLE_CXX_FLAGS CACHE)
ELSE()
  SET(D2K_HAVE_USABLE_CXX_FLAGS TRUE CACHE INTERNAL "")
ENDIF()
SET(CACHED_D2K_CXX_FLAGS_SAVED "${D2K_CXX_FLAGS_SAVED}" CACHE INTERNAL "" FORCE)
SET(CACHED_D2K_LINKER_FLAGS_SAVED "${D2K_LINKER_FLAGS_SAVED}" CACHE INTERNAL "" FORCE)

# Initialize all CMAKE_REQUIRED_* variables a this point:
RESET_CMAKE_REQUIRED()

CHECK_CXX_SOURCE_COMPILES(
  "int main(){ return 0; }"
  D2K_HAVE_USABLE_CXX_FLAGS)

IF(NOT D2K_HAVE_USABLE_CXX_FLAGS)
  UNSET(D2K_HAVE_USABLE_CXX_FLAGS CACHE)
  MESSAGE(FATAL_ERROR "
Configuration error: Cannot compile with the user supplied flags:
CXX flags: ${D2K_CXX_FLAGS_SAVED}
LD flags: ${D2K_LINKER_FLAGS_SAVED}
Please check the CMake variables D2K_CXX_FLAGS, D2K_LINKER_FLAGS
and the environment variables CXXFLAGS, LDFLAGS.\n\n"
    )
ENDIF()


########################################################################
#                                                                      #
#                           Compiler setup:                            #
#                                                                      #
########################################################################

IF(D2K_SETUP_DEFAULT_COMPILER_FLAGS)
  #
  # *Hooray* We are allowed to set compiler flags :-]
  #

  #
  # General setup for GCC and compilers sufficiently close to GCC:
  #
  IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
      CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
    VERBOSE_INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags_gnu.cmake)
    SET(D2K_KNOWN_COMPILER TRUE)
  ENDIF()

  #
  # Setup for ICC compiler (version >= 10):
  #
  IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    VERBOSE_INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags_intel.cmake)
    SET(D2K_KNOWN_COMPILER TRUE)
  ENDIF()

  #
  # Setup for MSVC compiler (version >= 2012):
  #
   IF(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    VERBOSE_INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags_msvc.cmake)
    SET(D2K_KNOWN_COMPILER TRUE)
  ENDIF()

  IF(NOT D2K_KNOWN_COMPILER)
    MESSAGE(FATAL_ERROR "\n"
      "Unknown compiler!\n"
      "If you're serious about it, set D2K_SETUP_DEFAULT_COMPILER_FLAGS=OFF "
      "and set the relevant compiler options by hand.\n\n"
      )
  ENDIF()

ELSE(D2K_SETUP_DEFAULT_COMPILER_FLAGS)

  MESSAGE(STATUS
    "Skipped setup of default compiler flags "
    "(D2K_SETUP_DEFAULT_COMPILER_FLAGS=OFF)"
    )
ENDIF(D2K_SETUP_DEFAULT_COMPILER_FLAGS)