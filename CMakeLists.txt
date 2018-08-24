cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0042 NEW)

set(CMAKE_USER_MAKE_RULES_OVERRIDE_Fortran "${CMAKE_SOURCE_DIR}/cmake/fortran_flags_override.cmake")

project(libmbd LANGUAGES Fortran)

option(ENABLE_SCALAPACK "Enable parallelisation with ScaLAPACK/MPI")
option(BUILD_SHARED_LIBS "Build shared rather than static library" ON)

set(DEFAULT_BUILD_TYPE "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(DEFAULT_BUILD_TYPE "Debug")
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to ${DEFAULT_BUILD_TYPE} as none was set")
    set(CMAKE_BUILD_TYPE ${DEFAULT_BUILD_TYPE})
endif()

include(GNUInstallDirs)

if(NOT LAPACK_LIBRARIES)
    find_package(LAPACK)
    if(NOT LAPACK_FOUND)
        message(SEND_ERROR "BLAS/LAPACK not found, specify custom installation with LAPACK_LIBRARIES variable")
    endif()
endif()

if(ENABLE_SCALAPACK)
    find_package(MPI)
    if(NOT MPI_Fortran_FOUND)
        message(SEND_ERROR "MPI not found, specify the MPI Fortran compiler with MPI_Fortran_COMPILER variable")
    endif()
    if(APPLE)
        # -flat_namespace (Darwin-specific) is causing a crash (seg fault) when the
        # Fortran library is called from Python and one writes into a character
        # variable, but only when the kind is the default one. It causes the
        # written to variable to appear as being four times shorter than it is.
        # Only mention of anything possibly related I could find is at
        #
        #   https://trac.mpich.org/projects/mpich/ticket/1590
        string(REPLACE "-Wl,-flat_namespace" "" MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
        string(STRIP "${MPI_Fortran_LINK_FLAGS}" MPI_Fortran_LINK_FLAGS)
    endif()
    if(NOT SCALAPACK_LIBRARIES)
        find_package(Scalapack QUIET)
        if(Scalapack_FOUND)
            message(STATUS "ScaLAPACK CMake package found in ${Scalapack_DIR}")
        else()
            message(SEND_ERROR "ScaLAPACK CMake package not found, specify custom installation with SCALAPACK_LIBRARIES variable")
        endif()
    endif()
endif()

enable_testing()

set(CMAKE_CTEST_COMMAND ctest -V)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

add_subdirectory(src)
