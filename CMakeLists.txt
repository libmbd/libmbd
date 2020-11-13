cmake_minimum_required(VERSION 3.1)
cmake_policy(SET CMP0042 NEW)
if(NOT CMAKE_VERSION VERSION_LESS 3.13)
    cmake_policy(SET CMP0077 NEW)
endif()

set(CMAKE_USER_MAKE_RULES_OVERRIDE_Fortran
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/fortran_flags_override.cmake")

project(Libmbd
    VERSION 0.10.3
    DESCRIPTION "Many-body dispersion library"
    LANGUAGES Fortran C)

include(GNUInstallDirs)
include(CMakeDependentOption)
include(CTest)

option(ENABLE_SCALAPACK_MPI "Enable parallelisation with ScaLAPACK/MPI")
CMAKE_DEPENDENT_OPTION(ENABLE_ELSI "Enable ELSI interface" OFF ENABLE_SCALAPACK_MPI OFF)
option(ENABLE_C_API "Enable C API" ON)

option(BUILD_SHARED_LIBS "Build shared rather than static library" ON)

set(DEFAULT_BUILD_TYPE "Release")
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    set(DEFAULT_BUILD_TYPE "Debug")
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to ${DEFAULT_BUILD_TYPE} as none was set")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}")
endif()

separate_arguments(MKL_LIBRARIES)
if(NOT LAPACK_LIBRARIES AND NOT MKL_LIBRARIES)
    find_package(LAPACK)
    if(NOT LAPACK_FOUND)
        message(SEND_ERROR "BLAS/LAPACK not found, specify custom installation with \
LAPACK_LIBRARIES or MKL_LIBRARIES variable")
    endif()
endif()

if(ENABLE_SCALAPACK_MPI)
    find_package(MPI)
    if(NOT MPI_Fortran_FOUND)
        message(SEND_ERROR "MPI not found, specify the MPI Fortran compiler with \
MPI_Fortran_COMPILER variable")
    endif()
    if(APPLE)
        # -flat_namespace (Darwin-specific) is causing a crash (seg fault) when the
        # Fortran library is called from Python and one writes into a character
        # variable, but only when the kind is the default one. It causes the
        # written to variable to appear as being four times shorter than it is.
        # Only mention of anything possibly related I could find is at
        #
        #   https://trac.mpich.org/projects/mpich/ticket/1590
        string(REPLACE "-Wl,-flat_namespace" ""
            MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
        string(STRIP "${MPI_Fortran_LINK_FLAGS}" MPI_Fortran_LINK_FLAGS)
    endif()
    if(NOT SCALAPACK_LIBRARIES AND NOT MKL_LIBRARIES)
        find_package(PkgConfig)
        if(PkgConfig_FOUND)
            pkg_search_module(scalapack scalapack scalapack-openmpi)
            if(scalapack_FOUND)
                message(STATUS "ScaLAPACK pkg-config package found, version ${scalapack_VERSION}")
                set(SCALAPACK_LIBRARIES "${scalapack_LDFLAGS}")
            endif()
        endif()
        if(NOT scalapack_FOUND)
            unset(scalapack_FOUND CACHE)
            find_package(scalapack QUIET)
            if(scalapack_FOUND)
                message(STATUS "ScaLAPACK CMake package found in ${scalapack_DIR}")
            endif()
        endif()
        if(NOT scalapack_FOUND)
            message(SEND_ERROR "ScaLAPACK pkg-config or Cmake package not found, \
specify custom installation with SCALAPACK_LIBRARIES or MKL_LIBRARIES variable")
        endif()
    endif()
endif()

if(ENABLE_ELSI)
    if(NOT ELSI_LIBRARIES)
        find_package(elsi 2.0 QUIET)
        if(elsi_FOUND)
            message(STATUS "ELSI CMake package found in ${elsi_DIR}")
        else()
            message(SEND_ERROR "ELSI CMake package not found, specify custom \
installation with ELSI_LIBRARIES variable")
        endif()
    endif()
endif()

add_subdirectory(src)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
