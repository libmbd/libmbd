cmake_minimum_required(VERSION 3.14)
cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0077 NEW)

set(CMAKE_USER_MAKE_RULES_OVERRIDE_Fortran
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/fortran_flags_override.cmake")

project(Libmbd
    DESCRIPTION "Many-body dispersion library"
    LANGUAGES Fortran C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(LibmbdVersion)
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

find_package(LAPACK REQUIRED)
if(CMAKE_VERSION VERSION_LESS 3.18)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    target_link_libraries(scalapack INTERFACE "${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}")
endif()

if(ENABLE_SCALAPACK_MPI)
    find_package(MPI REQUIRED)
    if(APPLE)
        # -flat_namespace (Darwin-specific) is causing a crash (seg fault) when the
        # Fortran library is called from Python and one writes into a character
        # variable, but only when the kind is the default one. It causes the
        # written to variable to appear as being four times shorter than it is.
        # Only mention of anything possibly related I could find is at
        #
        #   https://trac.mpich.org/projects/mpich/ticket/1590
        get_target_property(_VALUE MPI::MPI_Fortran INTERFACE_LINK_OPTIONS)
        string(REGEX REPLACE "-Wl\\$<COMMA>-flat_namespace ?" "" _VALUE "${_VALUE}")
        set_target_properties(MPI::MPI_Fortran PROPERTIES INTERFACE_LINK_OPTIONS "${_VALUE}")
    endif()
    if(NOT DEFINED SCALAPACK_LIBRARY)
        find_package(scalapack QUIET)
        if(scalapack_FOUND)
            message(STATUS "ScaLAPACK CMake package found in ${scalapack_DIR}")
        else()
            find_package(PkgConfig)
            if(PkgConfig_FOUND)
                unset(scalapack_FOUND CACHE)
                pkg_search_module(scalapack scalapack scalapack-openmpi)
                if(scalapack_FOUND)
                    message(STATUS "ScaLAPACK pkg-config package found, version ${scalapack_VERSION}")
                    set(SCALAPACK_LIBRARY "${scalapack_LDFLAGS}")
                endif()
            endif()
        endif()
        if(NOT scalapack_FOUND)
            message(SEND_ERROR "ScaLAPACK Cmake or pkg-config package not found, \
specify custom installation with SCALAPACK_LIBRARY")
        endif()
    endif()
    if(DEFINED SCALAPACK_LIBRARY)
        add_library(scalapack INTERFACE IMPORTED)
        target_link_libraries(scalapack INTERFACE "${SCALAPACK_LIBRARY}")
    endif()
endif()

if(ENABLE_ELSI)
    if(NOT DEFINED ELSI_LIBRARY)
        find_package(elsi 2.0 QUIET)
        if(elsi_FOUND)
            message(STATUS "ELSI CMake package found in ${elsi_DIR}")
        else()
            message(SEND_ERROR "ELSI CMake package not found, specify custom \
installation with ELSI_LIBRARY variable")
        endif()
    else()
        add_library(elsi::elsi INTERFACE IMPORTED)
        target_link_libraries(elsi::elsi INTERFACE "${ELSI_LIBRARY}")
    endif()
endif()

add_subdirectory(src)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
