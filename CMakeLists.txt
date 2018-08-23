cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0042 NEW)

set(CMAKE_USER_MAKE_RULES_OVERRIDE_Fortran "${CMAKE_SOURCE_DIR}/cmake/fortran_flags_override.cmake")

project(libmbd LANGUAGES Fortran)

option(ENABLE_SCALAPACK "Enable parallelisation with ScaLAPACK/MPI")

set(DEFAULT_BUILD_TYPE "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(DEFAULT_BUILD_TYPE "Debug")
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to ${DEFAULT_BUILD_TYPE} as none was set")
    set(CMAKE_BUILD_TYPE ${DEFAULT_BUILD_TYPE})
endif()

include(GNUInstallDirs)

set(SOURCES
    src/mbd.F90
    src/mbd_c_api.F90
    src/mbd_api.f90
    src/mbd_vdw_param.f90
    src/mbd_parallel.F90
    src/mbd_coulomb.f90
    # src/mbd_density.f90
    src/mbd_common.f90
    src/mbd_linalg.F90
    src/mbd_types.F90
    src/mbd_defaults.f90
)

add_library(mbd SHARED ${SOURCES})

if(NOT LAPACK_LIBRARIES)
    find_package(LAPACK)
    if(NOT LAPACK_FOUND)
        message(SEND_ERROR "BLAS/LAPACK not found, specify custom installation with LAPACK_LIBRARIES variable")
    endif()
endif()
target_link_libraries(mbd ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})

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
    target_link_libraries(mbd ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
    target_include_directories(mbd PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    if(SCALAPACK_LIBRARIES)
        target_link_libraries(mbd ${SCALAPACK_LIBRARIES})
    else()
        find_package(Scalapack QUIET)
        if(Scalapack_FOUND)
            message(STATUS "ScaLAPACK CMake package found in ${Scalapack_DIR}")
        else()
            message(SEND_ERROR "ScaLAPACK CMake package not found, specify custom installation with SCALAPACK_LIBRARIES variable")
        endif()
        target_link_libraries(mbd scalapack)
    endif()
    foreach(SRC ${SOURCES})
        if(SRC MATCHES ".*\\.F90")
            set_property(SOURCE ${SRC} APPEND PROPERTY COMPILE_DEFINITIONS WITH_SCALAPACK)
        endif()
    endforeach()
endif()

set_target_properties(mbd PROPERTIES PUBLIC_HEADER "src/mbd.h;${PROJECT_BINARY_DIR}/mbd_api.mod")
install(TARGETS mbd
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

enable_testing()
set(CMAKE_CTEST_COMMAND ctest -V)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

add_executable(mbd_tests EXCLUDE_FROM_ALL src/mbd_tests.F90)
target_link_libraries(mbd_tests mbd)
if(ENABLE_SCALAPACK)
    target_include_directories(mbd_tests PRIVATE ${MPI_Fortran_INCLUDE_PATH})
endif()
add_test(mbd_tests mbd_tests)

add_executable(mbd_api_tests EXCLUDE_FROM_ALL src/mbd_api_tests.F90)
target_link_libraries(mbd_api_tests mbd)
add_test(mbd_api_tests mbd_api_tests)

add_dependencies(check mbd_tests mbd_api_tests)
