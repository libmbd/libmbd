cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0042 NEW)

project(libmbd)
enable_language(Fortran)
find_package(LAPACK REQUIRED)

set(SRCS
    src/mbd.F90
    src/mbd_c_api.F90
    src/mbd_api.f90
    src/mbd_vdw_param.f90
    src/mbd_parallel.F90
    # src/mbd_math.f90
    # src/mbd_density.f90
    # src/mbd_repulsion.f90
    src/mbd_common.f90
    src/mbd_linalg.F90
    src/mbd_types.f90
    src/mbd_defaults.f90
)

add_library(mbd SHARED ${SRCS})

target_link_libraries(mbd ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})

install(TARGETS mbd
    LIBRARY DESTINATION lib
)

enable_testing()
set(CMAKE_CTEST_COMMAND ctest -V)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

add_executable(mbd_tests EXCLUDE_FROM_ALL src/mbd_tests.F90)
target_link_libraries(mbd_tests mbd)
add_test(mbd_tests mbd_tests)

add_executable(mbd_api_tests EXCLUDE_FROM_ALL src/mbd_api_tests.F90)
target_link_libraries(mbd_api_tests mbd)
add_test(mbd_api_tests mbd_api_tests)

add_dependencies(check mbd_tests mbd_api_tests)

find_package(MPI)
if(MPI_Fortran_FOUND)
    if(NOT SCALAPACK)
        message(FATAL_ERROR "MPI found, but no ScaLAPACK specified")
    endif()
    target_link_libraries(mbd ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES} ${SCALAPACK})
    target_include_directories(mbd PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    set_property(
        SOURCE src/mbd_parallel.F90 src/mbd_c_api.F90 src/mbd_linalg.F90
            src/mbd.F90 src/mbd_tests.F90 src/mbd_api_tests.F90
        APPEND PROPERTY COMPILE_DEFINITIONS WITH_SCALAPACK
    )
endif()
