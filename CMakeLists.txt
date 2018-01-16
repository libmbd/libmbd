cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0042 NEW)

project(pymbd)
enable_language(Fortran)
find_package(LAPACK REQUIRED)
find_package(MPI)

set(WITH_MPI .false.)

set(SRCS
    src/mbd.f90
    src/mbd_c_api.f90
    src/mbd_api.f90
    src/mbd_vdw_param.f90
    src/mbd_helper_dev.f90
    src/mbd_math.f90
    src/mbd_mpi.f90
    src/mbd_common.f90
    src/mbd_linalg.f90
    src/mbd_types.f90
    src/mbd_repulsion.f90
    mbd_build_flags.f90)
if(NOT ${MPI_Fortran_FOUND})
    list(APPEND SRCS
        src/mpi_stubs.f90)
endif()

add_library(mbd SHARED ${SRCS})

target_link_libraries(mbd ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})

if(${MPI_Fortran_FOUND})
    set_source_files_properties(src/mbd_mpi.f90
        PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
    target_link_libraries(mbd ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
    target_include_directories(mbd PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    set(WITH_MPI .true.)
endif()

configure_file(src/mbd_build_flags.in.f90 mbd_build_flags.f90)

install(TARGETS mbd
    LIBRARY DESTINATION lib)

enable_testing()
set(CMAKE_CTEST_COMMAND ctest -V)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

add_executable(tests EXCLUDE_FROM_ALL src/tests.f90)
target_link_libraries(tests mbd)
add_test(tests tests)
add_dependencies(check tests)
