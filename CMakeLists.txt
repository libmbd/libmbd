cmake_minimum_required(VERSION 3.1)
cmake_policy(SET CMP0042 NEW)

project(pymbd)
enable_language(Fortran)
find_package(LAPACK REQUIRED)
find_package(MPI)

add_library(mbd SHARED
    src/mbd.f90
    src/mbd_c_api.f90
    src/mbd_helper_dev.f90
    src/mbd_math.f90
    src/mbd_helper.f90
    src/mbd_interface.f90
    src/mbd_repulsion.f90)

target_link_libraries(mbd ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})

if(${MPI_Fortran_FOUND})
    set_source_files_properties(src/mbd_interface.f90
        PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
    target_link_libraries(mbd ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
    target_include_directories(mbd PRIVATE ${MPI_Fortran_INCLUDE_PATH})
else()
    target_sources(mbd PRIVATE
        src/mpi_stubs.f90)
endif()
