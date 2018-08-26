set(SRCS
    mbd.F90
    mbd_api.f90
    mbd_c_api.F90
    mbd_common.f90
    mbd_constants.f90
    mbd_coulomb.f90
    mbd_lapack.f90
    mbd_linalg.F90
    mbd_matrix_type.F90
    mbd_system_type.F90
    mbd_ts.f90
    mbd_vdw_param.f90
    # mbd_density.f90
)

if(ENABLE_SCALAPACK_MPI)
    list(APPEND SRCS mbd_blacs.f90 mbd_scalapack.f90)
endif()

add_library(mbd ${SRCS})

target_link_libraries(mbd PRIVATE ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})

if(ENABLE_SCALAPACK_MPI)
    target_link_libraries(mbd PRIVATE ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
    target_include_directories(mbd PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    if(SCALAPACK_LIBRARIES)
        target_link_libraries(mbd PRIVATE ${SCALAPACK_LIBRARIES})
    else()
        target_link_libraries(mbd PRIVATE scalapack)
    endif()
    set_target_properties(mbd PROPERTIES COMPILE_DEFINITIONS "WITH_MPI;WITH_SCALAPACK")
endif()

set_target_properties(mbd PROPERTIES PUBLIC_HEADER "mbd.h;${CMAKE_CURRENT_BINARY_DIR}/mbd_api.mod")
install(TARGETS mbd
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(CMAKE_TESTING_ENABLED)
    add_executable(mbd_tests EXCLUDE_FROM_ALL mbd_tests.F90)
    target_link_libraries(mbd_tests PRIVATE mbd)
    add_test(mbd_tests mbd_tests)

    add_executable(mbd_api_tests EXCLUDE_FROM_ALL mbd_api_tests.F90)
    target_link_libraries(mbd_api_tests PRIVATE mbd)
    add_test(mbd_api_tests mbd_api_tests)

    if(ENABLE_SCALAPACK_MPI)
        target_include_directories(mbd_tests PRIVATE ${MPI_Fortran_INCLUDE_PATH})
        target_include_directories(mbd_api_tests PRIVATE ${MPI_Fortran_INCLUDE_PATH})
        target_link_libraries(mbd_tests PRIVATE ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
        target_link_libraries(mbd_api_tests PRIVATE ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
        set_target_properties(mbd_tests mbd_api_tests PROPERTIES COMPILE_DEFINITIONS "WITH_MPI;WITH_SCALAPACK")
    endif()

    add_dependencies(check mbd_tests mbd_api_tests)
endif()
