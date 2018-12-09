if(CMAKE_TESTING_ENABLED)
    add_executable(mbd_unit_tests EXCLUDE_FROM_ALL mbd_unit_tests.F90)
    target_include_directories(mbd_unit_tests PRIVATE $<TARGET_FILE_DIR:mbd>)
    target_link_libraries(mbd_unit_tests PRIVATE mbd)
    add_test(mbd_unit_tests mbd_unit_tests)

    add_executable(mbd_api_tests EXCLUDE_FROM_ALL mbd_api_tests.F90)
    target_include_directories(mbd_api_tests PRIVATE $<TARGET_FILE_DIR:mbd>)
    target_link_libraries(mbd_api_tests PRIVATE mbd)
    add_test(mbd_api_tests mbd_api_tests)

    if(ENABLE_SCALAPACK_MPI)
        target_include_directories(mbd_unit_tests PRIVATE ${MPI_Fortran_INCLUDE_PATH})
        target_include_directories(mbd_api_tests PRIVATE ${MPI_Fortran_INCLUDE_PATH})
        target_link_libraries(mbd_unit_tests PRIVATE ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
        target_link_libraries(mbd_api_tests PRIVATE ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
        set_target_properties(mbd_unit_tests mbd_api_tests PROPERTIES COMPILE_DEFINITIONS "WITH_MPI;WITH_SCALAPACK")
    endif()

    add_dependencies(check mbd_unit_tests mbd_api_tests)
endif()

