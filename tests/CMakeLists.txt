add_executable(mbd_grad_tests mbd_grad_tests.F90 mbd_grad_test_cases.F90)
add_executable(mbd_api_tests mbd_api_tests.F90)
foreach(TEST_APP mbd_grad_tests mbd_api_tests)
    target_include_directories(${TEST_APP} PRIVATE $<TARGET_FILE_DIR:mbd>)
    target_link_libraries(${TEST_APP} PRIVATE mbd)
    if(ENABLE_SCALAPACK_MPI)
        target_link_libraries(${TEST_APP} PRIVATE MPI::MPI_Fortran)
        set_property(TARGET ${TEST_APP} APPEND PROPERTY COMPILE_DEFINITIONS WITH_MPI WITH_SCALAPACK)
    endif()
    if(ENABLE_ELSI)
        set_property(TARGET ${TEST_APP} APPEND PROPERTY COMPILE_DEFINITIONS WITH_ELSI)
    endif()
endforeach()

execute_process(
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/collect-mbd-tests.py
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE TESTS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message("${TESTS}")
if(ENABLE_SCALAPACK_MPI)
    set(MPI_NODES $ENV{MPI_NODES})
    if(NOT "${MPI_NODES}")
        set(MPI_NODES ${MPIEXEC_MAX_NUMPROCS})
    endif()
    set(MPIEXEC_EXTRA_FLAGS $ENV{MPIEXEC_EXTRA_FLAGS})
endif()
foreach(TEST ${TESTS})
    message("${TEST}")
    string(REPLACE "/" ";" TEST_APP_CASE ${TEST})
    message("${TEST}")
    list(GET TEST_APP_CASE 0 TEST_APP)
    list(GET TEST_APP_CASE 1 TEST_CASE)
    message("${TEST_APP} ${TEST_APP}")
    add_test(NAME "${TEST_APP}/${TEST_CASE}" COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_EXTRA_FLAGS}
        ${MPIEXEC_NUMPROC_FLAG}
        ${MPI_NODES}
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:mbd_${TEST_APP}_tests>
        ${TEST_CASE}
        ${MPIEXEC_POSTFLAGS}
    )
endforeach()
