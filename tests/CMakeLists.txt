set(TESTS
    grad/T_bare_deriv
    grad/T_GG_deriv_expl
    grad/T_GG_deriv_impl
    grad/T_erfc_deriv_expl
    grad/T_fermi_deriv_impl
    grad/mbd_deriv_expl
    grad/mbd_ewald_deriv_expl
    grad/mbd_ewald_deriv_stress
    grad/scs_deriv_expl
    grad/scs_ewald_deriv_expl
    grad/scs_ewald_deriv_stress
    grad/scs_deriv_impl_alpha
    grad/scs_ewald_deriv_impl_alpha
    grad/scs_deriv_impl_vdw
    grad/scs_ewald_deriv_impl_vdw
    grad/mbd_deriv_impl_alpha
    grad/mbd_ewald_deriv_impl_alpha
    grad/mbd_deriv_impl_omega
    grad/mbd_ewald_deriv_impl_omega
    grad/mbd_deriv_impl_vdw
    grad/mbd_ewald_deriv_impl_vdw
    grad/mbd_ewald_deriv_impl_q
    grad/mbd_rsscs_deriv_expl
    grad/mbd_rsscs_ewald_deriv_expl
    grad/mbd_rsscs_ewald_deriv_stress
    grad/mbd_rsscs_deriv_impl_alpha
    grad/mbd_rsscs_deriv_impl_C6
    grad/mbd_rsscs_deriv_impl_vdw
    grad/mbd_rsscs_ewald_deriv_impl_alpha
    grad/mbd_rsscs_ewald_deriv_impl_C6
    grad/mbd_rsscs_ewald_deriv_impl_vdw
    api/energy
    api/energy_from_ratios
    api/gradients
    api/exception
)

add_executable(mbd_grad_tests mbd_grad_tests.F90 mbd_grad_test_cases.f90)
add_executable(mbd_api_tests mbd_api_tests.F90)
foreach(TEST_APP mbd_grad_tests mbd_api_tests)
    target_include_directories(${TEST_APP} PRIVATE $<TARGET_FILE_DIR:mbd>)
    target_link_libraries(${TEST_APP} PRIVATE mbd)
    if(ENABLE_SCALAPACK_MPI)
        target_include_directories(${TEST_APP} PRIVATE ${MPI_Fortran_INCLUDE_PATH})
        target_link_libraries(${TEST_APP} PRIVATE ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
        set_property(TARGET ${TEST_APP} APPEND PROPERTY COMPILE_DEFINITIONS WITH_MPI WITH_SCALAPACK)
    endif()
    if(ENABLE_ELSI)
        set_property(TARGET ${TEST_APP} APPEND PROPERTY COMPILE_DEFINITIONS WITH_ELSI)
    endif()
endforeach()

foreach(TEST ${TESTS})
    string(REPLACE "/" ";" TEST_APP_CASE ${TEST})
    list(GET TEST_APP_CASE 0 TEST_APP)
    list(GET TEST_APP_CASE 1 TEST_CASE)
    add_test(NAME "${TEST_APP}/${TEST_CASE}" COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG}
        ${MPIEXEC_MAX_NUMPROCS}
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:mbd_${TEST_APP}_tests>
        ${TEST_CASE}
        ${MPIEXEC_POSTFLAGS}
    )
endforeach()
