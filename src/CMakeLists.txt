set(SRCS
    mbd.F90
    mbd_c_api.F90
    mbd_calc.F90
    mbd_constants.f90
    mbd_coulomb.f90
    mbd_damping.F90
    # mbd_density.f90
    mbd_dipole.F90
    mbd_formulas.f90
    mbd_geom.F90
    mbd_gradients.f90
    mbd_hamiltonian.F90
    mbd_lapack.f90
    mbd_linalg.F90
    mbd_matrix.F90
    mbd_methods.F90
    mbd_rpa.F90
    mbd_scs.f90
    mbd_ts.f90
    mbd_utils.F90
    mbd_vdw_param.f90
)

if(ENABLE_SCALAPACK_MPI)
    list(APPEND SRCS mbd_mpi.F90 mbd_blacs.f90 mbd_scalapack.f90)
endif()

if(ENABLE_ELSI)
    list(APPEND SRCS mbd_elsi.F90)
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
    set_property(TARGET mbd APPEND PROPERTY COMPILE_DEFINITIONS WITH_MPI WITH_SCALAPACK)
endif()

if(ENABLE_ELSI)
    if(ELSI_LIBRARIES)
        target_link_libraries(mbd PRIVATE ${ELSI_LIBRARIES})
    else()
        target_link_libraries(mbd PRIVATE elsi::elsi)
    endif()
    set_property(TARGET mbd APPEND PROPERTY COMPILE_DEFINITIONS WITH_ELSI)
endif()

set_property(TARGET mbd PROPERTY PUBLIC_HEADER mbd.h "${CMAKE_CURRENT_BINARY_DIR}/mbd.mod")
install(TARGETS mbd
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
