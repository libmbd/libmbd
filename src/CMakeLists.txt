configure_file(mbd_version.f90.in "${CMAKE_CURRENT_BINARY_DIR}/mbd_version.f90")

add_library(mbd
    mbd.F90
    mbd_constants.f90
    mbd_coulomb.f90
    mbd_damping.F90
    mbd_defaults.f90
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
    mbd_ts.F90
    mbd_utils.F90
    "${CMAKE_CURRENT_BINARY_DIR}/mbd_version.f90"
    mbd_vdw_param.f90
)

if(ENABLE_SCALAPACK_MPI)
    target_sources(mbd PRIVATE mbd_mpi.F90 mbd_blacs.f90 mbd_scalapack.f90)
    if(NOT MPI_Fortran_HAVE_F08_MODULE)
        if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
            set(mismatch_flag "-mismatch")
        elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU"
                AND "${CMAKE_Fortran_COMPILER_VERSION}" VERSION_GREATER_EQUAL "10")
            set(mismatch_flag "-fallow-argument-mismatch -Wno-pedantic")
        endif()
        if(DEFINED mismatch_flag)
            set_source_files_properties(
                SOURCE mbd_geom.F90 mbd_methods.F90 mbd_mpi.F90 mbd_ts.F90
                PROPERTY COMPILE_FLAGS ${mismatch_flag}
            )
        endif()
    endif()
endif()

if(ENABLE_ELSI)
    target_sources(mbd PRIVATE mbd_elsi.F90)
endif()

if(ENABLE_C_API)
    target_sources(mbd PRIVATE mbd_c_api.F90)
endif()

set(moduledir "${CMAKE_CURRENT_BINARY_DIR}/modules")
set(includedir "${CMAKE_INSTALL_INCLUDEDIR}/mbd")
set_target_properties(mbd PROPERTIES Fortran_MODULE_DIRECTORY "${moduledir}")

target_include_directories(mbd
    PRIVATE
        ${moduledir}
    INTERFACE
        $<BUILD_INTERFACE:${moduledir}>
        $<INSTALL_INTERFACE:${includedir}>
)

if(ENABLE_SCALAPACK_MPI)
    target_link_libraries(mbd PRIVATE MPI::MPI_Fortran scalapack)
    set_property(TARGET mbd APPEND PROPERTY COMPILE_DEFINITIONS WITH_MPI WITH_SCALAPACK)
    if(MPI_Fortran_HAVE_F08_MODULE)
        set_property(TARGET mbd APPEND PROPERTY COMPILE_DEFINITIONS WITH_MPIF08)
    endif()
endif()

target_link_libraries(mbd PRIVATE LAPACK::LAPACK)

if(ENABLE_ELSI)
    target_link_libraries(mbd PRIVATE elsi::elsi)
    set_property(TARGET mbd APPEND PROPERTY COMPILE_DEFINITIONS WITH_ELSI)
endif()

if(ENABLE_C_API)
    set_property(TARGET mbd PROPERTY PUBLIC_HEADER mbd.h)
endif()

if(CMAKE_INSTALL_LIBDIR)
    add_library(Mbd INTERFACE)
    target_link_libraries(Mbd INTERFACE mbd)
    install(TARGETS mbd Mbd EXPORT MbdConfig
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        PUBLIC_HEADER DESTINATION "${includedir}"
    )
    install(DIRECTORY "${moduledir}/" DESTINATION "${includedir}")
    install(EXPORT MbdConfig NAMESPACE Mbd:: DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/mbd")
endif()
