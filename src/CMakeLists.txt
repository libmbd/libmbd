set(SOURCES
    mbd.F90
    mbd_c_api.F90
    mbd_api.f90
    mbd_vdw_param.f90
    mbd_parallel.F90
    mbd_coulomb.f90
    # mbd_density.f90
    mbd_common.f90
    mbd_linalg.F90
    mbd_types.F90
    mbd_defaults.f90
)

add_library(mbd ${SOURCES})

target_link_libraries(mbd PRIVATE ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})

if(ENABLE_SCALAPACK_MPI)
    target_link_libraries(mbd PRIVATE ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
    target_include_directories(mbd PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    if(SCALAPACK_LIBRARIES)
        target_link_libraries(mbd PRIVATE ${SCALAPACK_LIBRARIES})
    else()
        target_link_libraries(mbd PRIVATE scalapack)
    endif()
    foreach(SRC ${SOURCES} mbd_tests.F90 mbd_api_tests.F90)
        if(SRC MATCHES ".*\\.F90")
            set_property(SOURCE ${SRC} APPEND PROPERTY COMPILE_DEFINITIONS WITH_MPI WITH_SCALAPACK)
        endif()
    endforeach()
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
    endif()

    add_dependencies(check mbd_tests mbd_api_tests)
endif()
