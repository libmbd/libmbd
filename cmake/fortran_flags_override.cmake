if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS_INIT
        "-fall-intrinsics -std=f2008ts -pedantic -Wall -Wcharacter-truncation -Wimplicit-procedure -Wextra -Wno-maybe-uninitialized")
    set(CMAKE_Fortran_FLAGS_DEBUG_INIT
        "-Og -fcheck=all -fno-check-array-temporaries")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
    # Covers both ifort ("Intel") and ifx ("IntelLLVM")
    set(CMAKE_Fortran_FLAGS_INIT "-warn all")
    set(CMAKE_Fortran_FLAGS_DEBUG_INIT
        "-O0 -g -traceback -check all,noarg_temp_created")
endif()
