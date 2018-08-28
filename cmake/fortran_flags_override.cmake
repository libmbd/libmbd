if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS_DEBUG_INIT
        "-Og -fcheck=all -Wall -Wcharacter-truncation -Wimplicit-procedure -Wextra -std=f2003 -pedantic -fall-intrinsics -Wno-maybe-uninitialized")
endif()
