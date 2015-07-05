c_warnings_off = incompatible-pointer-types parentheses-equality shorten-64-to-32 \#warnings
f_warnings_off = maybe-uninitialized
CFLAGS = $(addprefix -Wno-,${c_warnings_off})
FFLAGS = $(addprefix -Wno-,${f_warnings_off})
FVENDOR = gnu95
CVENDOR = unix
FC = mpifort
LDFLAGS = -lscalapack
