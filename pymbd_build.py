import cffi

ffibuilder = cffi.FFI()
ffibuilder.set_source(
    'pymbd._libmbd',
    f'#include "mbd.h"',
    libraries=['mbd'],
)
with open('src/mbd.h') as f:
    ffibuilder.cdef(f.read())
