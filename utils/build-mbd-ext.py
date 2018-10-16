#!/usr/bin/env python3
import cffi

ffibuilder = cffi.FFI()
ffibuilder.set_source(
    'pymbd._libmbd',
    '#include "mbd.h"',
    libraries=['mbd']
)
with open('src/mbd.h') as f:
    ffibuilder.cdef(f.read())

if __name__ == "__main__":
    ffibuilder.distutils_extension('.')
