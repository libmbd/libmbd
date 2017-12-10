import sys
import cffi

if sys.platform == 'darwin':
    kwargs = {'extra_link_args': ['-rpath', 'build']}
else:
    kwargs = {'runtime_library_dirs': ['build']}

ffibuilder = cffi.FFI()
ffibuilder.set_source(
    'pymbd._libmbd',
    '#include "mbd.h"',
    libraries=['mbd'],
    include_dirs=['src'],
    library_dirs=['build'],
    **kwargs
)
with open('src/mbd.h') as f:
    ffibuilder.cdef(f.read())
