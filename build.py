import os
import sys

import cffi

MBD_H = 'src/mbd.h'
LIBMBD_PREFIX = os.environ.get('LIBMBD_PREFIX')

if LIBMBD_PREFIX == '':

    def build(setup_kwargs):
        pass


else:
    if LIBMBD_PREFIX:
        ext_kwargs = {
            'include_dirs': ['{}/include'.format(LIBMBD_PREFIX)],
            'library_dirs': ['{}/lib'.format(LIBMBD_PREFIX)],
            'runtime_library_dirs': ['{}/lib'.format(LIBMBD_PREFIX)],
        }
    else:
        ext_kwargs = {}
    ffibuilder = cffi.FFI()
    ffibuilder.set_source(
        'pymbd._libmbd', '#include "mbd.h"', libraries=['mbd'], **ext_kwargs,
    )
    with open(MBD_H) as f:
        ffibuilder.cdef(f.read())

    def build(setup_kwargs):
        setup_kwargs.update({'cffi_modules': ['build.py:ffibuilder']})
        if sys.platform[:6] == 'darwin':
            from distutils.unixccompiler import UnixCCompiler

            UnixCCompiler.runtime_library_dir_option = lambda self, dir: ['-rpath', dir]
