import os
import sys

import cffi

MBD_H = 'src/mbd.h'
LIBMBD_PREFIX = os.environ.get('LIBMBD_PREFIX')

if LIBMBD_PREFIX == '':

    def build(setup_kwargs):
        pass


else:

    # some Conda environments do not add their include dir into default includes
    CONDA_PREFIX = os.environ.get('CONDA_PREFIX')
    if not LIBMBD_PREFIX and CONDA_PREFIX:
        LIBMBD_PREFIX = CONDA_PREFIX

    if LIBMBD_PREFIX:
        ext_kwargs = {
            'include_dirs': [f'{LIBMBD_PREFIX}/include'],
            'library_dirs': [f'{LIBMBD_PREFIX}/lib'],
            'runtime_library_dirs': [f'{LIBMBD_PREFIX}/lib'],
        }
    else:
        ext_kwargs = {}
    ffibuilder = cffi.FFI()
    ffibuilder.set_source(
        'pymbd._libmbd', '#include "mbd.h"', libraries=['mbd'], **ext_kwargs
    )
    with open(MBD_H) as f:
        ffibuilder.cdef(f.read())

    def build(setup_kwargs):
        setup_kwargs.update({'cffi_modules': ['build.py:ffibuilder']})
        if sys.platform[:6] == 'darwin':
            from distutils.unixccompiler import UnixCCompiler

            UnixCCompiler.runtime_library_dir_option = lambda self, dir: ['-rpath', dir]
