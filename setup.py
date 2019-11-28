import os
import sys

from setuptools import Extension, setup

if '--no-ext' in sys.argv:
    sys.argv.remove('--no-ext')
    setup()
    sys.exit()

MBD_BUILDER = 'utils/build-mbd-ext.py'
LIBMBDC = 'pymbd/_libmbd.c'

if os.path.exists(MBD_BUILDER):
    kwargs = {'cffi_modules': ['utils/build-mbd-ext.py:ffibuilder']}
else:
    assert os.path.exists(LIBMBDC)
    ext = Extension('pymbd._libmbd', sources=[LIBMBDC], libraries=['mbd'])
    kwargs = {'ext_modules': [ext]}

# fix incorrect rpath flags on darwin
if True:
    from distutils.unixccompiler import UnixCCompiler

    runtime_library_dir_option = UnixCCompiler.runtime_library_dir_option

    def _runtime_library_dir_option(self, dir):
        if sys.platform[:6] == 'darwin':
            return ['-rpath', dir]
        return runtime_library_dir_option(self, dir)

    UnixCCompiler.runtime_library_dir_option = _runtime_library_dir_option

setup(**kwargs)
