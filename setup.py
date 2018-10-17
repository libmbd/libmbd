import os
import sys

from setuptools import setup, Extension

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

# fix https://bugs.python.org/issue1703178
if sys.version_info[0] < 3:
    from distutils.command.build_ext import build_ext
    finalize_options = build_ext.finalize_options

    def _finalize_options(self):
        finalize_options(self)
        if isinstance(self.link_objects, str):
            self.link_objects = self.link_objects.split()
    build_ext.finalize_options = _finalize_options

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
