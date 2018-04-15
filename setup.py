import os
import sys
import tempfile
import cffi
import shutil
from setuptools import setup

library_dirs = ['build'] if os.path.exists('build') else []
sources = [
    'src/mbd_common.f90',
    'src/mbd_defaults.f90',
    'src/mbd_types.f90',
    'src/mbd_linalg.f90',
    'src/mbd.F90',
    'src/mbd_c_api.f90',
]


def libmbd_exists():
    import distutils.sysconfig
    import distutils.ccompiler
    from distutils.errors import LinkError

    tmpdir = tempfile.mkdtemp()
    src = os.path.join(tmpdir, 'test.c')
    with open(src, 'w') as f:
        f.write('int main() {}')
    compiler = distutils.ccompiler.new_compiler()
    distutils.sysconfig.customize_compiler(compiler)
    try:
        compiler.link_executable(
            compiler.compile([src], output_dir=tmpdir),
            'test',
            libraries=['mbd'],
            library_dirs=library_dirs,
            output_dir=tmpdir
        )
    except LinkError:
        return False
    else:
        return True
    finally:
        shutil.rmtree(tmpdir)


def setup_mpi():
    try:
        import mpi4py
    except ImportError:
        sources.insert(0, 'src/mpi_stubs.f90')
        return
    # patch find_executables to insert the MPI compiler before FC or
    # --f90exec is checked
    from numpy.distutils.fcompiler import FCompiler
    from functools import wraps
    _find_executables = FCompiler.find_executables
    @wraps(FCompiler.find_executables)  # noqa
    def find_executables(self):
        self.executables['compiler_f90'][0] = mpi4py.get_config()['mpif90']
        return _find_executables(self)
    FCompiler.find_executables = find_executables


def update_dict(dct, update):
    for key, val in update.items():
        if key in dct:
            dct[key].extend(val)
        else:
            dct[key] = val


if libmbd_exists():
    kwargs = {'libraries': ['mbd']}
    if sys.platform == 'darwin':
        kwargs['extra_link_args'] = ['-rpath', os.path.realpath('build')]
    else:
        kwargs['runtime_library_dirs'] = [os.path.realpath('build')]
else:
    from numpy.distutils.core import setup  # noqa
    from numpy.distutils.system_info import get_info
    setup_mpi()
    kwargs = {'libraries': [('mbd', {'sources': sources, 'language': 'f90'})]}
    update_dict(kwargs, get_info('lapack_opt', 2))


update_dict(kwargs, {
    'include_dirs': ['src'],
    'library_dirs': library_dirs
})
ffibuilder = cffi.FFI()
ffibuilder.set_source(
    'pymbd._libmbd',
    '#include "mbd.h"',
    **kwargs
)
with open('src/mbd.h') as f:
    ffibuilder.cdef(f.read())


setup(
    name='pymbd',
    version='0.3.2',
    description='Many-body dispersion method',
    author='Jan Hermann',
    author_email='dev@janhermann.cz',
    url='https://github.com/azag0/pymbd',
    packages=['pymbd'],
    package_data={'pymbd': ['vdw-params.csv']},
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Fortran',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    license='Mozilla Public License 2.0',
    install_requires=['cffi', 'numpy'],
    ext_modules=[ffibuilder.distutils_extension()],
)
