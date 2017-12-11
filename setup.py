import os
import sys
import tempfile
import cffi
import shutil
from setuptools import setup

library_dirs = ['build'] if os.path.exists('build') else []
sources = [
    'src/mpi_stubs.f90',
    'src/mbd_common.f90',
    'src/mbd_interface.f90',
    'src/mbd_linalg.f90',
    'src/mbd.f90',
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
        )
    except LinkError:
        return False
    else:
        return True
    finally:
        shutil.rmtree(tmpdir)


if libmbd_exists():
    kwargs = {'libraries': ['mbd']}
    if sys.platform == 'darwin':
        kwargs['extra_link_args'] = ['-rpath', 'build']
    else:
        kwargs['runtime_library_dirs'] = ['build']
else:
    from numpy.distutils.core import setup  # noqa
    from numpy.distutils.system_info import get_info
    kwargs = {'libraries': [('mbd', {'sources': sources, 'language': 'f90'})]}
    for arg, val in get_info('lapack_opt', 2).items():
        if arg == 'libraries':
            kwargs['libraries'].extend(val)
        else:
            kwargs[arg] = val


ffibuilder = cffi.FFI()
ffibuilder.set_source(
    'pymbd._libmbd',
    '#include "mbd.h"',
    include_dirs=['src'],
    library_dirs=library_dirs,
    **kwargs
)
with open('src/mbd.h') as f:
    ffibuilder.cdef(f.read())


setup(
    name='pymbd',
    version='0.3.1',
    description='Many-body dispersion method',
    author='Jan Hermann',
    author_email='dev@hermann.in',
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
