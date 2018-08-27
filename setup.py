import os
import sys
import tempfile
import cffi
import shutil
from setuptools import setup
if sys.version_info[0] > 2:
    from configparser import ConfigParser
else:
    from ConfigParser import ConfigParser  # noqa

blddir = os.environ.get('MBDBLDDIR', 'build/src')
library_dirs = [blddir] if os.path.exists(blddir) else []

conf = ConfigParser()
conf.read('setup.cfg')
sources = ['src/' + name for name in conf.get('info', 'sources').split()]


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


def update_dict(dct, update):
    for key, val in update.items():
        if key in dct:
            dct[key].extend(val)
        else:
            dct[key] = val


if libmbd_exists():
    ext_args = {'libraries': ['mbd']}
    if library_dirs:
        if sys.platform == 'darwin':
            ext_args['extra_link_args'] = ['-rpath', os.path.realpath(blddir)]
        else:
            ext_args['runtime_library_dirs'] = [os.path.realpath(blddir)]
else:
    from numpy.distutils.core import setup  # noqa
    from numpy.distutils.system_info import get_info
    ext_args = {'libraries': [('mbd', {'sources': sources, 'language': 'f90'})]}
    update_dict(ext_args, get_info('lapack_opt', 2))
    library_dirs = ['build']


if library_dirs:
    update_dict(ext_args, {
        'include_dirs': ['src'],
        'library_dirs': library_dirs
    })
ffibuilder = cffi.FFI()
ffibuilder.set_source(
    'pymbd._libmbd',
    '#include "mbd.h"',
    **ext_args
)
with open('src/mbd.h') as f:
    ffibuilder.cdef(f.read())


setup(
    name='pymbd',
    version='0.4.0a1',
    description='Many-body dispersion method',
    author='Jan Hermann',
    author_email='dev@janhermann.cz',
    url='https://github.com/azag0/libmbd',
    packages=['pymbd'],
    package_data={'pymbd': ['vdw-params.csv']},
    classifiers=[
        'Development Status :: 3 - Alpha',
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
    install_requires=['cffi', 'numpy', 'scipy'],
    ext_modules=[ffibuilder.distutils_extension()],
    extras_require={
        'MPI': ['mpi4py'],
    },
)
