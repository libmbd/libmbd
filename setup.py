import os
import sys
import shutil
import tempfile
if sys.version_info[0] > 2:
    from configparser import ConfigParser
else:
    from ConfigParser import ConfigParser  # noqa

from setuptools import setup, Extension  # noqa


def update_dict(dct, update):
    for key, val in update.items():
        if key in dct:
            dct[key].extend(val)
        else:
            dct[key] = val


def libmbd_exists():
    import distutils.sysconfig
    import distutils.ccompiler
    from distutils.errors import LinkError

    LIBMBD = os.environ.get('LIBMBD', 'build/src')
    library_dirs = [LIBMBD] if os.path.exists(LIBMBD) else []
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
        return None
    else:
        return library_dirs
    finally:
        shutil.rmtree(tmpdir)


conf = ConfigParser()
conf.read('setup.cfg')

library_dirs = libmbd_exists()
if library_dirs is not None:
    ext_kwargs = {'libraries': ['mbd']}
    if library_dirs:
        if sys.platform == 'darwin':
            ext_kwargs['extra_link_args'] = ['-rpath', os.path.realpath(library_dirs[0])]
        else:
            ext_kwargs['runtime_library_dirs'] = [os.path.realpath(library_dirs[0])]
else:
    from numpy.distutils.core import setup, Extension  # noqa
    from numpy.distutils.system_info import get_info

    sources = ['src/' + name for name in conf.get('libmbd:info', 'sources').split()]
    ext_kwargs = {'libraries': [('mbd', {'sources': sources, 'language': 'f90'})]}
    update_dict(ext_kwargs, get_info('lapack_opt', 2))
    library_dirs = ['build']

if library_dirs:
    update_dict(ext_kwargs, {
        'include_dirs': ['src'],
        'library_dirs': library_dirs
    })

LIBMBDC = 'pymbd/_libmbd.c'
EXTPATH = 'pymbd._libmbd'
try:
    import cffi
except ImportError:
    if not os.path.exists(LIBMBDC):
        raise
    ext = Extension(EXTPATH, sources=[LIBMBDC], **ext_kwargs)
else:
    ffibuilder = cffi.FFI()
    ffibuilder.set_source(EXTPATH, '#include "mbd.h"', **ext_kwargs)
    with open('src/mbd.h') as f:
        ffibuilder.cdef(f.read())
    ext = ffibuilder.distutils_extension(tmpdir='.')

with open('README.md') as f:
    long_description = f.read()

setup(
    name='pymbd',
    version=conf.get('libmbd:info', 'version'),
    description=conf.get('libmbd:info', 'description'),
    author=conf.get('libmbd:info', 'author'),
    author_email='dev@janhermann.cz',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/azag0/libmbd',
    packages=['pymbd'],
    package_data={'pymbd': ['vdw-params.csv']},
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Fortran',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    install_requires=[
        'cffi',
        'numpy',
        'scipy'
    ],
    ext_modules=[ext],
    extras_require={
        'mpi': ['mpi4py'],
    },
    python_requires='>=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*',
)
