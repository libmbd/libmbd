import sys
import cffi
from setuptools import setup

ffibuilder = cffi.FFI()
if sys.platform == 'darwin':
    kwargs = {'extra_link_args': ['-rpath', 'build']}
else:
    kwargs = {'runtime_library_dirs': ['build']}
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

setup(
    name='pymbd',
    version='0.3',
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
    ext_modules=[ffibuilder.distutils_extension()]
)
