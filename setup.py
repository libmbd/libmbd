from setuptools import setup

setup(
    name='pymbd',
    version='0.3',
    description='Many-body dispersion method',
    author='Jan Hermann',
    author_email='dev@hermann.in',
    url='https://github.com/azag0/pymbd',
    packages=['pymbd'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Fortran',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    license='Mozilla Public License 2.0',
    setup_requires=['cffi>=1.0.0'],
    install_requires=['cffi>=1.0.0'],
    cffi_modules=['pymbd_build.py:ffibuilder'],
)
