import setuptools  # noqa
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension


metadata = dict(
    name='pymbd',
    version='0.2',
    description='Many-body dispersion method',
    long_description='See README.md for details.',
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
        'Operating System :: Unix',
        'Programming Language :: Fortran',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    license='Mozilla Public License 2.0',
    install_requires=['numpy'],
    test_suite='pymbd.tests',
)


def get_setup():
    kwargs = metadata.copy()
    kwargs['ext_modules'] = [
        get_mbd_extension(),
    ]
    kwargs['cmdclass'] = {'build_ext': get_build_ext()}
    return kwargs


def get_mbd_extension():
    sources = ['src/mbd_interface.f90', 'src/mbd_helper.f90']
    if MPI_STUBS:
        sources.insert(0, MPI_STUBS)
    kwargs = dict(
        name='pymbd.lib',
        sources=['src/mbd.f90'],
        libraries=[
            ('mbdlib', dict(
                sources=sources,
                language='f90',
            ))
        ],
    )
    mod_lapack(kwargs)
    return Extension(**kwargs)


def get_build_ext():
    # patch build_ext to compile everything in the top-level temporary directory
    from numpy.distutils.command.build_ext import build_ext

    class my_build_ext(build_ext):
        def get_ext_filename(self, ext_name):
            filename = build_ext.get_ext_filename(self, ext_name)
            return filename.split('/')[-1]

    return my_build_ext


def mod_lapack(kwargs):
    from numpy.distutils.system_info import get_info

    for arg, val in get_info('lapack_opt', 2).items():
        if arg == 'libraries':
            kwargs['libraries'].extend(val)
        else:
            kwargs[arg] = val


MPI_STUBS = None


def setup_mpi():
    global MPI_STUBS
    try:
        import mpi4py
    except ImportError:
        MPI_STUBS = 'src/mpi_stubs.f90'
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


setup_mpi()
setup(**get_setup())
