# `pymbd` — Many-body dispersion method

Python 2/3 package for calculating [many-body dispersion](http://dx.doi.org/10.1063/1.4865104) energies.

Most functionality is implemented in the Fortran module in `src/mbd.f90`.

## Installation

Installation of Pymbd requires Numpy with properly configured Lapack libraries.

On macOS, this requirement is almost always satisfied with whatever Python with simple

```
python/python3 -m pip install numpy
```

This can be followed with

```
python/python3 -m pip install git+https://github.com/azag0/mbd.git
```

You can run tests with

```
python/python3 -m unittest pymbd.tests -v
```

On Linux, the pip-installed Numpy sometimes doesn't report the correct Lapack [1]. To avoid this issue, use the Anaconda Python with the no-MKL Numpy (if you have Anaconda already installed, [make sure](https://www.continuum.io/blog/developer-blog/anaconda-25-release-now-mkl-optimizations) that you don't use the MKL Numpy):

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
[...]
<path to anaconda>/bin/conda install nomkl numpy
```

To install Pymbd, run

```
<path to anaconda>/bin/pip install git+https://github.com/azag0/mbd.git
```

And finally test with

```
<path to anaconda>/bin/python -m unittest pymbd.tests -v
```

[1]: What really matters is whether `numpy.distutils.system_info.get_info('lapack_opt')` returns the correct information. On macOS, it always returns the Accelerate framework. With Anaconda and no-MKL Numpy, it always returns the internal openBLAS libraries. But as long as it returns a valid Lapack library, you can use whatever Python installation.

## Usage

Pymbd doesn't have any input files, it is called directly from Python scripts. 

For documentation, see the [basic examples](http://nbviewer.jupyter.org/github/azag0/mbd/blob/master/examples/basic.ipynb) or the [more advanced use](http://nbviewer.jupyter.org/github/azag0/mbd/blob/master/examples/advanced.ipynb).
