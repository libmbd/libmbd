# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.10.3] - 2020-11-13

## Fixed

- Compatibility with NAG compiler

## [0.10.2] - 2020-10-22

### Added

- Search for Scalapack with pkg-config

## [0.10.1] - 2020-10-08

### Fixed

- Pymbd installation

## [0.10.0] - 2020-10-07

### Added

- Version checking of Libmbd from Pymbd
- Support for IBM Fortran compiler by circumventing ICE

### Removed

- Python 3.5 support

## [0.9.3] - 2020-08-02

### Added

- Simple `python -m pymbd` installation check

### Fixed

- Missing `cffi` dependency of Pymbd

## [0.9.2] - 2020-08-01

### Added

- Numerical gradients for TS
- Support for Numpy>1.15

### Fixed

- Compiler error with GCC 10
- Compatibility with Cmake 3.1 when part of a superbuild
- Cmake crash on some platforms

## [0.9.1] - 2020-06-26

### Fixed

- Incorrect initialization of MBD-NL damping parameters from XC functional
- Support PGI 2019 compiler

## [0.9.0] - 2020-06-23

### Added

- MBD-NL damping parameters
- Export of Cmake packages
- ENABLE_C_API build option

- Improved default Scalapack block size

### Removed

- Python 2 support
- Support for Cmake <3.1

## [0.8.0] - 2019-10-30

Minor additions and fixes.

## [0.7.0] - 2019-10-06

### Added

- Optional rescaling of RPA eigenvalues as in [10.1021/acs.jctc.6b00925](http://dx.doi.org/10.1021/acs.jctc.6b00925)

## [0.6.0] - 2019-04-23

### Added

- C/Python API for getting eigenvalues and eigenvectors for crystals
- C/Python API for custom k-point grids

### Fixed

- Support Hessian evaluation with Tensorflow implementation

## [0.5.0] - 2019-03-01

### Changed

- Python/Fortran/C API changed

## [0.4.3] - 2019-02-28

### Added

- Fortran/Python API for RPA evaluation

### Fixed

- Evaluation of RPA orders

## [0.4.2] - 2019-01-20

### Added

- Optional parallelization over k-points

### Fixed

- Ifort compiler bug
- `WITH_MPIFH` build

## [0.4.1] - 2019-01-15

### Changed

- Numpy requirement restricted to <=1.15

## [0.4.0] - 2019-01-13

Completely reworked.

### Added

- Analytical gradients including lattice-vector derivatives.
- Scalapack parallelization of all calculations.

[unreleased]: https://github.com/jhrmnn/libmbd/compare/0.10.3...HEAD
[0.10.3]: https://github.com/jhrmnn/libmbd/compare/0.10.2...0.10.3
[0.10.2]: https://github.com/jhrmnn/libmbd/compare/0.10.1...0.10.2
[0.10.1]: https://github.com/jhrmnn/libmbd/compare/0.10.0...0.10.1
[0.10.0]: https://github.com/jhrmnn/libmbd/compare/0.9.3...0.10.0
[0.9.3]: https://github.com/jhrmnn/libmbd/compare/0.9.2...0.9.3
[0.9.2]: https://github.com/jhrmnn/libmbd/compare/0.9.1...0.9.2
[0.9.1]: https://github.com/jhrmnn/libmbd/compare/0.9.0...0.9.1
[0.9.0]: https://github.com/jhrmnn/libmbd/compare/0.8.0...0.9.0
[0.8.0]: https://github.com/jhrmnn/libmbd/compare/0.7.0...0.8.0
[0.7.0]: https://github.com/jhrmnn/libmbd/compare/0.6.0...0.7.0
[0.6.0]: https://github.com/jhrmnn/libmbd/compare/0.5.0...0.6.0
[0.5.0]: https://github.com/jhrmnn/libmbd/compare/0.4.3...0.5.0
[0.4.3]: https://github.com/jhrmnn/libmbd/compare/0.4.2...0.4.3
[0.4.2]: https://github.com/jhrmnn/libmbd/compare/0.4.1...0.4.2
[0.4.1]: https://github.com/jhrmnn/libmbd/compare/0.4.0...0.4.1
[0.4.0]: https://github.com/jhrmnn/libmbd/releases/tag/0.4.0
