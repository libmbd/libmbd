Libmbd — many-body dispersion library
=========================================

This project contains implementations of the `many-body dispersion <http://dx.doi.org/10.1063/1.4865104>`_ (MBD) method in several programming languages and frameworks:

- The Fortran implementation is the reference, most advanced implementation, with support for analytical gradients and distributed parallelism, and additional functionality beyond the MBD method itself. It provides a low-level and a high-level Fortran API, and a C API. Furthermore, Python bindings to the C API are provided.
- The Python/Numpy implementation is intended for prototyping, and as a high-level language reference.
- The Python/Tensorflow implemntation is an experiment that should enable rapid prototyping of machine learning applications with MBD.

Since all implementations provide Python bindings, the project is structured as a Python package called Pymbd, however, the implementations in languages other than Python can be used as standalone libraries in their respective languages.

.. toctree::

   fortran-c-api
