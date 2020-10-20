****************************************************
ScalapackFx - Modern Fortran Interface for ScaLAPACK
****************************************************

The open source library `ScalapackFx <https://github.com/dftbplus/scalapackfx>`_
provides convenient modern Fortran (Fortran 2003) wrappers for the routines of
the ScaLAPACK library. Currently mostly the routines related to diagonalization
are covered.

The documentation is included inside the repository, but is also available at
`dftbplus.github.io <https://dftbplus.github.io/>`_.


Installation
============

Prerequisites
-------------

* CMake (version >= 3.16)

* Fortran 2003 compatible Fortran compiler

* MPI-library and wrappers for your compiler

* ScaLAPACK and LAPACK libraries

* `Fypp preprocessor <https://github.com/aradi/fypp>`_


Building and installing the library
-----------------------------------

The library can be built and installed with the usual CMake-workflow::

  FC=gfortran cmake -B _build -DCMAKE_INSTALL_PREFIX=$HOME/opt/scalapackfx
  cmake --build _build
  cmake --install _build

You can influence the configuration via CMake-variables. A list of those
variables can be found in `config.cmake <config.cmake>`_. You can either modify
the values directly there or pass them as command line option at the
configuration phase, e.g.::

  FC=ifort cmake -B _build \
    -DSCALAPACK_LIBRARY=mkl_scalapack_lp64;mkl_blacs_intelmpi_lp64 \
    -DLAPACK_LIBRARY=mkl_intel_lp64;mkl_sequential;mkl_core

When customizing the external libraries as above, each of them will be searched
for and checked for existence. You can disable these checks by setting the
``SCALAPACK_DETECTION`` and ``LAPACK_DETECTION`` variables to ``False`` (not
recommended).


Testing
-------

A few tests (and usage examples) can be found in the `test/` subdirectory. The
compiled test examples must be invoked from this directory, e.g.::

  cd test
  mpirun -n 2 ../_build/test/test_det


Using the library
=================

CMake build
-----------

* Make sure to add the root folder of the installed library to the
  ``CMAKE_PREFIX_PATH`` environment variable.

* Use ``find_package()`` in `CMakeLists.txt` to locate the library and link
  ``ScalapackFx::ScalapackFx`` to every target which relies directly on the
  library ::

    cmake_minimum_required(VERSION 3.16)

    project(TestScalapackFxBuild LANGUAGES Fortran)

    find_package(ScalapackFx REQUIRED)

    add_executable(test_build test_build.f90)
    target_link_libraries(test_build ScalapackFx::ScalapackFx)

Note: When ScalapackFx is found by CMake, it will try to resolve its own
dependencies (ScaLAPACK and LAPACK) via the ``find_dependency()`` function. You
can use your own finders for those dependencies as long as they export the
targets ``Scalapack::Scalapack`` and ``LAPACK::LAPACK`` with the proper
configuration as ScalapackFx uses those as link targets through the
``target_link_libraries()`` function. If those targets do not exists when
ScalapackFx is being searched for, the custom finders shipped with ScalapackFx
will be invoked to search for the necessary libraries.


Pkg-config build
----------------

* Make sure to add the `lib/pkgconfig` folder of the installed library to the
  ``PKG_CONFIG_PATH`` environment variable.

* Query the include and library options needed for the build with the usual
  ``pkg-config`` commands::

    mpifort $(pkg-config --cflags mpifx) test_mpifxbuild.f90 $(pkg-config --libs mpifx)

  Note, that neither ``-cflags`` or ``--libs`` return any options related to
  your MPI-framework nor is the MPI-framework specified as dependency in the
  pkg-config file. Use the MPI-wrapper of your compiler to compile and link your
  executable or pass the additional include and library options by hand.


License
=======

ScalapackFx is licensed under the `BSD-2-Clause License <LICENSE>`_.
