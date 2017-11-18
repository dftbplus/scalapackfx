About SCALAPACKFX
=================

`SCALAPACKFX <https://github.com/dftbplus/scalapackfx>`_ is a library containing
modern Fortran (Fortran 2003) wrappers around SCALAPACK, PBLAS and BLACS
routines. The goal is to make the use of those libraries as simple as possible
in Fortran.

Consider for example a simple broadcast in BLACS. In order to broadcast an
integer array (aa) with 5x5 elements using the appropriate BLACS routine, you
have to issue::

    call igebs2d(ictxt, "All", " ", 5, 5, aa, 5)

Additional to the object to be broadcasted and the communicator, you also
*must* specify following arguments:

- type of the array as first letter of the subroutine name (although it is
  *known* at compile-time)

- number of row, number of columns and leading dimension of the array (although
  those are *known* at run-time)

- scope of the broadcast (could be optional as very often it is "All" any way)

- communcation pattern for the broadcast (could be optional as very often
  " " is used)

Using SCALAPACKFX the call above is as simple as::

    call blacsfx_gebs(mygrid, aa)

No redundant arguments, sensible defaults. Nevertheless the full functionality
still available via optional parameters if needed. E.g. if you wanted to the
scope, you could write::

    call blacsfx_gebs(mygrid, aa, scope="Row")

Also, the array aa does not have to be a rank two array, it can be also rank one
(vector) or rank zero (scalar).

A few routines are already covered (see :ref:`sec_routines`). If your desired
routine is not among them yet, you are cordially invited to extend SCALAPACKFX
and to share it in order to let others profit from your work (SCALAPACKFX is
licensed under the simplified BSD license). For more details see the `project
home page <https://github.com/dftbplus/scalapackfx>`_.
