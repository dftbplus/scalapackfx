!> \mainpage Modern Fortran wrappers for SCALAPACK, PBLAS, BLACS
!!
!! The open source library
!! [SCALAPACKFX](<https://github.com/dftbplus/scalapackfx) is an effort to
!! provide modern Fortran (Fortran 2003) wrappers around routines of the
!! SCALAPACK, PBLAS and BLACS libraries to make their use as simple as possible.
!!
!! For more information see the following sources:
!! * [Online documentation](https://github.com/dftbplus/scalapackfx)
!!   for installation and usage of the library
!! * [API documentation](annotated.html) for the reference manual.
!! * [Project home page](https://github.com/dftbplus/scalapackfx)
!!   for the source code, bug tracker and further information on the project.


!> Exports the functionality of the ScaLAPACKFX library as one module.
module libscalapackfx_module
  use blacsfx_module
  use pblasfx_module
  use scalapackfx_module
  use scalapackfx_tools_module
  implicit none

  public

end module libscalapackfx_module
