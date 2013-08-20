!> \mainpage Modern Fortran wrappers for SCALAPACK, PBLAS, BLACS
!!
!! \tableofcontents
!!
!! \section About About
!!
!! The ScaLAPACKFX library provides high level Fortran (Fortran 2003)
!! wrappers around the routines in the SCALAPACK, PBLAS and BLACS libraries and
!! some (hopefully) useful extensions to the libraries.
!! The subroutines can be found in the modules
!!
!!   * \ref blacsfx_module "blacsfx_module"
!!   * \ref pblasfx_module "pblasfx_module"
!!   * \ref scalapackfx_module "scalapackfx_module"
!!   * \ref scalapackfx_tools_module "scalapackfx_tools_module"
!!
!! You can also directly import the full functionality as one module
!! by importing the module
!!
!!   * \ref libscalapackfx_module "libscalapackfx_module"
!!
!! For the typical application you should use the modules above. If for some
!! reasons you need the original (low level) SCALAPACK, PBLAS or BLACS routines,
!! you find their Fortran 2003 interfaces in the modules
!!
!!   * \ref blacs_module "blacs_module"
!!   * \ref pblas_module "pblas_module"
!!   * \ref scalapack_module "scalapack_module"
!!
!! However, their use should not be really necesary as the high level
!! wrappers expose the full functionality of the original library routines.
!!
!!
!! \section Using Using the library
!!
!! You can access the wraped functionality by importing the appropriate
!! module(s) in your program. The following program snippet demonstrates this by
!! solving the generalized eigenvalue problem using the library. Further
!! examples you find in the \c test/ folder.
!!
!!     program test_diag
!!       use libscalapackfx_module
!!
!!       integer, parameter :: nprow = 4, npcol = 4   ! process rows/cols
!!       integer, parameter :: bsize = 64             ! block size
!!       integer, parameter :: nn = 1000              ! matrix size
!!
!!       type(blacsfx) :: mygrid  ! BLACS descriptor
!!       real(dp), allocatable :: aa(:,:), bb(:,:), eigvecs(:,:), eigvals(:)
!!       integer :: desc(DLEN_)    ! matrix descriptor
!!       integer :: mloc, nloc     ! nr. of local rows/columns of the matrices
!!
!!
!!       ! Initialize your BLACS grid
!!       call mygrid%init()
!!       call mygrid%gridinit(nprow, npcol)
!!
!!       ! Allocate the local part of the distributed matrices
!!       call scalafx_getdescriptor(mygrid, nn, nn, bsize, bsize, desc)
!!       call scalafx_getlocalshape(mygrid, desc, mloc, nloc)
!!       allocate(aa(mloc, nloc))
!!       allocate(bb(mloc, nloc))
!!       allocate(eigvecs(mloc, nloc))
!!       allocate(eigvals(nn))
!!       ...
!!       ! Here comes the code which distributes your matrix
!!       ...
!!       ! Get eigenvalues (on all nodes) and eigenvectors (distributed)
!!       call psygvd(aa, desc, bb, desc, eigvals, eigvecs, desc, jobz="V", uplo="L")
!!       ...
!!     end program test_diag
!!


!> Exports the functionality of the ScaLAPACKFX library as one module.
module libscalapackfx_module
  use blacsfx_module
  use pblasfx_module
  use scalapackfx_module
  use scalapackfx_tools_module
  implicit none

  public

end module libscalapackfx_module
