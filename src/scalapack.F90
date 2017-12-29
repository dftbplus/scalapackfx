include(scalapack.m4)

!> Wrapper functions for scalapack.
module scalapack_module
  use scalapackfx_common_module
  implicit none
  private

  public :: psygst, phegst, psyev, pheev, psyevd, pheevd, psyevr, pheevr
  public :: ptrsm, ppotrf, ppotri, ptrtri, pgesvd
  public :: sl_init, numroc, infog2l, indxl2g, descinit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCALAPACK CORE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Cholesky factorization of symmetric/Hermitian positive definite matrix.
  interface ppotrf
    _subroutine_interface_ppotrf(real, s, real(sp))
    _subroutine_interface_ppotrf(dreal, d, real(dp))
    _subroutine_interface_ppotrf(complex, c, complex(sp))
    _subroutine_interface_ppotrf(dcomplex, z, complex(dp))
  end interface ppotrf

  !> Inversion of a Cholesky decomposed symmetric/Hermitian matrix.
  interface ppotri
    _subroutine_interface_ppotri(real, s, real(sp))
    _subroutine_interface_ppotri(dreal, d, real(dp))
    _subroutine_interface_ppotri(complex, c, complex(sp))
    _subroutine_interface_ppotri(dcomplex, z, complex(dp))
  end interface ppotri

  !> Inversion of a triangular matrix.
  interface ptrtri
    _subroutine_interface_ptrtri(real, s, real(sp))
    _subroutine_interface_ptrtri(dreal, d, real(dp))
    _subroutine_interface_ptrtri(complex, c, complex(sp))
    _subroutine_interface_ptrtri(dcomplex, z, complex(dp))
  end interface ptrtri

  !> Reduces generalized symmetric eigenvalue problem to standard form.
  interface psygst
    _subroutine_interface_psygst(real, s, sp)
    _subroutine_interface_psygst(dreal, d, dp)
  end interface psygst

  !> Reduces generalized Hermitian eigenvalue problem to standard form.
  interface phegst
    _subroutine_interface_phegst(complex, c, sp)
    _subroutine_interface_phegst(dcomplex, z, dp)
  end interface phegst

  !> Solves the symmetric eigenvalue problem.
  interface psyev
    _subroutine_interface_psyev(real, s, sp)
    _subroutine_interface_psyev(dreal, d, dp)
  end interface psyev

  !> Solves the Hermitian eigenvalue problem.
  interface pheev
    _subroutine_interface_pheev(complex, c, sp)
    _subroutine_interface_pheev(complex, z, dp)
  end interface pheev

  !> Solves the symmetric eigenvalue problem by divide and conquer algorithm.
  interface psyevd
    _subroutine_interface_psyevd(real, s, sp)
    _subroutine_interface_psyevd(dreal, d, dp)
  end interface psyevd

  !> Solves the Hermitian eigenvalue problem by divide and conquer algorithm.
  interface pheevd
    _subroutine_interface_pheevd(complex, c, sp)
    _subroutine_interface_pheevd(dcomplex, z, dp)
  end interface pheevd
  
  !> Solves the symmetric eigenvalue problem by the MRRR algorithm.
  interface psyevr
    _subroutine_interface_psyevr(real, s, sp)
    _subroutine_interface_psyevr(dreal, d, dp)
  end interface psyevr
  
  !> Solves the Hermitian eigenvalue problem by the MRRR algorithm.
  interface pheevr
    _subroutine_interface_pheevr(complex, c, sp)
    _subroutine_interface_pheevr(dcomplex, z, dp)
  end interface pheevr
  
  !> Singular value decomposition of a matrix
  interface pgesvd
    _subroutine_interface_prgesvd(real, s, sp)
    _subroutine_interface_prgesvd(dreal, d, dp)
    _subroutine_interface_pcgesvd(complex, c, sp)
    _subroutine_interface_pcgesvd(dcomplex, z, dp)
  end interface pgesvd
  
  !> Linear system of equation for triangular matrix.
  interface ptrsm
    _subroutine_interface_ptrsm(real, s, real(sp))
    _subroutine_interface_ptrsm(dreal, d, real(dp))
    _subroutine_interface_ptrsm(complex, c, complex(sp))
    _subroutine_interface_ptrsm(dcomplex, z, complex(dp))
  end interface ptrsm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCALAPACK TOOLS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface

    !> Scalapack initialization routine.
    subroutine sl_init(ictxt, nprow, npcol)
      integer, intent(out) :: ictxt
      integer, intent(in) :: nprow, npcol
    end subroutine sl_init

    !> Number of rows or columns of distributed matrix owned by local process.
    function numroc(nn, nb, iproc, isrcproc, nproc)
      integer, intent(in) :: nn, nb, iproc, isrcproc, nproc
      integer :: numroc
    end function numroc
    
    !> Converts global matrix index into local.
    subroutine infog2l(grindx, gcindx, desc, nprow, npcol, myrow, mycol,&
        & lrindx, lcindx, rsrc, csrc)
      import DLEN_
      integer, intent(in) :: grindx, gcindx, desc(DLEN_)
      integer, intent(in) :: nprow, npcol, myrow, mycol
      integer, intent(out) :: lrindx, lcindx, rsrc, csrc
    end subroutine infog2l
    
    !> Converts local matrix index into global.
    function indxl2g(indxglob, nb, iproc, isrcproc, nprocs)
      integer :: indxl2g
      integer, intent(in) :: indxglob, nb, iproc, isrcproc, nprocs
    end function indxl2g
    
    !> Initializes a descriptor for a distributed array.
    subroutine descinit(desc, mm, nn, mb, nb, irsrc, icsrc, ictxt, lld, info)
      import DLEN_
      integer, intent(out) :: desc(DLEN_)
      integer, intent(in) :: mm, nn, mb, nb, irsrc, icsrc, ictxt, lld
      integer, intent(out) :: info
    end subroutine descinit

  end interface

  
end module scalapack_module
