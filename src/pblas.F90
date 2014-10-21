include(pblas.m4)

!> Interface to PBLAS routines
module pblas_module
  use scalapackfx_common_module
  implicit none
  private

  public :: psyr, pher
  public :: psyrk, pherk
  public :: psymv, phemv
  public :: pgemm
  public :: ptrmm

  !> Symmetric rank one update.
  interface psyr
    _subroutine_interface_psyr_pher(real, pssyr, real, sp)
    _subroutine_interface_psyr_pher(dreal, pdsyr, real, dp)
  end interface psyr

  !> Hermitian rank one update.
  interface pher
    _subroutine_interface_psyr_pher(complex, pcher, complex, sp)
    _subroutine_interface_psyr_pher(dcomplex, pzher, complex, dp)
  end interface pher

  !> Symmetric rank-k update.
  interface psyrk
    _subroutine_interface_psyrk_pherk(real, pssyrk, real, sp)
    _subroutine_interface_psyrk_pherk(dreal, pdsyrk, real, dp)
  end interface psyrk

  !> Hermitian rank-k update.
  interface pherk
    _subroutine_interface_psyrk_pherk(complex, pcherk, complex, sp)
    _subroutine_interface_psyrk_pherk(dcomplex, pzherk, complex, dp)
  end interface pherk

  !> Symmetric matrix vector product
  interface psymv
    _subroutine_interface_psymv_phemv(real, pssymv, real, sp)
    _subroutine_interface_psymv_phemv(dreal, pdsymv, real, dp)
  end interface psymv

  !> Symmetric matrix vector product
  interface phemv
    _subroutine_interface_psymv_phemv(complex, pchemv, complex, sp)
    _subroutine_interface_psymv_phemv(dcomplex, pzhemv, complex, dp)
  end interface phemv

  !> Triangular matrix matrix product
  interface ptrmm
    _subroutine_interface_ptrmm(real, pstrmm, real, sp)
    _subroutine_interface_ptrmm(dreal, pdtrmm, real, dp)
    _subroutine_interface_ptrmm(complex, pctrmm, complex, sp)
    _subroutine_interface_ptrmm(dcomplex, pztrmm, complex, dp)
  end interface ptrmm

  !> Genereal matrix matrix product
  interface pgemm
    _subroutine_interface_pgemm(real, psgemm, real, sp)
    _subroutine_interface_pgemm(dreal, pdgemm, real, dp)
    _subroutine_interface_pgemm(complex, pcgemm, complex, sp)
    _subroutine_interface_pgemm(dcomplex, pzgemm, complex, dp)
  end interface pgemm

end module pblas_module
