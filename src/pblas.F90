include(pblas.m4)

!> Interface to PBLAS routines
module pblas_module
  use scalapackfx_common_module
  implicit none
  private
  
  public :: psyr, pher
  public :: psyrk, pherk
  public :: psymv, phemv
  public :: psymm, phemm
  public :: pgemm
  public :: ptrmm
  public :: ptran, ptranu
  public :: ptranc

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

  !> Hermitian matrix vector product
  interface phemv
    _subroutine_interface_psymv_phemv(complex, pchemv, complex, sp)
    _subroutine_interface_psymv_phemv(dcomplex, pzhemv, complex, dp)
  end interface phemv
  
  !> Symmetric matrix general matrix product
  interface psymm
    _subroutine_interface_psymm_phemm(real, pssymm, real, sp)
    _subroutine_interface_psymm_phemm(dreal, pdsymm, real, dp)
  end interface psymm
  
  !> Hermitian matrix general matrix product
  interface phemm
    _subroutine_interface_psymm_phemm(complex, pchemm, complex, sp)
    _subroutine_interface_psymm_phemm(dcomplex, pzhemm, complex, dp)
  end interface phemm
  
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

  !> Real matrix transpose.
  interface ptran
    _subroutine_interface_ptranx(real, pstran, real, sp)
    _subroutine_interface_ptranx(dreal, pdtran, real, dp)
  end interface ptran

  !> Complex matrix transpose.
  interface ptranu
    _subroutine_interface_ptranx(complex, pctranu, complex, sp)
    _subroutine_interface_ptranx(dcomplex, pztranu, complex, dp)
  end interface ptranu

  !> Complex hermitian matrix transpose.
  interface ptranc
    _subroutine_interface_ptranx(complex, pctranc, complex, sp)
    _subroutine_interface_ptranx(dcomplex, pztranc, complex, dp)

  end interface ptranc

end module pblas_module
