include(pblasfx.m4)

!> High level Fortran wrappers for the PBLAS library.
module pblasfx_module
  use scalapackfx_common_module
  use pblas_module
  implicit none
  private
  
  public :: pblasfx_psyr, pblasfx_pher
  public :: pblasfx_psyrk, pblasfx_pherk
  public :: pblasfx_psymv, pblasfx_phemv
  public :: pblasfx_psymm, pblasfx_phemm
  public :: pblasfx_pgemm
  public :: pblasfx_ptrmm
  public :: pblasfx_ptran, pblasfx_ptranu
  public :: pblasfx_ptranc
  
  interface pblasfx_psyr
    module procedure pblasfx_psyr_real, pblasfx_psyr_dreal
  end interface pblasfx_psyr

  interface pblasfx_pher
    module procedure pblasfx_pher_complex, pblasfx_pher_dcomplex
  end interface pblasfx_pher
  
  interface pblasfx_psyrk
    module procedure pblasfx_psyrk_real, pblasfx_psyrk_dreal
  end interface pblasfx_psyrk
  
  interface pblasfx_pherk
    module procedure pblasfx_pherk_complex, pblasfx_pherk_dcomplex
  end interface pblasfx_pherk
  
  interface pblasfx_psymv
    module procedure pblasfx_psymv_real, pblasfx_psymv_dreal
  end interface pblasfx_psymv
  
  interface pblasfx_phemv
    module procedure pblasfx_phemv_complex, pblasfx_phemv_dcomplex
  end interface pblasfx_phemv
  
  interface pblasfx_psymm
    module procedure pblasfx_psymm_real, pblasfx_psymm_dreal
  end interface pblasfx_psymm
  
  interface pblasfx_phemm
    module procedure pblasfx_phemm_complex, pblasfx_phemm_dcomplex
  end interface pblasfx_phemm
  
  interface pblasfx_ptrmm
    module procedure pblasfx_ptrmm_real, pblasfx_ptrmm_dreal, &
        & pblasfx_ptrmm_complex, pblasfx_ptrmm_dcomplex
  end interface pblasfx_ptrmm

  interface pblasfx_pgemm
    module procedure pblasfx_pgemm_real, pblasfx_pgemm_dreal, &
        & pblasfx_pgemm_complex, pblasfx_pgemm_dcomplex
  end interface pblasfx_pgemm


  interface pblasfx_ptran
    module procedure pblasfx_ptran_real, pblasfx_ptran_dreal
  end interface pblasfx_ptran
  
  interface pblasfx_ptranu
    module procedure pblasfx_ptranu_complex, pblasfx_ptranu_dcomplex
  end interface pblasfx_ptranu
  
  interface pblasfx_ptranc
    module procedure pblasfx_ptranc_complex, pblasfx_ptranc_dcomplex
  end interface pblasfx_ptranc

contains

  _subroutine_pblasfx_psyr_pher(psyr_real, real, sp, real, psyr)
  _subroutine_pblasfx_psyr_pher(psyr_dreal, real, dp, real, psyr)
  _subroutine_pblasfx_psyr_pher(pher_complex, complex, sp, cmplx, pher)
  _subroutine_pblasfx_psyr_pher(pher_dcomplex, complex, dp, cmplx, pher)

  _subroutine_pblasfx_psyrk_pherk(psyrk_real, real, sp, real, psyrk)
  _subroutine_pblasfx_psyrk_pherk(psyrk_dreal, real, dp, real, psyrk)
  _subroutine_pblasfx_psyrk_pherk(pherk_complex, complex, sp, cmplx, pherk)
  _subroutine_pblasfx_psyrk_pherk(pherk_dcomplex, complex, dp, cmplx, pherk)

  _subroutine_pblasfx_psymv_phemv(psymv_real, real, sp, real, psymv)
  _subroutine_pblasfx_psymv_phemv(psymv_dreal, real, dp, real, psymv)
  _subroutine_pblasfx_psymv_phemv(phemv_complex, complex, sp, cmplx, phemv)
  _subroutine_pblasfx_psymv_phemv(phemv_dcomplex, complex, dp, cmplx, phemv)

  _subroutine_pblasfx_psymm_phemm(psymm_real, real, sp, real, psymm)
  _subroutine_pblasfx_psymm_phemm(psymm_dreal, real, dp, real, psymm)
  _subroutine_pblasfx_psymm_phemm(phemm_complex, complex, sp, cmplx, phemm)
  _subroutine_pblasfx_psymm_phemm(phemm_dcomplex, complex, dp, cmplx, phemm)

  _subroutine_pblasfx_pgemm(real, real, sp, real)
  _subroutine_pblasfx_pgemm(dreal, real, dp, real)
  _subroutine_pblasfx_pgemm(complex, complex, sp, cmplx)
  _subroutine_pblasfx_pgemm(dcomplex, complex, dp, cmplx)

  _subroutine_pblasfx_ptrmm(real, real, sp, real)
  _subroutine_pblasfx_ptrmm(dreal, real, dp, real)
  _subroutine_pblasfx_ptrmm(complex, complex, sp, cmplx)
  _subroutine_pblasfx_ptrmm(dcomplex, complex, dp, cmplx)

  _subroutine_pblasfx_ptran(ptran_real, real, sp, real, ptran)
  _subroutine_pblasfx_ptran(ptran_dreal, real, dp, real, ptran)
  _subroutine_pblasfx_ptranu(ptranu_complex, complex, sp, complex, ptranu)
  _subroutine_pblasfx_ptranu(ptranu_dcomplex, complex, dp, complex, ptranu)
  _subroutine_pblasfx_ptranc(ptranc_complex, complex, sp, complex, ptranc)
  _subroutine_pblasfx_ptranc(ptranc_dcomplex, complex, dp, complex, ptranc)

end module pblasfx_module
