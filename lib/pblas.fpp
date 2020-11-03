
! ************************************************************************
! *** psyr / pher
! ************************************************************************


#:def interface_psyr_pher_template(COMMENT, NAME, TYPE, KIND)

  !> Symmetric/hermitian rank one update (${COMMENT}$).
  subroutine ${NAME}$(uplo, nn, alpha, xx, ix, jx, descx, incx, aa, ia, ja, desca)
    import
    character, intent(in) :: uplo
    integer, intent(in) :: nn
    real(${KIND}$), intent(in) :: alpha
    integer, intent(in) :: descx(*)
    ${TYPE}$(${KIND}$), intent(in) :: xx(descx(LLD_), *)
    integer, intent(in) :: ix, jx
    integer, intent(in) :: incx
    integer, intent(in) :: desca(*)
    ${TYPE}$(${KIND}$), intent(inout) :: aa(desca(LLD_), *)
    integer, intent(in) :: ia, ja
  end subroutine ${NAME}$

#:enddef interface_psyr_pher_template



! ************************************************************************
! *** psyrk / pherk
! ************************************************************************

#:def interface_psyrk_pherk_template(COMMENT, NAME, TYPE, KIND)

  !> Symmetric/hermitian rank-k update (${COMMENT}$).
  subroutine ${NAME}$(uplo, trans, nn, kk, alpha, aa, ia, ja, desca, beta, cc,&
      & ic, jc, descc)
    import
    character, intent(in) :: uplo, trans
    integer, intent(in) :: nn, kk
    real(${KIND}$), intent(in) :: alpha
    integer, intent(in) :: desca(*)
    ${TYPE}$(${KIND}$), intent(in) :: aa(desca(LLD_), *)
    integer, intent(in) :: ia, ja
    real(${KIND}$), intent(in) :: beta
    integer, intent(in) :: descc(*)
    ${TYPE}$(${KIND}$), intent(inout) :: cc(descc(LLD_), *)
    integer, intent(in) :: ic, jc
  end subroutine ${NAME}$

#:enddef interface_psyrk_pherk_template



! ************************************************************************
! *** psymv / phemv
! ************************************************************************

#:def interface_psymv_phemv_template(COMMENT, NAME, TYPE, KIND)

  !> Symmetric/hermitian matrix vector product ($1).
  subroutine ${NAME}$(uplo, nn, alpha, aa, ia, ja, desca, xx, ix, jx, descx, incx, &
      & beta, yy, iy, jy, descy, incy)
    import
    character, intent(in) :: uplo
    integer, intent(in) :: nn
    ${TYPE}$(${KIND}$), intent(in) :: alpha
    integer, intent(in) :: desca(*)
    ${TYPE}$(${KIND}$), intent(in) :: aa(desca(LLD_), *)
    integer, intent(in) :: ia, ja
    integer, intent(in) :: descx(*)
    ${TYPE}$(${KIND}$), intent(in) :: xx(descx(LLD_), *)
    integer, intent(in) :: ix, jx, incx
    ${TYPE}$(${KIND}$), intent(in) :: beta
    integer, intent(in) :: descy(*)
    ${TYPE}$(${KIND}$), intent(inout) :: yy(descy(LLD_), *)
    integer, intent(in) :: iy, jy, incy
  end subroutine ${NAME}$

#:enddef interface_psymv_phemv_template



! ************************************************************************
! *** psymm / phemm
! ************************************************************************

#:def interface_psymm_phemm_template(COMMENT, NAME, TYPE, KIND)

  !> Symmetric/hermitian matrix with general matrix product ($1).
  subroutine ${NAME}$(side, uplo, mm, nn, alpha, aa, ia, ja, desca, &
      & bb, ib, jb, descb, beta, cc, ic, jc, descc)
    import
    character, intent(in) :: side
    character, intent(in) :: uplo
    integer, intent(in) :: mm
    integer, intent(in) :: nn
    ${TYPE}$(${KIND}$), intent(in) :: alpha
    integer, intent(in) :: desca(*)
    ${TYPE}$(${KIND}$), intent(in) :: aa(desca(LLD_), *)
    integer, intent(in) :: ia, ja
    ${TYPE}$(${KIND}$), intent(in) :: beta
    integer, intent(in) :: descb(*)
    ${TYPE}$(${KIND}$), intent(in) :: bb(descb(LLD_), *)
    integer, intent(in) :: ib, jb
    integer, intent(in) :: descc(*)
    ${TYPE}$(${KIND}$), intent(in) :: cc(descc(LLD_), *)
    integer, intent(in) :: ic, jc
  end subroutine ${NAME}$

#:enddef interface_psymm_phemm_template


! ************************************************************************
! *** ptrmm
! ************************************************************************

#:def interface_ptrmm_template(COMMENT, NAME, TYPE, KIND)

  !> Symmetric/hermitian matrix vector product (${COMMENT}$).
  subroutine ${NAME}$(side, uplo, transa, diag, mm, nn, alpha, aa, ia, ja, desca, &
      & bb, ib, jb, descb)
    import
    character, intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: mm, nn
    ${TYPE}$(${KIND}$), intent(in) :: alpha
    integer, intent(in) :: desca(*)
    ${TYPE}$(${KIND}$), intent(in) :: aa(desca(LLD_), *)
    integer, intent(in) :: ia, ja
    integer, intent(in) :: descb(*)
    ${TYPE}$(${KIND}$), intent(in) :: bb(descb(LLD_), *)
    integer, intent(in) :: ib, jb
  end subroutine ${NAME}$

#:enddef interface_ptrmm_template


! ************************************************************************
! *** pgemm
! ************************************************************************

#:def interface_pgemm_template(COMMENT, NAME, TYPE, KIND)

  !> General matrix matrix product (${COMMENT}$).
  subroutine ${NAME}$(transa, transb, mm, nn, kk, alpha, aa, ia, ja, desca, &
      & bb, ib, jb, descb, beta, cc, ic, jc, descc)
    import
    character, intent(in) :: transa, transb
    integer, intent(in) :: mm, nn, kk
    ${TYPE}$(${KIND}$), intent(in) :: alpha
    integer, intent(in) :: desca(*)
    ${TYPE}$(${KIND}$), intent(in) :: aa(desca(LLD_), *)
    integer, intent(in) :: ia, ja
    integer, intent(in) :: descb(*)
    ${TYPE}$(${KIND}$), intent(in) :: bb(descb(LLD_), *)
    integer, intent(in) :: ib, jb
    ${TYPE}$(${KIND}$), intent(in) :: beta
    integer, intent(in) :: descc(*)
    ${TYPE}$(${KIND}$), intent(inout) :: cc(descb(LLD_), *)
    integer, intent(in) :: ic, jc
  end subroutine ${NAME}$

#:enddef interface_pgemm_template


! ************************************************************************
! *** ptran / ptranu / ptranc
! ************************************************************************

#:def interface_ptranx_template(COMMENT, NAME, TYPE, KIND)

  !> Transpose of a distributed matrix ($COMMENT$).
  subroutine ${NAME}$(mm, nn, alpha, aa, ia, ja, desca, beta, cc, ic, jc, descc)
    import
    integer, intent(in)   :: mm, nn
    ${TYPE}$(${KIND}$), intent(in)    :: alpha
    integer, intent(in)   :: ia, ja
    integer, intent(in)   :: desca(*)
    ${TYPE}$(${KIND}$), intent(in)    :: aa(desca(LLD_), *)
    ${TYPE}$(${KIND}$), intent(in)    :: beta
    integer, intent(in)   :: ic, jc
    integer, intent(in)   :: descc(*)
    ${TYPE}$(${KIND}$), intent(inout) :: cc(descc(LLD_), *)
  end subroutine ${NAME}$

#:enddef interface_ptranx_template



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
    @:interface_psyr_pher_template(real, pssyr, real, sp)
    @:interface_psyr_pher_template(dreal, pdsyr, real, dp)
  end interface psyr

  !> Hermitian rank one update.
  interface pher
    @:interface_psyr_pher_template(complex, pcher, complex, sp)
    @:interface_psyr_pher_template(dcomplex, pzher, complex, dp)
  end interface pher

  !> Symmetric rank-k update.
  interface psyrk
    @:interface_psyrk_pherk_template(real, pssyrk, real, sp)
    @:interface_psyrk_pherk_template(dreal, pdsyrk, real, dp)
  end interface psyrk

  !> Hermitian rank-k update.
  interface pherk
    @:interface_psyrk_pherk_template(complex, pcherk, complex, sp)
    @:interface_psyrk_pherk_template(dcomplex, pzherk, complex, dp)
  end interface pherk

  !> Symmetric matrix vector product
  interface psymv
    @:interface_psymv_phemv_template(real, pssymv, real, sp)
    @:interface_psymv_phemv_template(dreal, pdsymv, real, dp)
  end interface psymv

  !> Hermitian matrix vector product
  interface phemv
    @:interface_psymv_phemv_template(complex, pchemv, complex, sp)
    @:interface_psymv_phemv_template(dcomplex, pzhemv, complex, dp)
  end interface phemv

  !> Symmetric matrix general matrix product
  interface psymm
    @:interface_psymm_phemm_template(real, pssymm, real, sp)
    @:interface_psymm_phemm_template(dreal, pdsymm, real, dp)
  end interface psymm

  !> Hermitian matrix general matrix product
  interface phemm
    @:interface_psymm_phemm_template(complex, pchemm, complex, sp)
    @:interface_psymm_phemm_template(dcomplex, pzhemm, complex, dp)
  end interface phemm

  !> Triangular matrix matrix product
  interface ptrmm
    @:interface_ptrmm_template(real, pstrmm, real, sp)
    @:interface_ptrmm_template(dreal, pdtrmm, real, dp)
    @:interface_ptrmm_template(complex, pctrmm, complex, sp)
    @:interface_ptrmm_template(dcomplex, pztrmm, complex, dp)
  end interface ptrmm

  !> General matrix matrix product
  interface pgemm
    @:interface_pgemm_template(real, psgemm, real, sp)
    @:interface_pgemm_template(dreal, pdgemm, real, dp)
    @:interface_pgemm_template(complex, pcgemm, complex, sp)
    @:interface_pgemm_template(dcomplex, pzgemm, complex, dp)
  end interface pgemm

  !> Real matrix transpose.
  interface ptran
    @:interface_ptranx_template(real, pstran, real, sp)
    @:interface_ptranx_template(dreal, pdtran, real, dp)
  end interface ptran

  !> Complex matrix transpose.
  interface ptranu
    @:interface_ptranx_template(complex, pctranu, complex, sp)
    @:interface_ptranx_template(dcomplex, pztranu, complex, dp)
  end interface ptranu

  !> Complex hermitian matrix transpose.
  interface ptranc
    @:interface_ptranx_template(complex, pctranc, complex, sp)
    @:interface_ptranx_template(dcomplex, pztranc, complex, dp)
  end interface ptranc

end module pblas_module
