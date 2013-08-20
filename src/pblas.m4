include(common.m4)

dnl ************************************************************************
dnl *** psyr / pher
dnl ************************************************************************ 

define(`_subroutine_interface_psyr_pher',`
dnl $1: comment
dnl $2: subroutine name
dnl $3: dummy arguments type
!> Symmetric/hermitian rank one update ($1).
subroutine $2(uplo, nn, alpha, xx, ix, jx, descx, incx, aa, ia, ja, desca)
  import
  character, intent(in) :: uplo
  integer, intent(in) :: nn
  real($4), intent(in) :: alpha
  integer, intent(in) :: descx(DLEN_)
  $3($4), intent(in) :: xx(descx(LLD_), *)
  integer, intent(in) :: ix, jx
  integer, intent(in) :: incx
  integer, intent(in) :: desca(DLEN_)
  $3($4), intent(inout) :: aa(desca(LLD_), *)
  integer, intent(in) :: ia, ja
end subroutine $2
')

dnl ************************************************************************
dnl *** psyrk / pherk
dnl ************************************************************************ 

define(`_subroutine_interface_psyrk_pherk',`
dnl $1: comment
dnl $2: subroutine name
dnl $3: dummy arguments type
dnl $4: dummy arguments kind
!> Symmetric/hermitian rank-k update ($1).
subroutine $2(uplo, trans, nn, kk, alpha, aa, ia, ja, desca, beta, cc,&
    & ic, jc, descc)
  import
  character, intent(in) :: uplo, trans
  integer, intent(in) :: nn, kk
  real($4), intent(in) :: alpha
  integer, intent(in) :: desca(DLEN_)
  $3($4), intent(in) :: aa(desca(LLD_), *)
  integer, intent(in) :: ia, ja
  real($4), intent(in) :: beta
  integer, intent(in) :: descc(DLEN_)
  $3($4), intent(inout) :: cc(descc(LLD_), *)
  integer, intent(in) :: ic, jc
end subroutine $2
')

dnl ************************************************************************
dnl *** psymv / phemv
dnl ************************************************************************ 

define(`_subroutine_interface_psymv_phemv',`
dnl $1: comment
dnl $2: subroutine name
dnl $3: dummy arguments type
!> Symmetric/hermitian matrix vector product ($1).
subroutine $2(uplo, nn, alpha, aa, ia, ja, desca, xx, ix, jx, descx, incx, &
    & beta, yy, iy, jy, descy, incy)
  import
  character, intent(in) :: uplo
  integer, intent(in) :: nn
  $3($4), intent(in) :: alpha
  integer, intent(in) :: desca(DLEN_)
  $3($4), intent(in) :: aa(desca(LLD_), *)
  integer, intent(in) :: ia, ja
  integer, intent(in) :: descx(DLEN_)
  $3($4), intent(in) :: xx(descx(LLD_), *)
  integer, intent(in) :: ix, jx, incx
  $3($4), intent(in) :: beta
  integer, intent(in) :: descy(DLEN_)
  $3($4), intent(inout) :: yy(descy(LLD_), *)
  integer, intent(in) :: iy, jy, incy
end subroutine $2
')
