include(common.m4)

dnl ************************************************************************
dnl *** psyr / pher
dnl ************************************************************************ 

define(`_subroutine_interface_psyr_pher',`
dnl $1: comment
dnl $2: subroutine name
dnl $3: dummy arguments type
dnl $4: dummy arguments kind
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
dnl $4: dummy arguments kind
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

dnl ************************************************************************
dnl *** psymm / phemm
dnl ************************************************************************ 

define(`_subroutine_interface_psymm_phemm',`
dnl $1: comment
dnl $2: subroutine name
dnl $3: dummy arguments type
dnl $4: dummy arguments kind
!> Symmetric/hermitian matrix with general matrix product ($1).
subroutine $2(side, uplo, mm, nn, alpha, aa, ia, ja, desca, &
    & bb, ib, jb, descb, beta, cc, ic, jc, descc)
  import
  character, intent(in) :: side
  character, intent(in) :: uplo
  integer, intent(in) :: mm
  integer, intent(in) :: nn
  $3($4), intent(in) :: alpha
  integer, intent(in) :: desca(DLEN_)
  $3($4), intent(in) :: aa(desca(LLD_), *)
  integer, intent(in) :: ia, ja
  $3($4), intent(in) :: beta
  integer, intent(in) :: descb(DLEN_)
  $3($4), intent(in) :: bb(descb(LLD_), *)
  integer, intent(in) :: ib, jb
  integer, intent(in) :: descc(DLEN_)
  $3($4), intent(in) :: cc(descc(LLD_), *)
  integer, intent(in) :: ic, jc
end subroutine $2
')

dnl ************************************************************************
dnl *** ptrmm
dnl ************************************************************************ 

define(`_subroutine_interface_ptrmm',`
dnl $1: comment
dnl $2: subroutine name
dnl $3: dummy arguments type
!> Symmetric/hermitian matrix vector product ($1).
subroutine $2(side, uplo, transa, diag, mm, nn, alpha, aa, ia, ja, desca, &
    & bb, ib, jb, descb)
  import
  character, intent(in) :: side, uplo, transa, diag
  integer, intent(in) :: mm, nn
  $3($4), intent(in) :: alpha
  integer, intent(in) :: desca(DLEN_)
  $3($4), intent(in) :: aa(desca(LLD_), *)
  integer, intent(in) :: ia, ja
  integer, intent(in) :: descb(DLEN_)
  $3($4), intent(in) :: bb(descb(LLD_), *)
  integer, intent(in) :: ib, jb
end subroutine $2
')

dnl ************************************************************************
dnl *** pgemm
dnl ************************************************************************ 

define(`_subroutine_interface_pgemm',`
dnl $1: comment
dnl $2: subroutine name
dnl $3: dummy arguments type
dnl $4: dummy arguments kind
!> Symmetric/hermitian matrix vector product ($1).
subroutine $2(transa, transb, mm, nn, kk, alpha, aa, ia, ja, desca, &
    & bb, ib, jb, descb, beta, cc, ic, jc, descc)
  import
  character, intent(in) :: transa, transb
  integer, intent(in) :: mm, nn, kk
  $3($4), intent(in) :: alpha
  integer, intent(in) :: desca(DLEN_)
  $3($4), intent(in) :: aa(desca(LLD_), *)
  integer, intent(in) :: ia, ja
  integer, intent(in) :: descb(DLEN_)
  $3($4), intent(in) :: bb(descb(LLD_), *)
  integer, intent(in) :: ib, jb
  $3($4), intent(in) :: beta
  integer, intent(in) :: descc(DLEN_)
  $3($4), intent(inout) :: cc(descb(LLD_), *)
  integer, intent(in) :: ic, jc
end subroutine $2
')


dnl ************************************************************************
dnl *** ptran
dnl ************************************************************************ 

define(`_subroutine_interface_ptran',`
dnl $1: comment
dnl $2: subroutine name
dnl $3: dummy arguments type
dnl $4: dummy arguments kind
!> Matrix transpose ($1).
subroutine $2(mm, nn, alpha, aa, ia, ja, desca, &
    & beta, cc, ic, jc, descc)
  import
  integer, intent(in) :: mm, nn
  $3($4), intent(in) :: alpha
  integer, intent(in) :: desca(DLEN_)
  $3($4), intent(in) :: aa(desca(LLD_), *)
  integer, intent(in) :: ia, ja
  $3($4), intent(in) :: beta
  integer, intent(in) :: descc(DLEN_)
  $3($4), intent(inout) :: cc(descc(LLD_), *)
  integer, intent(in) :: ic, jc
end subroutine $2
')


dnl ************************************************************************
dnl *** ptranc
dnl ************************************************************************ 

define(`_subroutine_interface_ptranc',`
dnl $1: comment
dnl $2: subroutine name
dnl $3: dummy arguments type
dnl $4: dummy arguments kind
!> Conjugated matrix transpose ($1).
subroutine $2(mm, nn, alpha, aa, ia, ja, desca, &
    & beta, cc, ic, jc, descc)
  import
  integer, intent(in) :: mm, nn
  $3($4), intent(in) :: alpha
  integer, intent(in) :: desca(DLEN_)
  $3($4), intent(in) :: aa(desca(LLD_), *)
  integer, intent(in) :: ia, ja
  $3($4), intent(in) :: beta
  integer, intent(in) :: descc(DLEN_)
  $3($4), intent(inout) :: cc(descc(LLD_), *)
  integer, intent(in) :: ic, jc
end subroutine $2
')
