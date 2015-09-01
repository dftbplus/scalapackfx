include(common.m4)

dnl ************************************************************************
dnl *** psymv/phemv
dnl ************************************************************************

define(`_subroutine_pblasfx_psymv_phemv',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl $3 dummy arguments kind
dnl $4 conversion function
dnl $5 pblas subroutine name
subroutine pblasfx_$1(aa, desca, xx, descx, yy, descy, uplo, alpha, beta, &
    & nn, ia, ja, ix, jx, incx, iy, jy, incy)
  $2($3), intent(in) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  $2($3), intent(in) :: xx(:,:)
  integer, intent(in) :: descx(DLEN_)
  $2($3), intent(inout) :: yy(:,:)
  integer, intent(in) :: descy(DLEN_)
  character, intent(in), optional :: uplo
  $2($3), intent(in), optional :: alpha, beta
  integer, intent(in), optional :: nn, ia, ja, ix, jx, incx, iy, jy, incy

  $2($3) :: alpha0, beta0
  character :: uplo0
  integer :: nn0, ia0, ja0, ix0, jx0, incx0, iy0, jy0, incy0

  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(alpha0, alpha, $4(1, kind=$3))
  _handle_inoptflag(beta0, beta, $4(0, kind=$3))
  _handle_inoptflag(nn0, nn, desca(M_))
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(ix0, ix, 1)
  _handle_inoptflag(jx0, jx, 1)
  _handle_inoptflag(incx0, incx, 1)
  _handle_inoptflag(iy0, iy, 1)
  _handle_inoptflag(jy0, jy, 1)
  _handle_inoptflag(incy0, incy, 1)
  call $5(uplo0, nn0, alpha0, aa, ia0, ja0, desca, xx, ix0, jx0, descx, &
      & incx0, beta0, yy, iy0, jy0, descy, incy0)

end subroutine pblasfx_$1
')


dnl ************************************************************************
dnl *** psyr/pher
dnl ************************************************************************

define(`_subroutine_pblasfx_psyr_pher',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl $3 dummy arguments kind
dnl $4 conversion function
dnl $5 pblas subroutine name
subroutine pblasfx_$1(xx, descx, aa, desca, uplo, nn, alpha, ix, jx, incx,&
    & ia, ja)
  $2($3), intent(in) :: xx(:,:)
  integer, intent(in) :: descx(DLEN_)
  $2($3), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  character, intent(in), optional :: uplo
  integer, intent(in), optional :: nn
  real($3), intent(in), optional :: alpha
  integer, intent(in), optional :: ix, jx, incx, ia, ja

  real($3) :: alpha0
  character :: uplo0
  integer :: nn0, ix0, jx0, incx0, ia0, ja0

  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(nn0, nn, desca(M_))
  _handle_inoptflag(alpha0, alpha, real(1, kind=$3))
  _handle_inoptflag(ix0, ix, 1)
  _handle_inoptflag(jx0, jx, 1)
  _handle_inoptflag(incx0, incx, 1)
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  call $5(uplo0, nn0, alpha0, xx, ix0, jx0, descx, incx0, aa, ia0, ja0, desca)
  
end subroutine pblasfx_$1
')


dnl ************************************************************************
dnl *** psyrk/pherk
dnl ************************************************************************

define(`_subroutine_pblasfx_psyrk_pherk',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl $3 dummy arguments kind
dnl $4 conversion function
dnl $5 pblas subroutine name
!> Symmetric/Hermitian rank-k update.
!! \param aa  Matrix to update with.
!! \param desca  Descriptor of aa.
!! \param cc  Matrix to be updated.
!! \param desccc Descriptor of cc.
!! \param uplo "U" for for upper, "L" for lower triangle matrix (default: "L").
!! \param trans  "N" for normal, "T" for transposed aa (default: "N").
!! \param alpha  Prefactor.
subroutine pblasfx_$1(aa, desca, cc, descc, uplo, trans, alpha, beta,&
    & nn, kk, ia, ja, ic, jc)
  $2($3), intent(in) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  $2($3), intent(inout) :: cc(:,:)
  integer, intent(in) :: descc(DLEN_)
  character, intent(in), optional :: uplo, trans
  real($3), intent(in), optional :: alpha, beta
  integer, intent(in), optional :: nn, kk
  integer, intent(in), optional :: ia, ja, ic, jc

  real($3) :: alpha0, beta0
  character :: uplo0, trans0
  integer :: nn0, kk0, ia0, ja0, ic0, jc0

  _handle_inoptflag(alpha0, alpha, real(1, kind=$3))
  _handle_inoptflag(beta0, beta, real(0, kind=$3))
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(trans0, trans, "N")
  if (trans0 == "N") then
    _handle_inoptflag(nn0, nn, desca(M_))
    _handle_inoptflag(kk0, kk, desca(N_))
  else
    _handle_inoptflag(nn0, nn, desca(N_))
    _handle_inoptflag(kk0, kk, desca(M_))
  end if
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(ic0, ic, 1)
  _handle_inoptflag(jc0, jc, 1)
  call $5(uplo0, trans0, nn0, kk0, alpha0, aa, ia0, ja0, desca, beta0,&
      & cc, ic0, jc0, descc)
  
end subroutine pblasfx_$1
')

dnl ************************************************************************
dnl *** ptran
dnl ************************************************************************

define(`_subroutine_pblasfx_ptran',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl $3 dummy arguments kind
dnl $4 conversion function
dnl $5 pblas subroutine name
!> Real matrix transpose.
!! \param aa  Matrix to update with.
!! \param desca  Descriptor of aa.
!! \param cc  Matrix to be updated.
!! \param desccc Descriptor of cc.
!! \param alpha  Prefactor.
!! \param beta  Prefactor.
subroutine pblasfx_$1(aa, desca, cc, descc, alpha, beta,&
    & mm, nn, ia, ja, ic, jc)
  $2($3), intent(in) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  $2($3), intent(inout) :: cc(:,:)
  integer, intent(in) :: descc(DLEN_)
  real($3), intent(in), optional :: alpha, beta
  integer, intent(in), optional :: mm, nn
  integer, intent(in), optional :: ia, ja, ic, jc

  real($3) :: alpha0, beta0
  integer :: mm0, nn0, ia0, ja0, ic0, jc0

  _handle_inoptflag(alpha0, alpha, real(1, kind=$3))
  _handle_inoptflag(beta0, beta, real(0, kind=$3))

  _handle_inoptflag(mm0, mm, desca(M_))
  _handle_inoptflag(nn0, nn, desca(N_))

  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(ic0, ic, 1)
  _handle_inoptflag(jc0, jc, 1)
  call $5(mm0, nn0, alpha0, aa, ia0, ja0, desca, beta0,&
      & cc, ic0, jc0, descc)

end subroutine pblasfx_$1
')

dnl ************************************************************************
dnl *** ptranu
dnl ************************************************************************

define(`_subroutine_pblasfx_ptranu',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl $3 dummy arguments kind
dnl $4 conversion function
dnl $5 pblas subroutine name
!> Complex matrix transpose.
!! \param aa  Matrix to update with.
!! \param desca  Descriptor of aa.
!! \param cc  Matrix to be updated.
!! \param desccc Descriptor of cc.
!! \param alpha  Prefactor.
!! \param beta  Prefactor.
subroutine pblasfx_$1(aa, desca, cc, descc, alpha, beta,&
    & mm, nn, ia, ja, ic, jc)
  $2($3), intent(in) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  $2($3), intent(inout) :: cc(:,:)
  integer, intent(in) :: descc(DLEN_)
  complex($3), intent(in), optional :: alpha, beta
  integer, intent(in), optional :: mm, nn
  integer, intent(in), optional :: ia, ja, ic, jc

  complex($3) :: alpha0, beta0
  integer :: mm0, nn0, ia0, ja0, ic0, jc0

  _handle_inoptflag(alpha0, alpha, cmplx(1, 0, kind=$3))
  _handle_inoptflag(beta0, beta, cmplx(0, 0, kind=$3))

  _handle_inoptflag(mm0, mm, desca(M_))
  _handle_inoptflag(nn0, nn, desca(N_))

  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(ic0, ic, 1)
  _handle_inoptflag(jc0, jc, 1)
  call $5(mm0, nn0, alpha0, aa, ia0, ja0, desca, beta0,&
      & cc, ic0, jc0, descc)

end subroutine pblasfx_$1
')

dnl ************************************************************************
dnl *** ptranc
dnl ************************************************************************

define(`_subroutine_pblasfx_ptranc',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl $3 dummy arguments kind
dnl $4 conversion function
dnl $5 pblas subroutine name
!> Complex matrix hermitian transpose.
!! \param aa  Matrix to update with.
!! \param desca  Descriptor of aa.
!! \param cc  Matrix to be updated.
!! \param desccc Descriptor of cc.
!! \param alpha  Prefactor.
!! \param beta  Prefactor.
subroutine pblasfx_$1(aa, desca, cc, descc, alpha, beta,&
    & mm, nn, ia, ja, ic, jc)
  $2($3), intent(in) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  $2($3), intent(inout) :: cc(:,:)
  integer, intent(in) :: descc(DLEN_)
  complex($3), intent(in), optional :: alpha, beta
  integer, intent(in), optional :: mm, nn
  integer, intent(in), optional :: ia, ja, ic, jc

  complex($3) :: alpha0, beta0
  integer :: mm0, nn0, ia0, ja0, ic0, jc0

  _handle_inoptflag(alpha0, alpha, cmplx(1, 0, kind=$3))
  _handle_inoptflag(beta0, beta, cmplx(0, 0, kind=$3))

  _handle_inoptflag(mm0, mm, desca(M_))
  _handle_inoptflag(nn0, nn, desca(N_))

  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(ic0, ic, 1)
  _handle_inoptflag(jc0, jc, 1)
  call $5(mm0, nn0, alpha0, aa, ia0, ja0, desca, beta0,&
      & cc, ic0, jc0, descc)

end subroutine pblasfx_$1
')
