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
dnl *** psymm/phemm
dnl ************************************************************************

define(`_subroutine_pblasfx_psymm_phemm',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl $3 dummy arguments kind
dnl $4 conversion function
dnl $5 pblas subroutine name
!> Symmetric/Hermitian matrix with general matrix product 
!! \param aa  Symmetric/Hermitian matrix.
!! \param desca  Descriptor of aa.
!! \param bb  general matrix.
!! \param descb  Descriptor of bb.
!! \param cc  Matrix to store result
!! \param descc  Descriptor of cc.
!! \param side "L" for for left, "R" for right (default: "L"), 
!!        if "L" C := alpha * A * B + beta*C
!!        if "R" C := alpha * B * A + beta*C
!! \param uplo "U" for for upper, "L" for lower triangle matrix (default: "L").
!! \param alpha  Prefactor.
!! \param beta  Prefactor.
subroutine pblasfx_$1(aa, desca, bb, descb, cc, descc, side, uplo, &
    & alpha, beta, mm, nn, ia, ja, ib, jb, ic, jc)
  $2($3), intent(in) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  $2($3), intent(in) :: bb(:,:)
  integer, intent(in) :: descb(DLEN_)
  $2($3), intent(inout) :: cc(:,:)
  integer, intent(in) :: descc(DLEN_)
  character, intent(in), optional :: side
  character, intent(in), optional :: uplo
  $2($3), intent(in), optional :: alpha, beta
  integer, intent(in), optional :: mm, nn, ia, ja, ib, jb, ic, jc
  
  $2($3) :: alpha0, beta0
  character :: side0, uplo0
  integer :: mm0, nn0, ia0, ja0, ib0, jb0, ic0, jc0
  
  _handle_inoptflag(side0, side, "L")
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(alpha0, alpha, $4(1, kind=$3))
  _handle_inoptflag(beta0, beta, $4(0, kind=$3))
  _handle_inoptflag(mm0, mm, desca(M_))
  _handle_inoptflag(nn0, nn, desca(N_))
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(ib0, ib, 1)
  _handle_inoptflag(jb0, jb, 1)
  _handle_inoptflag(ic0, ic, 1)
  _handle_inoptflag(jc0, jc, 1)
  call $5(side0, uplo0, mm0, nn0, alpha0, aa, ia0, ja0, desca, &
      & bb, ib0, jb0, descb, beta0, cc, ic0, jc0, descc)
end subroutine pblasfx_$1
')


dnl ************************************************************************
dnl *** pgemm
dnl ************************************************************************

define(`_subroutine_pblasfx_pgemm',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl $3 dummy arguments kind
dnl $4 conversion function
!> Matrix matrix product: alpha * A * B + beta * C.
!!
!! \see PBLAS documentation (p?gemm routines)
!!
subroutine pblasfx_pgemm_$1(aa, desca, bb, descb, cc, descc, alpha, beta, &
    & transa, transb, ia, ja, ib, jb, ic, jc, mm, nn, kk)

  !> Left operand matrix A.
  $2($3), intent(in) :: aa(:,:)

  !> Descriptor of A.
  integer, intent(in) :: desca(DLEN_)

  !> Right operand matrix B.
  $2($3), intent(in) :: bb(:,:)

  !> Descriptor of B.
  integer, intent(in) :: descb(DLEN_)

  !> Added matrix C.
  $2($3), intent(inout) :: cc(:,:)

  !> Descriptor of C.
  integer, intent(in) :: descc(DLEN_)

  !> Prefactor alpha (alpha * A * B). Default: 1.0
  $2($3), intent(in), optional :: alpha

  !> Prefactor beta (beta * C). Default: 0.0
  $2($3), intent(in), optional :: beta

  !> Whether A should be unchanged ("N"), transposed ("T") or transposed
  !! conjugated ("C"). Default: "N".
  character, intent(in), optional :: transa

  !> Whether B should be unchanged ("N"), transposed ("T") or transposed
  !! conjugated ("C"). Default: "N".
  character, intent(in), optional :: transb

  !> First row of submatrix of A. Default: 1
  integer, intent(in), optional :: ia

  !> First column of submatrix of A. Default: 1
  integer, intent(in), optional :: ja
  
  !> First row of submatrix of B. Default: 1
  integer, intent(in), optional :: ib

  !> First column of submatrix of B. Default: 1
  integer, intent(in), optional :: jb

  !> First row of submatrix of C. Default: 1
  integer, intent(in), optional :: ic

  !> First column of submatrix of C. Default: 1
  integer, intent(in), optional :: jc

  !> Number of rows in the submatrix of A and C. Default: desca(M_)
  integer, intent(in), optional :: mm

  !> Number of colums in the submatrix of B and C. Default: descb(N_)
  integer, intent(in), optional :: nn

  !> Number of columns/rows in the submatrix A, B. Default: desca(N_)
  integer, intent(in), optional :: kk

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  $2($3) :: alpha0, beta0
  integer :: ia0, ja0, ib0, jb0, ic0, jc0, mm0, nn0, kk0
  character :: transa0, transb0
  
  _handle_inoptflag(alpha0, alpha, $4(1, kind=$3))
  _handle_inoptflag(beta0, beta, $4(0, kind=$3))
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(ib0, ib, 1)
  _handle_inoptflag(jb0, jb, 1)
  _handle_inoptflag(ic0, ic, 1)
  _handle_inoptflag(jc0, jc, 1)
  _handle_inoptflag(mm0, mm, desca(M_))
  _handle_inoptflag(nn0, nn, descb(N_))
  _handle_inoptflag(kk0, kk, desca(N_))
  _handle_inoptflag(transa0, transa, "N")
  _handle_inoptflag(transb0, transb, "N")

  call pgemm(transa0, transb0, mm0, nn0, kk0, alpha0, aa, ia0, ja0, desca, &
      & bb, ib0, jb0, descb, beta0, cc, ic0, jc0, descc)

end subroutine pblasfx_pgemm_$1
')

dnl ************************************************************************
dnl *** ptrmm
dnl ************************************************************************

define(`_subroutine_pblasfx_ptrmm',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl $3 dummy arguments kind
dnl $4 conversion function
dnl $5 pblas subroutine name
!> Computes matrix-matrix product with one triangle matrix
!!
!! \see PBLAS documentation (routines p?trmm)
!!
subroutine pblasfx_ptrmm_$1(aa, desca, bb, descb, alpha, side, uplo, transa, &
    & diag, ia, ja, ib, jb, mm, nn)

  !> Unit or non-unit lower or upper triangular matrix A
  $2($3), intent(in) :: aa(:,:)

  !> Descriptor of A.
  integer, intent(in) :: desca(DLEN_)

  !> Second operand (general matrix) B on entry, result on exit.
  $2($3), intent(inout) :: bb(:,:)

  !> Descriptor of B.
  integer, intent(in) :: descb(DLEN_)

  !> Prefactor. Default: 1.0
  $2($3), intent(in), optional :: alpha

  !> From which side is B multiplied by A ("L"/"R"). Default: "L"
  character, intent(in), optional :: side

  !> Whether A is upper ("U") or lower("L") triangle. Default: "L".
  character, intent(in), optional :: uplo

  !> Whether A should be unchanged ("N"), transposed ("T") or transposed
  !! conjugated ("C"). Default: "N".
  character, intent(in), optional :: transa

  !> Whether A is unit triangular ("U") or not ("N"). Default: "N".
  character, intent(in), optional :: diag

  !> First row of matrix A to consider. Default: 1
  integer, intent(in), optional :: ia

  !> First column of matrix A to consider. Default: 1
  integer, intent(in), optional :: ja
  
  !> First row of matrix B to consider. Default: 1
  integer, intent(in), optional :: ib

  !> First column of matrix B to consider. Default: 1
  integer, intent(in), optional :: jb

  !> Number of rows for matrix A. Default: desca(M_)
  integer, intent(in), optional :: mm
  
  !> Number of columns for matrix A. Default: desca(N_)
  integer, intent(in), optional :: nn

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  integer :: ia0, ja0, ib0, jb0, mm0, nn0
  $2($3) :: alpha0
  character :: side0, uplo0, transa0, diag0

  _handle_inoptflag(alpha0, alpha, $4(1, kind=$3))
  _handle_inoptflag(side0, side, "L")
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(transa0, transa, "N")
  _handle_inoptflag(diag0, diag, "N")
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(ib0, ib, 1)
  _handle_inoptflag(jb0, jb, 1)
  _handle_inoptflag(mm0, mm, desca(M_))
  _handle_inoptflag(nn0, nn, desca(N_))

  call ptrmm(side0, uplo0, transa0, diag0, mm0, nn0, alpha0, aa, ia0, ja0, &
      & desca, bb, ib0, jb0, descb)

end subroutine pblasfx_ptrmm_$1
')

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
