#:include 'scalapackfx.fypp'
#:set TYPES = NUMERIC_TYPES
#:set RANKS = [2, 1, 0]


#! ************************************************************************
#! *** psymv/phemv
#! ************************************************************************

#:def pblasfx_psymv_phemv_template(SUFFIX, TYPE, KIND, FUNCTION, NAME)

  subroutine pblasfx_${SUFFIX}$(aa, desca, xx, descx, yy, descy, uplo, alpha, beta, &
      & nn, ia, ja, ix, jx, incx, iy, jy, incy)
    ${TYPE}$(${KIND}$), intent(in) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    ${TYPE}$(${KIND}$), intent(in) :: xx(:,:)
    integer, intent(in) :: descx(DLEN_)
    ${TYPE}$(${KIND}$), intent(inout) :: yy(:,:)
    integer, intent(in) :: descy(DLEN_)
    character, intent(in), optional :: uplo
    ${TYPE}$(${KIND}$), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: nn, ia, ja, ix, jx, incx, iy, jy, incy

    ${TYPE}$(${KIND}$) :: alpha0, beta0
    character :: uplo0
    integer :: nn0, ia0, ja0, ix0, jx0, incx0, iy0, jy0, incy0

    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(alpha0, alpha, ${FUNCTION}$(1, kind=${KIND}$))
    @:inoptflags(beta0, beta, ${FUNCTION}$(0, kind=${KIND}$))
    @:inoptflags(nn0, nn, desca(M_))
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(ix0, ix, 1)
    @:inoptflags(jx0, jx, 1)
    @:inoptflags(incx0, incx, 1)
    @:inoptflags(iy0, iy, 1)
    @:inoptflags(jy0, jy, 1)
    @:inoptflags(incy0, incy, 1)
    call ${NAME}$(uplo0, nn0, alpha0, aa, ia0, ja0, desca, xx, ix0, jx0, descx, &
        & incx0, beta0, yy, iy0, jy0, descy, incy0)

  end subroutine pblasfx_${SUFFIX}$

#:enddef pblasfx_psymv_phemv_template



#!************************************************************************
#!*** psyr/pher
#!************************************************************************

#:def pblasfx_psyr_pher_template(SUFFIX, TYPE, KIND, FUNCTION, NAME)

  subroutine pblasfx_${SUFFIX}$(xx, descx, aa, desca, uplo, nn, alpha, ix, jx, incx,&
      & ia, ja)
    ${TYPE}$(${KIND}$), intent(in) :: xx(:,:)
    integer, intent(in) :: descx(DLEN_)
    ${TYPE}$(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    character, intent(in), optional :: uplo
    integer, intent(in), optional :: nn
    real(${KIND}$), intent(in), optional :: alpha
    integer, intent(in), optional :: ix, jx, incx, ia, ja

    real(${KIND}$) :: alpha0
    character :: uplo0
    integer :: nn0, ix0, jx0, incx0, ia0, ja0

    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(nn0, nn, desca(M_))
    @:inoptflags(alpha0, alpha, real(1, kind=${KIND}$))
    @:inoptflags(ix0, ix, 1)
    @:inoptflags(jx0, jx, 1)
    @:inoptflags(incx0, incx, 1)
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    call ${NAME}$(uplo0, nn0, alpha0, xx, ix0, jx0, descx, incx0, aa, ia0, ja0, desca)

  end subroutine pblasfx_${SUFFIX}$

#:enddef pblasfx_psyr_pher_template



#!************************************************************************
#!*** psyrk/pherk
#!************************************************************************

#:def pblasfx_psyrk_pherk_template(SUFFIX, TYPE, KIND, FUNCTION, NAME)

  !> Symmetric/Hermitian rank-k update.
  !! \param aa  Matrix to update with.
  !! \param desca  Descriptor of aa.
  !! \param cc  Matrix to be updated.
  !! \param desccc Descriptor of cc.
  !! \param uplo "U" for for upper, "L" for lower triangle matrix (default: "L").
  !! \param trans  "N" for normal, "T" for transposed aa (default: "N").
  !! \param alpha  Prefactor.
  subroutine pblasfx_${SUFFIX}$(aa, desca, cc, descc, uplo, trans, alpha, beta,&
      & nn, kk, ia, ja, ic, jc)
    ${TYPE}$(${KIND}$), intent(in) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    ${TYPE}$(${KIND}$), intent(inout) :: cc(:,:)
    integer, intent(in) :: descc(DLEN_)
    character, intent(in), optional :: uplo, trans
    real(${KIND}$), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: nn, kk
    integer, intent(in), optional :: ia, ja, ic, jc

    real(${KIND}$) :: alpha0, beta0
    character :: uplo0, trans0
    integer :: nn0, kk0, ia0, ja0, ic0, jc0

    @:inoptflags(alpha0, alpha, real(1, kind=${KIND}$))
    @:inoptflags(beta0, beta, real(0, kind=${KIND}$))
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(trans0, trans, "N")
    if (trans0 == "N") then
      @:inoptflags(nn0, nn, desca(M_))
      @:inoptflags(kk0, kk, desca(N_))
    else
      @:inoptflags(nn0, nn, desca(N_))
      @:inoptflags(kk0, kk, desca(M_))
    end if
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(ic0, ic, 1)
    @:inoptflags(jc0, jc, 1)
    call ${NAME}$(uplo0, trans0, nn0, kk0, alpha0, aa, ia0, ja0, desca, beta0,&
        & cc, ic0, jc0, descc)

  end subroutine pblasfx_${SUFFIX}$

#:enddef pblasfx_psyrk_pherk_template




#! ************************************************************************
#! *** psymm/phemm
#! ************************************************************************

#:def pblasfx_psymm_phemm_template(SUFFIX, TYPE, KIND, FUNCTION, NAME)

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
  subroutine pblasfx_${SUFFIX}$(aa, desca, bb, descb, cc, descc, side, uplo, &
      & alpha, beta, mm, nn, ia, ja, ib, jb, ic, jc)
    ${TYPE}$(${KIND}$), intent(in) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    ${TYPE}$(${KIND}$), intent(in) :: bb(:,:)
    integer, intent(in) :: descb(DLEN_)
    ${TYPE}$(${KIND}$), intent(inout) :: cc(:,:)
    integer, intent(in) :: descc(DLEN_)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    ${TYPE}$(${KIND}$), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: mm, nn, ia, ja, ib, jb, ic, jc

    ${TYPE}$(${KIND}$) :: alpha0, beta0
    character :: side0, uplo0
    integer :: mm0, nn0, ia0, ja0, ib0, jb0, ic0, jc0

    @:inoptflags(side0, side, "L")
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(alpha0, alpha, ${FUNCTION}$(1, kind=${KIND}$))
    @:inoptflags(beta0, beta, ${FUNCTION}$(0, kind=${KIND}$))
    @:inoptflags(mm0, mm, desca(M_))
    @:inoptflags(nn0, nn, desca(N_))
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(ib0, ib, 1)
    @:inoptflags(jb0, jb, 1)
    @:inoptflags(ic0, ic, 1)
    @:inoptflags(jc0, jc, 1)
    call ${NAME}$(side0, uplo0, mm0, nn0, alpha0, aa, ia0, ja0, desca, &
        & bb, ib0, jb0, descb, beta0, cc, ic0, jc0, descc)
  end subroutine pblasfx_${SUFFIX}$

#:enddef pblasfx_psymm_phemm_template




#!************************************************************************
#!*** pgemm
#!************************************************************************

#:def pblasfx_pgemm_template(SUFFIX, TYPE, KIND, FUNCTION)

  !> Matrix matrix product: alpha * A * B + beta * C.
  !!
  !! \see PBLAS documentation (p?gemm routines)
  !!
  subroutine pblasfx_pgemm_${SUFFIX}$(aa, desca, bb, descb, cc, descc, alpha, beta, &
      & transa, transb, ia, ja, ib, jb, ic, jc, mm, nn, kk)

    !> Left operand matrix A.
    ${TYPE}$(${KIND}$), intent(in) :: aa(:,:)

    !> Descriptor of A.
    integer, intent(in) :: desca(DLEN_)

    !> Right operand matrix B.
    ${TYPE}$(${KIND}$), intent(in) :: bb(:,:)

    !> Descriptor of B.
    integer, intent(in) :: descb(DLEN_)

    !> Added matrix C.
    ${TYPE}$(${KIND}$), intent(inout) :: cc(:,:)

    !> Descriptor of C.
    integer, intent(in) :: descc(DLEN_)

    !> Prefactor alpha (alpha * A * B). Default: 1.0
    ${TYPE}$(${KIND}$), intent(in), optional :: alpha

    !> Prefactor beta (beta * C). Default: 0.0
    ${TYPE}$(${KIND}$), intent(in), optional :: beta

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

    ${TYPE}$(${KIND}$) :: alpha0, beta0
    integer :: ia0, ja0, ib0, jb0, ic0, jc0, mm0, nn0, kk0
    character :: transa0, transb0

    @:inoptflags(alpha0, alpha, ${FUNCTION}$(1, kind=${KIND}$))
    @:inoptflags(beta0, beta, ${FUNCTION}$(0, kind=${KIND}$))
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(ib0, ib, 1)
    @:inoptflags(jb0, jb, 1)
    @:inoptflags(ic0, ic, 1)
    @:inoptflags(jc0, jc, 1)
    @:inoptflags(mm0, mm, desca(M_))
    @:inoptflags(nn0, nn, descb(N_))
    @:inoptflags(kk0, kk, desca(N_))
    @:inoptflags(transa0, transa, "N")
    @:inoptflags(transb0, transb, "N")

    call pgemm(transa0, transb0, mm0, nn0, kk0, alpha0, aa, ia0, ja0, desca, &
        & bb, ib0, jb0, descb, beta0, cc, ic0, jc0, descc)

  end subroutine pblasfx_pgemm_${SUFFIX}$

#:enddef pblasfx_pgemm_template


#!************************************************************************
#!*** ptrmm
#!************************************************************************

#:def pblasfx_ptrmm_template(SUFFIX, TYPE, KIND, FUNCTION)

  !> Computes matrix-matrix product with one triangle matrix
  !!
  !! \see PBLAS documentation (routines p?trmm)
  !!
  subroutine pblasfx_ptrmm_${SUFFIX}$(aa, desca, bb, descb, alpha, side, uplo, transa, &
      & diag, ia, ja, ib, jb, mm, nn)

    !> Unit or non-unit lower or upper triangular matrix A
    ${TYPE}$(${KIND}$), intent(in) :: aa(:,:)

    !> Descriptor of A.
    integer, intent(in) :: desca(DLEN_)

    !> Second operand (general matrix) B on entry, result on exit.
    ${TYPE}$(${KIND}$), intent(inout) :: bb(:,:)

    !> Descriptor of B.
    integer, intent(in) :: descb(DLEN_)

    !> Prefactor. Default: 1.0
    ${TYPE}$(${KIND}$), intent(in), optional :: alpha

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
    ${TYPE}$(${KIND}$) :: alpha0
    character :: side0, uplo0, transa0, diag0

    @:inoptflags(alpha0, alpha, ${FUNCTION}$(1, kind=${KIND}$))
    @:inoptflags(side0, side, "L")
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(transa0, transa, "N")
    @:inoptflags(diag0, diag, "N")
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(ib0, ib, 1)
    @:inoptflags(jb0, jb, 1)
    @:inoptflags(mm0, mm, desca(M_))
    @:inoptflags(nn0, nn, desca(N_))

    call ptrmm(side0, uplo0, transa0, diag0, mm0, nn0, alpha0, aa, ia0, ja0, &
        & desca, bb, ib0, jb0, descb)

  end subroutine pblasfx_ptrmm_${SUFFIX}$

#:enddef pblasfx_ptrmm_template


#!************************************************************************
#!*** ptran
#!************************************************************************

#:def pblasfx_ptran_template(SUFFIX, TYPE, KIND, FUNCTION, NAME)

  !> Real matrix transpose.
  !! \param aa  Matrix to update with.
  !! \param desca  Descriptor of aa.
  !! \param cc  Matrix to be updated.
  !! \param desccc Descriptor of cc.
  !! \param alpha  Prefactor.
  !! \param beta  Prefactor.
  subroutine pblasfx_${SUFFIX}$(aa, desca, cc, descc, alpha, beta,&
      & mm, nn, ia, ja, ic, jc)
    ${TYPE}$(${KIND}$), intent(in) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    ${TYPE}$(${KIND}$), intent(inout) :: cc(:,:)
    integer, intent(in) :: descc(DLEN_)
    real(${KIND}$), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: mm, nn
    integer, intent(in), optional :: ia, ja, ic, jc

    real(${KIND}$) :: alpha0, beta0
    integer :: mm0, nn0, ia0, ja0, ic0, jc0

    @:inoptflags(alpha0, alpha, real(1, kind=${KIND}$))
    @:inoptflags(beta0, beta, real(0, kind=${KIND}$))

    @:inoptflags(mm0, mm, desca(M_))
    @:inoptflags(nn0, nn, desca(N_))

    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(ic0, ic, 1)
    @:inoptflags(jc0, jc, 1)
    call ${NAME}$(mm0, nn0, alpha0, aa, ia0, ja0, desca, beta0,&
        & cc, ic0, jc0, descc)
  end subroutine pblasfx_${SUFFIX}$

#:enddef pblasfx_ptran_template



#!************************************************************************
#!*** ptranu
#!************************************************************************

#:def pblasfx_ptranu_template(SUFFIX, TYPE, KIND, FUNCTION, NAME)

  !> Complex matrix transpose.
  !! \param aa  Matrix to update with.
  !! \param desca  Descriptor of aa.
  !! \param cc  Matrix to be updated.
  !! \param desccc Descriptor of cc.
  !! \param alpha  Prefactor.
  !! \param beta  Prefactor.
  subroutine pblasfx_${SUFFIX}$(aa, desca, cc, descc, alpha, beta,&
      & mm, nn, ia, ja, ic, jc)
    ${TYPE}$(${KIND}$), intent(in) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    ${TYPE}$(${KIND}$), intent(inout) :: cc(:,:)
    integer, intent(in) :: descc(DLEN_)
    complex(${KIND}$), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: mm, nn
    integer, intent(in), optional :: ia, ja, ic, jc

    complex(${KIND}$) :: alpha0, beta0
    integer :: mm0, nn0, ia0, ja0, ic0, jc0

    @:inoptflags(alpha0, alpha, cmplx(1, 0, kind=${KIND}$))
    @:inoptflags(beta0, beta, cmplx(0, 0, kind=${KIND}$))

    @:inoptflags(mm0, mm, desca(M_))
    @:inoptflags(nn0, nn, desca(N_))

    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(ic0, ic, 1)
    @:inoptflags(jc0, jc, 1)
    call ${NAME}$(mm0, nn0, alpha0, aa, ia0, ja0, desca, beta0,&
        & cc, ic0, jc0, descc)

  end subroutine pblasfx_${SUFFIX}$

#:enddef pblasfx_ptranu_template



#!************************************************************************
#!*** ptranc
#!************************************************************************

#:def pblasfx_ptranc_template(SUFFIX, TYPE, KIND, FUNCTION, NAME)

  !> Complex matrix hermitian transpose.
  !! \param aa  Matrix to update with.
  !! \param desca  Descriptor of aa.
  !! \param cc  Matrix to be updated.
  !! \param desccc Descriptor of cc.
  !! \param alpha  Prefactor.
  !! \param beta  Prefactor.
  subroutine pblasfx_${SUFFIX}$(aa, desca, cc, descc, alpha, beta,&
      & mm, nn, ia, ja, ic, jc)
    ${TYPE}$(${KIND}$), intent(in) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    ${TYPE}$(${KIND}$), intent(inout) :: cc(:,:)
    integer, intent(in) :: descc(DLEN_)
    complex(${KIND}$), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: mm, nn
    integer, intent(in), optional :: ia, ja, ic, jc

    complex(${KIND}$) :: alpha0, beta0
    integer :: mm0, nn0, ia0, ja0, ic0, jc0

    @:inoptflags(alpha0, alpha, cmplx(1, 0, kind=${KIND}$))
    @:inoptflags(beta0, beta, cmplx(0, 0, kind=${KIND}$))

    @:inoptflags(mm0, mm, desca(M_))
    @:inoptflags(nn0, nn, desca(N_))

    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(ic0, ic, 1)
    @:inoptflags(jc0, jc, 1)
    call ${NAME}$(mm0, nn0, alpha0, aa, ia0, ja0, desca, beta0,&
        & cc, ic0, jc0, descc)

  end subroutine pblasfx_${SUFFIX}$

#:enddef pblasfx_ptranc_template

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

  @:pblasfx_psyr_pher_template(psyr_real, real, sp, real, psyr)
  @:pblasfx_psyr_pher_template(psyr_dreal, real, dp, real, psyr)
  @:pblasfx_psyr_pher_template(pher_complex, complex, sp, cmplx, pher)
  @:pblasfx_psyr_pher_template(pher_dcomplex, complex, dp, cmplx, pher)

  @:pblasfx_psyrk_pherk_template(psyrk_real, real, sp, real, psyrk)
  @:pblasfx_psyrk_pherk_template(psyrk_dreal, real, dp, real, psyrk)
  @:pblasfx_psyrk_pherk_template(pherk_complex, complex, sp, cmplx, pherk)
  @:pblasfx_psyrk_pherk_template(pherk_dcomplex, complex, dp, cmplx, pherk)

  @:pblasfx_psymv_phemv_template(psymv_real, real, sp, real, psymv)
  @:pblasfx_psymv_phemv_template(psymv_dreal, real, dp, real, psymv)
  @:pblasfx_psymv_phemv_template(phemv_complex, complex, sp, cmplx, phemv)
  @:pblasfx_psymv_phemv_template(phemv_dcomplex, complex, dp, cmplx, phemv)

  @:pblasfx_psymm_phemm_template(psymm_real, real, sp, real, psymm)
  @:pblasfx_psymm_phemm_template(psymm_dreal, real, dp, real, psymm)
  @:pblasfx_psymm_phemm_template(phemm_complex, complex, sp, cmplx, phemm)
  @:pblasfx_psymm_phemm_template(phemm_dcomplex, complex, dp, cmplx, phemm)

  @:pblasfx_pgemm_template(real, real, sp, real)
  @:pblasfx_pgemm_template(dreal, real, dp, real)
  @:pblasfx_pgemm_template(complex, complex, sp, cmplx)
  @:pblasfx_pgemm_template(dcomplex, complex, dp, cmplx)

  @:pblasfx_ptrmm_template(real, real, sp, real)
  @:pblasfx_ptrmm_template(dreal, real, dp, real)
  @:pblasfx_ptrmm_template(complex, complex, sp, cmplx)
  @:pblasfx_ptrmm_template(dcomplex, complex, dp, cmplx)

  @:pblasfx_ptran_template(ptran_real, real, sp, real, ptran)
  @:pblasfx_ptran_template(ptran_dreal, real, dp, real, ptran)
  @:pblasfx_ptranu_template(ptranu_complex, complex, sp, complex, ptranu)
  @:pblasfx_ptranu_template(ptranu_dcomplex, complex, dp, complex, ptranu)
  @:pblasfx_ptranc_template(ptranc_complex, complex, sp, complex, ptranc)
  @:pblasfx_ptranc_template(ptranc_dcomplex, complex, dp, complex, ptranc)

end module pblasfx_module
