#:include 'scalapackfx.fypp'

!> High level Fortran wrappers for the SCALAPACK library.
module scalapackfx_module
  use scalapackfx_common_module
  use blacsfx_module
  use scalapack_module
  implicit none
  private

  public :: DLEN_, DT_, CTXT_, M_, N_, MB_, NB_, RSRC_, CSRC_, LLD_
  public :: scalafx_ppotrf
  public :: scalafx_ppotri
  public :: scalafx_ptrtri
  public :: scalafx_pgetrf
  public :: scalafx_psygst
  public :: scalafx_phegst
  public :: scalafx_psyev
  public :: scalafx_psygv
  public :: scalafx_pheev
  public :: scalafx_phegv
  public :: scalafx_psyevd
  public :: scalafx_psygvd
  public :: scalafx_pheevd
  public :: scalafx_phegvd
  public :: scalafx_psyevr
  public :: scalafx_pheevr
  public :: scalafx_psygvr
  public :: scalafx_phegvr
  public :: scalafx_ptrsm
  public :: scalafx_getdescriptor
  public :: scalafx_getlocalshape
  public :: scalafx_infog2l
  public :: scalafx_indxl2g
  public :: scalafx_localindices
  public :: scalafx_creatematrix
  public :: scalafx_pgesvd

  !> Cholesky factorization of a symmetric/Hermitian positive definite matrix.
  interface scalafx_ppotrf
    module procedure scalafx_ppotrf_real, scalafx_ppotrf_dreal
    module procedure scalafx_ppotrf_complex, scalafx_ppotrf_dcomplex
  end interface scalafx_ppotrf

  !> Inverse of a Cholesky decomposed symmetric/Hermitian matrix.
  interface scalafx_ppotri
    module procedure scalafx_ppotri_real, scalafx_ppotri_dreal
    module procedure scalafx_ppotri_complex, scalafx_ppotri_dcomplex
  end interface scalafx_ppotri

  !> Inverse of a triangular matrix
  interface scalafx_ptrtri
    module procedure scalafx_ptrtri_real, scalafx_ptrtri_dreal
    module procedure scalafx_ptrtri_complex, scalafx_ptrtri_dcomplex
  end interface scalafx_ptrtri

  !> LU decomposition of a general matrix with pivoting.
  interface scalafx_pgetrf
    module procedure scalafx_pgetrf_real, scalafx_pgetrf_dreal
    module procedure scalafx_pgetrf_complex, scalafx_pgetrf_dcomplex
  end interface scalafx_pgetrf

  !> Reduces symmetric definite generalized eigenvalue problem to standard form.
  interface scalafx_psygst
    module procedure scalafx_psygst_real, scalafx_psygst_dreal
  end interface scalafx_psygst

  !> Reduces Hermitian definite generalized eigenvalue problem to standard form.
  interface scalafx_phegst
    module procedure scalafx_phegst_complex, scalafx_phegst_dcomplex
  end interface scalafx_phegst

  !> Solves symmetric eigenvalue problem by the QR algorithm.
  interface scalafx_psyev
    module procedure scalafx_psyev_real, scalafx_psyev_dreal
  end interface scalafx_psyev

  !> Solves Hermitian eigenvalue problem by the QR algorithm.
  interface scalafx_pheev
    module procedure scalafx_pheev_complex, scalafx_pheev_dcomplex
  end interface scalafx_pheev

  !> Solves generalized symmetric eigenvalue problem by the QR algorithm.
  interface scalafx_psygv
    module procedure scalafx_psygv_real, scalafx_psygv_dreal
  end interface scalafx_psygv

  !> Solves generalized Hermitian eigenvalue problem by the QR algorithm.
  interface scalafx_phegv
    module procedure scalafx_phegv_complex, scalafx_phegv_dcomplex
  end interface scalafx_phegv

  !> Solves symmetric eigenvalue problem by the divide and conquer algorithm.
  interface scalafx_psyevd
    module procedure scalafx_psyevd_real, scalafx_psyevd_dreal
  end interface scalafx_psyevd

  !> Solves Hermitian eigenvalue problem by the divide and conquer algorithm.
  interface scalafx_pheevd
    module procedure scalafx_pheevd_complex, scalafx_pheevd_dcomplex
  end interface scalafx_pheevd

  !> Solves generalized symmetric eigenvalue problem by the divide and conquer
  !! algorithm.
  interface scalafx_psygvd
    module procedure scalafx_psygvd_real, scalafx_psygvd_dreal
  end interface scalafx_psygvd

  !> Solves generalized Hermitian eigenvalue problem by the divide and conquer
  !! algorithm.
  interface scalafx_phegvd
    module procedure scalafx_phegvd_complex, scalafx_phegvd_dcomplex
  end interface scalafx_phegvd

  !> Solves symmetric eigenvalue problem by the divide and conquer algorithm.
  interface scalafx_psyevr
    module procedure scalafx_psyevr_real, scalafx_psyevr_dreal
  end interface scalafx_psyevr

  !> Solves Hermitian eigenvalue problem by the divide and conquer algorithm.
  interface scalafx_pheevr
    module procedure scalafx_pheevr_complex, scalafx_pheevr_dcomplex
  end interface scalafx_pheevr

  !> Solves generalized symmetric eigenvalue problem by the divide and conquer
  !! algorithm.
  interface scalafx_psygvr
    module procedure scalafx_psygvr_real, scalafx_psygvr_dreal
  end interface scalafx_psygvr

  !> Solves generalized Hermitian eigenvalue problem by the divide and conquer
  !! algorithm.
  interface scalafx_phegvr
    module procedure scalafx_phegvr_complex, scalafx_phegvr_dcomplex
  end interface scalafx_phegvr

  !> Singular value decomposition
  interface scalafx_pgesvd
    module procedure scalafx_r_pgesvd_real, scalafx_r_pgesvd_dreal
    module procedure scalafx_c_pgesvd_complex, scalafx_c_pgesvd_dcomplex
  end interface scalafx_pgesvd

  !> Solves triangular matrix equation.
  interface scalafx_ptrsm
    module procedure scalafx_ptrsm_real, scalafx_ptrsm_dreal
    module procedure scalafx_ptrsm_complex, scalafx_ptrsm_dcomplex
  end interface scalafx_ptrsm

  !> Creates a distributed matrix and allocates local storage.
  interface scalafx_creatematrix
    module procedure scalafx_creatematrix_int
    module procedure scalafx_creatematrix_real, scalafx_creatematrix_dreal
    module procedure scalafx_creatematrix_complex, scalafx_creatematrix_dcomplex
  end interface scalafx_creatematrix

!************************************************************************
!*** ppotrf
!************************************************************************

#:def scalafx_ppotrf_template(TYPE, FTYPE)

  !> Computes the Cholesky factorization of a Hermitian positive definite matrix.
  !!
  !! \param desca  Descriptor of the matrix.
  !! \param aa  Matrix.
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param ia  First row of the submatrix (default: 1)
  !! \param ja  First column of the submatrix (default: 1)
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !!
  !! \see SCALAPACK documentation (routine p?potrf).
  !!
  subroutine scalafx_ppotrf_${TYPE}$(aa, desca, uplo, nn, ia, ja, info)
    ${FTYPE}$, intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    character, intent(in), optional :: uplo
    integer, intent(in), optional :: nn, ia, ja
    integer, intent(out), optional :: info

    character :: uplo0
    integer :: nn0, ia0, ja0, info0

    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(nn0, nn, desca(M_))
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    call ppotrf(uplo0, nn0, aa, ia0, ja0, desca, info0)
    call handle_infoflag(info0, "ppotrf in scalafx_ppotrf_${TYPE}$", info)

  end subroutine scalafx_ppotrf_${TYPE}$

#:enddef scalafx_ppotrf_template


!************************************************************************
!*** ppotri
!************************************************************************

#:def scalafx_ppotri_template(TYPE, FTYPE)

  !> Computes the inverse of a symmetric/Hermitian positive definite matrix.
  !!
  subroutine scalafx_ppotri_${TYPE}$(aa, desca, uplo, ia, ja, nn, info)

    !> Cholesky decomposed matrix A on entry, inverse on exit.
    ${FTYPE}$, intent(inout) :: aa(:,:)

    !> Descriptor of A.
    integer, intent(in) :: desca(DLEN_)

    !> Specifies whether lower ("L") or upper ("U") part of A contains the
    !! matrix. Default: "L".
    character, intent(in), optional :: uplo

    !> First row of the submatrix of A. Default: 1
    integer, intent(in), optional :: ia

    !> First column of the submatrix of A. Default: 1
    integer, intent(in), optional :: ja

    !> Number of columns in the submatrix of A. Default: desca(M_)
    integer, intent(in), optional :: nn

    !> Info flag. If not specified and error occurs, the subroutine stops.
    integer, intent(out), optional :: info

    !------------------------------------------------------------------------

    character :: uplo0
    integer :: ia0, ja0, nn0, info0

    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(nn0, nn, desca(M_))
    call ppotri(uplo0, nn0, aa, ia0, ja0, desca, info0)
    call handle_infoflag(info0, "ppotri in scalafx_ppotri_${TYPE}$", info)

  end subroutine scalafx_ppotri_${TYPE}$

#:enddef scalafx_ppotri_template


!************************************************************************
!*** ptrtri
!************************************************************************

#:def scalafx_ptrtri_template(TYPE, FTYPE)

  !> Computes the inverse of a symmetric/Hermitian positive definite matrix.
  !!
  subroutine scalafx_ptrtri_${TYPE}$(aa, desca, uplo, diag, ia, ja, nn, info)

    !> Cholesky decomposed matrix A on entry, inverse on exit.
    ${FTYPE}$, intent(inout) :: aa(:,:)

    !> Descriptor of A.
    integer, intent(in) :: desca(DLEN_)

    !> Specifies whether lower ("L") or upper ("U") part of A contains the
    !! matrix. Default: "L".
    character, intent(in), optional :: uplo

    !> Specifies whether A is unit triangular ("U") or not ("N"). Default: "N".
    character, intent(in), optional :: diag

    !> First row of the submatrix of A. Default: 1
    integer, intent(in), optional :: ia

    !> First column of the submatrix of A. Default: 1
    integer, intent(in), optional :: ja

    !> Number of columns in the submatrix of A. Default: desca(M_).
    integer, intent(in), optional :: nn

    !> Info flag. If not specified and error occurs, the subroutine stops.
    integer, intent(out), optional :: info

    !------------------------------------------------------------------------

    character :: uplo0, diag0
    integer :: ia0, ja0, nn0, info0

    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(diag0, diag, "N")
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(nn0, nn, desca(M_))
    call ptrtri(uplo0, diag0, nn0, aa, ia0, ja0, desca, info0)
    call handle_infoflag(info0, "ptrtri in scalafx_ptrtri_${TYPE}$", info)

  end subroutine scalafx_ptrtri_${TYPE}$

#:enddef scalafx_ptrtri_template


!************************************************************************
!*** pgetrf
!************************************************************************

#:def scalafx_pgetrf_template(TYPE, FTYPE)

  !> LU factorization of a general matrix with pivoting
  !!
  subroutine scalafx_pgetrf_${TYPE}$(aa, desca, ipiv, ia, ja, mm, nn, info)

    !> LU decomposition on exit, pivoted by ipiv
    ${FTYPE}$, intent(inout) :: aa(:,:)

    !> Descriptor of A.
    integer, intent(in) :: desca(DLEN_)

    !> Pivot matrix
    integer, intent(out) :: ipiv(:)

    !> First row of the submatrix of A. Default: 1
    integer, intent(in), optional :: ia

    !> First column of the submatrix of A. Default: 1
    integer, intent(in), optional :: ja

    !> Number of columns in the submatrix of A. Default: desca(M_)
    integer, intent(in), optional :: mm

    !> Number of rows in the submatrix of A. Default: desca(N_)
    integer, intent(in), optional :: nn

    !> Info flag. If not specified and error occurs, the subroutine stops.
    integer, intent(out), optional :: info

    !------------------------------------------------------------------------

    integer :: ia0, ja0, mm0, nn0, info0

    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(mm0, mm, desca(M_))
    @:inoptflags(nn0, nn, desca(N_))
    call pgetrf(mm0, nn0, aa, ia0, ja0, desca, ipiv, info0)
    call handle_infoflag(info0, "pgetrf in scalafx_pgetrf_${TYPE}$", info)

  end subroutine scalafx_pgetrf_${TYPE}$

#:enddef scalafx_pgetrf_template


!************************************************************************
!*** psygst / phegst
!************************************************************************

#:def scalafx_psygst_phegst_template(TYPE, FTYPE)

  #:set PSYGST_PHEGST_HANDLE = {'real': 'psygst', 'dreal': 'psygst',&
      & 'complex': 'phegst', 'dcomplex': 'phegst'}
  #:set NAME = PSYGST_PHEGST_HANDLE[TYPE]
      
  !> Reduces Hermitian-definite generalized eigenvalue problem to standard form.
  !!
  !! \param ibtype  Type of the problem (1, 2, 3).
  !! \param aa  Matrix A.
  !! \param desca  Descriptor of matrix A.
  !! \param bb  Right hand side of the eigenvalue equation (B).
  !! \param desb  Descriptor of matrix B.
  !! \param scale  Scaling factors on return.
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param nn  Number of rows and columns of the submatrices A and B
  !!     (default: desca(M_)).
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param ib  First row of the submatrix B (default: 1)
  !! \param jb  First column of the submatrix B (default: 1)
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !!
  !! \see SCALAPACK documentation (routines p?sygst/p?hegst).
  !!
  subroutine scalafx_${NAME}$_${TYPE}$(ibtype, aa, desca, bb, descb, scale, uplo, &
      & nn, ia, ja, ib, jb, info)
    integer, intent(in) :: ibtype
    ${FTYPE}$, intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    ${FTYPE}$, intent(in) :: bb(:,:)
    integer, intent(in) :: descb(DLEN_)
    real(${KIND}$), intent(out) :: scale
    character, intent(in), optional :: uplo
    integer, intent(in), optional :: nn, ia, ja, ib, jb
    integer, intent(out), optional :: info

    integer :: nn0, ia0, ja0, ib0, jb0, info0
    character :: uplo0

    @:inoptflags(nn0, nn, desca(M_))
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(ib0, ib, 1)
    @:inoptflags(jb0, jb, 1)
    @:inoptflags(uplo0, uplo, "L")
    call ${NAME}$(ibtype, uplo0, nn0, aa, ia0, ja0, desca, bb,&
        & ib0, jb0, descb, scale, info0)
    call handle_infoflag(info0, "$4 in scalafx_$4_${TYPE}$", info)

  end subroutine scalafx_${NAME}$_${TYPE}$

#:enddef scalafx_psygst_phegst_template

!************************************************************************
!*** psyev
!************************************************************************

#:def scalafx_psyev_template(TYPE, FTYPE)

  !> Solves real eigenvalue problem.
  !!
  !! \param aa  Matrix to diagonalize (A).
  !! \param desca  Descriptor of matrix A.
  !! \param ww  Eigenvalues on exit.
  !! \param zz  Eigenvectors on exit (Z).
  !! \param descz  Descriptor of the eigenvector matrix.
  !! \param jobz  Job type (default: "V")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param iz  First row of the submatrix Z (default: 1)
  !! \param jz  First column of the submatrix Z (default: 1)
  !! \param work  Working array (if not specified, allocated automatically)
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !!
  !! \see SCALAPACK documentation (routine p?syev).
  !!
  subroutine scalafx_psyev_${TYPE}$(aa, desca, ww, zz, descz, jobz, uplo, ia, ja,&
      & iz, jz, work, info)
    real(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    real(${KIND}$), intent(out) :: ww(:)
    real(${KIND}$), intent(out) :: zz(:,:)
    integer, intent(in) :: descz(DLEN_)
    character, intent(in), optional :: jobz, uplo
    integer, intent(in), optional :: ia, ja, iz, jz
    real(${KIND}$), intent(inout), allocatable, optional :: work(:)
    integer, intent(out), optional :: info

    integer :: nn, lwork, info0, ia0, ja0, iz0, jz0
    character :: uplo0, jobz0
    real(${KIND}$) :: rtmp(1)
    real(${KIND}$), allocatable :: work0(:)

    ! Handle optional flags
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(jobz0, jobz, "V")
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(iz0, iz, 1)
    @:inoptflags(jz0, jz, 1)

    ! Allocate real workspace
    nn = desca(M_)
    call psyev(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
        & rtmp, -1, info0)
    call handle_infoflag(info0, "psev in scalafx_psyev", info)
    @:move_minoptalloc(work0, int(rtmp(1)), lwork, work)

    ! Diagonalization
    ! Initializing workspace as SCALAPACK apparently accesses uninitialized
    ! elements in it (nagfors -nan flag causes *sometimes* arithmetic exception)
    work0(:) = 0.0_${KIND}$
    call psyev(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
        & work0, lwork, info0)
    call handle_infoflag(info0, "psyev in scalafx_psyev_${TYPE}$", info)

    ! Save work space allocations, if dummy arguments present
    @:optmovealloc(work0, work)

  end subroutine scalafx_psyev_${TYPE}$

#:enddef scalafx_psyev_template


!************************************************************************
!*** pheev
!************************************************************************

#:def scalafx_pheev_template(TYPE, KIND)

  !> Solves complex eigenvalue problem.
  !!
  !! \param aa  Matrix to diagonalize (A).
  !! \param desca  Descriptor of matrix A.
  !! \param ww  Eigenvalues on exit.
  !! \param zz  Eigenvectors on exit (Z).
  !! \param descz  Descriptor of the eigenvector matrix.
  !! \param jobz  Job type (default: "V")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param iz  First row of the submatrix Z (default: 1)
  !! \param jz  First column of the submatrix Z (default: 1)
  !! \param work  Complex working array (if not specified, allocated automatically)
  !! \param rwork Real working array (if not specified, allocated automatically)
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !!
  !! \see SCALAPACK documentation (routine p?heev).
  !!
  subroutine scalafx_pheev_${TYPE}$(aa, desca, ww, zz, descz, jobz, uplo, ia, ja,&
      & iz, jz, work, rwork, info)
    complex(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    real(${KIND}$), intent(out) :: ww(:)
    complex(${KIND}$), intent(out) :: zz(:,:)
    integer, intent(in) :: descz(DLEN_)
    character, intent(in), optional :: jobz, uplo
    integer, intent(in), optional :: ia, ja, iz, jz
    complex(${KIND}$), intent(inout), allocatable, optional :: work(:)
    real(${KIND}$), intent(inout), allocatable, optional :: rwork(:)
    integer, intent(out), optional :: info

    integer :: nn, lwork, lrwork, lrwork_tmp, info0, ia0, ja0, iz0, jz0
    character :: uplo0, jobz0
    real(${KIND}$) :: rtmp(1)
    complex(${KIND}$) :: ctmp(1)
    complex(${KIND}$), allocatable :: work0(:)
    real(${KIND}$), allocatable :: rwork0(:)

    ! Handle optional flags
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(jobz0, jobz, "V")
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(iz0, iz, 1)
    @:inoptflags(jz0, jz, 1)

    ! Allocate  workspace
    nn = desca(M_)
    call pheev(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
        & ctmp, -1, rtmp, -1, info0)
    call handle_infoflag(info0, "pheev in scalafx_pheev", info)
    @:move_minoptalloc(work0, int(ctmp(1)), lwork, work)

    ! Apparently pheev sometimes returns considerably smaller work space size as
    ! required in the documentation. We take the bigger value to be safe.
    if (jobz0 == "N") then
      lrwork_tmp = max(2 * nn, int(rtmp(1)))
    else
      lrwork_tmp = max(2 * nn + 2 * nn - 2, int(rtmp(1)))
    end if
    @:move_minoptalloc(rwork0, lrwork_tmp, lrwork, rwork)

    ! Diagonalization
    ! Initializing workspace as SCALAPACK apparently accesses uninitialized
    ! elements in it (nagfors -nan flag causes *sometimes* arithmetic exception)
    work0(:) = cmplx(0, kind=${KIND}$)
    rwork0(:) = 0.0_${KIND}$
    call pheev(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
        & work0, lwork, rwork0, lrwork, info0)
    call handle_infoflag(info0, "pheev in scalafx_pheev_${TYPE}$", info)

    ! Save work space allocations, if dummy arguments present
    @:optmovealloc(rwork0, rwork)
    @:optmovealloc(work0, work)


  end subroutine scalafx_pheev_${TYPE}$

#:enddef scalafx_pheev_template



#!************************************************************************
#!*** psyevd
#!************************************************************************

#:def scalafx_psyevd_template(TYPE, KIND)

  !> Solves symmetric eigenvalue problem by the divide and conquer algorithm.
  !!
  !! \param aa  Matrix to diagonalize (A).
  !! \param desca  Descriptor of matrix A.
  !! \param ww  Eigenvalues on exit.
  !! \param zz  Eigenvectors on exit (Z).
  !! \param descz  Descriptor of the eigenvector matrix.
  !! \param jobz  Job type (default: "V")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param iz  First row of the submatrix Z (default: 1)
  !! \param jz  First column of the submatrix Z (default: 1)
  !! \param work  Working array (if not specified, allocated automatically)
  !! \param iwork Integer working array (if not specified, allocated
  !!     automatically)
  !! \param allocfix  If yes, the routine tries to enlarge the workspace size
  !!     as returned by the appropriate p?syevd() routine by some empirical
  !!     values.
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !!
  !! \note
  !!     Unfortunately, SCALAPACK seems to return the wrong real work space
  !!     sizes for many cases. This routine (when allocfix had been set to true)
  !!     tries to improve on it by two ways:
  !! \note
  !!     * It queries also the QR routine (p?syev) for workspace size and takes
  !!       this, if bigger than returned by p?syevd. That should ensure
  !!       that the pdormtr() routine does not encounter any difficulties.
  !! \note
  !!     * It additionally enlarges the real workspace size by the amount of
  !!       memory needed by the pdlasrt() routine, to make sure this would not
  !!       fail either.
  !! \note
  !!     Those fixes are empirical, may lead to oversized workspace allocations
  !!     and probably would not even fix the allocation problem, but are the best
  !!     we could find so far.
  !!
  !! \see SCALAPACK documentation (routine p?syevd).
  !!
  subroutine scalafx_psyevd_${TYPE}$(aa, desca, ww, zz, descz, jobz, uplo, ia,&
      & ja, iz, jz, work, iwork, allocfix, info)
    real(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    real(${KIND}$), intent(out) :: ww(:)
    real(${KIND}$), intent(out) :: zz(:,:)
    integer, intent(in) :: descz(DLEN_)
    character, intent(in), optional :: jobz, uplo
    integer, intent(in), optional :: ia, ja, iz, jz
    real(${KIND}$), intent(inout), allocatable, optional :: work(:)
    integer, intent(inout), allocatable, optional :: iwork(:)
    logical, intent(in), optional :: allocfix
    integer, intent(out), optional :: info

    integer :: nn, liwork, lwork, lwmin, liwmin, info0, ia0, ja0, iz0, jz0
    character :: uplo0, jobz0
    integer :: itmp(1)
    real(${KIND}$) :: rtmp(1), rtmp2(1)
    real(${KIND}$), allocatable :: work0(:)
    integer, allocatable :: iwork0(:)
    logical :: allocfix0

    ! Handle optional flags
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(jobz0, jobz, "V")
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(iz0, iz, 1)
    @:inoptflags(jz0, jz, 1)
    @:inoptflags(allocfix0, allocfix, .false.)

    ! Allocate workspaces
    nn = desca(M_)
    call psyevd(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
        & rtmp, -1, itmp, 1, info0)
    call handle_infoflag(info0, "psevd in scalafx_psyevd", info)

    ! Check workspace size of psyev() and take that one if bigger in the hope
    ! pdormtr() would work. Additionally extend workspace in the hope pdlasrt()
    ! work as well.
    if (allocfix0) then
      call psyev(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz, &
        & rtmp2, -1, info0)
      call handle_infoflag(info0, "psev in scalafx_psyevd", info)
      lwmin = max(int(rtmp(1)), int(rtmp2(1)))
      lwmin = lwmin + MAX(nn, size(aa, dim=1) * (desca(NB_) + size(aa, dim=2)))
      liwmin = itmp(1)
    else
      lwmin = int(rtmp(1))
      liwmin = itmp(1)
    end if

    @:move_minoptalloc(work0, lwmin, lwork, work)
    @:move_minoptalloc(iwork0, liwmin, liwork, iwork)

    ! Diagonalization
    ! Initializing workspace as SCALAPACK apparently accesses uninitialized
    ! elements in it (nagfors -nan flag causes *sometimes* arithmetic exception)
    work0(:) = 0.0_${KIND}$
    iwork0(:) = 0
    call psyevd(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
        & work0, lwork, iwork0, liwork, info0)
    call handle_infoflag(info0, "psyevd in scalafx_psyevd_${TYPE}$", info)

    ! Save work space allocations, if dummy arguments present
    @:optmovealloc(iwork0, iwork)
    @:optmovealloc(work0, work)

  end subroutine scalafx_psyevd_${TYPE}$

#:enddef scalafx_psyevd_template



#! ************************************************************************
#! *** pheevd
#! ************************************************************************

#:def scalafx_pheevd_template(TYPE, KIND)

  !> Solves Hermitian eigenvalue problem by the divide and conquer algorithm.
  !!
  !! \param aa  Matrix to diagonalize (A).
  !! \param desca  Descriptor of matrix A.
  !! \param ww  Eigenvalues on exit.
  !! \param zz  Eigenvectors on exit (Z).
  !! \param descz  Descriptor of the eigenvector matrix.
  !! \param jobz  Job type (default: "V")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param iz  First row of the submatrix Z (default: 1)
  !! \param jz  First column of the submatrix Z (default: 1)
  !! \param work  Complex working array (if not specified, allocated
  !!     automatically)
  !! \param rwork Real working array (if not specified, allocated automatically)
  !! \param iwork Integer working array (if not specified, allocated
  !!     automatically)
  !! \param allocfix  If yes, the routine tries to enlarge the workspace size
  !!     as returned by the appropriate p?syevd() routine by some empirical
  !!     values.
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !!
  !! \note
  !!     Unfortunately, SCALAPACK seems to return the wrong real work space
  !!     sizes for many cases. This routine (when allocfix had been set to true)
  !!     tries to improve on it by two ways:
  !! \note
  !!     * It queries also the QR routine (p?heev) for workspace size and takes
  !!       this, if bigger than returned by p?heevd. That should ensure
  !!       that the pzunmtr() routine does not encounter any difficulties.
  !! \note
  !!     * It additionally enlarges the real workspace size by the amount of
  !!       memory needed by the pdlasrt() routine, to make sure this would not
  !!       fail either.
  !! \note
  !!     Those fixes are empirical, may lead to oversized workspace allocations
  !!     and probably would not even fix the allocation problem, but are the best
  !!     we could find so far.
  !!
  !! \see SCALAPACK documentation (routine p?heevd).
  !!
  subroutine scalafx_pheevd_${TYPE}$(aa, desca, ww, zz, descz, jobz, uplo, ia, ja, iz,&
      & jz, work, rwork, iwork, allocfix, info)
    complex(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    real(${KIND}$), intent(out) :: ww(:)
    complex(${KIND}$), intent(out) :: zz(:,:)
    integer, intent(in) :: descz(DLEN_)
    character, intent(in), optional :: jobz, uplo
    integer, intent(in), optional :: ia, ja, iz, jz
    complex(${KIND}$), intent(inout), allocatable, optional :: work(:)
    real(${KIND}$), intent(inout), allocatable, optional :: rwork(:)
    integer, intent(inout), allocatable, optional :: iwork(:)
    logical, intent(in), optional :: allocfix
    integer, intent(out), optional :: info

    integer :: nn, liwork, lwork, lrwork, info0, ia0, ja0, iz0, jz0
    integer :: lwmin, lrwmin, liwmin
    character :: uplo0, jobz0
    real(${KIND}$) :: rtmp(1), rtmp2(1)
    complex(${KIND}$) :: ctmp(1), ctmp2(1)
    integer :: itmp(1)
    complex(${KIND}$), allocatable :: work0(:)
    real(${KIND}$), allocatable :: rwork0(:)
    integer, allocatable :: iwork0(:)
    logical :: allocfix0

    ! Handle optional flags
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(jobz0, jobz, "V")
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(iz0, iz, 1)
    @:inoptflags(jz0, jz, 1)
    @:inoptflags(allocfix0, allocfix, .false.)

    ! Allocate workspaces
    nn = desca(M_)
    call pheevd(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
        & ctmp, -1, rtmp, -1, itmp, 1, info0)
    call handle_infoflag(info0, "pheevd in scalafx_pheevd", info)

    ! Check workspace size of pheev() and take that one if bigger in the hope
    ! pdormtr() would work. Additionally extend workspace in the hope pdlasrt()
    ! work as well.
    if (allocfix0) then
      call pheev(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
        & ctmp2, -1, rtmp2, -1, info0)
      call handle_infoflag(info0, "pheev in scalafx_pheevd", info)
      if (jobz0 == "N") then
        lrwmin = max(2 * nn, int(rtmp2(1)))
      else
        lrwmin = max(2 * nn + 2 * nn - 2, int(rtmp2(1)))
      end if
      lwmin = max(int(ctmp(1)), int(ctmp2(1)))
      lrwmin = max(int(rtmp(1)), lrwmin)
      lrwmin = lrwmin + MAX(nn, size(aa, dim=1) * (desca(NB_) + size(aa, dim=2)))
      liwmin = itmp(1)
    else
      lwmin = int(ctmp(1))
      lrwmin = int(rtmp(1))
      liwmin = itmp(1)
    end if

    @:move_minoptalloc(work0, lwmin, lwork, work)
    @:move_minoptalloc(rwork0, lrwmin, lrwork, rwork)
    @:move_minoptalloc(iwork0, liwmin, liwork, iwork)

    ! Diagonalization
    ! Initializing workspace as SCALAPACK apparently accesses uninitialized
    ! elements in it (nagfors -nan flag causes *sometimes* arithmetic exception)
    work0(:) = cmplx(0, kind=${KIND}$)
    rwork0(:) = 0.0_${KIND}$
    iwork0(:) = 0
    call pheevd(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
        & work0, lwork, rwork0, lrwork, iwork0, liwork, info0)
    call handle_infoflag(info0, "pheevd in scalafx_pheevd_${TYPE}$", info)

    ! Save work space allocations, if dummy arguments present
    @:optmovealloc(work0, work)
    @:optmovealloc(rwork0, rwork)
    @:optmovealloc(iwork0, iwork)

  end subroutine scalafx_pheevd_${TYPE}$

#:enddef scalafx_pheevd_template


#! ************************************************************************
#! *** psyevr
#! ************************************************************************

#:def scalafx_psyevr_template(TYPE, KIND)

  !> Solves symmetric eigenvalue problem by the MRRR algorithm.
  !!
  !! \param aa  Matrix to diagonalize (A).
  !! \param desca  Descriptor of matrix A.
  !! \param ww  Eigenvalues on exit.
  !! \param zz  Eigenvectors on exit (Z).
  !! \param descz  Descriptor of the eigenvector matrix.
  !! \param jobz  Job type (default: "V")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param nn  Number of rows and colums of the submatrix A (default: desca(M_))
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param iz  First row of the submatrix Z (default: 1)
  !! \param jz  First column of the submatrix Z (default: 1)
  !! \param work  Working array (if not specified, allocated automatically)
  !! \param iwork Integer working array (if not specified, allocated
  !!     automatically)
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !!
  !! \see SCALAPACK documentation (routine p?syevr).
  !!
  subroutine scalafx_psyevr_${TYPE}$(aa, desca, ww, zz, descz, vl, vu, il, iu, &
      & jobz, uplo, nn, ia, ja, iz, jz, work, iwork, mm, nz, info)
    real(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    real(${KIND}$), intent(out) :: ww(:)
    real(${KIND}$), intent(out) :: zz(:,:)
    integer, intent(in) :: descz(DLEN_)
    real(${KIND}$), intent(in), optional :: vl, vu
    integer, intent(in), optional :: il, iu
    character, intent(in), optional :: jobz, uplo
    integer, intent(in), optional :: nn, ia, ja, iz, jz
    real(${KIND}$), intent(inout), allocatable, optional :: work(:)
    integer, intent(inout), allocatable, optional :: iwork(:)
    integer, intent(out), optional :: mm, nz
    integer, intent(out), optional :: info

    integer :: liwork, lwork, lwmin, liwmin, info0, nn0, ia0, ja0, iz0, jz0
    integer :: mm0, nz0, il0, iu0
    real(${KIND}$) :: vl0, vu0
    character :: uplo0, jobz0, range
    integer :: itmp(1)
    real(${KIND}$) :: rtmp(1)
    real(${KIND}$), allocatable :: work0(:)
    integer, allocatable :: iwork0(:)

    ! Handle optional flags
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(jobz0, jobz, "V")
    @:inoptflags(nn0, nn, desca(M_))
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(iz0, iz, 1)
    @:inoptflags(jz0, jz, 1)

    if (present(il) .or. present(iu)) then
      @:inoptflags(il0, il, 1)
      @:inoptflags(iu0, iu, nn0)
      range = "I"
    elseif (present(vl) .or. present(vu)) then
      if (range == "I") then
        call error("Eigenvalue subset must be specified by index or by&
            & energy range, but not by both.")
      end if
      if (.not. (present(vl) .and. present(vu))) then
        call error("When eigenvalue subset is defined by energy range, both, &
            & upper and lower limit must be present.")
      end if
      vl0 = vl
      vu0 = vu
      range ="V"
    else
      range = "A"
    end if

    if (range /= "V" .and. present(nz)) then
      call error("Optional argument nz only valid if subset is selected by &
          &energy range")
    end if

    ! Allocate workspaces
    call psyevr(jobz0, range, uplo0, nn0, aa, ia0, ja0, desca, vl0, vu0, &
        & il0, iu0, mm0, nz0, ww, zz, iz0, jz0, descz, rtmp, -1, itmp, -1, info0)
    call handle_infoflag(info0, "psyevr (workspace query) in scalafx_psyevr_${TYPE}$", &
        & info)

    ! Check workspace size of psyevr() and take that.
    lwmin = int(rtmp(1))
    liwmin = itmp(1)

    @:move_minoptalloc(work0, lwmin, lwork, work)
    @:move_minoptalloc(iwork0, liwmin, liwork, iwork)

    ! Diagonalization
    ! Initializing workspace as SCALAPACK apparently accesses uninitialized
    ! elements in it (nagfors -nan flag causes *sometimes* arithmetic exception)
    work0(:) = 0.0_${KIND}$
    iwork0(:) = 0
    call psyevr(jobz0, range, uplo0, nn0, aa, ia0, ja0, desca, vl0, vu0, &
        & il0, iu0, mm0, nz0, ww, zz, iz0, jz0, descz, work0, lwork, iwork0, &
        & liwork, info0)
    call handle_infoflag(info0, "psyevr (diagonalization) in scalafx_psyevr_${TYPE}$", &
        & info)

    ! Save work space allocations, if dummy arguments present
    @:optmovealloc(iwork0, iwork)
    @:optmovealloc(work0, work)

    @:handle_outoptflag(mm, mm0)
    if (range == "V") then
      @:handle_outoptflag(nz, nz0)
    end if

  end subroutine scalafx_psyevr_${TYPE}$

#:enddef scalafx_psyevr_template



#! ************************************************************************
#! *** pheevr
#! ************************************************************************

#:def scalafx_pheevr_template(TYPE, KIND)

  !> Solves Hermitian eigenvalue problem by the divide and conquer algorithm.
  !!
  !! \param aa  Matrix to diagonalize (A).
  !! \param desca  Descriptor of matrix A.
  !! \param ww  Eigenvalues on exit.
  !! \param zz  Eigenvectors on exit (Z).
  !! \param descz  Descriptor of the eigenvector matrix.
  !! \param jobz  Job type (default: "V")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param iz  First row of the submatrix Z (default: 1)
  !! \param jz  First column of the submatrix Z (default: 1)
  !! \param work  Complex working array (if not specified, allocated
  !!     automatically)
  !! \param rwork Real working array (if not specified, allocated automatically)
  !! \param iwork Integer working array (if not specified, allocated
  !!     automatically)
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !!
  !! \see SCALAPACK documentation (routine p?heevr).
  !!
  subroutine scalafx_pheevr_${TYPE}$(aa, desca, ww, zz, descz, range, vl, vu,&
      & il, iu, m, nz, jobz, uplo, ia, ja, iz, jz, work, rwork, iwork, info)
    complex(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    real(${KIND}$), intent(out) :: ww(:)
    complex(${KIND}$), intent(out) :: zz(:,:)
    integer, intent(in) :: descz(DLEN_)
    character, intent(in) :: range
    real(${KIND}$), intent(in) :: vl, vu
    integer, intent(in) :: il, iu
    integer, intent(out) :: m, nz
    character, intent(in), optional :: jobz, uplo
    integer, intent(in), optional :: ia, ja, iz, jz
    complex(${KIND}$), intent(inout), allocatable, optional :: work(:)
    real(${KIND}$), intent(inout), allocatable, optional :: rwork(:)
    integer, intent(inout), allocatable, optional :: iwork(:)
    integer, intent(out), optional :: info

    integer :: nn, liwork, lwork, lrwork, info0, ia0, ja0, iz0, jz0
    integer :: lwmin, lrwmin, liwmin
    character :: uplo0, jobz0
    real(${KIND}$) :: rtmp(1)
    complex(${KIND}$) :: ctmp(1)
    integer :: itmp(1)
    complex(${KIND}$), allocatable :: work0(:)
    real(${KIND}$), allocatable :: rwork0(:)
    integer, allocatable :: iwork0(:)

    ! Handle optional flags
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(jobz0, jobz, "V")
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(iz0, iz, 1)
    @:inoptflags(jz0, jz, 1)

    ! Allocate workspaces
    nn = desca(M_)
    call pheevr(jobz0, range, uplo0, nn, aa, ia0, ja0, desca, vl, vu, il, iu, &
        & m, nz, ww, zz, iz0, jz0, descz, ctmp, -1, rtmp, -1, itmp, -1, info0)
    call handle_infoflag(info0, "pheevr in scalafx_pheevr", info)

    lwmin = int(ctmp(1))
    lrwmin = int(rtmp(1))
    liwmin = itmp(1)

    @:move_minoptalloc(work0, lwmin, lwork, work)
    @:move_minoptalloc(rwork0, lrwmin, lrwork, rwork)
    @:move_minoptalloc(iwork0, liwmin, liwork, iwork)

    ! Diagonalization
    ! Initializing workspace as SCALAPACK apparently accesses uninitialized
    ! elements in it (nagfors -nan flag causes *sometimes* arithmetic exception)
    work0(:) = cmplx(0, kind=${KIND}$)
    rwork0(:) = 0.0_${KIND}$
    iwork0(:) = 0
    call pheevr(jobz0, range, uplo0, nn, aa, ia0, ja0, desca, vl, vu, il, iu, &
        & m, nz, ww, zz, iz0, jz0, descz, &
        & work0, lwork, rwork0, lrwork, iwork0, liwork, info0)
    call handle_infoflag(info0, "pheevr in scalafx_pheevr_${TYPE}$", info)

    ! Save work space allocations, if dummy arguments present
    @:optmovealloc(work0, work)
    @:optmovealloc(rwork0, rwork)
    @:optmovealloc(iwork0, iwork)

  end subroutine scalafx_pheevr_${TYPE}$

#:enddef scalafx_pheevr_template



#! ************************************************************************
#! *** psygv
#! ************************************************************************

#:def scalafx_psygv_template(TYPE, KIND)
  !> Solves symmetric generalized eigenvalue problem by the QR algorithm.
  !!
  !! \details Invokes SCALAPACK routines p?potrf, p?sygst, p?syev, p?trsm in
  !! order to transform the general eigenvalue problem to the standard form
  !! and transform the eigenvectors back.
  !!
  !! \param aa  Matrix to diagonalize (A), transformed matrix on exit.
  !! \param desca  Descriptor of matrix A.
  !! \param bb  Matrix on the right hand side (B), transformed matrix on exit.
  !! \param descb  Descriptor of matrix B.
  !! \param ww  Eigenvalues on exit.
  !! \param zz  Eigenvectors on exit (Z).
  !! \param descz  Descriptor of the eigenvector matrix.
  !! \param jobz  Job type (default: "V")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param ib  First row of the submatrix B (default: 1)
  !! \param jb  First column of the submatrix B (default: 1)
  !! \param iz  First row of the submatrix Z (default: 1)
  !! \param jz  First column of the submatrix Z (default: 1)
  !! \param work  Working array (if not specified, allocated automatically)
  !! \param skipchol  If true, the Cholesky transformation will be skipped.
  !!     Array bb must have the Cholesky transformed form.
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !! \see SCALAPACK documentation (routines p?potrf, p?hegst, p?heev, p?trsm).
  subroutine scalafx_psygv_${TYPE}$(aa, desca, bb, descb, ww, zz, descz, jobz, uplo,&
      & ia, ja, ib, jb, iz, jz, work, skipchol, info)
    real(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    real(${KIND}$), intent(inout) :: bb(:,:)
    integer, intent(in) :: descb(DLEN_)
    real(${KIND}$), intent(out) :: ww(:)
    real(${KIND}$), intent(out) :: zz(:,:)
    integer, intent(in) :: descz(DLEN_)
    character, intent(in), optional :: jobz, uplo
    integer, intent(in), optional :: ia, ja, ib, jb, iz, jz
    real(${KIND}$), intent(inout), allocatable, optional :: work(:)
    logical, intent(in), optional :: skipchol
    integer, intent(out), optional :: info

    real(${KIND}$) :: scale
    character :: jobz0, transa, uplo0
    logical :: skipchol0

    @:inoptflags(jobz0, jobz, "V")
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(skipchol0, skipchol, .false.)

    if (.not. skipchol0) then
      call scalafx_ppotrf(bb, descb, uplo=uplo0, ia=ib, ja=jb, info=info)
    end if
    call scalafx_psygst(1, aa, desca, bb, descb, scale, uplo=uplo0, ia=ia, ja=ja,&
        & ib=ib, jb=jb, info=info)
    call scalafx_psyev(aa, desca, ww, zz, descz, jobz=jobz0, uplo=uplo0, ia=ia,&
        & ja=ja, iz=iz, jz=jz, work=work, info=info)
    if (jobz0 == "V") then
      if (uplo0 == "L" .or. uplo0 == "l") then
         transa = "T"
      else
         transa = "N"
      end if
      call scalafx_ptrsm(bb, descb, zz, descz, side="L", uplo=uplo0,&
          & transa=transa, diag="N", ia=ib, ja=jb, ib=iz, jb=jz)
    end if

  end subroutine scalafx_psygv_${TYPE}$

#:enddef scalafx_psygv_template


#! ************************************************************************
#! *** phegv
#! ************************************************************************

#:def scalafx_phegv_template(TYPE, KIND)

  !> Solves Hermitian generalized eigenvalue problem by the QR algorithm.
  !!
  !! \details Invokes SCALAPACK routines p?potrf, p?hegst, p?heev, p?trsm in
  !! order to transform the general eigenvalue problem to the standard form
  !! and transform the eigenvectors back.
  !!
  !! \param aa  Matrix to diagonalize (A), transformed matrix on exit.
  !! \param desca  Descriptor of matrix A.
  !! \param bb  Matrix on the right hand side (B), transformed matrix on exit.
  !! \param descb  Descriptor of matrix B.
  !! \param ww  Eigenvalues on exit.
  !! \param zz  Eigenvectors on exit (Z).
  !! \param descz  Descriptor of the eigenvector matrix.
  !! \param jobz  Job type (default: "V")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param ib  First row of the submatrix B (default: 1)
  !! \param jb  First column of the submatrix B (default: 1)
  !! \param iz  First row of the submatrix Z (default: 1)
  !! \param jz  First column of the submatrix Z (default: 1)
  !! \param work  Complex working array (if not specified, allocated
  !!     (automatically)
  !! \param rwork Real working array (if not specified, allocated automatically)
  !! \param skipchol  If true, the Cholesky transformation will be skipped.
  !!     Array bb must have the Cholesky transformed form.
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !! \see SCALAPACK documentation (routines p?potrf, p?hegst, p?heev, p?trsm).
  subroutine scalafx_phegv_${TYPE}$(aa, desca, bb, descb, ww, zz, descz, jobz,&
      & uplo, ia, ja, ib, jb, iz, jz, work, rwork, skipchol, info)
    complex(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    complex(${KIND}$), intent(inout) :: bb(:,:)
    integer, intent(in) :: descb(DLEN_)
    real(${KIND}$), intent(out) :: ww(:)
    complex(${KIND}$), intent(out) :: zz(:,:)
    integer, intent(in) :: descz(DLEN_)
    character, intent(in), optional :: jobz, uplo
    integer, intent(in), optional :: ia, ja, ib, jb, iz, jz
    complex(${KIND}$), intent(inout), allocatable, optional :: work(:)
    real(${KIND}$), intent(inout), allocatable, optional :: rwork(:)
    logical, intent(in), optional :: skipchol
    integer, intent(out), optional :: info

    real(${KIND}$) :: scale
    character :: jobz0, transa, uplo0
    logical :: skipchol0

    @:inoptflags(jobz0, jobz, "V")
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(skipchol0, skipchol, .false.)

    if (.not. skipchol0) then
      call scalafx_ppotrf(bb, descb, uplo=uplo0, ia=ib, ja=jb, info=info)
    end if
    call scalafx_phegst(1, aa, desca, bb, descb, scale, uplo=uplo0, ia=ia, ja=ja,&
        & ib=ib, jb=jb, info=info)
    call scalafx_pheev(aa, desca, ww, zz, descz, jobz=jobz0, uplo=uplo0, ia=ia,&
        & ja=ja, iz=iz, jz=jz, work=work, rwork=rwork, info=info)
    if (jobz0 == "V") then
      if (uplo0 == "L" .or. uplo0 == "l") then
         transa = "C"
      else
         transa = "N"
      end if
      call scalafx_ptrsm(bb, descb, zz, descz, side="L", uplo=uplo0,&
          & transa=transa, diag="N", ia=ib, ja=jb, ib=iz, jb=jz)
    end if

  end subroutine scalafx_phegv_${TYPE}$

#:enddef scalafx_phegv_template



#! ************************************************************************
#! *** psygvd
#! ************************************************************************

#:def scalafx_psygvd_template(TYPE, KIND)

  !> Solves real generalized eigenvalue problem by the divide and conquer
  !! algorithm.
  !!
  !! \details Invokes SCALAPACK routines p?potrf, p?sygst, p?syevd, p?trsm in
  !! order to transform the general eigenvalue problem to the standard form
  !! and transform the eigenvectors back.
  !!
  !! \param aa  Matrix to diagonalize (A), transformed matrix on exit.
  !! \param desca  Descriptor of matrix A.
  !! \param bb  Matrix on the right hand side (B), transformed matrix on exit.
  !! \param descb  Descriptor of matrix B.
  !! \param ww  Eigenvalues on exit.
  !! \param zz  Eigenvectors on exit (Z).
  !! \param descz  Descriptor of the eigenvector matrix.
  !! \param jobz  Job type (default: "V")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param ib  First row of the submatrix B (default: 1)
  !! \param jb  First column of the submatrix B (default: 1)
  !! \param iz  First row of the submatrix Z (default: 1)
  !! \param jz  First column of the submatrix Z (default: 1)
  !! \param work  Working array (if not specified, allocated automatically)
  !! \param iwork Integer working array (if not specified, allocated
  !!     automatically)
  !! \param allocfix  If yes, the routine tries to enlarge the workspace size
  !!     as returned by the appropriate p?syevd() routine by some empirical
  !!     values. See the scalafx_psyevd_${TYPE}$() routine for details.
  !! \param skipchol  If true, the Cholesky transformation will be skipped.
  !!     Array bb must have the Cholesky transformed form.
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !! \see SCALAPACK documentation (routines p?potrf, p?sygst, p?syevd, p?trsm).
  subroutine scalafx_psygvd_${TYPE}$(aa, desca, bb, descb, ww, zz, descz, jobz, uplo,&
      & ia, ja, ib, jb, iz, jz, work, iwork, allocfix, skipchol, info)
    real(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    real(${KIND}$), intent(inout) :: bb(:,:)
    integer, intent(in) :: descb(DLEN_)
    real(${KIND}$), intent(out) :: ww(:)
    real(${KIND}$), intent(out) :: zz(:,:)
    integer, intent(in) :: descz(DLEN_)
    character, intent(in), optional :: jobz, uplo
    integer, intent(in), optional :: ia, ja, ib, jb, iz, jz
    real(${KIND}$), intent(inout), allocatable, optional :: work(:)
    integer, intent(inout), allocatable, optional :: iwork(:)
    logical, intent(in), optional :: allocfix, skipchol
    integer, intent(out), optional :: info

    real(${KIND}$) :: scale
    character :: jobz0, transa, uplo0
    logical :: skipchol0

    @:inoptflags(jobz0, jobz, "V")
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(skipchol0, skipchol, .false.)

    ! Cholesky transformation of B
    if (.not. skipchol0) then
      call scalafx_ppotrf(bb, descb, uplo=uplo0, ia=ib, ja=jb, info=info)
    end if
    ! Reducing to standard form
    call scalafx_psygst(1, aa, desca, bb, descb, scale, uplo=uplo0, ia=ia, ja=ja,&
        & ib=ib, jb=jb, info=info)
    ! Solving eigenvalue problem.
    call scalafx_psyevd(aa, desca, ww, zz, descz, jobz=jobz0, uplo=uplo0, ia=ia,&
        & ja=ja, iz=iz, jz=jz, work=work, iwork=iwork, allocfix=allocfix, &
        & info=info)
    ! Transforming eigenvectors back
    if (jobz0 == "V") then
      if (uplo0 == "L" .or. uplo0 == "l") then
         transa = "T"
      else
         transa = "N"
      end if
      call scalafx_ptrsm(bb, descb, zz, descz, side="L", uplo=uplo0,&
          & transa=transa, diag="N", ia=ib, ja=jb, ib=iz, jb=jz)
    end if

  end subroutine scalafx_psygvd_${TYPE}$

#:enddef scalafx_psygvd_template


#! ************************************************************************
#! *** phegvd
#! ************************************************************************

#:def scalafx_phegvd_template(TYPE, KIND)

  !> Solves Hermitian generalized eigenvalue problem by the divide and conquer
  !! algorithm.
  !!
  !! \details Invokes SCALAPACK routines p?potrf, p?hegst, p?heevd, p?trsm in
  !! order to transform the general eigenvalue problem to the standard form
  !! and transform the eigenvectors back.
  !!
  !! \param aa  Matrix to diagonalize (A), transformed matrix on exit.
  !! \param desca  Descriptor of matrix A.
  !! \param bb  Matrix on the right hand side (B), transformed matrix on exit.
  !! \param descb  Descriptor of matrix B.
  !! \param ww  Eigenvalues on exit.
  !! \param zz  Eigenvectors on exit (Z).
  !! \param descz  Descriptor of the eigenvector matrix.
  !! \param jobz  Job type (default: "V")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param ib  First row of the submatrix B (default: 1)
  !! \param jb  First column of the submatrix B (default: 1)
  !! \param iz  First row of the submatrix Z (default: 1)
  !! \param jz  First column of the submatrix Z (default: 1)
  !! \param work  Complex working array (if not specified, allocated
  !!     automatically)
  !! \param rwork Real working array (if not specified, allocated automatically)
  !! \param iwork Integer working array (if not specified, allocated
  !!     automatically)
  !! \param allocfix  If yes, the routine tries to enlarge the workspace size
  !!     as returned by the appropriate p?syevd() routine by some empirical
  !!     values. See the scalafx_pheevd_${TYPE}$() routine for details.
  !! \param skipchol  If true, the Cholesky transformation will be skipped.
  !!     Array bb must have the Cholesky transformed form.
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !! \see SCALAPACK documentation (routines p?potrf, p?hegst, p?heevd, p?trsm).
  subroutine scalafx_phegvd_${TYPE}$(aa, desca, bb, descb, ww, zz, descz, jobz,&
      & uplo, ia, ja, ib, jb, iz, jz, work, rwork, iwork, allocfix, skipchol, info)
    complex(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    complex(${KIND}$), intent(inout) :: bb(:,:)
    integer, intent(in) :: descb(DLEN_)
    real(${KIND}$), intent(out) :: ww(:)
    complex(${KIND}$), intent(out) :: zz(:,:)
    integer, intent(in) :: descz(DLEN_)
    character, intent(in), optional :: jobz, uplo
    integer, intent(in), optional :: ia, ja, ib, jb, iz, jz
    complex(${KIND}$), intent(inout), allocatable, optional :: work(:)
    real(${KIND}$), intent(inout), allocatable, optional :: rwork(:)
    integer, intent(inout), allocatable, optional :: iwork(:)
    logical, intent(in), optional :: allocfix, skipchol
    integer, intent(out), optional :: info

    real(${KIND}$) :: scale
    character :: jobz0, transa, uplo0
    logical :: skipchol0

    @:inoptflags(jobz0, jobz, "V")
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(skipchol0, skipchol, .false.)

    ! Cholesky transformation of B.
    if (.not. skipchol0) then
      call scalafx_ppotrf(bb, descb, uplo=uplo0, ia=ib, ja=jb, info=info)
    end if
    ! Reducing to standard form
    call scalafx_phegst(1, aa, desca, bb, descb, scale, uplo=uplo0, ia=ia, ja=ja,&
        & ib=ib, jb=jb, info=info)
    ! Solving eigenvalue problem.
    call scalafx_pheevd(aa, desca, ww, zz, descz, jobz=jobz0, uplo=uplo0, ia=ia,&
        & ja=ja, iz=iz, jz=jz, work=work, rwork=rwork, iwork=iwork, &
        & allocfix=allocfix, info=info)
    ! Transforming eigenvectors back
    if (jobz0 == "V") then
      if (uplo0 == "L" .or. uplo0 == "l") then
         transa = "C"
      else
         transa = "N"
      end if
      call scalafx_ptrsm(bb, descb, zz, descz, side="L", uplo=uplo0,&
          & transa=transa, diag="N", ia=ib, ja=jb, ib=iz, jb=jz)
    end if

  end subroutine scalafx_phegvd_${TYPE}$

#:enddef scalafx_phegvd_template


#! ************************************************************************
#! *** psygvr
#! ************************************************************************

#:def scalafx_psygvr_template(TYPE, KIND)

  !> Solves real generalized eigenvalue problem by the MRRR algorithm.
  !!
  !! \details Invokes SCALAPACK routines p?potrf, p?sygst, p?syevr, p?trsm in
  !! order to transform the general eigenvalue problem to the standard form
  !! and transform the eigenvectors back. Currently all eigenvalues calculated.
  !! \param aa  Matrix to diagonalize (A), transformed matrix on exit.
  !! \param desca  Descriptor of matrix A.
  !! \param bb  Matrix on the right hand side (B), transformed matrix on exit.
  !! \param descb  Descriptor of matrix B.
  !! \param ww  Eigenvalues on exit.
  !! \param zz  Eigenvectors on exit (Z).
  !! \param descz  Descriptor of the eigenvector matrix.
  !! \param jobz  Job type (default: "V")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param nn  Number of rows and columns of the matrix A (default: desca(M_))
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param ib  First row of the submatrix B (default: 1)
  !! \param jb  First column of the submatrix B (default: 1)
  !! \param iz  First row of the submatrix Z (default: 1)
  !! \param jz  First column of the submatrix Z (default: 1)
  !! \param work  Working array (if not specified, allocated automatically)
  !! \param iwork Integer working array (if not specified, allocated
  !!     automatically)
  !! \param skipchol  If true, the Cholesky transformation will be skipped.
  !!     Array bb must have the Cholesky transformed form.
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !!
  !! \see SCALAPACK documentation (routines p?potrf, p?sygst, p?syevr, p?trsm).
  !!
  subroutine scalafx_psygvr_${TYPE}$(aa, desca, bb, descb, ww, zz, descz, &
      & vl, vu, il, iu, jobz, uplo, nn, ia, ja, ib, jb, iz, jz, work, iwork, &
      & skipchol, mm, nz, info)
    real(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    real(${KIND}$), intent(inout) :: bb(:,:)
    integer, intent(in) :: descb(DLEN_)
    real(${KIND}$), intent(out) :: ww(:)
    real(${KIND}$), intent(out) :: zz(:,:)
    integer, intent(in) :: descz(DLEN_)
    real(${KIND}$), intent(in), optional :: vl, vu
    integer, intent(in), optional :: il, iu
    character, intent(in), optional :: jobz, uplo
    integer, intent(in), optional :: nn, ia, ja, ib, jb, iz, jz
    real(${KIND}$), intent(inout), allocatable, optional :: work(:)
    integer, intent(inout), allocatable, optional :: iwork(:)
    logical, intent(in), optional :: skipchol
    integer, intent(out), optional :: mm, nz
    integer, intent(out), optional :: info

    real(${KIND}$) :: scale
    character :: jobz0, transa, uplo0
    logical :: skipchol0
    integer :: mm0, ncol

    @:inoptflags(jobz0, jobz, "V")
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(skipchol0, skipchol, .false.)

    ! Cholesky transformation of B
    if (.not. skipchol0) then
      call scalafx_ppotrf(bb, descb, uplo=uplo0, nn=nn, ia=ib, ja=jb, info=info)
    end if
    ! Reducing to standard form
    call scalafx_psygst(1, aa, desca, bb, descb, scale, uplo=uplo0, nn=nn, &
        & ia=ia, ja=ja, ib=ib, jb=jb, info=info)
    ! Solving eigenvalue problem.
    call scalafx_psyevr(aa, desca, ww, zz, descz, vl=vl, vu=vu, il=il, iu=iu, &
        & jobz=jobz0, uplo=uplo0, nn=nn, ia=ia, ja=ja, iz=iz, jz=jz, work=work, &
        & iwork=iwork, mm=mm0, nz=nz, info=info)
    ! Transforming eigenvectors back
    if (jobz0 == "V") then
      if (uplo0 == "L" .or. uplo0 == "l") then
         transa = "T"
      else
         transa = "N"
      end if
      if (present(nz)) then
        ncol = nz
      else
        ncol = mm0
      end if
      call scalafx_ptrsm(bb, descb, zz, descz, side="L", uplo=uplo0,&
          & transa=transa, diag="N", mm=nn, nn=ncol, ia=ib, ja=jb, ib=iz, jb=jz)
    end if
    @:handle_outoptflag(mm, mm0)

  end subroutine scalafx_psygvr_${TYPE}$

#:enddef scalafx_psygvr_template



#! ************************************************************************
#! *** phegvr
#! ************************************************************************

#:def scalafx_phegvr_template(TYPE, KIND)

  !> Solves Hermitian generalized eigenvalue problem by the MRRR algorithm.
  !!
  !! \details Invokes SCALAPACK routines p?potrf, p?hegst, p?heevr, p?trsm in
  !! order to transform the general eigenvalue problem to the standard form
  !! and transform the eigenvectors back. Currently all eigenvalues calculated.
  !!
  !! \param aa  Matrix to diagonalize (A), transformed matrix on exit.
  !! \param desca  Descriptor of matrix A.
  !! \param bb  Matrix on the right hand side (B), transformed matrix on exit.
  !! \param descb  Descriptor of matrix B.
  !! \param ww  Eigenvalues on exit.
  !! \param zz  Eigenvectors on exit (Z).
  !! \param descz  Descriptor of the eigenvector matrix.
  !! \param jobz  Job type (default: "V")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param ib  First row of the submatrix B (default: 1)
  !! \param jb  First column of the submatrix B (default: 1)
  !! \param iz  First row of the submatrix Z (default: 1)
  !! \param jz  First column of the submatrix Z (default: 1)
  !! \param work  Complex working array (if not specified, allocated
  !!     automatically)
  !! \param rwork Real working array (if not specified, allocated automatically)
  !! \param iwork Integer working array (if not specified, allocated
  !!     automatically)
  !! \param skipchol  If true, the Cholesky transformation will be skipped.
  !!     Array bb must have the Cholesky transformed form.
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !!
  !! \see SCALAPACK documentation (routines p?potrf, p?hegst, p?heevr, p?trsm).
  !!
  subroutine scalafx_phegvr_${TYPE}$(aa, desca, bb, descb, ww, zz, descz, jobz, uplo,&
      & ia, ja, ib, jb, iz, jz, work, rwork, iwork, skipchol, info)
    complex(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    complex(${KIND}$), intent(inout) :: bb(:,:)
    integer, intent(in) :: descb(DLEN_)
    real(${KIND}$), intent(out) :: ww(:)
    complex(${KIND}$), intent(out) :: zz(:,:)
    integer, intent(in) :: descz(DLEN_)
    character, intent(in), optional :: jobz, uplo
    integer, intent(in), optional :: ia, ja, ib, jb, iz, jz
    complex(${KIND}$), intent(inout), allocatable, optional :: work(:)
    real(${KIND}$), intent(inout), allocatable, optional :: rwork(:)
    integer, intent(inout), allocatable, optional :: iwork(:)
    logical, intent(in), optional :: skipchol
    integer, intent(out), optional :: info

    real(${KIND}$) :: scale
    character :: jobz0, transa, uplo0, range0
    logical :: skipchol0
    ! would be arguments for range /= "A" cases :
    real(${KIND}$) :: vl, vu
    integer :: il, iu
    integer :: m
    integer :: nz

    @:inoptflags(jobz0, jobz, "V")
    @:inoptflags(uplo0, uplo, "L")
    range0 = "A" ! all eigenvalue/vectors - need to generalise this
    @:inoptflags(skipchol0, skipchol, .false.)

    ! Cholesky transformation of B.
    if (.not. skipchol0) then
      call scalafx_ppotrf(bb, descb, uplo=uplo0, ia=ib, ja=jb, info=info)
    end if
    ! Reducing to standard form
    call scalafx_phegst(1, aa, desca, bb, descb, scale, uplo=uplo0, ia=ia, ja=ja,&
        & ib=ib, jb=jb, info=info)
    ! Solving eigenvalue problem.
    call scalafx_pheevr(aa, desca, ww, zz, descz, range0, vl, vu, il, iu, m, nz,&
        & jobz=jobz0, uplo=uplo0, ia=ia, ja=ja, iz=iz, jz=jz, work=work, &
        & rwork=rwork, iwork=iwork, info=info)
    ! Transforming eigenvectors back
    if (jobz0 == "V") then
      if (uplo0 == "L" .or. uplo0 == "l") then
         transa = "C"
      else
         transa = "N"
      end if
      call scalafx_ptrsm(bb, descb, zz, descz, side="L", uplo=uplo0,&
          & transa=transa, diag="N", ia=ib, ja=jb, ib=iz, jb=jz)
    end if

  end subroutine scalafx_phegvr_${TYPE}$

#:enddef scalafx_phegvr_template



#! ************************************************************************
#! *** pgesvd
#! ************************************************************************

#:def scalafx_r_pgesvd_template(TYPE, KIND)

  !> Singular value decomposition
  !!
  !! \param aa  Matrix to decompose (A).
  !! \param desca  Descriptor of matrix A.
  !! \param uu  Left singular vectors (U).
  !! \param descu  Descriptor of the left singular vectors.
  !! \param sigma  Singular values on exit.
  !! \param vt  Right singular vectors, transposed (Vt).
  !! \param descvt  Descriptor of the right singular vectors.
  !! \param jobu  Job type for U matrix (default: "V")
  !! \param jobu  Job type for vt matrix (default: "V")
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param iu  First row of the submatrix U (default: 1)
  !! \param ju  First column of the submatrix U (default: 1)
  !! \param ivt  First row of the submatrix vt (default: 1)
  !! \param jvt  First column of the submatrix vt (default: 1)
  !! \param mm  Number of columns of the matrix A (default: desca(M_))
  !! \param nn  Number of rows of the matrix A (default: desca(N_))
  !! \param work  Working array (if not specified, allocated automatically)
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !!
  !! \see SCALAPACK documentation (real routines p[sd]gesvd).
  !!
  subroutine scalafx_r_pgesvd_${TYPE}$(aa, desca, uu, descu, sigma, vt, descvt,&
      & jobu, jobvt, ia, ja, iu, ju, ivt, jvt, mm, nn, work, info)
    real(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    real(${KIND}$), intent(inout) :: uu(:,:)
    integer, intent(in) :: descu(DLEN_)
    real(${KIND}$), intent(out) :: sigma(:)
    real(${KIND}$), intent(out) :: vt(:,:)
    integer, intent(in) :: descvt(DLEN_)
    character, intent(in), optional :: jobu, jobvt
    integer, intent(in), optional :: ia, ja, iu, ju, ivt, jvt, mm, nn
    real(${KIND}$), intent(inout), allocatable, optional :: work(:)
    integer, intent(out), optional :: info

    integer :: mm0, nn0, lwork, lwmin, info0, ia0, ja0, iu0, ju0, ivt0, jvt0
    character :: jobu0, jobvt0
    real(${KIND}$) :: rtmp(1)
    real(${KIND}$), allocatable :: work0(:)

    ! Handle optional flags
    @:inoptflags(jobu0, jobu, "V")
    @:inoptflags(jobvt0, jobvt, "V")
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(iu0, iu, 1)
    @:inoptflags(ju0, ju, 1)
    @:inoptflags(ivt0, ivt, 1)
    @:inoptflags(jvt0, jvt, 1)
    @:inoptflags(mm0, mm, desca(M_))
    @:inoptflags(nn0, nn, desca(N_))


    ! Allocate  workspace
    call pgesvd(jobu0, jobvt0, mm0, nn0, aa, ia0, ja0, desca, sigma, uu, iu0, ju0, descu, &
        & vt, ivt0, jvt0, descvt, rtmp, -1, info0)
    call handle_infoflag(info0, "pgesvd in scalafx_r_pgesvd_${TYPE}$", info)
    lwmin = int(rtmp(1))
    @:move_minoptalloc(work0, lwmin, lwork, work)
    work0(:) = 0.0_${KIND}$

    ! SVD
    call pgesvd(jobu0, jobvt0, mm0, nn0, aa, ia0, ja0, desca, sigma, uu, iu0, ju0, descu, &
        & vt, ivt0, jvt0, descvt, work0, lwork, info0)
    call handle_infoflag(info0, "pgesvd in scalafx_r_pgesvd_${TYPE}$", info)

    ! Save work space allocations, if dummy arguments present
    @:optmovealloc(work0, work)

    end subroutine scalafx_r_pgesvd_${TYPE}$

#:enddef scalafx_r_pgesvd_template


#:def scalafx_c_pgesvd_template(TYPE, KIND)

  !> Singular value decomposition
  !!
  !! \param aa  Matrix to decompose (A).
  !! \param desca  Descriptor of matrix A.
  !! \param uu  Left singular vectors (U).
  !! \param descu  Descriptor of the left singular vectors.
  !! \param sigma  Singular values on exit.
  !! \param vt  Right singular vectors, transposed (Vt).
  !! \param descvt  Descriptor of the right singular vectors.
  !! \param jobu  Job type for U matrix (default: "V")
  !! \param jobu  Job type for vt matrix (default: "V")
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param iu  First row of the submatrix U (default: 1)
  !! \param ju  First column of the submatrix U (default: 1)
  !! \param ivt  First row of the submatrix vt (default: 1)
  !! \param jvt  First column of the submatrix vt (default: 1)
  !! \param mm  Number of columns of the matrix A (default: desca(M_))
  !! \param nn  Number of rows of the matrix A (default: desca(N_))
  !! \param work  Working array (if not specified, allocated automatically)
  !! \param rwork Real working array (if not specified, allocated automatically)
  !! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
  !!     subroutine stops.
  !!
  !! \see SCALAPACK documentation (complex routines p[cz]gesvd).
  !!
  subroutine scalafx_c_pgesvd_${TYPE}$(aa, desca, uu, descu, sigma, vt, descvt,&
      & jobu, jobvt, ia, ja, iu, ju, ivt, jvt, mm, nn, work, rwork, info)
    complex(${KIND}$), intent(inout) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    complex(${KIND}$), intent(inout) :: uu(:,:)
    integer, intent(in) :: descu(DLEN_)
    real(${KIND}$), intent(out) :: sigma(:)
    complex(${KIND}$), intent(out) :: vt(:,:)
    integer, intent(in) :: descvt(DLEN_)
    character, intent(in), optional :: jobu, jobvt
    integer, intent(in), optional :: ia, ja, iu, ju, ivt, jvt, mm, nn
    complex(${KIND}$), intent(inout), allocatable, optional :: work(:)
    real(${KIND}$), intent(inout), allocatable, optional :: rwork(:)
    integer, intent(out), optional :: info

    integer :: mm0, nn0, lwork, lrwmin, lrwork, lwmin, info0, ia0, ja0, iu0, ju0, ivt0, jvt0
    character :: jobu0, jobvt0
    complex(${KIND}$) :: ctmp(1)
    complex(${KIND}$), allocatable :: work0(:)
    real(${KIND}$) :: rtmp(1)
    real(${KIND}$), allocatable :: rwork0(:)

    ! Handle optional flags
    @:inoptflags(jobu0, jobu, "V")
    @:inoptflags(jobvt0, jobvt, "V")
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(iu0, iu, 1)
    @:inoptflags(ju0, ju, 1)
    @:inoptflags(ivt0, ivt, 1)
    @:inoptflags(jvt0, jvt, 1)
    @:inoptflags(mm0, mm, desca(M_))
    @:inoptflags(nn0, nn, desca(N_))

    ! Allocate  workspace
    call pgesvd(jobu0, jobvt0, mm0, nn0, aa, ia0, ja0, desca, sigma, uu, iu0,&
        & ju0, descu, vt, ivt0, jvt0, descvt, ctmp, -1, rtmp, info0)
    call handle_infoflag(info0, "pgesvd in scalafx_r_pgesvd_${TYPE}$", info)
    lwmin = int(real(ctmp(1)))
    lrwmin = int(rtmp(1))
    @:move_minoptalloc(work0, lwmin, lwork, work)
    @:move_minoptalloc(rwork0, lrwmin, lrwork, rwork)

    ! SVD
    call pgesvd(jobu0, jobvt0, mm0, nn0, aa, ia0, ja0, desca, sigma, uu, iu0,&
        & ju0, descu, vt, ivt0, jvt0, descvt, work0, lwmin, rwork0, info0)
    call handle_infoflag(info0, "pgesvd in scalafx_c_pgesvd_${TYPE}$", info)

    ! Save work space allocations, if dummy arguments present
    @:optmovealloc(work0, work)
    @:optmovealloc(rwork0, rwork)

  end subroutine scalafx_c_pgesvd_${TYPE}$

#:enddef scalafx_c_pgesvd_template



#! ************************************************************************
#! *** ptrsm
#! ************************************************************************

#:def scalafx_ptrsm_template(TYPE, FTYPE, HANDLE)

  !> Solves triangular matrix equation.
  !!
  !! \param aa  Left hand side of equation (A)
  !! \param desca  Descriptor of matrix A.
  !! \param bb  Right hand side (B).
  !! \param descb  Descriptor of matrix B.
  !! \param side  Side of A (default: "L")
  !! \param uplo  Upper or lower diagonal matrix (default: "L")
  !! \param transa  Transposition flag (default "N")
  !! \param diag  Specifieds whether matrix A is unit triangular (default: "N")
  !! \param alpha  Prefactor of B (default: 1.0)
  !! \param mm  Number of rows of the submatrix B (default: descb(M_))
  !! \param nn  Number of columns of the submatrix B (default: descb(N_))
  !! \param ia  First row of the submatrix A (default: 1)
  !! \param ja  First column of the submatrix A (default: 1)
  !! \param ib  First row of the submatrix B (default: 1)
  !! \param jb  First column of the submatrix B (default: 1)
  !!
  !! \see SCALAPACK documentation (routine p?trsm).
  !!
  subroutine scalafx_ptrsm_${TYPE}$(aa, desca, bb, descb, side, uplo, transa, diag,&
      & alpha, mm, nn, ia, ja, ib, jb)
    ${FTYPE}$, intent(in) :: aa(:,:)
    integer, intent(in) :: desca(DLEN_)
    ${FTYPE}$, intent(inout) :: bb(:,:)
    integer, intent(in) :: descb(DLEN_)
    character, intent(in), optional :: side, uplo, transa, diag
    ${FTYPE}$, intent(in), optional :: alpha
    integer, intent(in), optional :: mm, nn, ia, ja, ib, jb

    integer :: ia0, ja0, ib0, jb0, mm0, nn0
    character :: side0, uplo0, transa0, diag0
    ${FTYPE}$ :: alpha0

    @:inoptflags(side0, side, "L")
    @:inoptflags(uplo0, uplo, "L")
    @:inoptflags(transa0, transa, "N")
    @:inoptflags(diag0, diag, "N")
    @:inoptflags(alpha0, alpha, ${HANDLE}$)
    @:inoptflags(ia0, ia, 1)
    @:inoptflags(ja0, ja, 1)
    @:inoptflags(ib0, ib, 1)
    @:inoptflags(jb0, jb, 1)
    @:inoptflags(mm0, mm, descb(M_))
    @:inoptflags(nn0, nn, descb(N_))
    call ptrsm(side0, uplo0, transa0, diag0, mm0, nn0, alpha0, aa, ia0, ja0, &
        & desca, bb, ib0, jb0, descb)

  end subroutine scalafx_ptrsm_${TYPE}$

#:enddef scalafx_ptrsm_template


#! ************************************************************************
#! *** creatematrix
#! ************************************************************************

#:def scalafx_creatematrix_template(TYPE, FTYPE)

  !> Creates a distributed matrix and allocates local storage.
  !!
  !! \param mygrid  BLACS descriptor.
  !! \param mm  Number of rows of global matrix.
  !! \param nn  Number of columns of global matrix.
  !! \param mb  Row block size.
  !! \param nb  Column block size.
  !! \param desc  Matrix descriptor on exit.
  !! \param mtxloc  Allocated local matrix on exit.
  !! \param rsrc  Process row, over which first row is distributed
  !!     (default: lead row).
  !! \param csrc  Process column, over which first column is distributed
  !!     (default: lead column).
  !! \param info  Info flag.
  !!
  subroutine scalafx_creatematrix_${TYPE}$(mygrid, mm, nn, mb, nb, mtxloc,&
      & desc, rsrc, csrc, info)
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: mm, nn, mb, nb
    ${FTYPE}$, allocatable, intent(out) :: mtxloc(:,:)
    integer, intent(out) :: desc(DLEN_)
    integer, intent(in), optional :: rsrc, csrc
    integer, intent(out), optional :: info

    integer :: nrowloc, ncolloc

    call scalafx_getdescriptor(mygrid, mm, nn, mb, nb, desc, rsrc, csrc, info)
    call scalafx_getlocalshape(mygrid, desc, nrowloc, ncolloc)
    allocate(mtxloc(nrowloc, ncolloc))

  end subroutine scalafx_creatematrix_${TYPE}$

#:enddef scalafx_creatematrix_template

contains

  #:for TYPE in FLOAT_TYPES
    #:set FTYPE = FORTRAN_TYPES[TYPE]
    #:set KIND = FORTRAN_KINDS[TYPE]

    $:scalafx_ppotrf_template(TYPE, FTYPE)
    $:scalafx_ppotri_template(TYPE, FTYPE)
    $:scalafx_ptrtri_template(TYPE, FTYPE)
    $:scalafx_pgetrf_template(TYPE, FTYPE)
    $:scalafx_psygst_phegst_template(TYPE, FTYPE)

    #:if TYPE in REAL_TYPES
      $:scalafx_psyevr_template(TYPE, KIND)
      $:scalafx_psygvd_template(TYPE, KIND)
      $:scalafx_psyevd_template(TYPE, KIND)
      $:scalafx_psygv_template(TYPE, KIND)
      $:scalafx_psygvr_template(TYPE, KIND)
      $:scalafx_psyev_template(TYPE, KIND)
      $:scalafx_r_pgesvd_template(TYPE, KIND)
      $:scalafx_ptrsm_template(TYPE, FTYPE, "real(1.0, " + str(KIND) + ")")
    #:else
      $:scalafx_pheevr_template(TYPE, KIND)
      $:scalafx_phegvd_template(TYPE, KIND)
      $:scalafx_pheevd_template(TYPE, KIND)
      $:scalafx_phegv_template(TYPE, KIND)
      $:scalafx_phegvr_template(TYPE, KIND)
      $:scalafx_pheev_template(TYPE, KIND)
      $:scalafx_c_pgesvd_template(TYPE, KIND)
      $:scalafx_ptrsm_template(TYPE, FTYPE, "cmplx(1, 0, " + str(KIND) + ")")
    #:endif

    $:scalafx_creatematrix_template(TYPE, FTYPE)
  #:endfor

  $:scalafx_creatematrix_template("int", "integer")

  !> Returns descriptor and size for the local part of a distributed matrix.
  !!
  !! \param mygrid  BLACS descriptor.
  !! \param mm  Number of rows of global matrix.
  !! \param nn  Number of columns of global matrix.
  !! \param mb  Row block size.
  !! \param nb  Column block size.
  !! \param desc  Matrix descriptor on exit.
  !! \param nrowloc  Number of rows for local matrix.
  !! \param ncolloc  Number of columns for local matrix.
  !! \param rsrc  Process row, over which first row is distributed
  !!     (default: lead row).
  !! \param csrc  Process column, over which first column is distributed
  !!     (default: lead column).
  !! \param info  Info flag.
  !!
  subroutine scalafx_getdescriptor(mygrid, mm, nn, mb, nb, desc, rsrc, csrc, &
      & info)
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: mm, nn, mb, nb
    integer, intent(out) :: desc(DLEN_)
    integer, intent(in), optional :: rsrc, csrc
    integer, intent(out), optional :: info

    integer :: rsrc0, csrc0, info0
    integer :: nrowloc

    @:inoptflags(rsrc0, rsrc, mygrid%leadrow)
    @:inoptflags(csrc0, csrc, mygrid%leadcol)
    ! At least one local row, as per definition LLD >= 1 in the descriptor.
    nrowloc = max(1, numroc(mm, mb, mygrid%myrow, rsrc0, mygrid%nrow))
    call descinit(desc, mm, nn, mb, nb, rsrc0, csrc0, mygrid%ctxt, nrowloc, &
        & info0)
    call handle_infoflag(info0, "descinit in scalafx_getdescriptor", info)

  end subroutine scalafx_getdescriptor


  !> Returns the shape of the local part of a distributed array.
  !!
  !! \param mygrid  BLACS grid descriptor.
  !! \param desc  Global matrix descriptor.
  !! \param nrowloc  Nr. of local rows.
  !! \param ncolloc  Nr. of local columns.
  !!
  subroutine scalafx_getlocalshape(mygrid, desc, nrowloc, ncolloc)
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: desc(DLEN_)
    integer, intent(out) :: nrowloc, ncolloc

    nrowloc = max(1, numroc(desc(M_), desc(MB_), mygrid%myrow, desc(RSRC_), &
        & mygrid%nrow))
    ncolloc = numroc(desc(N_), desc(NB_), mygrid%mycol, desc(CSRC_), &
        & mygrid%ncol)

  end subroutine scalafx_getlocalshape

  !> Maps global position in a distributed matrix to local one.
  !!
  !! \param mygrid  BLACS descriptor.
  !! \param desc  Descriptor of the distributed matrix.
  !! \param grow  Global row index.
  !! \param gcol  Global column index.
  !! \param lrow  Local row index on output.
  !! \param lcol  Local column index on output.
  !! \param rsrc  Row of the process owning the local matrix.
  !! \param csrc  Column of the process owning the local matrix.
  !!
  subroutine scalafx_infog2l(mygrid, desc, grow, gcol, lrow, lcol, rsrc, csrc)
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: desc(DLEN_)
    integer, intent(in) :: grow, gcol
    integer, intent(out) :: lrow, lcol
    integer, intent(out) :: rsrc, csrc

    call infog2l(grow, gcol, desc, mygrid%nrow, mygrid%ncol, mygrid%myrow,&
        & mygrid%mycol, lrow, lcol, rsrc, csrc)

  end subroutine scalafx_infog2l


  !> Maps local row or column index onto global matrix position.
  !!
  !! \param indxloc  Local index on input.
  !! \param mygrid  BLACS descriptor.
  !! \param blocksize  Block size for direction (row or column)
  !!
  function scalafx_indxl2g(indxloc, crB, mycr, crsrc, ncr)
    integer :: scalafx_indxl2g
    integer, intent(in) :: indxloc, crB, mycr, crsrc, ncr

    scalafx_indxl2g = indxl2g(indxloc, crB, mycr, crsrc, ncr)

  end function scalafx_indxl2g

  !> Maps a global position in a distributed matrix to local one.
  !!
  subroutine scalafx_localindices(mygrid, desc, grow, gcol, local, lrow, lcol)

    !> BLACS descriptor.
    type(blacsgrid), intent(in) :: mygrid

    !> Descriptor of the distributed matrix.
    integer, intent(in) :: desc(DLEN_)

    !> Global row index.
    integer, intent(in) :: grow

    !> Global column index
    integer, intent(in) :: gcol

    !> Indicates whether given global index is local for the process.
    logical, intent(out) :: local

    !> Row index in the local matrix (or 0 if global index is not local)
    integer, intent(out) :: lrow

    !> Column index in the local matrix (or 0 if global index is not local)
    integer, intent(out) :: lcol

    !------------------------------------------------------------------------

    integer :: rsrc, csrc

    call infog2l(grow, gcol, desc, mygrid%nrow, mygrid%ncol, mygrid%myrow,&
        & mygrid%mycol, lrow, lcol, rsrc, csrc)
    local = (rsrc == mygrid%myrow .and. csrc == mygrid%mycol)
    if (.not. local) then
      lrow = 0
      lcol = 0
    end if

  end subroutine scalafx_localindices


end module scalapackfx_module
