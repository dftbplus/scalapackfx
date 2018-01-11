include(common.m4)

dnl ************************************************************************
dnl *** ppotrf
dnl ************************************************************************

define(`_subroutine_scalafx_ppotrf',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl
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
subroutine scalafx_ppotrf_$1(aa, desca, uplo, nn, ia, ja, info)
  $2, intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  character, intent(in), optional :: uplo
  integer, intent(in), optional :: nn, ia, ja
  integer, intent(out), optional :: info

  character :: uplo0
  integer :: nn0, ia0, ja0, info0

  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(nn0, nn, desca(M_))
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  call ppotrf(uplo0, nn0, aa, ia0, ja0, desca, info0)
  call handle_infoflag(info0, "ppotrf in scalafx_ppotrf_$1", info)

end subroutine scalafx_ppotrf_$1
')

dnl ************************************************************************
dnl *** ppotri
dnl ************************************************************************

define(`_subroutine_scalafx_ppotri',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl
!> Computes the inverse of a symmetric/Hermitian positive definite matrix.
!!
subroutine scalafx_ppotri_$1(aa, desca, uplo, ia, ja, nn, info)

  !> Cholesky decomposed matrix A on entry, inverse on exit.
  $2, intent(inout) :: aa(:,:)

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

  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(nn0, nn, desca(M_))
  call ppotri(uplo0, nn0, aa, ia0, ja0, desca, info0)
  call handle_infoflag(info0, "ppotri in scalafx_ppotri_$1", info)

end subroutine scalafx_ppotri_$1
')

dnl ************************************************************************
dnl *** ptrtri
dnl ************************************************************************

define(`_subroutine_scalafx_ptrtri',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl
!> Computes the inverse of a symmetric/Hermitian positive definite matrix.
!!
subroutine scalafx_ptrtri_$1(aa, desca, uplo, diag, ia, ja, nn, info)

  !> Cholesky decomposed matrix A on entry, inverse on exit.
  $2, intent(inout) :: aa(:,:)

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

  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(diag0, diag, "N")
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(nn0, nn, desca(M_))
  call ptrtri(uplo0, diag0, nn0, aa, ia0, ja0, desca, info0)
  call handle_infoflag(info0, "ptrtri in scalafx_ptrtri_$1", info)

end subroutine scalafx_ptrtri_$1
')

dnl ************************************************************************
dnl *** pgetrf
dnl ************************************************************************

define(`_subroutine_scalafx_pgetrf',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl
!> LU factorization of a general matrix with pivoting
!!
subroutine scalafx_pgetrf_$1(aa, desca, ipiv, ia, ja, mm, nn, info)

  !> LU decomposition on exit, pivoted by ipiv
  $2, intent(inout) :: aa(:,:)

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

  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(mm0, mm, desca(M_))
  _handle_inoptflag(nn0, nn, desca(N_))
  call pgetrf(mm0, nn0, aa, ia0, ja0, desca, ipiv, info0)
  call handle_infoflag(info0, "pgetrf in scalafx_pgetrf_$1", info)

end subroutine scalafx_pgetrf_$1
')

dnl ************************************************************************
dnl *** psygst / phegst
dnl ************************************************************************

define(`_subroutine_scalafx_psygst_phegst',`
dnl
dnl $1 subroutine suffix
dnl $2 real/complex dummy argument type
dnl $3 real dummy argument type
dnl $4 scalapack routine name
dnl
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
subroutine scalafx_$4_$1(ibtype, aa, desca, bb, descb, scale, uplo, &
    & nn, ia, ja, ib, jb, info)
  integer, intent(in) :: ibtype
  $2, intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  $2, intent(in) :: bb(:,:)
  integer, intent(in) :: descb(DLEN_)
  $3, intent(out) :: scale
  character, intent(in), optional :: uplo
  integer, intent(in), optional :: nn, ia, ja, ib, jb
  integer, intent(out), optional :: info

  integer :: nn0, ia0, ja0, ib0, jb0, info0
  character :: uplo0

  _handle_inoptflag(nn0, nn, desca(M_))
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(ib0, ib, 1)
  _handle_inoptflag(jb0, jb, 1)
  _handle_inoptflag(uplo0, uplo, "L")
  call $4(ibtype, uplo0, nn0, aa, ia0, ja0, desca, bb,&
      & ib0, jb0, descb, scale, info0)
  call handle_infoflag(info0, "$4 in scalafx_$4_$1", info)
  
end subroutine scalafx_$4_$1
')

dnl ************************************************************************
dnl *** psyev
dnl ************************************************************************

define(`_subroutine_scalafx_psyev',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument kind
dnl
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
subroutine scalafx_psyev_$1(aa, desca, ww, zz, descz, jobz, uplo, ia, ja, iz,&
    & jz, work, info)
  real($2), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  real($2), intent(out) :: ww(:)
  real($2), intent(out) :: zz(:,:)
  integer, intent(in) :: descz(DLEN_)
  character, intent(in), optional :: jobz, uplo
  integer, intent(in), optional :: ia, ja, iz, jz
  real($2), intent(inout), allocatable, optional :: work(:)
  integer, intent(out), optional :: info

  integer :: nn, lwork, info0, ia0, ja0, iz0, jz0
  character :: uplo0, jobz0
  real($2) :: rtmp(1)
  real($2), allocatable :: work0(:)

  ! Handle optional flags
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(jobz0, jobz, "V")
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(iz0, iz, 1)
  _handle_inoptflag(jz0, jz, 1)

  ! Allocate real workspace
  nn = desca(M_)
  call psyev(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
      & rtmp, -1, info0)
  call handle_infoflag(info0, "psev in scalafx_psyev", info)
  _move_minoptalloc(work0, int(rtmp(1)), lwork, work)

  ! Diagonalization
  ! Initializing workspace as SCALAPACK apparently accesses uninitialized
  ! elements in it (nagfors -nan flag causes *sometimes* arithmetic exception)
  work0(:) = 0.0_$2  
  call psyev(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
      & work0, lwork, info0)
  call handle_infoflag(info0, "psyev in scalafx_psyev_$1", info)

  ! Save work space allocations, if dummy arguments present
  _optmovealloc(work0, work)

end subroutine scalafx_psyev_$1
')


dnl ************************************************************************
dnl *** pheev
dnl ************************************************************************

define(`_subroutine_scalafx_pheev',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument kind
dnl
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
subroutine scalafx_pheev_$1(aa, desca, ww, zz, descz, jobz, uplo, ia, ja, iz,&
    & jz, work, rwork, info)
  complex($2), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  real($2), intent(out) :: ww(:)
  complex($2), intent(out) :: zz(:,:)
  integer, intent(in) :: descz(DLEN_)
  character, intent(in), optional :: jobz, uplo
  integer, intent(in), optional :: ia, ja, iz, jz
  complex($2), intent(inout), allocatable, optional :: work(:)
  real($2), intent(inout), allocatable, optional :: rwork(:)
  integer, intent(out), optional :: info

  integer :: nn, lwork, lrwork, lrwork_tmp, info0, ia0, ja0, iz0, jz0
  character :: uplo0, jobz0
  real($2) :: rtmp(1)
  complex($2) :: ctmp(1)
  complex($2), allocatable :: work0(:)
  real($2), allocatable :: rwork0(:)

  ! Handle optional flags
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(jobz0, jobz, "V")
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(iz0, iz, 1)
  _handle_inoptflag(jz0, jz, 1)

  ! Allocate  workspace
  nn = desca(M_)
  call pheev(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
      & ctmp, -1, rtmp, -1, info0)
  call handle_infoflag(info0, "pheev in scalafx_pheev", info)
  _move_minoptalloc(work0, int(ctmp(1)), lwork, work)

  ! Apparently pheev sometimes returns considerably smaller work space size as
  ! required in the documentation. We take the bigger value to be safe.
  if (jobz0 == "N") then
    lrwork_tmp = max(2 * nn, int(rtmp(1)))
  else
    lrwork_tmp = max(2 * nn + 2 * nn - 2, int(rtmp(1)))
  end if
  _move_minoptalloc(rwork0, lrwork_tmp, lrwork, rwork)

  ! Diagonalization
  ! Initializing workspace as SCALAPACK apparently accesses uninitialized
  ! elements in it (nagfors -nan flag causes *sometimes* arithmetic exception)
  work0(:) = cmplx(0, kind=$2)
  rwork0(:) = 0.0_$2
  call pheev(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
      & work0, lwork, rwork0, lrwork, info0)
  call handle_infoflag(info0, "pheev in scalafx_pheev_$1", info)

  ! Save work space allocations, if dummy arguments present
  _optmovealloc(rwork0, rwork)
  _optmovealloc(work0, work)


end subroutine scalafx_pheev_$1
')


dnl ************************************************************************
dnl *** psyevd
dnl ************************************************************************

define(`_subroutine_scalafx_psyevd',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument kind
dnl
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
subroutine scalafx_psyevd_$1(aa, desca, ww, zz, descz, jobz, uplo, ia, ja, iz,&
    & jz, work, iwork, allocfix, info)
  real($2), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  real($2), intent(out) :: ww(:)
  real($2), intent(out) :: zz(:,:)
  integer, intent(in) :: descz(DLEN_)
  character, intent(in), optional :: jobz, uplo
  integer, intent(in), optional :: ia, ja, iz, jz
  real($2), intent(inout), allocatable, optional :: work(:)
  integer, intent(inout), allocatable, optional :: iwork(:)
  logical, intent(in), optional :: allocfix
  integer, intent(out), optional :: info

  integer :: nn, liwork, lwork, lwmin, liwmin, info0, ia0, ja0, iz0, jz0
  character :: uplo0, jobz0
  integer :: itmp(1)
  real($2) :: rtmp(1), rtmp2(1)
  real($2), allocatable :: work0(:)
  integer, allocatable :: iwork0(:)
  logical :: allocfix0

  ! Handle optional flags
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(jobz0, jobz, "V")
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(iz0, iz, 1)
  _handle_inoptflag(jz0, jz, 1)
  _handle_inoptflag(allocfix0, allocfix, .false.)

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

  _move_minoptalloc(work0, lwmin, lwork, work)
  _move_minoptalloc(iwork0, liwmin, liwork, iwork)

  ! Diagonalization
  ! Initializing workspace as SCALAPACK apparently accesses uninitialized
  ! elements in it (nagfors -nan flag causes *sometimes* arithmetic exception)
  work0(:) = 0.0_$2
  iwork0(:) = 0
  call psyevd(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
      & work0, lwork, iwork0, liwork, info0)
  call handle_infoflag(info0, "psyevd in scalafx_psyevd_$1", info)

  ! Save work space allocations, if dummy arguments present
  _optmovealloc(iwork0, iwork)
  _optmovealloc(work0, work)

end subroutine scalafx_psyevd_$1
')

dnl ************************************************************************
dnl *** pheevd
dnl ************************************************************************

define(`_subroutine_scalafx_pheevd',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument kind
dnl
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
subroutine scalafx_pheevd_$1(aa, desca, ww, zz, descz, jobz, uplo, ia, ja, iz,&
    & jz, work, rwork, iwork, allocfix, info)
  complex($2), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  real($2), intent(out) :: ww(:)
  complex($2), intent(out) :: zz(:,:)
  integer, intent(in) :: descz(DLEN_)
  character, intent(in), optional :: jobz, uplo
  integer, intent(in), optional :: ia, ja, iz, jz
  complex($2), intent(inout), allocatable, optional :: work(:)
  real($2), intent(inout), allocatable, optional :: rwork(:)
  integer, intent(inout), allocatable, optional :: iwork(:)
  logical, intent(in), optional :: allocfix
  integer, intent(out), optional :: info

  integer :: nn, liwork, lwork, lrwork, info0, ia0, ja0, iz0, jz0
  integer :: lwmin, lrwmin, liwmin
  character :: uplo0, jobz0
  real($2) :: rtmp(1), rtmp2(1)
  complex($2) :: ctmp(1), ctmp2(1)
  integer :: itmp(1)
  complex($2), allocatable :: work0(:)
  real($2), allocatable :: rwork0(:)
  integer, allocatable :: iwork0(:)
  logical :: allocfix0

  ! Handle optional flags
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(jobz0, jobz, "V")
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(iz0, iz, 1)
  _handle_inoptflag(jz0, jz, 1)
  _handle_inoptflag(allocfix0, allocfix, .false.)

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

  _move_minoptalloc(work0, lwmin, lwork, work)
  _move_minoptalloc(rwork0, lrwmin, lrwork, rwork)
  _move_minoptalloc(iwork0, liwmin, liwork, iwork)

  ! Diagonalization
  ! Initializing workspace as SCALAPACK apparently accesses uninitialized
  ! elements in it (nagfors -nan flag causes *sometimes* arithmetic exception)
  work0(:) = cmplx(0, kind=$2)
  rwork0(:) = 0.0_$2
  iwork0(:) = 0
  call pheevd(jobz0, uplo0, nn, aa, ia0, ja0, desca, ww, zz, iz0, jz0, descz,&
      & work0, lwork, rwork0, lrwork, iwork0, liwork, info0)
  call handle_infoflag(info0, "pheevd in scalafx_pheevd_$1", info)

  ! Save work space allocations, if dummy arguments present
  _optmovealloc(work0, work)
  _optmovealloc(rwork0, rwork)
  _optmovealloc(iwork0, iwork)

end subroutine scalafx_pheevd_$1
')

dnl ************************************************************************
dnl *** psyevr
dnl ************************************************************************

define(`_subroutine_scalafx_psyevr',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument kind
dnl
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
subroutine scalafx_psyevr_$1(aa, desca, ww, zz, descz, vl, vu, il, iu, &
    & jobz, uplo, nn, ia, ja, iz, jz, work, iwork, mm, nz, info)
  real($2), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  real($2), intent(out) :: ww(:)
  real($2), intent(out) :: zz(:,:)
  integer, intent(in) :: descz(DLEN_)
  real($2), intent(in), optional :: vl, vu
  integer, intent(in), optional :: il, iu
  character, intent(in), optional :: jobz, uplo
  integer, intent(in), optional :: nn, ia, ja, iz, jz
  real($2), intent(inout), allocatable, optional :: work(:)
  integer, intent(inout), allocatable, optional :: iwork(:)
  integer, intent(out), optional :: mm, nz
  integer, intent(out), optional :: info
  
  integer :: liwork, lwork, lwmin, liwmin, info0, nn0, ia0, ja0, iz0, jz0
  integer :: mm0, nz0, il0, iu0
  real($2) :: vl0, vu0
  character :: uplo0, jobz0, range
  integer :: itmp(1)
  real($2) :: rtmp(1)
  real($2), allocatable :: work0(:)
  integer, allocatable :: iwork0(:)

  ! Handle optional flags
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(jobz0, jobz, "V")
  _handle_inoptflag(nn0, nn, desca(M_))
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(iz0, iz, 1)
  _handle_inoptflag(jz0, jz, 1)

  if (present(il) .or. present(iu)) then
    _handle_inoptflag(il0, il, 1)
    _handle_inoptflag(iu0, iu, nn0)
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
  call handle_infoflag(info0, "psyevr (workspace query) in scalafx_psyevr_$1", &
      & info)
  
  ! Check workspace size of psyevr() and take that.
  lwmin = int(rtmp(1))
  liwmin = itmp(1)
  
  _move_minoptalloc(work0, lwmin, lwork, work)
  _move_minoptalloc(iwork0, liwmin, liwork, iwork)
  
  ! Diagonalization
  ! Initializing workspace as SCALAPACK apparently accesses uninitialized
  ! elements in it (nagfors -nan flag causes *sometimes* arithmetic exception)
  work0(:) = 0.0_$2
  iwork0(:) = 0
  call psyevr(jobz0, range, uplo0, nn0, aa, ia0, ja0, desca, vl0, vu0, &
      & il0, iu0, mm0, nz0, ww, zz, iz0, jz0, descz, work0, lwork, iwork0, &
      & liwork, info0)
  call handle_infoflag(info0, "psyevr (diagonalization) in scalafx_psyevr_$1", &
      & info)
  
  ! Save work space allocations, if dummy arguments present
  _optmovealloc(iwork0, iwork)
  _optmovealloc(work0, work)

  _handle_outoptflag(mm, mm0)
  if (range == "V") then
    _handle_outoptflag(nz, nz0)
  end if
  
end subroutine scalafx_psyevr_$1
')

dnl ************************************************************************
dnl *** pheevr
dnl ************************************************************************

define(`_subroutine_scalafx_pheevr',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument kind
dnl
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
subroutine scalafx_pheevr_$1(aa, desca, ww, zz, descz, range, vl, vu, il, iu, &
    & m, nz, jobz, uplo, ia, ja, iz, jz, work, rwork, iwork, info)
  complex($2), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  real($2), intent(out) :: ww(:)
  complex($2), intent(out) :: zz(:,:)
  integer, intent(in) :: descz(DLEN_)
  character, intent(in) :: range
  real($2), intent(in) :: vl, vu
  integer, intent(in) :: il, iu
  integer, intent(out) :: m, nz
  character, intent(in), optional :: jobz, uplo
  integer, intent(in), optional :: ia, ja, iz, jz
  complex($2), intent(inout), allocatable, optional :: work(:)
  real($2), intent(inout), allocatable, optional :: rwork(:)
  integer, intent(inout), allocatable, optional :: iwork(:)
  integer, intent(out), optional :: info
  
  integer :: nn, liwork, lwork, lrwork, info0, ia0, ja0, iz0, jz0
  integer :: lwmin, lrwmin, liwmin
  character :: uplo0, jobz0
  real($2) :: rtmp(1)
  complex($2) :: ctmp(1)
  integer :: itmp(1)
  complex($2), allocatable :: work0(:)
  real($2), allocatable :: rwork0(:)
  integer, allocatable :: iwork0(:)
  
  ! Handle optional flags
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(jobz0, jobz, "V")
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(iz0, iz, 1)
  _handle_inoptflag(jz0, jz, 1)
  
  ! Allocate workspaces
  nn = desca(M_)
  call pheevr(jobz0, range, uplo0, nn, aa, ia0, ja0, desca, vl, vu, il, iu, &
      & m, nz, ww, zz, iz0, jz0, descz, ctmp, -1, rtmp, -1, itmp, -1, info0)
  call handle_infoflag(info0, "pheevr in scalafx_pheevr", info)
    
  lwmin = int(ctmp(1))
  lrwmin = int(rtmp(1))
  liwmin = itmp(1)
  
  _move_minoptalloc(work0, lwmin, lwork, work)
  _move_minoptalloc(rwork0, lrwmin, lrwork, rwork)
  _move_minoptalloc(iwork0, liwmin, liwork, iwork)
  
  ! Diagonalization
  ! Initializing workspace as SCALAPACK apparently accesses uninitialized
  ! elements in it (nagfors -nan flag causes *sometimes* arithmetic exception)
  work0(:) = cmplx(0, kind=$2)
  rwork0(:) = 0.0_$2
  iwork0(:) = 0
  call pheevr(jobz0, range, uplo0, nn, aa, ia0, ja0, desca, vl, vu, il, iu, &
      & m, nz, ww, zz, iz0, jz0, descz, &
      & work0, lwork, rwork0, lrwork, iwork0, liwork, info0)
  call handle_infoflag(info0, "pheevr in scalafx_pheevr_$1", info)
  
  ! Save work space allocations, if dummy arguments present
  _optmovealloc(work0, work)
  _optmovealloc(rwork0, rwork)
  _optmovealloc(iwork0, iwork)
  
end subroutine scalafx_pheevr_$1
')

dnl ************************************************************************
dnl *** psygv
dnl ************************************************************************

define(`_subroutine_scalafx_psygv',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument kind
dnl
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
subroutine scalafx_psygv_$1(aa, desca, bb, descb, ww, zz, descz, jobz, uplo,&
    & ia, ja, ib, jb, iz, jz, work, skipchol, info)
  real($2), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  real($2), intent(inout) :: bb(:,:)
  integer, intent(in) :: descb(DLEN_)
  real($2), intent(out) :: ww(:)
  real($2), intent(out) :: zz(:,:)
  integer, intent(in) :: descz(DLEN_)
  character, intent(in), optional :: jobz, uplo
  integer, intent(in), optional :: ia, ja, ib, jb, iz, jz
  real($2), intent(inout), allocatable, optional :: work(:)
  logical, intent(in), optional :: skipchol
  integer, intent(out), optional :: info

  real($2) :: scale
  character :: jobz0, transa, uplo0
  logical :: skipchol0

  _handle_inoptflag(jobz0, jobz, "V")
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(skipchol0, skipchol, .false.)

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

end subroutine scalafx_psygv_$1
')

dnl ************************************************************************
dnl *** phegv
dnl ************************************************************************

define(`_subroutine_scalafx_phegv',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument kind
dnl
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
subroutine scalafx_phegv_$1(aa, desca, bb, descb, ww, zz, descz, jobz, uplo,&
    & ia, ja, ib, jb, iz, jz, work, rwork, skipchol, info)
  complex($2), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  complex($2), intent(inout) :: bb(:,:)
  integer, intent(in) :: descb(DLEN_)
  real($2), intent(out) :: ww(:)
  complex($2), intent(out) :: zz(:,:)
  integer, intent(in) :: descz(DLEN_)
  character, intent(in), optional :: jobz, uplo
  integer, intent(in), optional :: ia, ja, ib, jb, iz, jz
  complex($2), intent(inout), allocatable, optional :: work(:)
  real($2), intent(inout), allocatable, optional :: rwork(:)
  logical, intent(in), optional :: skipchol
  integer, intent(out), optional :: info

  real($2) :: scale
  character :: jobz0, transa, uplo0
  logical :: skipchol0

  _handle_inoptflag(jobz0, jobz, "V")
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(skipchol0, skipchol, .false.)

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

end subroutine scalafx_phegv_$1
')

dnl ************************************************************************
dnl *** psygvd
dnl ************************************************************************

define(`_subroutine_scalafx_psygvd',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument kind
dnl
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
!!     values. See the scalafx_psyevd_$1() routine for details.
!! \param skipchol  If true, the Cholesky transformation will be skipped. 
!!     Array bb must have the Cholesky transformed form.
!! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
!!     subroutine stops.
!! \see SCALAPACK documentation (routines p?potrf, p?sygst, p?syevd, p?trsm).
subroutine scalafx_psygvd_$1(aa, desca, bb, descb, ww, zz, descz, jobz, uplo,&
    & ia, ja, ib, jb, iz, jz, work, iwork, allocfix, skipchol, info)
  real($2), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  real($2), intent(inout) :: bb(:,:)
  integer, intent(in) :: descb(DLEN_)
  real($2), intent(out) :: ww(:)
  real($2), intent(out) :: zz(:,:)
  integer, intent(in) :: descz(DLEN_)
  character, intent(in), optional :: jobz, uplo
  integer, intent(in), optional :: ia, ja, ib, jb, iz, jz
  real($2), intent(inout), allocatable, optional :: work(:)
  integer, intent(inout), allocatable, optional :: iwork(:)
  logical, intent(in), optional :: allocfix, skipchol
  integer, intent(out), optional :: info

  real($2) :: scale
  character :: jobz0, transa, uplo0
  logical :: skipchol0

  _handle_inoptflag(jobz0, jobz, "V")
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(skipchol0, skipchol, .false.)

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

end subroutine scalafx_psygvd_$1
')

dnl ************************************************************************
dnl *** phegvd
dnl ************************************************************************

define(`_subroutine_scalafx_phegvd',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument kind
dnl
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
!!     values. See the scalafx_pheevd_$1() routine for details.
!! \param skipchol  If true, the Cholesky transformation will be skipped. 
!!     Array bb must have the Cholesky transformed form.
!! \param info  Info flag. If not specified and SCALAPACK calls returns nozero,
!!     subroutine stops.
!! \see SCALAPACK documentation (routines p?potrf, p?hegst, p?heevd, p?trsm).
subroutine scalafx_phegvd_$1(aa, desca, bb, descb, ww, zz, descz, jobz, uplo,&
    & ia, ja, ib, jb, iz, jz, work, rwork, iwork, allocfix, skipchol, info)
  complex($2), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  complex($2), intent(inout) :: bb(:,:)
  integer, intent(in) :: descb(DLEN_)
  real($2), intent(out) :: ww(:)
  complex($2), intent(out) :: zz(:,:)
  integer, intent(in) :: descz(DLEN_)
  character, intent(in), optional :: jobz, uplo
  integer, intent(in), optional :: ia, ja, ib, jb, iz, jz
  complex($2), intent(inout), allocatable, optional :: work(:)
  real($2), intent(inout), allocatable, optional :: rwork(:)
  integer, intent(inout), allocatable, optional :: iwork(:)
  logical, intent(in), optional :: allocfix, skipchol
  integer, intent(out), optional :: info

  real($2) :: scale
  character :: jobz0, transa, uplo0
  logical :: skipchol0

  _handle_inoptflag(jobz0, jobz, "V")
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(skipchol0, skipchol, .false.)

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

end subroutine scalafx_phegvd_$1
')

dnl ************************************************************************
dnl *** psygvr
dnl ************************************************************************

define(`_subroutine_scalafx_psygvr',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument kind
dnl
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
subroutine scalafx_psygvr_$1(aa, desca, bb, descb, ww, zz, descz, &
    & vl, vu, il, iu, jobz, uplo, nn, ia, ja, ib, jb, iz, jz, work, iwork, &
    & skipchol, mm, nz, info)
  real($2), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  real($2), intent(inout) :: bb(:,:)
  integer, intent(in) :: descb(DLEN_)
  real($2), intent(out) :: ww(:)
  real($2), intent(out) :: zz(:,:)
  integer, intent(in) :: descz(DLEN_)
  real($2), intent(in), optional :: vl, vu
  integer, intent(in), optional :: il, iu
  character, intent(in), optional :: jobz, uplo
  integer, intent(in), optional :: nn, ia, ja, ib, jb, iz, jz
  real($2), intent(inout), allocatable, optional :: work(:)
  integer, intent(inout), allocatable, optional :: iwork(:)
  logical, intent(in), optional :: skipchol
  integer, intent(out), optional :: mm, nz
  integer, intent(out), optional :: info
  
  real($2) :: scale
  character :: jobz0, transa, uplo0
  logical :: skipchol0
  integer :: mm0, ncol
  
  _handle_inoptflag(jobz0, jobz, "V")
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(skipchol0, skipchol, .false.)
  
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
  _handle_outoptflag(mm, mm0)

end subroutine scalafx_psygvr_$1
')

dnl ************************************************************************
dnl *** phegvr
dnl ************************************************************************

define(`_subroutine_scalafx_phegvr',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument kind
dnl
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
subroutine scalafx_phegvr_$1(aa, desca, bb, descb, ww, zz, descz, jobz, uplo,&
    & ia, ja, ib, jb, iz, jz, work, rwork, iwork, skipchol, info)
  complex($2), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  complex($2), intent(inout) :: bb(:,:)
  integer, intent(in) :: descb(DLEN_)
  real($2), intent(out) :: ww(:)
  complex($2), intent(out) :: zz(:,:)
  integer, intent(in) :: descz(DLEN_)
  character, intent(in), optional :: jobz, uplo
  integer, intent(in), optional :: ia, ja, ib, jb, iz, jz
  complex($2), intent(inout), allocatable, optional :: work(:)
  real($2), intent(inout), allocatable, optional :: rwork(:)
  integer, intent(inout), allocatable, optional :: iwork(:)
  logical, intent(in), optional :: skipchol
  integer, intent(out), optional :: info
  
  real($2) :: scale
  character :: jobz0, transa, uplo0, range0
  logical :: skipchol0
  ! would be arguments for range /= "A" cases :
  real($2) :: vl, vu
  integer :: il, iu
  integer :: m
  integer :: nz
    
  _handle_inoptflag(jobz0, jobz, "V")
  _handle_inoptflag(uplo0, uplo, "L")  
  range0 = "A" ! all eigenvalue/vectors - need to generalise this
  _handle_inoptflag(skipchol0, skipchol, .false.)
  
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

end subroutine scalafx_phegvr_$1
')

dnl ************************************************************************
dnl *** pgesvd
dnl ************************************************************************

define(`_subroutine_scalafx_r_pgesvd',`
dnl
dnl $1 subroutine suffix
dnl $2: type letter
dnl $3 dummy argument kind
dnl
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
subroutine scalafx_r_pgesvd_$1(aa, desca, uu, descu, sigma, vt, descvt, jobu, jobvt, &
      & ia, ja, iu, ju, ivt, jvt, mm, nn, work, info)
  real($3), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  real($3), intent(inout) :: uu(:,:)
  integer, intent(in) :: descu(DLEN_)
  real($3), intent(out) :: sigma(:)
  real($3), intent(out) :: vt(:,:)
  integer, intent(in) :: descvt(DLEN_)
  character, intent(in), optional :: jobu, jobvt
  integer, intent(in), optional :: ia, ja, iu, ju, ivt, jvt, mm, nn
  real($3), intent(inout), allocatable, optional :: work(:)
  integer, intent(out), optional :: info

  integer :: mm0, nn0, lwork, lwmin, info0, ia0, ja0, iu0, ju0, ivt0, jvt0
  character :: jobu0, jobvt0
  real($3) :: rtmp(1)
  real($3), allocatable :: work0(:)

  ! Handle optional flags
  _handle_inoptflag(jobu0, jobu, "V")
  _handle_inoptflag(jobvt0, jobvt, "V")
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(iu0, iu, 1)
  _handle_inoptflag(ju0, ju, 1)
  _handle_inoptflag(ivt0, ivt, 1)
  _handle_inoptflag(jvt0, jvt, 1)
  _handle_inoptflag(mm0, mm, desca(M_))
  _handle_inoptflag(nn0, nn, desca(N_))


  ! Allocate  workspace
  call pgesvd(jobu0, jobvt0, mm0, nn0, aa, ia0, ja0, desca, sigma, uu, iu0, ju0, descu, &
      & vt, ivt0, jvt0, descvt, rtmp, -1, info0)
  call handle_infoflag(info0, "pgesvd in scalafx_r_pgesvd_$1", info)
  lwmin = int(rtmp(1))
  _move_minoptalloc(work0, lwmin, lwork, work)
  work0(:) = 0.0_$3

  ! SVD
  call pgesvd(jobu0, jobvt0, mm0, nn0, aa, ia0, ja0, desca, sigma, uu, iu0, ju0, descu, &
      & vt, ivt0, jvt0, descvt, work0, lwork, info0)
  call handle_infoflag(info0, "pgesvd in scalafx_r_pgesvd_$1", info)

  ! Save work space allocations, if dummy arguments present
  _optmovealloc(work0, work)

  end subroutine scalafx_r_pgesvd_$1
')

define(`_subroutine_scalafx_c_pgesvd',`
dnl
dnl $1 subroutine suffix
dnl $2: type letter
dnl $3 dummy argument kind
dnl
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
subroutine scalafx_c_pgesvd_$1(aa, desca, uu, descu, sigma, vt, descvt, jobu, jobvt, &
      & ia, ja, iu, ju, ivt, jvt, mm, nn, work, rwork, info)
  complex($3), intent(inout) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  complex($3), intent(inout) :: uu(:,:)
  integer, intent(in) :: descu(DLEN_)
  real($3), intent(out) :: sigma(:)
  complex($3), intent(out) :: vt(:,:)
  integer, intent(in) :: descvt(DLEN_)
  character, intent(in), optional :: jobu, jobvt
  integer, intent(in), optional :: ia, ja, iu, ju, ivt, jvt, mm, nn
  complex($3), intent(inout), allocatable, optional :: work(:)
  real($3), intent(inout), allocatable, optional :: rwork(:)
  integer, intent(out), optional :: info

  integer :: mm0, nn0, lwork, lrwmin, lrwork, lwmin, info0, ia0, ja0, iu0, ju0, ivt0, jvt0
  character :: jobu0, jobvt0
  complex($3) :: ctmp(1)
  complex($3), allocatable :: work0(:)
  real($3) :: rtmp(1)
  real($3), allocatable :: rwork0(:)

  ! Handle optional flags
  _handle_inoptflag(jobu0, jobu, "V")
  _handle_inoptflag(jobvt0, jobvt, "V")
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(iu0, iu, 1)
  _handle_inoptflag(ju0, ju, 1)
  _handle_inoptflag(ivt0, ivt, 1)
  _handle_inoptflag(jvt0, jvt, 1)
  _handle_inoptflag(mm0, mm, desca(M_))
  _handle_inoptflag(nn0, nn, desca(N_))

  ! Allocate  workspace
  call pgesvd(jobu0, jobvt0, mm0, nn0, aa, ia0, ja0, desca, sigma, uu, iu0, ju0, descu, &
      & vt, ivt0, jvt0, descvt, ctmp, -1, rtmp, info0)
  call handle_infoflag(info0, "pgesvd in scalafx_r_pgesvd_$1", info)
  lwmin = int(real(ctmp(1)))
  lrwmin = int(rtmp(1))
  _move_minoptalloc(work0, lwmin, lwork, work)
  _move_minoptalloc(rwork0, lrwmin, lrwork, rwork)

  ! SVD
  call pgesvd(jobu0, jobvt0, mm0, nn0, aa, ia0, ja0, desca, sigma, uu, iu0, ju0, descu, &
      & vt, ivt0, jvt0, descvt, work0, lwmin, rwork0, info0)
  call handle_infoflag(info0, "pgesvd in scalafx_c_pgesvd_$1", info)

  ! Save work space allocations, if dummy arguments present
  _optmovealloc(work0, work)
  _optmovealloc(rwork0, rwork)

  end subroutine scalafx_c_pgesvd_$1
')


dnl ************************************************************************
dnl *** ptrsm
dnl ************************************************************************

define(`_subroutine_scalafx_ptrsm',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl
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
subroutine scalafx_ptrsm_$1(aa, desca, bb, descb, side, uplo, transa, diag,&
    & alpha, mm, nn, ia, ja, ib, jb)
  $2, intent(in) :: aa(:,:)
  integer, intent(in) :: desca(DLEN_)
  $2, intent(inout) :: bb(:,:)
  integer, intent(in) :: descb(DLEN_)
  character, intent(in), optional :: side, uplo, transa, diag
  $2, intent(in), optional :: alpha
  integer, intent(in), optional :: mm, nn, ia, ja, ib, jb

  integer :: ia0, ja0, ib0, jb0, mm0, nn0
  character :: side0, uplo0, transa0, diag0
  $2 :: alpha0

  _handle_inoptflag(side0, side, "L")
  _handle_inoptflag(uplo0, uplo, "L")
  _handle_inoptflag(transa0, transa, "N")
  _handle_inoptflag(diag0, diag, "N")
  _handle_inoptflag(alpha0, alpha, $3)
  _handle_inoptflag(ia0, ia, 1)
  _handle_inoptflag(ja0, ja, 1)
  _handle_inoptflag(ib0, ib, 1)
  _handle_inoptflag(jb0, jb, 1)
  _handle_inoptflag(mm0, mm, descb(M_))
  _handle_inoptflag(nn0, nn, descb(N_))
  call ptrsm(side0, uplo0, transa0, diag0, mm0, nn0, alpha0, aa, ia0, ja0, &
      & desca, bb, ib0, jb0, descb)

end subroutine scalafx_ptrsm_$1
')

dnl ************************************************************************
dnl *** creatematrix
dnl ************************************************************************

define(`_subroutine_scalafx_creatematrix',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl
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
!!     (default: master row).
!! \param csrc  Process column, over which first column is distributed
!!     (default: master column).
!! \param info  Info flag.
!!
subroutine scalafx_creatematrix_$1(mygrid, mm, nn, mb, nb, mtxloc, desc,&
    & rsrc, csrc, info)
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: mm, nn, mb, nb
  $2, allocatable, intent(out) :: mtxloc(:,:)
  integer, intent(out) :: desc(DLEN_)
  integer, intent(in), optional :: rsrc, csrc
  integer, intent(out), optional :: info

  integer :: nrowloc, ncolloc
  
  call scalafx_getdescriptor(mygrid, mm, nn, mb, nb, desc, rsrc, csrc, info)
  call scalafx_getlocalshape(mygrid, desc, nrowloc, ncolloc)
  allocate(mtxloc(nrowloc, ncolloc))
  
end subroutine scalafx_creatematrix_$1
')
