include(common.m4)

dnl ************************************************************************
dnl *** ppotrf
dnl ************************************************************************ 

define(`_subroutine_interface_ppotrf',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments type
!> Cholesky factorization of symmetric/Hermitian pos.def. matrix ($1).
subroutine p$2potrf(uplo, nn, aa, ia, ja, desca, info)
  import
  character, intent(in) :: uplo
  integer, intent(in) :: nn
  integer, intent(in) :: ia, ja, desca(DLEN_)
  $3, intent(inout) :: aa(desca(LLD_), *)
  integer, intent(out) :: info
end subroutine p$2potrf
')

dnl ************************************************************************
dnl *** ppotri
dnl ************************************************************************ 

define(`_subroutine_interface_ppotri',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments type
!> Inversion of a Cholesky decomposed symmetric/Hermitian matrix ($1).
subroutine p$2potri(uplo, nn, aa, ia, ja, desca, info)
  import
  character, intent(in) :: uplo
  integer, intent(in) :: nn
  integer, intent(in) :: ia, ja, desca(DLEN_)
  $3, intent(inout) :: aa(desca(LLD_), *)
  integer, intent(out) :: info
end subroutine p$2potri
')

dnl ************************************************************************
dnl *** ptrtri
dnl ************************************************************************ 

define(`_subroutine_interface_ptrtri',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments type
!> Inversion of a Cholesky decomposed symmetric/Hermitian matrix ($1).
subroutine p$2trtri(uplo, diag, nn, aa, ia, ja, desca, info)
  import
  character, intent(in) :: uplo, diag
  integer, intent(in) :: nn
  integer, intent(in) :: ia, ja, desca(DLEN_)
  $3, intent(inout) :: aa(desca(LLD_), *)
  integer, intent(out) :: info
end subroutine p$2trtri
')

dnl ************************************************************************
dnl *** psygst
dnl ************************************************************************ 

define(`_subroutine_interface_psygst',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments kind
!> Reduces generalized symmetric eigenvalue problem to standard form ($1).
subroutine p$2sygst(ibtype, uplo, nn, aa, ia, ja, desca, bb, ib, jb, descb,&
    & scale, info)
  import
  integer, intent(in) :: ibtype
  character, intent(in) :: uplo
  integer, intent(in) :: nn
  integer, intent(in) :: ia, ja, desca(DLEN_)
  real($3), intent(inout) :: aa(desca(LLD_), *)
  integer, intent(in) :: ib, jb, descb(DLEN_)
  real($3), intent(in) :: bb(descb(LLD_), *)
  real($3), intent(out) :: scale
  integer, intent(out) :: info
end subroutine p$2sygst
')

dnl ************************************************************************
dnl *** phegst
dnl ************************************************************************ 

define(`_subroutine_interface_phegst',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments kind
!> Reduces generalized Hermitian eigenvalue problem to standard form ($1).
subroutine p$2hegst(ibtype, uplo, nn, aa, ia, ja, desca, bb, ib, jb, descb,&
    & scale, info)
  import
  integer, intent(in) :: ibtype
  character, intent(in) :: uplo
  integer, intent(in) :: nn
  integer, intent(in) :: ia, ja, desca(DLEN_)
  complex($3), intent(inout) :: aa(desca(LLD_), *)
  integer, intent(in) :: ib, jb, descb(DLEN_)
  complex($3), intent(in) :: bb(descb(LLD_), *)
  real($3), intent(out) :: scale
  integer, intent(out) :: info
end subroutine p$2hegst
')

dnl ************************************************************************
dnl *** psyev
dnl ************************************************************************ 

define(`_subroutine_interface_psyev',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments kind
!> Eigenvalues and eigenvectors by divide and conquer algorithm ($1)
subroutine p$2syev(jobz, uplo, nn, aa, ia, ja, desca, ww, zz, iz, jz,&
    & descz, work, lwork, info)
  import
  character, intent(in) :: jobz, uplo
  integer, intent(in) :: nn
  integer, intent(in) :: ia, ja, desca(DLEN_)
  real($3), intent(inout) :: aa(desca(LLD_), *)
  integer, intent(in) :: iz, jz, descz(DLEN_)
  real($3), intent(out) :: ww(nn), zz(descz(LLD_),*)
  real($3), intent(inout) :: work(*)
  integer, intent(in) :: lwork
  integer, intent(out) :: info
end subroutine p$2syev
')

dnl ************************************************************************
dnl *** pheev
dnl ************************************************************************ 

define(`_subroutine_interface_pheev',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments kind
!> Eigenvalues and eigenvectors by divide and conquer algorithm ($1)
subroutine p$2heev(jobz, uplo, nn, aa, ia, ja, desca, ww, zz, iz, jz,&
    & descz, work, lwork, rwork, lrwork, info)
  import
  character, intent(in) :: jobz, uplo
  integer, intent(in) :: nn
  integer, intent(in) :: ia, ja, desca(DLEN_)
  complex($3), intent(inout) :: aa(desca(LLD_), *)
  integer, intent(in) :: iz, jz, descz(DLEN_)
  real($3), intent(out) :: ww(nn)
  complex($3), intent(out) ::  zz(descz(LLD_),*)
  complex($3), intent(inout) :: work(*)
  integer, intent(in) :: lwork
  real($3), intent(inout) :: rwork(*)
  integer, intent(in) :: lrwork
  integer, intent(out) :: info
end subroutine p$2heev
')

dnl ************************************************************************
dnl *** psyevd
dnl ************************************************************************ 

define(`_subroutine_interface_psyevd',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments kind
!> Eigenvalues and eigenvectors by divide and conquer algorithm ($1)
subroutine p$2syevd(jobz, uplo, nn, aa, ia, ja, desca, ww, zz, iz, jz,&
    & descz, work, lwork, iwork, liwork, info)
  import
  character, intent(in) :: jobz, uplo
  integer, intent(in) :: nn
  integer, intent(in) :: ia, ja, desca(DLEN_)
  real($3), intent(inout) :: aa(desca(LLD_), *)
  integer, intent(in) :: iz, jz, descz(DLEN_)
  real($3), intent(out) :: ww(nn), zz(descz(LLD_),*)
  real($3), intent(inout) :: work(*)
  integer, intent(in) :: lwork
  integer, intent(inout) :: iwork(*)
  integer, intent(in) :: liwork
  integer, intent(out) :: info
end subroutine p$2syevd
')

dnl ************************************************************************
dnl *** pheevd
dnl ************************************************************************ 

define(`_subroutine_interface_pheevd',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments kind
!> Eigenvalues and eigenvectors by divide and conquer algorithm ($1)
subroutine p$2heevd(jobz, uplo, nn, aa, ia, ja, desca, ww, zz, iz, jz,&
    & descz, work, lwork, rwork, lrwork, iwork, liwork, info)
  import
  character, intent(in) :: jobz, uplo
  integer, intent(in) :: nn
  integer, intent(in) :: ia, ja, desca(DLEN_)
  complex($3), intent(inout) :: aa(desca(LLD_), *)
  integer, intent(in) :: iz, jz, descz(DLEN_)
  real($3), intent(out) :: ww(nn)
  complex($3), intent(out) ::  zz(descz(LLD_),*)
  complex($3), intent(inout) :: work(*)
  integer, intent(in) :: lwork
  real($3), intent(inout) :: rwork(*)
  integer, intent(in) :: lrwork
  integer, intent(inout) :: iwork(*)
  integer, intent(in) :: liwork
  integer, intent(out) :: info
end subroutine p$2heevd
')

dnl ************************************************************************
dnl *** psyevr
dnl ************************************************************************ 

define(`_subroutine_interface_psyevr',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments kind
!> Eigenvalues and eigenvectors by MRRR algorithm ($1)
subroutine p$2syevr(jobz, range, uplo, nn, aa, ia, ja, desca, vl, vu, il, iu, &
    & mm, nz, ww, zz, iz, jz, descz, work, lwork, iwork, liwork, info)
  import
  character, intent(in) :: jobz, range, uplo
  integer, intent(in) :: nn
  integer, intent(in) :: desca(DLEN_)
  real($3), intent(inout) :: aa(desca(LLD_), *)
  integer, intent(in) :: ia, ja
  real($3), intent(in) :: vl, vu
  integer, intent(in) :: il, iu
  integer, intent(out) :: mm, nz
  real($3), intent(out) :: ww(nn)
  integer, intent(in) :: descz(DLEN_)
  real($3), intent(out) :: zz(descz(LLD_),*)
  integer, intent(in) :: iz, jz
  real($3), intent(inout) :: work(*)
  integer, intent(in) :: lwork
  integer, intent(inout) :: iwork(*)
  integer, intent(in) :: liwork
  integer, intent(out) :: info
end subroutine p$2syevr
')

dnl ************************************************************************
dnl *** pheevr
dnl ************************************************************************ 

define(`_subroutine_interface_pheevr',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments kind
!> Eigenvalues and eigenvectors by MRRR algorithm ($1)
subroutine p$2heevr(jobz, range, uplo, nn, aa, ia, ja, desca, vl, vu, il, iu, &
    & mm, nz, ww, zz, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork,&
    & info)
  import
  character, intent(in) :: jobz, range, uplo
  integer, intent(in) :: nn
  integer, intent(in) :: desca(DLEN_)
  complex($3), intent(inout) :: aa(desca(LLD_), *)
  integer, intent(in) :: ia, ja
  real($3), intent(in) :: vl, vu
  integer, intent(in) :: il, iu
  integer, intent(out) :: mm, nz
  real($3), intent(out) :: ww(nn)
  integer, intent(in) :: descz(DLEN_)
  complex($3), intent(out) ::  zz(descz(LLD_),*)
  integer, intent(in) :: iz, jz
  complex($3), intent(inout) :: work(*)
  integer, intent(in) :: lwork
  real($3), intent(inout) :: rwork(*)
  integer, intent(in) :: lrwork
  integer, intent(inout) :: iwork(*)
  integer, intent(in) :: liwork
  integer, intent(out) :: info
end subroutine p$2heevr
')

dnl ************************************************************************
dnl *** prgesvd
dnl ************************************************************************

define(`_subroutine_interface_prgesvd',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments kind
!> Singular values and vectors ($1)
subroutine p$2gesvd(jobu, jobvt, mm, nn, aa, ia, ja, desca, sigma, uu, iu, ju, descu, &
    & vt, ivt, jvt, descvt, work, lwork, info)
  import
  character, intent(in) :: jobu, jobvt
  integer, intent(in) :: mm, nn
  integer, intent(in) :: ia, ja, desca(DLEN_)
  real($3), intent(inout) :: aa(desca(LLD_), *)
  real($3), intent(out) :: sigma(*)
  integer, intent(in) :: iu, ju, descu(DLEN_)
  real($3), intent(out) :: uu(descu(LLD_), *)
  integer, intent(in) :: ivt, jvt, descvt(DLEN_)
  real($3), intent(out) :: vt(descvt(LLD_), *)
  real($3), intent(inout) :: work(*)
  integer, intent(in) :: lwork
  integer, intent(out) :: info
end subroutine p$2gesvd
')

dnl ************************************************************************
dnl *** pcgesvd
dnl ************************************************************************

define(`_subroutine_interface_pcgesvd',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments kind
!> Singular values and vectors ($1)
subroutine p$2gesvd(jobu, jobvt, mm, nn, aa, ia, ja, desca, sigma, uu, iu, ju, descu, &
    & vt, ivt, jvt, descvt, work, lwork, rwork, info)
  import
  character, intent(in) :: jobu, jobvt
  integer, intent(in) :: mm, nn
  integer, intent(in) :: ia, ja, desca(DLEN_)
  complex($3), intent(inout) :: aa(desca(LLD_), *)
  real($3), intent(out) :: sigma(*)
  integer, intent(in) :: iu, ju, descu(DLEN_)
  complex($3), intent(out) :: uu(descu(LLD_), *)
  integer, intent(in) :: ivt, jvt, descvt(DLEN_)
  complex($3), intent(out) :: vt(descvt(LLD_), *)
  complex($3), intent(inout) :: work(*)
  integer, intent(in) :: lwork
  real($3), intent(inout) :: rwork(*)
  integer, intent(out) :: info
end subroutine p$2gesvd
')

dnl ************************************************************************
dnl *** ptrsm
dnl ************************************************************************ 

define(`_subroutine_interface_ptrsm',`
dnl $1: comment
dnl $2: type letter
dnl $3: dummy arguments type
!> Solves a triangular matrix equation ($1).
subroutine p$2trsm(side, uplo, transa, diag, mm, nn, alpha, aa, ia, ja,&
    & desca, bb, ib, jb, descb)
  import
  character, intent(in) :: side, uplo, transa, diag
  integer, intent(in) :: mm, nn
  $3, intent(in) :: alpha
  integer, intent(in) :: desca(DLEN_)
  $3, intent(in) :: aa(desca(LLD_), *)
  integer, intent(in) :: ia, ja
  integer, intent(in) :: descb(DLEN_)
  $3, intent(inout) :: bb(descb(LLD_), *)
  integer, intent(in) :: ib, jb
end subroutine p$2trsm
')
