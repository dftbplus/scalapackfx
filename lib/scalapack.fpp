#:include 'scalapackfx.fypp'
#:set TYPES = FLOAT_TYPES

!> Wrapper functions for scalapack.
module scalapack_module
  use scalapackfx_common_module
  implicit none
  private

  public :: psygst, phegst, psyev, pheev, psyevd, pheevd, psyevr, pheevr
  public :: ptrsm, ppotrf, ppotri, ptrtri, pgetrf, pgesvd
  public :: sl_init, numroc, infog2l, indxl2g, descinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ppotrf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_ppotrf_template(TYPEABBREV, FTYPES)
  !> Cholesky factorization of symmetric/Hermitian pos.def. matrix (${TYPE}$).
  subroutine p${TYPEABBREV}$potrf(uplo, nn, aa, ia, ja, desca, info)
    import
    character, intent(in) :: uplo
    integer, intent(in) :: nn
    integer, intent(in) :: ia, ja, desca(DLEN_)
    ${FTYPES}$, intent(inout) :: aa(desca(LLD_), *)
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$potrf
#:enddef interface_ppotrf_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ppotri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_ppotri_template(TYPEABBREV, FTYPES)
  !> Inversion of a Cholesky decomposed symmetric/Hermitian matrix (${TYPE}$).
  subroutine p${TYPEABBREV}$potri(uplo, nn, aa, ia, ja, desca, info)
    import
    character, intent(in) :: uplo
    integer, intent(in) :: nn
    integer, intent(in) :: ia, ja, desca(DLEN_)
    ${FTYPES}$, intent(inout) :: aa(desca(LLD_), *)
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$potri
#:enddef interface_ppotri_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ptrtri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_ptrtri_template(TYPEABBREV, FTYPES)
  !> Inversion of a Cholesky decomposed symmetric/Hermitian matrix (${TYPE}$).
  subroutine p${TYPEABBREV}$trtri(uplo, diag, nn, aa, ia, ja, desca, info)
    import
    character, intent(in) :: uplo, diag
    integer, intent(in) :: nn
    integer, intent(in) :: ia, ja, desca(DLEN_)
    ${FTYPE}$, intent(inout) :: aa(desca(LLD_), *)
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$trtri
#:enddef interface_ptrtri_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! pgetrf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_pgetrf_template(TYPEABBREV, FTYPES)
  !> LU factorization of a general matrix with pivoting (${TYPE}$).
  subroutine p${TYPEABBREV}$getrf(mm, nn, aa, ia, ja, desca, ipiv, info)
    import
    integer, intent(in) :: mm
    integer, intent(in) :: nn
    integer, intent(in) :: ia, ja, desca(DLEN_)
    ${FTYPES}$, intent(inout) :: aa(desca(LLD_), *)
    integer, intent(out) :: ipiv(*)
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$getrf
#:enddef interface_pgetrf_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! psygst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




#:def interface_psygst_template(TYPEABBREV, KIND)
  !> Reduces generalized symmetric eigenvalue problem to standard form (${TYPE}$).
  subroutine p${TYPEABBREV}$sygst(ibtype, uplo, nn, aa, ia, ja, desca, bb, ib,&
    & jb, descb, scale, info)
    import
    integer, intent(in) :: ibtype
    character, intent(in) :: uplo
    integer, intent(in) :: nn
    integer, intent(in) :: ia, ja, desca(DLEN_)
    real(${KIND}$), intent(inout) :: aa(desca(LLD_), *)
    integer, intent(in) :: ib, jb, descb(DLEN_)
    real(${KIND}$), intent(in) :: bb(descb(LLD_), *)
    real(${KIND}$), intent(out) :: scale
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$sygst
#:enddef interface_psygst_template


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! phegst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_phegst_template(TYPEABBREV, KIND)
  !> Reduces generalized Hermitian eigenvalue problem to standard form (${TYPE}$).
  subroutine p${TYPEABBREV}$hegst(ibtype, uplo, nn, aa, ia, ja, desca, bb, ib,&
    & jb, descb, scale, info)
    import
    integer, intent(in) :: ibtype
    character, intent(in) :: uplo
    integer, intent(in) :: nn
    integer, intent(in) :: ia, ja, desca(DLEN_)
    complex(${KIND}$), intent(inout) :: aa(desca(LLD_), *)
    integer, intent(in) :: ib, jb, descb(DLEN_)
    complex(${KIND}$), intent(in) :: bb(descb(LLD_), *)
    real(${KIND}$), intent(out) :: scale
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$hegst
#:enddef interface_phegst_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! psyev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_psyev_template(TYPEABBREV, KIND)
  !> Eigenvalues and eigenvectors by divide and conquer algorithm (${TYPE}$)
  subroutine p${TYPEABBREV}$syev(jobz, uplo, nn, aa, ia, ja, desca, ww, zz,&
    & iz, jz, descz, work, lwork, info)
    import
    character, intent(in) :: jobz, uplo
    integer, intent(in) :: nn
    integer, intent(in) :: ia, ja, desca(DLEN_)
    real(${KIND}$), intent(inout) :: aa(desca(LLD_), *)
    integer, intent(in) :: iz, jz, descz(DLEN_)
    real(${KIND}$), intent(out) :: ww(nn), zz(descz(LLD_),*)
    real(${KIND}$), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$syev
#:enddef interface_psyev_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! pheev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_pheev_template(TYPEABBREV, KIND)
  !> Eigenvalues and eigenvectors by divide and conquer algorithm (${TYPE}$)
  subroutine p${TYPEABBREV}$heev(jobz, uplo, nn, aa, ia, ja, desca, ww, zz, iz, jz,&
      & descz, work, lwork, rwork, lrwork, info)
    import
    character, intent(in) :: jobz, uplo
    integer, intent(in) :: nn
    integer, intent(in) :: ia, ja, desca(DLEN_)
    complex(${KIND}$), intent(inout) :: aa(desca(LLD_), *)
    integer, intent(in) :: iz, jz, descz(DLEN_)
    real(${KIND}$), intent(out) :: ww(nn)
    complex(${KIND}$), intent(out) ::  zz(descz(LLD_),*)
    complex(${KIND}$), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    real(${KIND}$), intent(inout) :: rwork(*)
    integer, intent(in) :: lrwork
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$heev
#:enddef interface_pheev_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! psyevd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_psyevd_template(TYPEABBREV, KIND)
  !> Eigenvalues and eigenvectors by divide and conquer algorithm (${TYPE}$)
  subroutine p${TYPEABBREV}$syevd(jobz, uplo, nn, aa, ia, ja, desca, ww, zz, iz, jz,&
      & descz, work, lwork, iwork, liwork, info)
    import
    character, intent(in) :: jobz, uplo
    integer, intent(in) :: nn
    integer, intent(in) :: ia, ja, desca(DLEN_)
    real(${KIND}$), intent(inout) :: aa(desca(LLD_), *)
    integer, intent(in) :: iz, jz, descz(DLEN_)
    real(${KIND}$), intent(out) :: ww(nn), zz(descz(LLD_),*)
    real(${KIND}$), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    integer, intent(inout) :: iwork(*)
    integer, intent(in) :: liwork
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$syevd
#:enddef interface_psyevd_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! pheevd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_pheevd_template(TYPEABBREV, KIND)
  !> Eigenvalues and eigenvectors by divide and conquer algorithm (${TYPE}$)
  subroutine p${TYPEABBREV}$heevd(jobz, uplo, nn, aa, ia, ja, desca, ww, zz, iz, jz,&
      & descz, work, lwork, rwork, lrwork, iwork, liwork, info)
    import
    character, intent(in) :: jobz, uplo
    integer, intent(in) :: nn
    integer, intent(in) :: ia, ja, desca(DLEN_)
    complex(${KIND}$), intent(inout) :: aa(desca(LLD_), *)
    integer, intent(in) :: iz, jz, descz(DLEN_)
    real(${KIND}$), intent(out) :: ww(nn)
    complex(${KIND}$), intent(out) ::  zz(descz(LLD_),*)
    complex(${KIND}$), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    real(${KIND}$), intent(inout) :: rwork(*)
    integer, intent(in) :: lrwork
    integer, intent(inout) :: iwork(*)
    integer, intent(in) :: liwork
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$heevd
#:enddef interface_pheevd_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! psyevr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_psyevr_template(TYPEABBREV, KIND)
  !> Eigenvalues and eigenvectors by MRRR algorithm (${TYPE}$)
  subroutine p${TYPEABBREV}$syevr(jobz, range, uplo, nn, aa, ia, ja, desca, vl, vu,&
    & il, iu, mm, nz, ww, zz, iz, jz, descz, work, lwork, iwork, liwork, info)
    import
    character, intent(in) :: jobz, range, uplo
    integer, intent(in) :: nn
    integer, intent(in) :: desca(DLEN_)
    real(${KIND}$), intent(inout) :: aa(desca(LLD_), *)
    integer, intent(in) :: ia, ja
    real(${KIND}$), intent(in) :: vl, vu
    integer, intent(in) :: il, iu
    integer, intent(out) :: mm, nz
    real(${KIND}$), intent(out) :: ww(nn)
    integer, intent(in) :: descz(DLEN_)
    real(${KIND}$), intent(out) :: zz(descz(LLD_),*)
    integer, intent(in) :: iz, jz
    real(${KIND}$), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    integer, intent(inout) :: iwork(*)
    integer, intent(in) :: liwork
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$syevr
#:enddef interface_psyevr_template


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! pheevr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_pheevr_template(TYPEABBREV, KIND)
  !> Eigenvalues and eigenvectors by MRRR algorithm (${TYPE}$)
  subroutine p${TYPEABBREV}$heevr(jobz, range, uplo, nn, aa, ia, ja, desca, vl,&
    & vu, il, iu, mm, nz, ww, zz, iz, jz, descz, work, lwork, rwork, lrwork,&
    & iwork, liwork, info)
    import
    character, intent(in) :: jobz, range, uplo
    integer, intent(in) :: nn
    integer, intent(in) :: desca(DLEN_)
    complex(${KIND}$), intent(inout) :: aa(desca(LLD_), *)
    integer, intent(in) :: ia, ja
    real(${KIND}$), intent(in) :: vl, vu
    integer, intent(in) :: il, iu
    integer, intent(out) :: mm, nz
    real(${KIND}$), intent(out) :: ww(nn)
    integer, intent(in) :: descz(DLEN_)
    complex(${KIND}$), intent(out) ::  zz(descz(LLD_),*)
    integer, intent(in) :: iz, jz
    complex(${KIND}$), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    real(${KIND}$), intent(inout) :: rwork(*)
    integer, intent(in) :: lrwork
    integer, intent(inout) :: iwork(*)
    integer, intent(in) :: liwork
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$heevr
#:enddef interface_pheevr_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! prgesvd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_prgesvd_template(TYPEABBREV, KIND)
  !> Singular values and vectors (${TYPE}$)
  subroutine p${TYPEABBREV}$gesvd(jobu, jobvt, mm, nn, aa, ia, ja, desca, sigma,&
    & uu, iu, ju, descu, vt, ivt, jvt, descvt, work, lwork, info)
    import
    character, intent(in) :: jobu, jobvt
    integer, intent(in) :: mm, nn
    integer, intent(in) :: ia, ja, desca(DLEN_)
    real(${KIND}$), intent(inout) :: aa(desca(LLD_), *)
    real(${KIND}$), intent(out) :: sigma(*)
    integer, intent(in) :: iu, ju, descu(DLEN_)
    real(${KIND}$), intent(out) :: uu(descu(LLD_), *)
    integer, intent(in) :: ivt, jvt, descvt(DLEN_)
    real(${KIND}$), intent(out) :: vt(descvt(LLD_), *)
    real(${KIND}$), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$gesvd
#:enddef interface_prgesvd_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! pcgesvd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_pcgesvd_template(TYPEABBREV, KIND)
  !> Singular values and vectors (${TYPE}$)
  subroutine p${TYPEABBREV}$gesvd(jobu, jobvt, mm, nn, aa, ia, ja, desca, sigma,&
    & uu, iu, ju, descu, vt, ivt, jvt, descvt, work, lwork, rwork, info)
    import
    character, intent(in) :: jobu, jobvt
    integer, intent(in) :: mm, nn
    integer, intent(in) :: ia, ja, desca(DLEN_)
    complex(${KIND}$), intent(inout) :: aa(desca(LLD_), *)
    real(${KIND}$), intent(out) :: sigma(*)
    integer, intent(in) :: iu, ju, descu(DLEN_)
    complex(${KIND}$), intent(out) :: uu(descu(LLD_), *)
    integer, intent(in) :: ivt, jvt, descvt(DLEN_)
    complex(${KIND}$), intent(out) :: vt(descvt(LLD_), *)
    complex(${KIND}$), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    real(${KIND}$), intent(inout) :: rwork(*)
    integer, intent(out) :: info
  end subroutine p${TYPEABBREV}$gesvd
#:enddef interface_pcgesvd_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ptrsm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def interface_ptrsm_template()
  !> Solves a triangular matrix equation (${TYPE}$).
  subroutine p${TYPEABBREV}$trsm(side, uplo, transa, diag, mm, nn, alpha, aa, ia, ja,&
      & desca, bb, ib, jb, descb)
    import
    character, intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: mm, nn
    ${FTYPE}$, intent(in) :: alpha
    integer, intent(in) :: desca(DLEN_)
    ${FTYPE}$, intent(in) :: aa(desca(LLD_), *)
    integer, intent(in) :: ia, ja
    integer, intent(in) :: descb(DLEN_)
    ${FTYPE}$, intent(inout) :: bb(descb(LLD_), *)
    integer, intent(in) :: ib, jb
  end subroutine p${TYPEABBREV}$trsm
#:enddef interface_ptrsm_template

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCALAPACK CORE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Cholesky factorization of symmetric/Hermitian positive definite matrix.
  interface ppotrf
    #:for TYPE in TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set FTYPE = FORTRAN_TYPES[TYPE]
      $:interface_ppotrf_template(TYPEABBREV, FTYPE)
    #:endfor
  end interface ppotrf

  !> Inversion of a Cholesky decomposed symmetric/Hermitian matrix.
  interface ppotri
    #:for TYPE in TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set FTYPE = FORTRAN_TYPES[TYPE]
      $:interface_ppotri_template(TYPEABBREV, FTYPE)
    #:endfor
  end interface ppotri

  !> Inversion of a triangular matrix.
  interface ptrtri
    #:for TYPE in TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set FTYPE = FORTRAN_TYPES[TYPE]
      $:interface_ptrtri_template(TYPEABBREV, FTYPE)
    #:endfor
  end interface ptrtri

  !> LU decomposition of a general matrix with pivoting
  interface pgetrf
    #:for TYPE in TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set FTYPE = FORTRAN_TYPES[TYPE]
      $:interface_pgetrf_template(TYPEABBREV, FTYPE)
    #:endfor
  end interface pgetrf

  !> Reduces generalized symmetric eigenvalue problem to standard form.
  interface psygst
    #:for TYPE in REAL_TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set KIND = FORTRAN_KINDS[TYPE]
      $:interface_psygst_template(TYPEABBREV, KIND)
    #:endfor
  end interface psygst

  !> Reduces generalized Hermitian eigenvalue problem to standard form.
  interface phegst
    #:for TYPE in COMPLEX_TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set KIND = FORTRAN_KINDS[TYPE]
      $:interface_phegst_template(TYPEABBREV, KIND)
    #:endfor
  end interface phegst

  !> Solves the symmetric eigenvalue problem.
  interface psyev
    #:for TYPE in REAL_TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set KIND = FORTRAN_KINDS[TYPE]
      $:interface_psyev_template(TYPEABBREV, KIND)
    #:endfor
  end interface psyev

  !> Solves the Hermitian eigenvalue problem.
  interface pheev
    #:for TYPE in COMPLEX_TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set KIND = FORTRAN_KINDS[TYPE]
      $:interface_pheev_template(TYPEABBREV, KIND)
    #:endfor
  end interface pheev

  !> Solves the symmetric eigenvalue problem by divide and conquer algorithm.
  interface psyevd
    #:for TYPE in REAL_TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set KIND = FORTRAN_KINDS[TYPE]
      $:interface_psyevd_template(TYPEABBREV, KIND)
    #:endfor
  end interface psyevd

  !> Solves the Hermitian eigenvalue problem by divide and conquer algorithm.
  interface pheevd
    #:for TYPE in COMPLEX_TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set KIND = FORTRAN_KINDS[TYPE]
      $:interface_pheevd_template(TYPEABBREV, KIND)
    #:endfor
  end interface pheevd

  !> Solves the symmetric eigenvalue problem by the MRRR algorithm.
  interface psyevr
    #:for TYPE in REAL_TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set KIND = FORTRAN_KINDS[TYPE]
      $:interface_psyevr_template(TYPEABBREV, KIND)
    #:endfor
  end interface psyevr

  !> Solves the Hermitian eigenvalue problem by the MRRR algorithm.
  interface pheevr
    #:for TYPE in COMPLEX_TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set KIND = FORTRAN_KINDS[TYPE]
      $:interface_pheevr_template(TYPEABBREV, KIND)
    #:endfor
  end interface pheevr

  !> Singular value decomposition of a matrix
  interface pgesvd
    #:for TYPE in REAL_TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set KIND = FORTRAN_KINDS[TYPE]
      $:interface_prgesvd_template(TYPEABBREV, KIND)
    #:endfor
    #:for TYPE in COMPLEX_TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set KIND = FORTRAN_KINDS[TYPE]
      $:interface_pcgesvd_template(TYPEABBREV, KIND)
    #:endfor
  end interface pgesvd

  !> Linear system of equation for triangular matrix.
  interface ptrsm
    #:for TYPE in TYPES
      #:set TYPEABBREV = TYPE_ABBREVS[TYPE]
      #:set FTYPE = FORTRAN_TYPES[TYPE]
      $:interface_ptrsm_template()
    #:endfor
  end interface ptrsm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCALAPACK TOOLS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface

    !> Scalapack initialization routine.
    subroutine sl_init(ictxt, nprow, npcol)
      integer, intent(out) :: ictxt
      integer, intent(in) :: nprow, npcol
    end subroutine sl_init

    !> Number of rows or columns of distributed matrix owned by local process.
    function numroc(nn, nb, iproc, isrcproc, nproc)
      integer, intent(in) :: nn, nb, iproc, isrcproc, nproc
      integer :: numroc
    end function numroc

    !> Converts global matrix index into local.
    subroutine infog2l(grindx, gcindx, desc, nprow, npcol, myrow, mycol,&
        & lrindx, lcindx, rsrc, csrc)
      import DLEN_
      integer, intent(in) :: grindx, gcindx, desc(DLEN_)
      integer, intent(in) :: nprow, npcol, myrow, mycol
      integer, intent(out) :: lrindx, lcindx, rsrc, csrc
    end subroutine infog2l

    !> Converts local matrix index into global.
    function indxl2g(indxglob, nb, iproc, isrcproc, nprocs)
      integer :: indxl2g
      integer, intent(in) :: indxglob, nb, iproc, isrcproc, nprocs
    end function indxl2g

    !> Initializes a descriptor for a distributed array.
    subroutine descinit(desc, mm, nn, mb, nb, irsrc, icsrc, ictxt, lld, info)
      import DLEN_
      integer, intent(out) :: desc(DLEN_)
      integer, intent(in) :: mm, nn, mb, nb, irsrc, icsrc, ictxt, lld
      integer, intent(out) :: info
    end subroutine descinit

  end interface


end module scalapack_module

