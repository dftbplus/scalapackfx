include(common.m4)

dnl ************************************************************************
dnl *** getblock
dnl ************************************************************************

define(`_subroutine_getblock_master',`
dnl
dnl $1 Subroutine suffix
dnl $2 Matrix element type.
dnl
!> Returns the given block of the distributed matrix (master, $1).
!!
!! \param self  Instance.
!! \param mygrid  BLACS descriptor.
!! \param ii  Row/Column index.
!! \param ib  Block index within given row/column.
!! \param locmtx  Local part of the global matrix.
!! \param buffer  Contains the given piece of the distributed matrix on exit.
!!    Its size should be greater than or equal to the BLACS block size
!!    along that dimension. The actual number of elements for a given block
!!    can be queried via the getblocksize() call.
!!
subroutine getblock_master_$1(self, mygrid, ii, ib, locmtx, &
    & buffer)
  class(linecomm), intent(in) :: self
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: ii, ib
  $2, intent(in) :: locmtx(:,:)
  $2, target, intent(out) :: buffer(:)

  integer :: prow, pcol, lrow, lcol, nrow, ncol
  $2, pointer :: work(:,:)

  call self%getpositions(mygrid, ii, ib, prow, pcol, lrow, lcol, nrow, &
      & ncol)
  work(1:nrow,1:ncol) => buffer(1:nrow*ncol)
  if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
    work = locmtx(lrow:lrow+nrow-1,lcol:lcol+ncol-1)
  else
    call blacsfx_gerv(mygrid, work, prow, pcol)
  end if
  
end subroutine getblock_master_$1
')

define(`_subroutine_getblock_slave',`
dnl
dnl $1 Subroutine suffix
dnl $2 Matrix element type.
dnl
!> Returns the given block of the distributed matrix (slave, $1).
!!
!! \param self  Instance.
!! \param mygrid  BLACS descriptor.
!! \param ii  Row/Column index.
!! \param ib  Block index within given row/column.
!! \param locmtx  Local part of the global matrix.
!!
subroutine getblock_slave_$1(self, mygrid, ii, ib, locmtx)
  class(linecomm), intent(in) :: self
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: ii, ib
  $2, target, intent(in) :: locmtx(:,:)

  integer :: prow, pcol, lrow, lcol, nrow, ncol
  $2, pointer :: work(:,:)

  call self%getpositions(mygrid, ii, ib, prow, pcol, lrow, lcol, nrow, &
      & ncol)
  if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
    work => locmtx(lrow:lrow+nrow-1, lcol:lcol+ncol-1)
    call blacsfx_gesd(mygrid, work, self%iorow, self%iocol)
  end if
  
end subroutine getblock_slave_$1
')


dnl ************************************************************************
dnl *** getline
dnl ************************************************************************

define(`_subroutine_getline_master',`
dnl
dnl $1 Subroutine suffix
dnl $2 Matrix element type.
dnl
!> Returns an entire row/column of a distributed matrix (master, $1)
!!
!! \param self  Instance.
!! \param mygrid  BLACS descriptor
!! \param ii  Number of the line (row or column) to collect.
!! \param locmtx  Local part of the global matrix.
!! \param buffer  Contains the collected line on exit. Its size should be
!!     big enough to contain the result (greater or equal to the size of
!!     the distributed matrix along that direction).
!!
subroutine getline_master_$1(self, mygrid, ii, locmtx, buffer)
  class(linecomm), intent(in) :: self
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: ii
  $2, intent(in) :: locmtx(:,:)
  $2, intent(out) :: buffer(:)

  integer :: ib, istart, iend

  iend = 0
  do ib = 1, self%nblock
    istart = iend + 1
    iend = istart + self%getblocksize(ib) - 1
    call self%getblock_master(mygrid, ii, ib, locmtx, buffer(istart:iend))
  end do
    
end subroutine getline_master_$1
')

define(`_subroutine_getline_slave',`
dnl
dnl $1 Subroutine suffix
dnl $2 Matrix element type.
dnl
!> Returns the entire row/column of a distributed matrix (slave)
!!
!! \param mygrid  BLACS descriptor
!! \param ii  Number of the line (row or column) to collect.
!! \param locmtx  Local part of the global matrix.
!!
subroutine getline_slave_$1(self, mygrid, ii, locmtx)
  class(linecomm), intent(in) :: self
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: ii
  $2, intent(in) :: locmtx(:,:)

  integer :: ib

  do ib = 1, self%nblock
    call self%getblock_slave(mygrid, ii, ib, locmtx)
  end do
    
end subroutine getline_slave_$1
')


dnl ************************************************************************
dnl *** setblock
dnl ************************************************************************

define(`_subroutine_setblock_master',`
dnl
dnl $1 Subroutine suffix
dnl $2 Matrix element type.
dnl
!> Sets the given block of the distributed matrix (master, $1).
!!
!! \param self  Instance.
!! \param mygrid  BLACS descriptor.
!! \param ii  Row/Column index.
!! \param ib  Block index within given row/column.
!! \param buffer  Contains the given piece to be distributed. It should contain
!!    the appropriate number of elements for the given block, as returned by
!!    the getblocksize() call. Its size can be bigger than that.
!! \param locmtx  Local part of the global matrix.
!!
subroutine setblock_master_$1(self, mygrid, ii, ib, buffer, locmtx)
  class(linecomm), intent(in) :: self
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: ii, ib
  $2, target, intent(in) :: buffer(:)
  $2, intent(inout) :: locmtx(:,:)

  integer :: prow, pcol, lrow, lcol, nrow, ncol
  $2, pointer :: work(:,:)

  call self%getpositions(mygrid, ii, ib, prow, pcol, lrow, lcol, nrow, &
      & ncol)
  work(1:nrow,1:ncol) => buffer(1:nrow*ncol)
  if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
    locmtx(lrow:lrow+nrow-1,lcol:lcol+ncol-1) = work
  else
    call blacsfx_gesd(mygrid, work, prow, pcol)
  end if
  
end subroutine setblock_master_$1
')

define(`_subroutine_setblock_slave',`
dnl
dnl $1 Subroutine suffix
dnl $2 Matrix element type.
dnl
!> Sets the given block of the distributed matrix (slave, $1).
!!
!! \param self  Instance.
!! \param mygrid  BLACS descriptor.
!! \param ii  Row/Column index.
!! \param ib  Block index within given row/column.
!! \param locmtx  Local part of the global matrix.
!!
subroutine setblock_slave_$1(self, mygrid, ii, ib, locmtx)
  use iso_fortran_env
  class(linecomm), intent(in) :: self
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: ii, ib
  $2, intent(inout), target :: locmtx(:,:)

  integer :: prow, pcol, lrow, lcol, nrow, ncol
  $2, pointer :: work(:,:)

  call self%getpositions(mygrid, ii, ib, prow, pcol, lrow, lcol, nrow, &
      & ncol)
  if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
    work => locmtx(lrow:lrow+nrow-1, lcol:lcol+ncol-1)
    call blacsfx_gerv(mygrid, work, self%iorow, self%iocol)
  end if
  
end subroutine setblock_slave_$1
')


dnl ************************************************************************
dnl *** setline
dnl ************************************************************************

define(`_subroutine_setline_master',`
dnl
dnl $1 Subroutine suffix
dnl $2 Matrix element type.
dnl
!> Sets an entire row/column of a distributed matrix (master, $1)
!!
!! \param self  Instance.
!! \param mygrid  BLACS descriptor
!! \param ii  Number of the line (row or column) to set.
!! \param buffer  Contains the line to distribute. It should contain all
!!     elements along the line. Its size can be bigger.
!! \param locmtx  Local part of the global matrix.
!!
subroutine setline_master_$1(self, mygrid, ii, buffer, locmtx)
  class(linecomm), intent(in) :: self
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: ii
  $2, intent(in) :: buffer(:)
  $2, intent(inout) :: locmtx(:,:)

  integer :: ib, istart, iend

  iend = 0
  do ib = 1, self%nblock
    istart = iend + 1
    iend = istart + self%getblocksize(ib) - 1
    call self%setblock_master(mygrid, ii, ib, buffer(istart:iend), locmtx)
  end do
    
end subroutine setline_master_$1
')

define(`_subroutine_setline_slave',`
dnl
dnl $1 Subroutine suffix
dnl $2 Matrix element type.
dnl
!> Sets the entire row/column of a distributed matrix (slave)
!!
!! \param mygrid  BLACS descriptor
!! \param ii  Number of the line (row or column) to set.
!! \param locmtx  Local part of the global matrix.
!!
subroutine setline_slave_$1(self, mygrid, ii, locmtx)
  class(linecomm), intent(in) :: self
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: ii
  $2, intent(inout) :: locmtx(:,:)

  integer :: ib

  do ib = 1, self%nblock
    call self%setblock_slave(mygrid, ii, ib, locmtx)
  end do
    
end subroutine setline_slave_$1
')
