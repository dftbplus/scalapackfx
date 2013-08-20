include(common.m4)

dnl ************************************************************************
dnl *** cpl2g
dnl ************************************************************************

define(`_subroutine_cpl2g',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
!> Copies the content of a local matrix to a global one ($1).
!!
!! \param mygrid BLACS descriptor
!! \param loc  Local matrix.
!! \param desc  Descriptor of the global matrix.
!! \param ii  Starting row in the global matrix.
!! \param jj  Starting column in the global matrix
!! \param glob  Local part of the global matrix.
!!
subroutine cpl2g_$1(mygrid, loc, desc, ii, jj, glob)
  type(blacsgrid), intent(in) :: mygrid
  $2, intent(in) :: loc(:,:)
  integer, intent(in) :: desc(DLEN_)
  integer, intent(in) :: ii, jj
  $2, intent(inout) :: glob(:,:)

  integer :: i2, j2, iloc, jloc, prow, pcol

  do j2 = 1, size(loc, dim=2)
    do i2 = 1, size(loc, dim=1)
      call scalafx_infog2l(mygrid, desc, i2 + ii - 1, j2 + jj - 1, &
          & iloc, jloc, prow, pcol)
      if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
        glob(iloc, jloc) = loc(i2, j2)
      end if
    end do
  end do
  
end subroutine cpl2g_$1
')

dnl ************************************************************************
dnl *** addl2g
dnl ************************************************************************

define(`_subroutine_addl2g',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
!> Adds the content of a local matrix to a global one ($1).
!!
!! \param mygrid BLACS descriptor
!! \param loc  Local matrix.
!! \param desc  Descriptor of the global matrix.
!! \param ii  Starting row in the global matrix.
!! \param jj  Starting column in the global matrix
!! \param glob  Local part of the global matrix.
!!
subroutine addl2g_$1(mygrid, loc, desc, ii, jj, glob)
  type(blacsgrid), intent(in) :: mygrid
  $2, intent(in) :: loc(:,:)
  integer, intent(in) :: desc(DLEN_)
  integer, intent(in) :: ii, jj
  $2, intent(inout) :: glob(:,:)

  integer :: i2, j2, iloc, jloc, prow, pcol

  do j2 = 1, size(loc, dim=2)
    do i2 = 1, size(loc, dim=1)
      call scalafx_infog2l(mygrid, desc, i2 + ii - 1, j2 + jj - 1, &
          & iloc, jloc, prow, pcol)
      if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
        glob(iloc, jloc) = glob(iloc, jloc) + loc(i2, j2)
      end if
    end do
  end do
  
end subroutine addl2g_$1
')

dnl ************************************************************************
dnl *** cpg2l
dnl ************************************************************************

define(`_subroutine_cpg2l',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
!> Copies the content from the global matrix into a local one.
!!
!! \param mygrid BLACS descriptor
!! \param desc  Descriptor of the global matrix.
!! \param ii  Starting row in the global matrix.
!! \param jj  Starting column in the global matrix
!! \param glob  Local part of the global matrix.
!! \param loc  Local matrix.
!!
subroutine cpg2l_$1(mygrid, desc, ii, jj, glob, loc)
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: desc(DLEN_)
  integer, intent(in) :: ii, jj
  $2, intent(in) :: glob(:,:)
  $2, intent(out) :: loc(:,:)
  
  integer :: i2, j2, iloc, jloc, prow, pcol

  loc(:,:) = 0.0_dp
  do j2 = 1, size(loc, dim=2)
    do i2 = 1, size(loc, dim=1)
      call scalafx_infog2l(mygrid, desc, i2 + ii - 1, j2 + jj - 1, &
          & iloc, jloc, prow, pcol)
      if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
        loc(i2, j2) = glob(iloc, jloc)
      end if
    end do
  end do
  
end subroutine cpg2l_$1
')

dnl ************************************************************************
dnl *** addg2l
dnl ************************************************************************

define(`_subroutine_addg2l',`
dnl $1 subroutine suffix
dnl $2 dummy argument type
!> Copies the content from the global matrix into a local one.
!!
!! \param mygrid BLACS descriptor
!! \param desc  Descriptor of the global matrix.
!! \param ii  Starting row in the global matrix.
!! \param jj  Starting column in the global matrix
!! \param glob  Local part of the global matrix.
!! \param loc  Local matrix.
!!
subroutine addg2l_$1(mygrid, desc, ii, jj, glob, loc)
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: desc(DLEN_)
  integer, intent(in) :: ii, jj
  $2, intent(in) :: glob(:,:)
  $2, intent(out) :: loc(:,:)
  
  integer :: i2, j2, iloc, jloc, prow, pcol

  loc(:,:) = 0.0_dp
  do j2 = 1, size(loc, dim=2)
    do i2 = 1, size(loc, dim=1)
      call scalafx_infog2l(mygrid, desc, i2 + ii - 1, j2 + jj - 1, &
          & iloc, jloc, prow, pcol)
      if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
        loc(i2, j2) = loc(i2, j2) + glob(iloc, jloc)
      end if
    end do
  end do
  
end subroutine addg2l_$1
')

dnl ************************************************************************
dnl *** writearray_master
dnl ************************************************************************

define(`_subroutine_writearray_master',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl
!> Writes a distributed array to a file (master, $1).
!!
!! \param mygrid  BLACS descriptor
!! \param fd  File descriptor of an opened file.
!! \param desc  Descriptor of the distributed matrix
!! \param mtxloc  Local portion of the distributed matrix
!! \param rowwise  If .true. matrix is dumped rowwise otherwise columnwise
!! \param elemformat  Formatting of one element (e.g. "(E23.15)"). If present
!!     matrix will be written formatted, otherwise the matrix is written
!!     unformatted. The file descriptor must accordingly belong to a formatted
!!     or an unformatted file! The formatting string must contain the
!!     delimiting parantheses.
!!
subroutine writearray_master_$1(mygrid, fd, desc, mtxloc, rowwise, elemformat)
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: fd, desc(DLEN_)
  $2, intent(in) :: mtxloc(:,:)
  logical, intent(in), optional :: rowwise
  character(*), intent(in), optional :: elemformat

  type(linecomm) :: distributor
  $2, allocatable :: buffer(:)
  character(:), allocatable :: lineformat
  logical :: rowwise0, formatted
  integer :: nline, linelen, ii, ndigit, nn

  _handle_inoptflag(rowwise0, rowwise, .false.)
  if (rowwise0) then
    call distributor%init(mygrid, desc, "r")
    nline = desc(N_)
    linelen = desc(M_)
  else
    call distributor%init(mygrid, desc, "c")
    nline = desc(M_)
    linelen = desc(N_)
  end if
  formatted = present(elemformat)
  if (formatted) then
    ndigit = floor(log(real(linelen, dp)) / log(10.0_dp)) + 1
    nn = ndigit + len_trim(elemformat) + 2
    allocate(character(nn) :: lineformat)
    write(lineformat, "(A,I0,A,A)") "(", linelen, trim(elemformat), ")"
  end if

  allocate(buffer(linelen))
  do ii = 1, nline
    call distributor%getline_master(mygrid, ii, mtxloc, buffer)
    if (formatted) then
      write(fd, lineformat) buffer
    else
      write(fd) buffer
    end if
  end do
  
end subroutine writearray_master_$1
')


dnl ************************************************************************
dnl *** writearray_slave
dnl ************************************************************************

define(`_subroutine_writearray_slave',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl
!> Writes a distributed array to a file (slave, $1).
!!
!! \param mygrid  BLACS descriptor
!! \param desc  Descriptor of the distributed matrix
!! \param mtxloc  Local portion of the distributed matrix
!! \param rowwise  If .true. matrix is dumped rowwise otherwise columnwise
!!
subroutine writearray_slave_$1(mygrid, desc, mtxloc, rowwise)
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: desc(DLEN_)
  $2, intent(in) :: mtxloc(:,:)
  logical, intent(in), optional :: rowwise

  type(linecomm) :: distributor
  integer :: ii, nline
  logical :: rowwise0

  _handle_inoptflag(rowwise0, rowwise, .false.)
  if (rowwise0) then
    call distributor%init(mygrid, desc, "r")
    nline = desc(N_)
  else
    call distributor%init(mygrid, desc, "c")
    nline = desc(M_)
  end if
  do ii = 1, nline
    call distributor%getline_slave(mygrid, ii, mtxloc)
  end do
  
end subroutine writearray_slave_$1
')

dnl ************************************************************************
dnl *** readarray_master
dnl ************************************************************************

define(`_subroutine_readarray_master',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl
!> Reads a distributed array from a file (master, $1).
!!
!! \param mygrid  BLACS descriptor
!! \param fd  File descriptor of an opened file.
!! \param desc  Descriptor of the distributed matrix
!! \param mtxloc  Local portion of the distributed matrix
!! \param rowwise  If .true. matrix is assumed to be stored rowwise otherwise
!!     columnwise (default: .false.)
!! \param formatted  If .true. matrix will be read formatted otherwise
!!     unformatted. The file descriptor must accordingly belong to a formatted
!!     or an unformatted file! (default: .false.)
!!
subroutine readarray_master_$1(mygrid, fd, desc, mtxloc, rowwise, formatted)
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: fd, desc(DLEN_)
  $2, intent(inout) :: mtxloc(:,:)
  logical, intent(in), optional :: rowwise, formatted

  type(linecomm) :: collector
  $2, allocatable :: buffer(:)
  logical :: rowwise0, formatted0
  integer :: nline, linelen, ii

  _handle_inoptflag(rowwise0, rowwise, .false.)
  _handle_inoptflag(formatted0, formatted, .false.)
  if (rowwise0) then
    call collector%init(mygrid, desc, "r")
    nline = desc(N_)
    linelen = desc(M_)
  else
    call collector%init(mygrid, desc, "c")
    nline = desc(M_)
    linelen = desc(N_)
  end if
  allocate(buffer(linelen))
  do ii = 1, nline
    if (formatted0) then
      read(fd, *) buffer
    else
      read(fd) buffer
    end if
    call collector%setline_master(mygrid, ii, buffer, mtxloc)
  end do
  
end subroutine readarray_master_$1
')

dnl ************************************************************************
dnl *** readarray_slave
dnl ************************************************************************

define(`_subroutine_readarray_slave',`
dnl
dnl $1 subroutine suffix
dnl $2 dummy argument type
dnl
!> Reads a distributed array from a file (slave, $1).
!!
!! \param mygrid  BLACS descriptor
!! \param desc  Descriptor of the distributed matrix
!! \param mtxloc  Local portion of the distributed matrix
!! \param rowwise  If .true. matrix is assumed to be stored rowwise otherwise
!!     columnwise (default: .false.)
!!
subroutine readarray_slave_$1(mygrid, desc, mtxloc, rowwise)
  type(blacsgrid), intent(in) :: mygrid
  integer, intent(in) :: desc(DLEN_)
  $2, intent(inout) :: mtxloc(:,:)
  logical, intent(in), optional :: rowwise

  type(linecomm) :: collector
  integer :: ii, nline
  logical :: rowwise0

  _handle_inoptflag(rowwise0, rowwise, .false.)
  if (rowwise0) then
    call collector%init(mygrid, desc, "r")
    nline = desc(N_)
  else
    call collector%init(mygrid, desc, "c")
    nline = desc(M_)
  end if
  do ii = 1, nline
    call collector%setline_slave(mygrid, ii, mtxloc)
  end do
  
end subroutine readarray_slave_$1
')
