#:include 'scalapackfx.fypp'
#:set TYPES = NUMERIC_TYPES
#:set CONTROLS = CONTROL_TYPES

!> Contains the linecomm type.
module linecomm_module
  use blacsfx_module
  use scalapackfx_common_module
  use scalapackfx_module
  implicit none
  private

  public :: linecomm


  !> Type for communicating a row or a column of a distributed matrix.
  !!
  !! \details The type linecomm collects/distributes a line (row or column)
  !! of a distributed matrix into/from a buffer on the master node.
  !! It communicate the entire line at once or blockwise, with the blocks
  !! having the size of the BLACS block size.
  !!
  !! The code below demonstrates how to write out a distributed matrix
  !! columnwise to disc with the help of linecomm:
  !!
  !!     type(linecomm) :: collector
  !!     real(dp), allocatable :: iobuffer(:)
  !!
  !!     allocate(iobuffer(desc(M_))
  !!     call collector%init(mygrid, desc, "c")
  !!     do icol = 1, desc(N_)
  !!       if (mygrid%master) then
  !!         call collector%getline_master(mygrid, icol, mtxloc, iobuffer)
  !!         write(fd, formatstr) iobuffer(:)
  !!       else
  !!         call collector%getline_slave(mygrid, icol, mtxloc)
  !!       end if
  !!     end do
  !!
  !! Similarly, to read from a file columnwise, you could do the following:
  !!
  !!     type(linecomm) :: distributor
  !!     real(dp), allocatable :: iobuffer(:)
  !!
  !!     allocate(iobuffer(desc(M_))
  !!     call distributor%init(desc, "c")
  !!     do icol = 1, ncol
  !!       if (mygrid%master) then
  !!         read(fd, *) iobuffer(:)
  !!         call distributor%setline_master(mygrid, icol, iobuffer, mtxloc)
  !!       else
  !!         call distributor%setline_slave(mygrid, icol, mtxloc)
  !!       end if
  !!     end do
  !!
  type :: linecomm
    private
    integer :: nn, nblock, blocksize, iorow, iocol
    logical :: rowcollect
    integer :: desc(DLEN_)
  contains
    procedure :: init
    procedure :: getnrblocks
    procedure :: getblocksize

    #:for CONTROL in CONTROLS
      #:for TYPE in TYPES
        procedure :: getblock_${CONTROL}$_${TYPE}$
      #:endfor
    #:endfor

    generic :: getblock_master => getblock_master_int, getblock_master_real, &
        & getblock_master_dreal, getblock_master_complex, &
        & getblock_master_dcomplex
    generic :: getblock_slave => getblock_slave_int, getblock_slave_real, &
        & getblock_slave_dreal, getblock_slave_complex, getblock_slave_dcomplex

    #:for CONTROL in CONTROLS
      #:for TYPE in TYPES
        procedure :: getline_${CONTROL}$_${TYPE}$
      #:endfor
    #:endfor

    generic :: getline_master => getline_master_int, getline_master_real, &
        & getline_master_dreal, getline_master_complex, getline_master_dcomplex
    generic :: getline_slave => getline_slave_int, getline_slave_real, &
        & getline_slave_dreal, getline_slave_complex, getline_slave_dcomplex

    #:for CONTROL in CONTROLS
      #:for TYPE in TYPES
        procedure :: setblock_${CONTROL}$_${TYPE}$
      #:endfor
    #:endfor

    generic :: setblock_master => setblock_master_int, setblock_master_real, &
        & setblock_master_dreal, setblock_master_complex, &
        & setblock_master_dcomplex
    generic :: setblock_slave => setblock_slave_int, setblock_slave_real, &
        & setblock_slave_dreal, setblock_slave_complex, setblock_slave_dcomplex

    #:for CONTROL in CONTROLS
      #:for TYPE in TYPES
        procedure :: setline_${CONTROL}$_${TYPE}$
      #:endfor
    #:endfor

    generic :: setline_master => setline_master_int, setline_master_real, &
        & setline_master_dreal, setline_master_complex, setline_master_dcomplex
    generic :: setline_slave => setline_slave_int, setline_slave_real, &
        & setline_slave_dreal, setline_slave_complex, setline_slave_dcomplex

    procedure, private :: getpositions
  end type linecomm


contains

  !> Initializes a linecomm instance.
  !!
  !! \param self  Initialized instance on exit.
  !! \param desc  Descriptor of distributed matrix.
  !! \param rowcol  If "r" or "R", a given row of the matrix is collected,
  !!     otherwise a given column.
  !! \param iorow  Row of process doing the io (default: row of master).
  !! \param iocol  Column of process doing the io (default: column of master).
  !!
  subroutine init(self, mygrid, desc, rowcol, iorow, iocol)
    class(linecomm), intent(out) :: self
    class(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: desc(DLEN_)
    character, intent(in) :: rowcol
    integer, intent(in), optional :: iorow, iocol

    @:inoptflags(self%iorow, iorow, mygrid%masterrow)
    @:inoptflags(self%iocol, iocol, mygrid%mastercol)
    self%rowcollect = (rowcol == "R" .or. rowcol == "r")
    self%desc(:) = desc
    if (self%rowcollect) then
      self%nn = desc(M_)
      self%blocksize = desc(MB_)
    else
      self%nn = desc(N_)
      self%blocksize = desc(NB_)
    end if
    self%nblock = self%nn / self%blocksize
    if (mod(self%nn, self%blocksize) /= 0) then
      self%nblock = self%nblock + 1
    end if

  end subroutine init


  !> Returns the nr. of blocks along the given row or column.
  !!
  !! \param self  Instance.
  !!
  function getnrblocks(self) result(res)
    class(linecomm), intent(in) :: self
    integer :: res

    res = self%nblock

  end function getnrblocks


  !> Returns the size of a block with the given index.
  !!
  !! \param self Instance.
  !! \param ib  Block index.
  !! \return  Size of the given block.
  !!
  function getblocksize(self, ib) result(res)
    class(linecomm), intent(in) :: self
    integer, intent(in) :: ib
    integer :: res

    if (ib < self%nblock) then
      res = self%blocksize
    else
      res = mod(self%nn, self%blocksize)
      if (res == 0) then
        res = self%blocksize
      end if
    end if

  end function getblocksize


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!! GETBLOCK
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def getblock_master_template(SUFFIX,TYPE)

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
  subroutine getblock_master_${SUFFIX}$(self, mygrid, ii, ib, locmtx, &
      & buffer)
    class(linecomm), intent(in) :: self
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: ii, ib
    ${TYPE}$, intent(in) :: locmtx(:,:)
    ${TYPE}$, target, intent(out) :: buffer(:)

    integer :: prow, pcol, lrow, lcol, nrow, ncol
    ${TYPE}$, pointer :: work(:,:)

    call self%getpositions(mygrid, ii, ib, prow, pcol, lrow, lcol, nrow, &
        & ncol)
    work(1:nrow,1:ncol) => buffer(1:nrow*ncol)
    if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
      work = locmtx(lrow:lrow+nrow-1,lcol:lcol+ncol-1)
    else
      call blacsfx_gerv(mygrid, work, prow, pcol)
    end if

  end subroutine getblock_master_${SUFFIX}$

#:enddef getblock_master_template


#:def getblock_slave_template(SUFFIX,TYPE)

  !> Returns the given block of the distributed matrix (slave, $1).
  !!
  !! \param self  Instance.
  !! \param mygrid  BLACS descriptor.
  !! \param ii  Row/Column index.
  !! \param ib  Block index within given row/column.
  !! \param locmtx  Local part of the global matrix.
  !!
  subroutine getblock_slave_${SUFFIX}$(self, mygrid, ii, ib, locmtx)
    class(linecomm), intent(in) :: self
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: ii, ib
    ${TYPE}$, target, intent(in) :: locmtx(:,:)

    integer :: prow, pcol, lrow, lcol, nrow, ncol
    ${TYPE}$, pointer :: work(:,:)

    call self%getpositions(mygrid, ii, ib, prow, pcol, lrow, lcol, nrow, &
        & ncol)
    if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
      work => locmtx(lrow:lrow+nrow-1, lcol:lcol+ncol-1)
      call blacsfx_gesd(mygrid, work, self%iorow, self%iocol)
    end if

  end subroutine getblock_slave_${SUFFIX}$

#:enddef getblock_slave_template


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!! GETLINE
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def getline_master_template(SUFFIX,TYPE)

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
  subroutine getline_master_${SUFFIX}$(self, mygrid, ii, locmtx, buffer)
    class(linecomm), intent(in) :: self
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: ii
    ${TYPE}$, intent(in) :: locmtx(:,:)
    ${TYPE}$, intent(out) :: buffer(:)

    integer :: ib, istart, iend

    iend = 0
    do ib = 1, self%nblock
      istart = iend + 1
      iend = istart + self%getblocksize(ib) - 1
      call self%getblock_master(mygrid, ii, ib, locmtx, buffer(istart:iend))
    end do

  end subroutine getline_master_${SUFFIX}$

#:enddef getline_master_template

#:def getline_slave_template(SUFFIX,TYPE)

  !> Returns the entire row/column of a distributed matrix (slave)
  !!
  !! \param mygrid  BLACS descriptor
  !! \param ii  Number of the line (row or column) to collect.
  !! \param locmtx  Local part of the global matrix.
  !!
  subroutine getline_slave_${SUFFIX}$(self, mygrid, ii, locmtx)
    class(linecomm), intent(in) :: self
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: ii
    ${TYPE}$, intent(in) :: locmtx(:,:)

    integer :: ib

    do ib = 1, self%nblock
      call self%getblock_slave(mygrid, ii, ib, locmtx)
    end do

  end subroutine getline_slave_${SUFFIX}$

#:enddef getline_slave_template


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!! SETBLOCK
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def setblock_master_template(SUFFIX,TYPE)

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
  subroutine setblock_master_${SUFFIX}$(self, mygrid, ii, ib, buffer, locmtx)
    class(linecomm), intent(in) :: self
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: ii, ib
    ${TYPE}$, target, intent(in) :: buffer(:)
    ${TYPE}$, intent(inout) :: locmtx(:,:)

    integer :: prow, pcol, lrow, lcol, nrow, ncol
    ${TYPE}$, pointer :: work(:,:)

    call self%getpositions(mygrid, ii, ib, prow, pcol, lrow, lcol, nrow, &
        & ncol)
    work(1:nrow,1:ncol) => buffer(1:nrow*ncol)
    if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
      locmtx(lrow:lrow+nrow-1,lcol:lcol+ncol-1) = work
    else
      call blacsfx_gesd(mygrid, work, prow, pcol)
    end if

  end subroutine setblock_master_${SUFFIX}$

#:enddef setblock_master_template

#:def setblock_slave_template(SUFFIX,TYPE)

  !> Sets the given block of the distributed matrix (slave, $1).
  !!
  !! \param self  Instance.
  !! \param mygrid  BLACS descriptor.
  !! \param ii  Row/Column index.
  !! \param ib  Block index within given row/column.
  !! \param locmtx  Local part of the global matrix.
  !!
  subroutine setblock_slave_${SUFFIX}$(self, mygrid, ii, ib, locmtx)
    use iso_fortran_env
    class(linecomm), intent(in) :: self
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: ii, ib
    ${TYPE}$, intent(inout), target :: locmtx(:,:)

    integer :: prow, pcol, lrow, lcol, nrow, ncol
    ${TYPE}$, pointer :: work(:,:)

    call self%getpositions(mygrid, ii, ib, prow, pcol, lrow, lcol, nrow, &
        & ncol)
    if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
      work => locmtx(lrow:lrow+nrow-1, lcol:lcol+ncol-1)
      call blacsfx_gerv(mygrid, work, self%iorow, self%iocol)
    end if

  end subroutine setblock_slave_${SUFFIX}$

#:enddef setblock_slave_template


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!! SETLINE
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:def setline_master_template(SUFFIX,TYPE)

  !> Sets an entire row/column of a distributed matrix (master, $1)
  !!
  !! \param self  Instance.
  !! \param mygrid  BLACS descriptor
  !! \param ii  Number of the line (row or column) to set.
  !! \param buffer  Contains the line to distribute. It should contain all
  !!     elements along the line. Its size can be bigger.
  !! \param locmtx  Local part of the global matrix.
  !!
  subroutine setline_master_${SUFFIX}$(self, mygrid, ii, buffer, locmtx)
    class(linecomm), intent(in) :: self
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: ii
    ${TYPE}$, intent(in) :: buffer(:)
    ${TYPE}$, intent(inout) :: locmtx(:,:)

    integer :: ib, istart, iend

    iend = 0
    do ib = 1, self%nblock
      istart = iend + 1
      iend = istart + self%getblocksize(ib) - 1
      call self%setblock_master(mygrid, ii, ib, buffer(istart:iend), locmtx)
    end do

  end subroutine setline_master_${SUFFIX}$

#:enddef setline_master_template


#:def setline_slave_template(SUFFIX,TYPE)

  !> Sets the entire row/column of a distributed matrix (slave)
  !!
  !! \param mygrid  BLACS descriptor
  !! \param ii  Number of the line (row or column) to set.
  !! \param locmtx  Local part of the global matrix.
  !!
  subroutine setline_slave_${SUFFIX}$(self, mygrid, ii, locmtx)
    class(linecomm), intent(in) :: self
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: ii
    ${TYPE}$, intent(inout) :: locmtx(:,:)

    integer :: ib

    do ib = 1, self%nblock
      call self%setblock_slave(mygrid, ii, ib, locmtx)
    end do

  end subroutine setline_slave_${SUFFIX}$

#:enddef setline_slave_template


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  #:for TYPE in TYPES
    #:set SUFFIX=TYPE
    #:set FTYPE=FORTRAN_TYPES[TYPE]
    $:getblock_master_template(SUFFIX, FTYPE)
    $:getblock_slave_template(SUFFIX, FTYPE)
    $:getline_master_template(SUFFIX,FTYPE)
    $:getline_slave_template(SUFFIX,FTYPE)
    $:setblock_master_template(SUFFIX,FTYPE)
    $:setblock_slave_template(SUFFIX,FTYPE)
    $:setline_master_template(SUFFIX,FTYPE)
    $:setline_slave_template(SUFFIX,FTYPE)
  #:endfor

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Helper routine to get global positions for given block.
  !!
  !! \param self  Instance
  !! \param mygrid BLACS descriptor
  !! \param ii  Index of the line (row/column).
  !! \param ib  Index of the block within the line.
  !! \param prow  Row of the processor owning the block.
  !! \param pcol  Column of the processor ownig the block.
  !! \param lrow  Starting row in the local matrix for the block.
  !! \param lcol  Starting column in the local matrix for the block.
  !! \param nrow  Nr. of rows corresponding to the block.
  !! \param ncol  Nr. of columns corresponding to the block.
  !!
  subroutine getpositions(self, mygrid, ii, ib, prow, pcol, lrow, lcol, &
      & nrow, ncol)
    class(linecomm), intent(in) :: self
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: ii, ib
    integer, intent(out) :: prow, pcol, lrow, lcol, nrow, ncol

    integer :: grow, gcol

    if (self%rowcollect) then
      grow = ii
      gcol = (ib - 1) * self%desc(NB_) + 1
      nrow = 1
      ncol = self%getblocksize(ib)
    else
      grow = (ib - 1) * self%desc(MB_) + 1
      gcol = ii
      nrow = self%getblocksize(ib)
      ncol = 1
    end if
    call scalafx_infog2l(mygrid, self%desc, grow, gcol, lrow, lcol, prow, pcol)

  end subroutine getpositions


end module linecomm_module
