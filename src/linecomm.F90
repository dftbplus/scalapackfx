include(linecomm.m4)

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
    
    procedure :: getblock_master_int
    procedure :: getblock_master_real
    procedure :: getblock_master_dreal
    procedure :: getblock_master_complex
    procedure :: getblock_master_dcomplex
    procedure :: getblock_slave_int
    procedure :: getblock_slave_real
    procedure :: getblock_slave_dreal
    procedure :: getblock_slave_complex
    procedure :: getblock_slave_dcomplex
    generic :: getblock_master => getblock_master_int, getblock_master_real, &
        & getblock_master_dreal, getblock_master_complex, &
        & getblock_master_dcomplex
    generic :: getblock_slave => getblock_slave_int, getblock_slave_real, &
        & getblock_slave_dreal, getblock_slave_complex, getblock_slave_dcomplex

    procedure :: getline_master_int
    procedure :: getline_master_real
    procedure :: getline_master_dreal
    procedure :: getline_master_complex
    procedure :: getline_master_dcomplex
    procedure :: getline_slave_int
    procedure :: getline_slave_real
    procedure :: getline_slave_dreal
    procedure :: getline_slave_complex
    procedure :: getline_slave_dcomplex
    generic :: getline_master => getline_master_int, getline_master_real, &
        & getline_master_dreal, getline_master_complex, getline_master_dcomplex
    generic :: getline_slave => getline_slave_int, getline_slave_real, &
        & getline_slave_dreal, getline_slave_complex, getline_slave_dcomplex
    
    procedure :: setblock_master_int
    procedure :: setblock_master_real
    procedure :: setblock_master_dreal
    procedure :: setblock_master_complex
    procedure :: setblock_master_dcomplex
    procedure :: setblock_slave_int
    procedure :: setblock_slave_real
    procedure :: setblock_slave_dreal
    procedure :: setblock_slave_complex
    procedure :: setblock_slave_dcomplex
    generic :: setblock_master => setblock_master_int, setblock_master_real, &
        & setblock_master_dreal, setblock_master_complex, &
        & setblock_master_dcomplex
    generic :: setblock_slave => setblock_slave_int, setblock_slave_real, &
        & setblock_slave_dreal, setblock_slave_complex, setblock_slave_dcomplex

    procedure :: setline_master_int
    procedure :: setline_master_real
    procedure :: setline_master_dreal
    procedure :: setline_master_complex
    procedure :: setline_master_dcomplex
    procedure :: setline_slave_int
    procedure :: setline_slave_real
    procedure :: setline_slave_dreal
    procedure :: setline_slave_complex
    procedure :: setline_slave_dcomplex
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

    _handle_inoptflag(self%iorow, iorow, mygrid%masterrow)
    _handle_inoptflag(self%iocol, iocol, mygrid%mastercol)
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


  _subroutine_getblock_master(int, integer)
  _subroutine_getblock_master(real, real(sp))
  _subroutine_getblock_master(dreal, real(dp))
  _subroutine_getblock_master(complex, complex(sp))
  _subroutine_getblock_master(dcomplex, complex(dp))

  _subroutine_getblock_slave(int, integer)
  _subroutine_getblock_slave(real, real(sp))
  _subroutine_getblock_slave(dreal, real(dp))
  _subroutine_getblock_slave(complex, complex(sp))
  _subroutine_getblock_slave(dcomplex, complex(dp))

  _subroutine_getline_master(int, integer)
  _subroutine_getline_master(real, real(sp))
  _subroutine_getline_master(dreal, real(dp))
  _subroutine_getline_master(complex, complex(sp))
  _subroutine_getline_master(dcomplex, complex(dp))

  _subroutine_getline_slave(int, integer)
  _subroutine_getline_slave(real, real(sp))
  _subroutine_getline_slave(dreal, real(dp))
  _subroutine_getline_slave(complex, complex(sp))
  _subroutine_getline_slave(dcomplex, complex(dp))

  _subroutine_setblock_master(int, integer)
  _subroutine_setblock_master(real, real(sp))
  _subroutine_setblock_master(dreal, real(dp))
  _subroutine_setblock_master(complex, complex(sp))
  _subroutine_setblock_master(dcomplex, complex(dp))

  _subroutine_setblock_slave(int, integer)
  _subroutine_setblock_slave(real, real(sp))
  _subroutine_setblock_slave(dreal, real(dp))
  _subroutine_setblock_slave(complex, complex(sp))
  _subroutine_setblock_slave(dcomplex, complex(dp))

  _subroutine_setline_master(int, integer)
  _subroutine_setline_master(real, real(sp))
  _subroutine_setline_master(dreal, real(dp))
  _subroutine_setline_master(complex, complex(sp))
  _subroutine_setline_master(dcomplex, complex(dp))

  _subroutine_setline_slave(int, integer)
  _subroutine_setline_slave(real, real(sp))
  _subroutine_setline_slave(dreal, real(dp))
  _subroutine_setline_slave(complex, complex(sp))
  _subroutine_setline_slave(dcomplex, complex(dp))


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
