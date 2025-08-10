#:include 'scalapackfx.fypp'
#:set TYPES = NUMERIC_TYPES


#! ************************************************************************
#! *** cpl2g
#! ************************************************************************

  
#:def cpl2g_template(SUFFIX, FTYPE)

  #! ${SUFFIX}$ subroutine suffix
  #! ${FTYPE}$ dummy argument type
  !> Copies the content of a local matrix to a global one (${SUFFIX}$).
  !!
  !! \param mygrid BLACS descriptor
  !! \param loc  Local matrix.
  !! \param desc  Descriptor of the global matrix.
  !! \param ii  Starting row in the global matrix.
  !! \param jj  Starting column in the global matrix
  !! \param glob  Local part of the global matrix.
  !!
  subroutine cpl2g_${SUFFIX}$(mygrid, loc, desc, ii, jj, glob)
    type(blacsgrid), intent(in) :: mygrid
    ${FTYPE}$, intent(in) :: loc(:,:)
    integer, intent(in) :: desc(DLEN_)
    integer, intent(in) :: ii, jj
    ${FTYPE}$, intent(inout) :: glob(:,:)
  
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
    
  end subroutine cpl2g_${SUFFIX}$

#:enddef

#! ************************************************************************
#! *** addl2g
#! ************************************************************************

  
#:def addl2g_template(SUFFIX, FTYPE)

  #! ${SUFFIX}$ subroutine suffix
  #! ${FTYPE}$ dummy argument type
  !> Adds the content of a local matrix to a global one (${SUFFIX}$).
  !!
  !! \param mygrid BLACS descriptor
  !! \param loc  Local matrix.
  !! \param desc  Descriptor of the global matrix.
  !! \param ii  Starting row in the global matrix.
  !! \param jj  Starting column in the global matrix
  !! \param glob  Local part of the global matrix.
  !!
  subroutine addl2g_${SUFFIX}$(mygrid, loc, desc, ii, jj, glob)
    type(blacsgrid), intent(in) :: mygrid
    ${FTYPE}$, intent(in) :: loc(:,:)
    integer, intent(in) :: desc(DLEN_)
    integer, intent(in) :: ii, jj
    ${FTYPE}$, intent(inout) :: glob(:,:)
  
    integer :: i2, j2, nr, nc
    integer, dimension(size(loc, dim=1)) :: irows, iloc, prow
    integer, dimension(size(loc, dim=2)) :: icols, jloc, pcol

    nr = size(loc, dim=1)
    nc = size(loc, dim=2)

    do i2 = 1, nr
      irows(i2) = i2 + ii - 1
    end do

    do j2 = 1, nc
      icols(j2) = j2 + jj - 1
    end do

    call scalafx_infog2l(mygrid, desc, irows, icols, iloc, jloc,&
        & prow, pcol, .false.)
  
    do j2 = 1, nc
      if (pcol(j2) /= mygrid%mycol) then
        cycle
      end if
      do i2 = 1, nr
        if (prow(i2) /= mygrid%myrow) then
          cycle
        end if
        glob(iloc(i2), jloc(j2)) = glob(iloc(i2), jloc(j2)) + loc(i2, j2)
      end do
    end do
    
  end subroutine addl2g_${SUFFIX}$

#:enddef

#! ************************************************************************
#! *** cpg2l
#! ************************************************************************

  
#:def cpg2l_template(SUFFIX, FTYPE)

  #! ${SUFFIX}$ subroutine suffix
  #! ${FTYPE}$ dummy argument type
  !> Copies the content from the global matrix into a local one.
  !!
  !! \param mygrid BLACS descriptor
  !! \param desc  Descriptor of the global matrix.
  !! \param ii  Starting row in the global matrix.
  !! \param jj  Starting column in the global matrix
  !! \param glob  Local part of the global matrix.
  !! \param loc  Local matrix.
  !!
  subroutine cpg2l_${SUFFIX}$(mygrid, desc, ii, jj, glob, loc)
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: desc(DLEN_)
    integer, intent(in) :: ii, jj
    ${FTYPE}$, intent(in) :: glob(:,:)
    ${FTYPE}$, intent(out) :: loc(:,:)
    
    integer :: i2, j2, nr, nc
    integer, dimension(size(loc, dim=1)) :: irows, iloc, prow
    integer, dimension(size(loc, dim=2)) :: icols, jloc, pcol

    nr = size(loc, dim=1)
    nc = size(loc, dim=2)

    do i2 = 1, nr
      irows(i2) = i2 + ii - 1
    end do

    do j2 = 1, nc
      icols(j2) = j2 + jj - 1
    end do

    call scalafx_infog2l(mygrid, desc, irows, icols, iloc, jloc,&
        & prow, pcol, .false.)
  
    do j2 = 1, nc
      if (pcol(j2) == mygrid%mycol) then
        do i2 = 1, nr
          if (prow(i2) == mygrid%myrow) then
            loc(i2, j2) = glob(iloc(i2), jloc(j2))
          else
            loc(i2, j2) = 0.0_dp
          end if
        end do
      else
        loc(:, j2) = 0.0_dp
      end if
    end do
    
  end subroutine cpg2l_${SUFFIX}$

#:enddef

#! ************************************************************************
#! *** addg2l
#! ************************************************************************

  
#:def addg2l_template(SUFFIX, FTYPE)

  #! ${SUFFIX}$ subroutine suffix
  #! ${FTYPE}$ dummy argument type
  !> Copies the content from the global matrix into a local one.
  !!
  !! \param mygrid BLACS descriptor
  !! \param desc  Descriptor of the global matrix.
  !! \param ii  Starting row in the global matrix.
  !! \param jj  Starting column in the global matrix
  !! \param glob  Local part of the global matrix.
  !! \param loc  Local matrix.
  !!
  subroutine addg2l_${SUFFIX}$(mygrid, desc, ii, jj, glob, loc)
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: desc(DLEN_)
    integer, intent(in) :: ii, jj
    ${FTYPE}$, intent(in) :: glob(:,:)
    ${FTYPE}$, intent(out) :: loc(:,:)
    
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
    
  end subroutine addg2l_${SUFFIX}$

#:enddef

#! ************************************************************************
#! *** writearray_lead
#! ************************************************************************

  
#:def writearray_lead_template(SUFFIX, FTYPE)

  #!
  #! ${SUFFIX}$ subroutine suffix
  #! ${FTYPE}$ dummy argument type
  #!
  !> Writes a distributed array to a file (lead, ${SUFFIX}$).
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
  subroutine writearray_lead_${SUFFIX}$(mygrid, fd, desc, mtxloc, rowwise, elemformat)
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: fd, desc(DLEN_)
    ${FTYPE}$, intent(in) :: mtxloc(:,:)
    logical, intent(in), optional :: rowwise
    character(*), intent(in), optional :: elemformat
  
    type(linecomm) :: distributor
    ${FTYPE}$, allocatable :: buffer(:)
    character(:), allocatable :: lineformat
    logical :: rowwise0, formatted
    integer :: nline, linelen, ii, ndigit, nn
  
    @:inoptflags(rowwise0, rowwise, .false.)
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
      call distributor%getline_lead(mygrid, ii, mtxloc, buffer)
      if (formatted) then
        write(fd, lineformat) buffer
      else
        write(fd) buffer
      end if
    end do
    
  end subroutine writearray_lead_${SUFFIX}$

#:enddef


#! ************************************************************************
#! *** writearray_follow
#! ************************************************************************

  
#:def writearray_follow_template(SUFFIX, FTYPE)

  #!
  #! ${SUFFIX}$ subroutine suffix
  #! ${FTYPE}$ dummy argument type
  #!
  !> Writes a distributed array to a file (follow, ${SUFFIX}$).
  !!
  !! \param mygrid  BLACS descriptor
  !! \param desc  Descriptor of the distributed matrix
  !! \param mtxloc  Local portion of the distributed matrix
  !! \param rowwise  If .true. matrix is dumped rowwise otherwise columnwise
  !!
  subroutine writearray_follow_${SUFFIX}$(mygrid, desc, mtxloc, rowwise)
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: desc(DLEN_)
    ${FTYPE}$, intent(in) :: mtxloc(:,:)
    logical, intent(in), optional :: rowwise
  
    type(linecomm) :: distributor
    integer :: ii, nline
    logical :: rowwise0
  
    @:inoptflags(rowwise0, rowwise, .false.)
    if (rowwise0) then
      call distributor%init(mygrid, desc, "r")
      nline = desc(N_)
    else
      call distributor%init(mygrid, desc, "c")
      nline = desc(M_)
    end if
    do ii = 1, nline
      call distributor%getline_follow(mygrid, ii, mtxloc)
    end do
    
  end subroutine writearray_follow_${SUFFIX}$

#:enddef

#! ************************************************************************
#! *** readarray_lead
#! ************************************************************************

  
#:def readarray_lead_template(SUFFIX, FTYPE)

  #!
  #! ${SUFFIX}$ subroutine suffix
  #! ${FTYPE}$ dummy argument type
  #!
  !> Reads a distributed array from a file (lead, ${SUFFIX}$).
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
  subroutine readarray_lead_${SUFFIX}$(mygrid, fd, desc, mtxloc, rowwise, formatted)
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: fd, desc(DLEN_)
    ${FTYPE}$, intent(inout) :: mtxloc(:,:)
    logical, intent(in), optional :: rowwise, formatted
  
    type(linecomm) :: collector
    ${FTYPE}$, allocatable :: buffer(:)
    logical :: rowwise0, formatted0
    integer :: nline, linelen, ii
  
    @:inoptflags(rowwise0, rowwise, .false.)
    @:inoptflags(formatted0, formatted, .false.)
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
      call collector%setline_lead(mygrid, ii, buffer, mtxloc)
    end do
    
  end subroutine readarray_lead_${SUFFIX}$

#:enddef

#! ************************************************************************
#! *** readarray_follow
#! ************************************************************************

  
#:def readarray_follow_template(SUFFIX, FTYPE)

  #!
  #! ${SUFFIX}$ subroutine suffix
  #! ${FTYPE}$ dummy argument type
  #!
  !> Reads a distributed array from a file (follow, ${SUFFIX}$).
  !!
  !! \param mygrid  BLACS descriptor
  !! \param desc  Descriptor of the distributed matrix
  !! \param mtxloc  Local portion of the distributed matrix
  !! \param rowwise  If .true. matrix is assumed to be stored rowwise otherwise
  !!     columnwise (default: .false.)
  !!
  subroutine readarray_follow_${SUFFIX}$(mygrid, desc, mtxloc, rowwise)
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: desc(DLEN_)
    ${FTYPE}$, intent(inout) :: mtxloc(:,:)
    logical, intent(in), optional :: rowwise
  
    type(linecomm) :: collector
    integer :: ii, nline
    logical :: rowwise0
  
    @:inoptflags(rowwise0, rowwise, .false.)
    if (rowwise0) then
      call collector%init(mygrid, desc, "r")
      nline = desc(N_)
    else
      call collector%init(mygrid, desc, "c")
      nline = desc(M_)
    end if
    do ii = 1, nline
      call collector%setline_follow(mygrid, ii, mtxloc)
    end do
    
  end subroutine readarray_follow_${SUFFIX}$

#:enddef
  
!> Some extension routines to the scalapack library making it more usable.
!!
!! \details Some of the types exported by this module are defined in other
!! modules. Therefore, see also the documentation of the following modules:
!!
!!   * \ref linecomm_module "linecomm_module"
!!
module scalapackfx_tools_module
  use scalapackfx_common_module
  use blacsfx_module
  use scalapackfx_module
  use linecomm_module
  implicit none
  private

  public :: scalafx_cpl2g, scalafx_cpg2l, scalafx_addl2g, scalafx_addg2l
  public :: scalafx_islocal
  public :: writearray_lead, writearray_follow
  public :: readarray_lead, readarray_follow
  public :: blocklist, size
  public :: linecomm

  interface scalafx_cpl2g
    module procedure cpl2g_real, cpl2g_dreal, cpl2g_complex, cpl2g_dcomplex
    module procedure cpl2g_int
  end interface scalafx_cpl2g

  interface scalafx_cpg2l
    module procedure cpg2l_real, cpg2l_dreal, cpg2l_complex, cpg2l_dcomplex
    module procedure cpg2l_int
  end interface scalafx_cpg2l

  interface scalafx_addl2g
    module procedure addl2g_real, addl2g_dreal, addl2g_complex, addl2g_dcomplex
    module procedure addl2g_int
  end interface scalafx_addl2g

  interface scalafx_addg2l
    module procedure addg2l_real, addg2l_dreal, addg2l_complex, addg2l_dcomplex
    module procedure addg2l_int
  end interface scalafx_addg2l

  interface writearray_lead
    module procedure writearray_lead_int
    module procedure writearray_lead_real, writearray_lead_dreal
    module procedure writearray_lead_complex, writearray_lead_dcomplex
  end interface writearray_lead

  interface writearray_follow
    module procedure writearray_follow_int
    module procedure writearray_follow_real, writearray_follow_dreal
    module procedure writearray_follow_complex, writearray_follow_dcomplex
  end interface writearray_follow

  interface readarray_lead
    module procedure readarray_lead_int
    module procedure readarray_lead_real, readarray_lead_dreal
    module procedure readarray_lead_complex, readarray_lead_dcomplex
  end interface readarray_lead

  interface readarray_follow
    module procedure readarray_follow_int
    module procedure readarray_follow_real, readarray_follow_dreal
    module procedure readarray_follow_complex, readarray_follow_dcomplex
  end interface readarray_follow
  

  !> List of the local blocks of a distributed matrix.
  !!
  !! \details This structure can be helpful when modifying a distributed
  !! matrix directly on the local nodes. For example, in order to multiply
  !! every column of a distributed matrix by a column dependent factor, 
  !! you could use the blocklist the following way:
  !!
  !!     type(blocklist) :: blocks
  !!     integer :: ii, jj, jglob, jloc, bsize
  !!     :
  !!     call blocks%init(mygrid, descaa, "c")
  !!     do ii = 1, size(blocks)
  !!       call blocks%getblock(ii, jglob, jloc, bsize)
  !!       do jj = 0, bsize  - 1
  !!         aa(:,jloc + jj) =  aa(:,jloc + jj) * ff(jglob + jj)
  !!       end do
  !!     end do
  !!
  type :: blocklist
    private
    integer :: nn, nb, nproc, myproc, srcproc, nblock
  contains
    !> Initializes the instance.
    procedure :: init => blocklist_init

    !> Returns the size of the blocklist.
    procedure :: getsize => blocklist_getsize

    !> Returns the indices (local and global) of a local block.
    procedure :: getblock => blocklist_getblock
  end type blocklist

  interface size
    module procedure blocklist_getsize
  end interface size

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Blocklist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initializes a blocklist instance.
  !! \param self  Initialized instance on exit.
  !! \param mygrid  BLACS descriptor
  !! \param desc  Descriptor of the distributed matrix.
  !! \param rowcol  "C" for column, "R" for row blocks.
  subroutine blocklist_init(self, mygrid, desc, rowcol)
    class(blocklist), intent(inout) :: self
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: desc(DLEN_)
    character, intent(in) :: rowcol

    integer :: nblockall, nextrablock, mydist

    if (rowcol == "c" .or. rowcol == "C") then
      self%nn = desc(N_)
      self%nb = desc(NB_)
      self%nproc = mygrid%ncol
      self%myproc = mygrid%mycol
      self%srcproc = desc(CSRC_)
    else
      self%nn = desc(M_)
      self%nb = desc(MB_)
      self%nproc = mygrid%nrow
      self%myproc = mygrid%myrow
      self%srcproc = desc(RSRC_)
    end if
    nblockall = self%nn / self%nb
    self%nblock = nblockall / self%nproc
    nextrablock = mod(nblockall, self%nproc)
    mydist = mod(self%nproc + self%myproc - self%srcproc, self%nproc)
    if (mydist < nextrablock) then
      self%nblock = self%nblock + 1
    elseif (mydist == nextrablock .and. mod(self%nn, self%nb) > 0) then
      self%nblock = self%nblock +1
    end if

  end subroutine blocklist_init


  !> Returns the size of the blocklist.
  !! \param self  Instance.
  !! \returns Number of local blocks of the distributed matrix.
  function blocklist_getsize(self) result(res)
    class(blocklist), intent(in) :: self
    integer :: res

    res = self%nblock

  end function blocklist_getsize


  !> Returns the indices (local and global) of a local block.
  !! \param self  Blocklist instance.
  !! \param iblock  Number of local block.
  !! \param iglob  Index of the first element of the block in the global matrix.
  !! \param iloc  Index of the first element of the block in the local matirx.
  !! \param bsize  Size of the block (number of elements in the block).
  subroutine blocklist_getblock(self, iblock, iglob, iloc, bsize)
    class(blocklist), intent(in) :: self
    integer, intent(in) :: iblock
    integer, intent(out) :: iglob, iloc, bsize

    integer :: mydist

    if (iblock >= 1 .and. iblock <= self%nblock) then
      mydist = mod(self%nproc + self%myproc - self%srcproc, self%nproc)
      iglob = ((iblock - 1) * self%nproc + mydist) * self%nb + 1
      iloc = (iblock - 1) * self%nb + mod(iglob - 1, self%nb) + 1
      bsize = min(self%nb, self%nn - iglob + 1)
    else
      iglob = 0
      iloc = 0
      bsize = 0
    end if

  end subroutine blocklist_getblock


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Block copy/adding routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  @:cpl2g_template(real, real(sp))
  @:cpl2g_template(dreal, real(dp))
  @:cpl2g_template(complex, complex(sp))
  @:cpl2g_template(dcomplex, complex(dp))
  @:cpl2g_template(int, integer)

  @:cpg2l_template(real, real(sp))
  @:cpg2l_template(dreal, real(dp))
  @:cpg2l_template(complex, complex(sp))
  @:cpg2l_template(dcomplex, complex(dp))
  @:cpg2l_template(int, integer)

  @:addl2g_template(real, real(sp))
  @:addl2g_template(dreal, real(dp))
  @:addl2g_template(complex, complex(sp))
  @:addl2g_template(dcomplex, complex(dp))
  @:addl2g_template(int, integer)

  @:addg2l_template(real, real(sp))
  @:addg2l_template(dreal, real(dp))
  @:addg2l_template(complex, complex(sp))
  @:addg2l_template(dcomplex, complex(dp))
  @:addg2l_template(int, integer)
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! writearray/readarray
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  @:writearray_lead_template(int, integer)
  @:writearray_lead_template(real, real(sp))
  @:writearray_lead_template(dreal, real(dp))
  @:writearray_lead_template(complex, complex(sp))
  @:writearray_lead_template(dcomplex, complex(dp))

  @:writearray_follow_template(int, integer)
  @:writearray_follow_template(real, real(sp))
  @:writearray_follow_template(dreal, real(dp))
  @:writearray_follow_template(complex, complex(sp))
  @:writearray_follow_template(dcomplex, complex(dp))

  @:readarray_lead_template(int, integer)
  @:readarray_lead_template(real, real(sp))
  @:readarray_lead_template(dreal, real(dp))
  @:readarray_lead_template(complex, complex(sp))
  @:readarray_lead_template(dcomplex, complex(dp))

  @:readarray_follow_template(int, integer)
  @:readarray_follow_template(real, real(sp))
  @:readarray_follow_template(dreal, real(dp))
  @:readarray_follow_template(complex, complex(sp))
  @:readarray_follow_template(dcomplex, complex(dp))

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Other routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !> Returns whether a given global position can be found in the local part
  !! of the distributed matrix.
  !!
  !! \param grid  Grid on which the matrix is distributed.
  !! \param desc  Matrix descriptor
  !! \param ii  Row of the global position
  !! \param jj  Column of the global position
  !! \param local  Whether the position is local
  !! \param iloc  Row of the position in the local matrix  (only meaningful,
  !!     if local is .true.)
  !! \param jloc  Columnt of the position in the local matrix  (only
  !!     meaningful, if local is .true.)
  !!
  subroutine scalafx_islocal(grid, desc, ii, jj, local, iloc, jloc)
    type(blacsgrid), intent(in) :: grid
    integer, intent(in) :: desc(DLEN_)
    integer, intent(in) :: ii, jj
    logical, intent(out) :: local
    integer, intent(out) :: iloc, jloc

    integer :: prow, pcol

    call scalafx_infog2l(grid, desc, ii, jj, iloc, jloc, prow, pcol)
    local = (prow == grid%myrow .and. pcol == grid%mycol)

  end subroutine scalafx_islocal

  
end module scalapackfx_tools_module
