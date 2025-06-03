#:include 'scalapackfx.fypp'

!> BLACS grid related routines.
module blacsgrid_module
  use scalapackfx_common_module
  use blacs_module
  implicit none
  private

  public :: blacsgrid

  !> Descriptor for BLACS grids.
  type :: blacsgrid
    integer :: ctxt = -1   !< Grid context.
    integer :: nrow = -1   !< Nr. of process rows in the grid.
    integer :: ncol = -1   !< Nr. of process columns in the grid.
    integer :: iproc = -1  !< Id of the current process in the grid.
    integer :: nproc = -1  !< Nr. of processes in the grid.
    integer :: myrow = -1  !< Row of the current process.
    integer :: mycol = -1  !< Column of the current process.
    integer :: leadproc = -1  !< Id of the lead process.
    integer :: leadrow = -1   !< Row of the lead process.
    integer :: leadcol = -1   !< Column of the lead process.
    logical :: lead = .false. !< Whether the current process is the lead.
  contains
    !> Creates process grid.
    procedure :: initgrid

    !> Creates separated equivalent process grids by equally splitting the current grid
    procedure :: initsplitgrids

    !> Created separated grid using explicit grid mapping
    procedure :: initmappedgrids

    !> Destructs grid.
    procedure :: destruct

    !> Common initialization routine obtaining BLACS system context.
    procedure, private :: initcontext

    !> Reset structure to initial value.
    procedure, private :: reset

  end type blacsgrid


contains


  !! Initializes the BLACS context for further operations.
  !!
  !! \param self  BLACS grid descriptor instance.
  !! \param context  System context (default: BLACS default context)
  !!
  subroutine initcontext(self, context)
    class(blacsgrid), intent(inout) :: self
    integer, intent(in), optional :: context

    if (present(context)) then
      self%ctxt = context
    else
      call blacs_get(-1, 0, self%ctxt)
    end if

  end subroutine initcontext


  !> Creates processor grid.
  !!
  !! Sets up a rectangular processor grid. All processes must call this routine
  !! collectively. If size of the grid is smaller than the number of available
  !! processes, those processes not fitting into the grid will obtain an
  !! uninitialized grid descriptor (with iproc = -1) at return.
  !!
  !! \param self  BLACS grid descriptor.
  !! \param nrow  Number of processor rows.
  !! \param ncol  Number of processor columns.
  !! \param colmajor  If .true., processes will be aligned in column major order
  !!     otherwise in row major order. (Default: .false.)
  !! \param leadrow  Lead row in each subgrid.
  !! \param leadcol  Lead column in each grid.
  !! \param context  BLACS system context (default: default system context)
  !! \param repeatable if present and T, forces topologies to be repeatable.
  !!     May degrade performance in this case.
  !!
  subroutine initgrid(self, nrow, ncol, colmajor, leadrow, &
      & leadcol, context, repeatable)
    class(blacsgrid), intent(inout) :: self
    integer, intent(in) :: nrow, ncol
    logical, intent(in), optional :: colmajor
    integer, intent(in), optional :: leadrow, leadcol, context
    logical, intent(in), optional :: repeatable

    logical :: colmajor0
    integer :: leadrow0, leadcol0
    integer :: nproc, iproc, gridsize
    integer :: val(1)

    ! Stop, if we do not have enough processess
    call blacs_pinfo(iproc, nproc)
    gridsize = nrow * ncol
    if (gridsize > nproc) then
      write(*, "(A,I0,A,I0,A)") "Nr. of processors (", nproc, &
          & ") less than size grid (", gridsize, ")"
      call blacs_exit(0)
      stop
    end if

    call self%initcontext(context)
    @:inoptflags(colmajor0, colmajor, .false.)
    @:inoptflags(leadrow0, leadrow, 0)
    @:inoptflags(leadcol0, leadcol, 0)
    if (colmajor0) then
      call blacs_gridinit(self%ctxt, "C", nrow, ncol)
    else
      call blacs_gridinit(self%ctxt, "R", nrow, ncol)
    end if

    if (present(repeatable)) then
      val = 0
      if (repeatable) then
        val = 1
      end if
      call blacs_set(self%ctxt,15,val)
    end if

    !! Processes not participating should obtain a cleared structure.
    if (iproc >= gridsize) then
      call self%reset()
    else
      call blacs_gridinfo(self%ctxt, self%nrow, self%ncol, self%myrow, &
          & self%mycol)
      self%nproc = gridsize
      if (colmajor0) then
        self%iproc = self%mycol * self%nrow + self%myrow
      else
        self%iproc = self%myrow * self%ncol + self%mycol
      end if
      self%leadrow = leadrow0
      self%leadcol = leadcol0
      self%leadproc = blacs_pnum(self%ctxt, self%leadrow, self%leadcol)
      self%lead = (self%myrow == self%leadrow .and. self%mycol == self%leadcol)
    end if

  end subroutine initgrid


  !> Creates equivalent independent subgrids by splitting the current one up.
  !!
  !! Sets up independent rectangular processor grids. All processes must call
  !! this routine collectively. If size of all grids is smaller than the
  !! number of available processes, those processes not fitting into any grid
  !! will obtain an uninitialized grid descriptor (with iproc = -1) at
  !! return. All other processes will obtain a grid descriptor for the subgrid
  !! they belong to.
  !!
  !! \param self  BLACS grid descriptor
  !! \param ngrid  Nr. of grids
  !! \param nrow  Number of rows in every grid
  !! \param ncol  Number of columns in every grid
  !! \param colmajor  If .true., processes will be aligned in column major order
  !!     otherwise in row major order. (Default: .false.)
  !! \param leadrow  Lead row in each subgrid.
  !! \param leadcol  Lead column in each grid.
  !! \param context  BLACS system context (default: default system context)
  !! \param leadgrid  If present, an additional (1, ngrid) shaped grid is
  !!     created which contains only the lead nodes from all grids.
  !! \param repeatable if present and T, forces topologies to be repeatable.
  !!     May degrade performance in this case.
  !!
  subroutine initsplitgrids(self, ngrid, nrow, ncol, colmajor, &
      & leadrow, leadcol, context, leadgrid, repeatable)
    class(blacsgrid), intent(inout) :: self
    integer, intent(in) :: nrow, ncol
    logical, intent(in), optional :: colmajor
    integer, intent(in), optional :: leadrow, leadcol, context
    class(blacsgrid), intent(out), optional :: leadgrid
    logical, intent(in), optional :: repeatable

    integer :: ngrid, ctxt0, gridsize, nproc, iproc
    integer :: ind, irow, icol, shift
    integer, allocatable :: imap(:,:), imap2(:,:)
    logical :: colmajor0
    integer :: leadrow0, leadcol0
    integer :: val(1)

    call blacs_pinfo(iproc, nproc)
    gridsize = nrow * ncol
    ! Stop, if we do not have enough processess
    if (ngrid * gridsize > nproc) then
      write(*, "(A,I0,A,I0,A)") "Nr. of processors (", nproc, &
          & ") not equal to size of all subgrids (", ngrid * gridsize, ")"
      call blacs_exit(0)
      stop
    end if

    call self%initcontext(context)
    ctxt0 = self%ctxt
    @:inoptflags(colmajor0, colmajor, .false.)
    @:inoptflags(leadrow0, leadrow, 0)
    @:inoptflags(leadcol0, leadcol, 0)

    ! Processes outside of the grid see the gridmap of the last process
    ! within the grid, so they realize they are not participating.
    allocate(imap(nrow, ncol))
    ind = min(iproc / gridsize, ngrid - 1) * gridsize
    if (colmajor0) then
      do icol = 1, ncol
        do irow = 1, nrow
          imap(irow, icol) = ind
          ind = ind + 1
        end do
      end do
    else
      do irow = 1, nrow
        do icol = 1, ncol
          imap(irow, icol) = ind
          ind = ind + 1
        end do
      end do
    end if

    ! Create grid
    self%ctxt = ctxt0
    call blacs_gridmap(self%ctxt, imap, size(imap, dim=1), nrow, ncol)

    if (present(repeatable)) then
      val = 0
      if (repeatable) then
        val = 1
      end if
      call blacs_set(self%ctxt,15,val)
    end if

    if (iproc >= ngrid * gridsize) then
      call self%reset()
    else
      call blacs_gridinfo(self%ctxt, self%nrow, self%ncol, &
          & self%myrow, self%mycol)
      self%nproc = gridsize
      if (colmajor0) then
        self%iproc = self%mycol * self%nrow + self%myrow
      else
        self%iproc = self%myrow * self%ncol + self%mycol
      end if
      self%leadrow = leadrow0
      self%leadcol = leadcol0
      self%lead = (self%myrow == self%leadrow &
          & .and. self%mycol == self%leadcol)
    end if

    ! Return here if leadgrid is not asked for.
    if (.not. present(leadgrid)) then
      return
    end if
    if (colmajor0) then
      shift = nrow * leadcol0 + leadrow0
    else
      shift = ncol * leadrow0 + leadcol0
    end if
    allocate(imap2(1, ngrid))
    do icol = 1, ngrid
      imap2(1, icol) = (icol - 1) * gridsize + shift
    end do
    leadgrid%ctxt = ctxt0
    call blacs_gridmap(leadgrid%ctxt, imap2, size(imap2, dim=1), 1, ngrid)
    if (mod(iproc, gridsize) /= shift) then
      call leadgrid%reset()
    else
      call blacs_gridinfo(leadgrid%ctxt, leadgrid%nrow, &
          & leadgrid%ncol, leadgrid%myrow, leadgrid%mycol)
      leadgrid%nproc = ngrid
      leadgrid%iproc = iproc / gridsize
      leadgrid%leadrow = 0
      leadgrid%leadcol = 0
      leadgrid%lead = (leadgrid%iproc == 0)
    end if

  end subroutine initsplitgrids



  !> Creates subgrids by explicitly mapping specified processes into a subgrid.
  !!
  !! Sets up a independent rectangular processor grids. All processes must call
  !! this routine at collectively. If size of all grids is smaller than the
  !! number of available processes, those processes not fitting into any grid
  !! will obtain an uninitialized grid descriptor (with iproc = -1) at
  !! return. All other processes will obtain a grid descriptor for the subgrid
  !! they belong to.
  !!
  !! \param self  BLACS grid descriptor
  !! \param ngrid  Nr. of grids
  !! \param nrow  Number of rows in every grid
  !! \param ncol  Number of columns in every grid
  !! \param colmajor  If .true., processes will be aligned in column major order
  !!     otherwise in row major order. (Default: .false.)
  !! \param leadrow  Lead row in each subgrid.
  !! \param leadcol  Lead column in each grid.
  !! \param context  BLACS system context (default: default system context)
  !! \param leadgrid  If present, an additional (1, ngrid) shaped grid is
  !!     created which contains only the lead nodes from all grids.
  !! \param repeatable if present and T, forces topologies to be repeatable.
  !!     May degrade performance in this case.
  !!
  subroutine initmappedgrids(self, gridmap, context, repeatable)
    class(blacsgrid), intent(inout) :: self
    integer, intent(in) :: gridmap(:,:)
    integer, intent(in), optional :: context
    logical, intent(in), optional :: repeatable

    integer :: ncol, nrow, gridsize, nproc, iproc, val(1)

    call blacs_pinfo(iproc, nproc)
    nrow = size(gridmap, dim=1)
    ncol = size(gridmap, dim=2)
    gridsize = size(gridmap)

    call self%initcontext(context)

    call blacs_gridmap(self%ctxt, gridmap, size(gridmap, dim=1), nrow, ncol)

    if (present(repeatable)) then
      val = 0
      if (repeatable) then
        val = 1
      end if
      call blacs_set(self%ctxt,15,val)
    end if

    if (.not. any(gridmap == iproc)) then
      call self%reset()
    else
      call blacs_gridinfo(self%ctxt, self%nrow, self%ncol, self%myrow, self%mycol)
      self%nproc = gridsize
      self%iproc = self%myrow * self%ncol + self%mycol
      self%leadrow = 0
      self%leadcol = 0
      self%lead = (self%myrow == self%leadrow .and. self%mycol == self%leadcol)
    end if

  end subroutine initmappedgrids


  !> Destructs processor grid.
  !!
  !! \param mygrid BLACS grid descriptor.
  !!
  subroutine destruct(self)
    class(blacsgrid), intent(inout) :: self

    if (self%iproc /= -1) then
      call blacs_gridexit(self%ctxt)
      call self%reset()
    end if

  end subroutine destruct


  !> Resets structure to its uninitialized value.
  !!
  !! \param mygrid BLACS grid descriptor.
  !!
  subroutine reset(self)
    class(blacsgrid), intent(inout) :: self

    self%ctxt = -1
    self%nrow = -1
    self%ncol = -1
    self%iproc = -1
    self%nproc = -1
    self%myrow = -1
    self%mycol = -1
    self%leadrow = -1
    self%leadcol = -1
    self%lead = .false.

  end subroutine reset


end module blacsgrid_module
