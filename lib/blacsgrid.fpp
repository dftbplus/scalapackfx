#:include 'scalapackfx.fypp'

!> BLACS grid related routines.
!!
!! Note: The grid creating routines in this module accept an optional BLACS system context. If it is
!! not present, the default system context (containing all processes) is used. Otherwise, the
!! context can be any valid BLACS system context, e.g. one corresponding to an MPI communicator
!! containing only a subset of all processes (for MPI-based BLACS implementations, this is usually
!! the Fortran handle of the communicator itself). In latter case, all processes of the system
!! context (and only those) must call the grid creating routines collectively. Process ids in the
!! grid maps are always interpreted relative to the system context. As BLACS can not be queried for
!! the number of processes in a system context, the consistency checks on the grid sizes are only
!! carried out when the default system context is used.
!!
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

    !> Reset structure to initial value.
    procedure, private :: reset

    !> Fills up the descriptor after the BLACS grid has been created.
    procedure, private :: setup_after_gridinit

  end type blacsgrid


contains


  !> Creates processor grid.
  !!
  !! Sets up a rectangular processor grid. All processes of the system context must call this
  !! routine collectively. If the size of the grid is smaller than the number of available
  !! processes, those processes not fitting into the grid will obtain an uninitialized grid
  !! descriptor (with iproc = -1) at return.
  !!
  !! \param self  BLACS grid descriptor.
  !! \param nrow  Number of processor rows.
  !! \param ncol  Number of processor columns.
  !! \param colmajor  If .true., processes will be aligned in column major order, otherwise in row
  !!     major order. (Default: .false.)
  !! \param leadrow  Lead row in each subgrid.
  !! \param leadcol  Lead column in each grid.
  !! \param context  BLACS system context (default: default system context)
  !! \param repeatable if present and T, forces topologies to be repeatable. May degrade performance
  !!     in this case.
  !! \param info  If present, it is set to 0 on success and to a non-zero value on error (and
  !!     the grid is left uninitialized). Otherwise if absent the code stops on error.
  !!
  subroutine initgrid(self, nrow, ncol, colmajor, leadrow, &
      & leadcol, context, repeatable, info)
    class(blacsgrid), intent(inout) :: self
    integer, intent(in) :: nrow, ncol
    logical, intent(in), optional :: colmajor
    integer, intent(in), optional :: leadrow, leadcol, context
    logical, intent(in), optional :: repeatable
    integer, intent(out), optional :: info

    logical :: colmajor0
    integer :: leadrow0, leadcol0
    integer :: ctxt0, nproc0
    integer :: gridsize
    character(200) :: msg

    if (present(info)) info = 0
    call self%reset()

    gridsize = nrow * ncol

    call get_syscontext(context, ctxt0, nproc0)
    if (nproc0 > 0 .and. gridsize > nproc0) then
      write(msg, "(A,I0,A,I0,A)") "Nr. of processors (", nproc0, &
          & ") less than size grid (", gridsize, ")"
      call handle_error(trim(msg), info)
      return
    end if

    self%ctxt = ctxt0
    @:inoptflags(colmajor0, colmajor, .false.)
    @:inoptflags(leadrow0, leadrow, 0)
    @:inoptflags(leadcol0, leadcol, 0)
    if (colmajor0) then
      call blacs_gridinit(self%ctxt, "C", nrow, ncol)
    else
      call blacs_gridinit(self%ctxt, "R", nrow, ncol)
    end if
    call self%setup_after_gridinit(colmajor0, leadrow0, leadcol0, repeatable)

  end subroutine initgrid


  !> Creates equivalent independent subgrids by splitting the current one up.
  !!
  !! Sets up independent rectangular processor grids. All processes of the system context must call
  !! this routine collectively. If size of all grids is smaller than the number of available
  !! processes, those processes not fitting into any grid will obtain an uninitialized grid
  !! descriptor (with iproc = -1) at return. All other processes will obtain a grid descriptor for
  !! the subgrid they belong to.
  !!
  !! Note: The subgrids are created one after the other by separate collective calls to
  !! blacs_gridmap (one for each subgrid), so that all processes pass the same map in each call.
  !! Every process keeps the context of the subgrid it belongs to.
  !!
  !! \param self  BLACS grid descriptor
  !! \param ngrid  Nr. of grids
  !! \param nrow  Number of rows in every grid
  !! \param ncol  Number of columns in every grid
  !! \param colmajor  If .true., processes will be aligned in column major order otherwise, in row
  !!     major order. (Default: .false.)
  !! \param leadrow  Lead row in each subgrid.
  !! \param leadcol  Lead column in each grid.
  !! \param context  BLACS system context (default: default system context)
  !! \param leadgrid  If present, an additional (1, ngrid) shaped grid is created which contains
  !!     only the lead nodes from all grids.
  !! \param repeatable if present and T, forces topologies to be repeatable.
  !!     May degrade performance in this case.
  !! \param info  If present, it is set to 0 on success and to a non-zero value on error (and
  !!     the grids are left uninitialized). Otherwise if absent the code stops on error.
  !!
  subroutine initsplitgrids(self, ngrid, nrow, ncol, colmajor, &
      & leadrow, leadcol, context, leadgrid, repeatable, info)
    class(blacsgrid), intent(inout) :: self
    integer, intent(in) :: ngrid
    integer, intent(in) :: nrow, ncol
    logical, intent(in), optional :: colmajor
    integer, intent(in), optional :: leadrow, leadcol, context
    class(blacsgrid), intent(out), optional :: leadgrid
    logical, intent(in), optional :: repeatable
    integer, intent(out), optional :: info

    integer :: ctxt0, nproc0
    integer :: gridsize, ind, irow, icol, shift, igrid
    integer :: ctxt, nrowtmp, ncoltmp, myrow, mycol
    integer, allocatable :: imap(:,:), imap2(:,:)
    logical :: colmajor0
    integer :: leadrow0, leadcol0
    character(200) :: msg

    if (present(info)) info = 0
    call self%reset()
    if (present(leadgrid)) call leadgrid%reset()

    gridsize = nrow * ncol

    call get_syscontext(context, ctxt0, nproc0)
    if (nproc0 > 0 .and. ngrid * gridsize > nproc0) then
      write(msg, "(A,I0,A,I0,A)") "Nr. of processors (", nproc0, &
          & ") less than size of all subgrids (", ngrid * gridsize, ")"
      call handle_error(trim(msg), info)
      return
    end if

    @:inoptflags(colmajor0, colmajor, .false.)
    @:inoptflags(leadrow0, leadrow, 0)
    @:inoptflags(leadcol0, leadcol, 0)

    ! Create subgrids by ngrid collective calls. Every process in the call passes the same gridmap
    ! (the one describing the subgrid to be created) to blacs_gridmap. After creation of each
    ! subgrid, each process checks whether it belongs to that subgrid and stores its context if so.
    ! As the subgrids are disjunct, every process belongs to maximally one subgrid. Processes not
    ! belonging to any subgrids (e.g. if there are more than ngrid * nrow * ncol processes) will
    ! obtain an uninitialized grid descriptor (with iproc = -1) at return.
    allocate(imap(nrow, ncol))
    do igrid = 0, ngrid - 1
      ind = igrid * gridsize
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
      ctxt = ctxt0
      call blacs_gridmap(ctxt, imap, size(imap, dim=1), nrow, ncol)
      call blacs_gridinfo(ctxt, nrowtmp, ncoltmp, myrow, mycol)
      if (myrow >= 0 .and. mycol >= 0) self%ctxt = ctxt
    end do
    if (self%ctxt /= -1) then
      call self%setup_after_gridinit(colmajor0, leadrow0, leadcol0, repeatable)
    end if

    ! Return here if leadgrid is not asked for.
    if (.not. present(leadgrid)) return
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
    call leadgrid%setup_after_gridinit(.false., 0, 0, repeatable)

  end subroutine initsplitgrids


  !> Creates subgrids by explicitly mapping specified processes into a subgrid.
  !!
  !! Sets up a independent rectangular processor grids. All processes of the system context must
  !! call this routine collectively, passing the same gridmap. Processes not contained in the
  !! gridmap will obtain an uninitialized grid descriptor (with iproc = -1) at return.
  !!
  !! \param self  BLACS grid descriptor
  !! \param gridmap  Ids of the processes (relative to the system context) forming the grid. The
  !!     shape of the array determines the shape of the grid.
  !! \param context  BLACS system context (default: default system context)
  !! \param repeatable  If present and T, forces topologies to be repeatable.
  !!     May degrade performance in this case.
  !! \param info  If present, it is set to 0 on success and to a non-zero value on error (and
  !!     the grid is left uninitialized). Otherwise if absent the code stops on error.
  !!
  subroutine initmappedgrids(self, gridmap, context, repeatable, info)
    class(blacsgrid), intent(inout) :: self
    integer, intent(in) :: gridmap(:,:)
    integer, intent(in), optional :: context
    logical, intent(in), optional :: repeatable
    integer, intent(out), optional :: info

    integer :: ctxt0, nproc0
    integer :: ncol, nrow
    character(200) :: msg

    if (present(info)) then
      info = 0
    end if
    call self%reset()

    nrow = size(gridmap, dim=1)
    ncol = size(gridmap, dim=2)

    call get_syscontext(context, ctxt0, nproc0)
    if (nproc0 > 0 .and. any(gridmap >= nproc0)) then
      write(msg, "(A,I0,A)") "Process id in grid map exceeds number of processes (",&
          & nproc0, ")"
      call handle_error(trim(msg), info)
      return
    end if

    self%ctxt = ctxt0
    call blacs_gridmap(self%ctxt, gridmap, size(gridmap, dim=1), nrow, ncol)
    call self%setup_after_gridinit(.false., 0, 0, repeatable)

  end subroutine initmappedgrids


  !> Destroys processor grid.
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
    self%leadproc = -1
    self%leadrow = -1
    self%leadcol = -1
    self%lead = .false.

  end subroutine reset


  !> Fills up the descriptor after the grid context had been created by BLACS.
  !!
  !! Membership in the grid is determined via blacs_gridinfo, which works for any system context.
  !! Processes not being part of the grid obtain a reset descriptor.
  !!
  !! \param self  BLACS grid descriptor with the context set by blacs_gridinit/blacs_gridmap.
  !! \param colmajor  Whether process ids should be enumerated in column major order.
  !! \param leadrow  Row of the lead process.
  !! \param leadcol  Column of the lead process.
  !! \param repeatable  If present and T, forces topologies to be repeatable.
  !!
  subroutine setup_after_gridinit(self, colmajor, leadrow, leadcol, repeatable)
    class(blacsgrid), intent(inout) :: self
    logical, intent(in) :: colmajor
    integer, intent(in) :: leadrow, leadcol
    logical, intent(in), optional :: repeatable

    integer :: val(1)

    call blacs_gridinfo(self%ctxt, self%nrow, self%ncol, self%myrow, self%mycol)
    if (self%myrow < 0 .or. self%mycol < 0) then
      call self%reset()
      return
    end if

    if (present(repeatable)) then
      val = 0
      if (repeatable) val = 1
      call blacs_set(self%ctxt, 15, val)
    end if

    self%nproc = self%nrow * self%ncol
    if (colmajor) then
      self%iproc = self%mycol * self%nrow + self%myrow
    else
      self%iproc = self%myrow * self%ncol + self%mycol
    end if
    self%leadrow = leadrow
    self%leadcol = leadcol
    self%leadproc = blacs_pnum(self%ctxt, self%leadrow, self%leadcol)
    self%lead = (self%myrow == self%leadrow .and. self%mycol == self%leadcol)

  end subroutine setup_after_gridinit


  !> Returns the system context to use.
  !!
  !! \param context  Optional system context passed by the caller.
  !! \param ctxt0  System context to use (context if present, the default system context otherwise).
  !! \param nproc0  Nr. of processes in the system context, if known (default system context),
  !!     -1 otherwise.
  !!
  subroutine get_syscontext(context, ctxt0, nproc0)
    integer, intent(in), optional :: context
    integer, intent(out) :: ctxt0, nproc0

    integer :: iproc

    if (present(context)) then
      ctxt0 = context
      nproc0 = -1
    else
      call blacs_get(-1, 0, ctxt0)
      call blacs_pinfo(iproc, nproc0)
    end if

  end subroutine get_syscontext


  !> Handles an error in the grid creating routines.
  !!
  !! \param msg  Error message.
  !! \param info If present, it is set to a non-zero value, otherwise if absent the code is stopped.
  !!
  subroutine handle_error(msg, info)
    character(*), intent(in) :: msg
    integer, intent(out), optional :: info

    if (present(info)) then
      info = -1
    else
      call error(msg)
    end if

  end subroutine handle_error


end module blacsgrid_module
