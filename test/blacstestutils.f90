module blacstestutils
  use libscalapackfx_module, only : blacsgrid, blacsfx_exit, blacsfx_pinfo
  use fortuno_mpi, only : test_item, mpi_case, check => mpi_check, failed => mpi_failed, num_ranks
  implicit none

  private
  public :: this_proc, num_procs
  public :: blacs_test
  public :: blacs_grid_env, get_grid_or_fail


  !> Implements a test class with BLACS initialization and destruction
  type, extends(mpi_case) :: blacs_case
  contains
    procedure :: run => blacs_case_run
  end type blacs_case


  abstract interface
    !> Interface of the test procedure
    subroutine blacs_test_procedure()
    end subroutine blacs_test_procedure
  end interface


  !> Implements a BLACS grid wrapper enforcing grid finalization
  type, extends(blacsgrid) :: blacs_grid_env

    !> Whether the grid contains all BLACS processes
    logical :: has_all_procs = .false.
  contains
    final :: final_blacs_grid_env
  end type blacs_grid_env


  ! Number of processes available in the BLACS framework
  integer :: num_procs_ = -1

  ! The id of the current process in the BLACS framework
  integer :: this_proc_ = -1

contains


  !> Returns the id of the current process in the BLACS framework
  function this_proc()
    integer :: this_proc
    this_proc = this_proc_
  end function this_proc


  !> Returns the number for processes available in the BLACS framework
  function num_procs()
    integer :: num_procs
    num_procs = num_procs_
  end function num_procs


  !> Wraps a blacs_case instance as test_item suitable for array constructors.
  function blacs_test(name, proc) result(testitem)
    character(*), intent(in) :: name
    procedure(blacs_test_procedure) :: proc

    type(test_item) :: testitem

    testitem = test_item(blacs_case(name=name, proc=proc))

  end function blacs_test


  !> Run procedure of the tempfile_case type.
  subroutine blacs_case_run(this)
    class(blacs_case), intent(in) :: this

    call blacsfx_pinfo(this_proc_, num_procs_)
    call check(num_procs_ == num_ranks(),&
        & "Number of BLACS processes differ from number of MPI ranks")
    if (failed()) return
    call this%proc()
    call blacsfx_exit(keepmpi=.true.)
    this_proc_ = -1
    num_procs_ = -1

  end subroutine blacs_case_run


  !> Returns a grid environment or sets the calling test to failed if not possible
  !!
  !! Note: This routine must be called from within fortuno MPI test procedures.
  !!
  subroutine get_grid_or_fail(this, nrow, ncol, includeall)

    !> Instance
    type(blacs_grid_env), intent(out) :: this

    !> Number of process rows
    integer, optional, intent(in) :: nrow

    !> Number of process columns
    integer, optional, intent(in) :: ncol

    !> Whether it should be ensured that all processes are included in the grid (default: .true.)
    logical, optional, intent(in) :: includeall

    integer :: nprocs
    integer :: nrow_, ncol_
    logical :: includeall_, hasall
    type(blacsgrid) :: grid

    includeall_ = .true.
    if (present(includeall)) includeall_ = includeall

    nprocs = num_procs()
    if (present(nrow) .and. present(ncol)) then
      nrow_ = nrow
      ncol_ = ncol
    else if (present(nrow)) then
      nrow_ = nrow
      ncol_ = nprocs / nrow_
    else if (present(ncol)) then
      ncol_ = ncol
      nrow_ = nprocs / ncol_
    else if (includeall_) then
      do nrow_ = floor(sqrt(real(nprocs))), 1, -1
        ncol_ = nprocs / nrow_
        if (ncol_ * nrow_ == nprocs) exit
      end do
    else
      nrow_ = floor(sqrt(real(nprocs)))
      ncol_ = nprocs / nrow_
    end if
    hasall = ncol_ * nrow_ == nprocs

    call check(nrow_ * ncol_ <= nprocs, msg="Required grid needs more processes than available")
    if (failed()) return
    call check(nrow_ > 0, msg="Could not set up grid with at least one process row")
    if (failed()) return
    call check(ncol_ > 0, msg="Could not set up grid with at least one process column")
    if (failed()) return
    call check(.not. includeall_ .or. hasall,&
        & msg="Could not include all processes in the required grid")
    if (failed()) return

    call this%blacsgrid%initgrid(nrow_, ncol_)
    this%has_all_procs = hasall

  end subroutine get_grid_or_fail


  subroutine final_blacs_grid_env(this)
    type(blacs_grid_env), intent(inout) :: this

    call this%blacsgrid%destruct()

  end subroutine final_blacs_grid_env

end module blacstestutils
