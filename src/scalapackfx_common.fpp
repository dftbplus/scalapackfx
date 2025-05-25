!> Helper routines for SCALAPACKFX.
!! \cond HIDDEN
module scalapackfx_common_module
  implicit none
  private

  public :: DLEN_
  public :: DT_, CTXT_, M_, N_, MB_, NB_, RSRC_, CSRC_, LLD_
  public :: sp, dp
  public :: error, handle_infoflag

  !> Descriptor length for distributed matrices
  integer, parameter :: DLEN_ = 9

  !> Index for data type.
  integer, parameter :: DT_ = 1

  !> Index for context.
  integer, parameter :: CTXT_ = 2

  !> Index for number of rows.
  integer, parameter :: M_ = 3

  !> Index for number of columns.
  integer, parameter :: N_ = 4

  !> Index for row block size.
  integer, parameter :: MB_ = 5

  !> Index for column block size.
  integer, parameter :: NB_ = 6

  !> Index for row of the source process.
  integer, parameter :: RSRC_ = 7

  !> Index for column of the source process.
  integer, parameter :: CSRC_ = 8

  !> Index for the leading dimension.
  integer, parameter :: LLD_ = 9

  !> Single precision kind.
  integer, parameter :: sp = kind(1.0)

  !> Double precision kind
  integer, parameter :: dp = kind(1.0d0)

contains


  !> Issues error and stops code execution.
  !!
  !! \param msg  Error message to print.
  !!
  subroutine error(msg)
    character(*), intent(in) :: msg

    write(*, "(A)") "!!! ERROR !!!"
    write(*, "(A)") msg
    stop

  end subroutine error


  !> Handles optional info flag.
  !!
  !! \param info0  Info flag as returned by some routine.
  !! \param msg  Msg to print out, if program is stopped.
  !! \param info  Optional info flag. If present, info0 is passed to it,
  !!     otherwise if info0 was not zero, the error message in msg is printed
  !!     and the program is stopped.
  !!
  subroutine handle_infoflag(info0, msg, info)
    integer, intent(in) :: info0
    character(*), intent(in) :: msg
    integer, intent(out), optional :: info

    if (present(info)) then
      info = info0
    elseif (info0 /= 0) then
      write(*, "(A)") "Scalapack operation failed!"
      write(*, "(A)") msg
      write(*, "(A,I0)") "Info: ", info0
      stop
    end if

  end subroutine handle_infoflag


end module scalapackfx_common_module

!> \endcond
