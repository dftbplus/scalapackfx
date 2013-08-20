!> Testing diagonalization.
!! \details 
!! Invoke it with the following parameters:
!!
!!   # number of processor rows
!!   # number of processor columns
!!   # block size
!!   # which diagonalizer to use (1 - QR, 2 - Divide and conquer)
!!   # file with the matrix to diagonalize
!!   # (Optional:) file with the matrix on the right hand side, if it is a
!!     generalized eigenvalue problem.
!!
!! The filess with the matrices should contain an integer number (n) in their
!! first line, indicating the dimension of the problem. Then n lines with
!! n entries should follow.
program test_diag
  use, intrinsic :: iso_fortran_env, stdout => output_unit
  use test_common_module
  use libscalapackfx_module
  implicit none

  call main()
  

contains

  subroutine main()
    type(blacsgrid) :: mygrid
    real(dp), allocatable :: hammtx(:,:), overmtx(:,:), eigvecs(:,:), eigvals(:)
    integer :: hamdesc(DLEN_), overdesc(DLEN_), eigvdesc(DLEN_)
    integer :: nn
    integer :: tsys1, tsys2, tcount
    real(dp) :: tcpu1, tcpu2
    character(100) :: hamfile, overfile
    logical :: overlap
    integer :: nprow, npcol, bsize, idiag
    integer :: iproc, nproc

    call getarguments(nprow, npcol, bsize, idiag, hamfile, overfile)
    overlap = len_trim(overfile) > 0

    call blacsfx_pinfo(iproc, nproc)
    if (npcol * nprow /= nproc) then
      write(stdout, *) "Incorrect number of processors:", npcol, nprow,&
          & nproc
      call blacsfx_exit(keepmpi=.false.)
      stop
    end if
    call mygrid%initgrid(nprow, npcol)
    if (mygrid%master) then
      write(stdout, "(A,2(1X,I0))") "# processor grid:", nprow, npcol
      write(stdout, "(A,1X,I0)") "# block size:", bsize
    end if

    call readfromfile(mygrid, hamfile, bsize, bsize, hammtx, hamdesc)
    if (overlap) then
      call readfromfile(mygrid, overfile, bsize, bsize, overmtx, overdesc)
    end if

    nn = hamdesc(M_)
    if (mygrid%master) then
      write(stdout, "(A,2(1X,I0))") "# matrix size:", nn, nn
      if (overlap) then
        write(stdout, "(A,1X,A)") "# eigenvalue problem:", "generalized"
      else
        write(stdout, "(A,1X,A)") "# eigenvalue problem:", "standard"
      end if
      write(stdout, "(A,2(1X,I0))") "# matrix size on master:",&
          & size(hammtx, dim=1), size(hammtx, dim=2)
      select case (idiag)
      case(1)
        write(stdout, "(A)") "# diagonalization: QR"
      case(2)
        write(stdout, "(A)") "# diagonalization: DAC"
      end select
    end if

    call system_clock(tsys1, tcount)
    call cpu_time(tcpu1)

    ! Allocate matrix for eigenvectors and eigenvals
    call scalafx_creatematrix(mygrid, nn, nn, bsize, bsize, eigvecs, eigvdesc)
    allocate(eigvals(nn))

    if (overlap) then
      select case (idiag)
      case(1)
        call scalafx_psygv(hammtx, hamdesc, overmtx, overdesc, eigvals,&
            & eigvecs, eigvdesc, jobz="V", uplo="L")
      case(2)
        call scalafx_psygvd(hammtx, hamdesc, overmtx, overdesc,&
            & eigvals, eigvecs, eigvdesc, jobz="V", uplo="L", &
            & allocfix=.true.)
      end select
    else
      select case (idiag)
      case(1)
        call scalafx_psyev(hammtx, hamdesc, eigvals, eigvecs, eigvdesc,&
            & jobz="V", uplo="L")
      case(2)
        call scalafx_psyevd(hammtx, hamdesc, eigvals, eigvecs,&
            & eigvdesc, jobz="V", uplo="L", allocfix=.true.)
      end select
    end if

    call cpu_time(tcpu2)
    call system_clock(tsys2)

    if (mygrid%master) then
      write(stdout, "(A,F8.1)") "# cpu time on master:", tcpu2 - tcpu1
      write(stdout, "(A,F8.1)") "# system time on master:",&
          & real(tsys2 - tsys1, dp) / real(tcount, dp)
    end if

    ! Write eigenvalues and eigenvectors
    if (mygrid%master) then
      open(12, file="eigvals.dat", form="formatted", status="replace")
      write(12, "(ES23.15)") eigvals
      close(12)
      write(stdout, "(A,A,A)") "# Eigenvalues written to '", "eigvals.dat", "'"
    end if
    call writetofile(mygrid, "eigvecs.dat", eigvecs, eigvdesc)
    if (mygrid%master) then
      write(stdout, "(A,A,A)") "# Eigenvectors written to '", "eigvecs.dat", "'"
    end if

    ! Destroy blacs communication layer.
    call blacsfx_exit()

  end subroutine main


  !> Process command line arguments.
  subroutine getarguments(nprow, npcol, bsize, idiag, hamfile, overfile)
    integer, intent(out) :: npcol, nprow, bsize, idiag
    character(*), intent(out) :: hamfile, overfile

    integer :: narg
    character(128) :: buffer

    narg = command_argument_count()
    if (narg /= 6 .and. narg /= 5) then
      write(stdout, *) "Incorrect number of arguments"
      write(stdout, *) "Need: nprow npcol bsize idiag hamiltonian [ overlap ]"
      stop
    end if
    call get_command_argument(1, buffer)
    read(buffer, *) nprow
    call get_command_argument(2, buffer)
    read(buffer, *) npcol
    call get_command_argument(3, buffer)
    read(buffer, *) bsize
    call get_command_argument(4, buffer)
    read(buffer, *) idiag
    if (idiag /= 1 .and. idiag /= 2) then
      stop "idiag must be one or two"
    end if
    call get_command_argument(5, hamfile)
    if (narg == 6) then
      call get_command_argument(6, overfile)
    else
      overfile = ""
    end if

  end subroutine getarguments


end program test_diag

