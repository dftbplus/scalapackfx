!> Testing communication of a line from a distributed matrix.
!! \details 
!! Invoke it with the following parameters:
!!
!!   # number of processor rows
!!   # number of processor columns
!!   # block size
!!   # file with the matrix
!!
!! The file with the matrix should contain an integer number (n) in the
!! first line, indicating the dimension of the problem. Then n lines with
!! n entries should follow.
program test_linecomm
  use, intrinsic :: iso_fortran_env, stdout => output_unit
  use test_common_module
  use libscalapackfx_module
  use linecomm_module
  implicit none

  call main()
  

contains

  subroutine main()
    type(blacsgrid) :: mygrid
    real(dp), allocatable :: matrix(:,:)
    integer :: matrixdesc(DLEN_)
    integer :: nn
    character(100) :: matrixfile
    integer :: nprow, npcol, icol, bsize
    integer :: iproc, nproc
    type(linecomm) :: collector
    real(dp), allocatable :: iobuffer(:)

    call getarguments(nprow, npcol, bsize, matrixfile)

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

    call readfromfile(mygrid, matrixfile, bsize, bsize, matrix, matrixdesc)

    nn = matrixdesc(M_)
    if (mygrid%master) then
      write(stdout, "(A,2(1X,I0))") "# matrix size:", nn, nn
      write(stdout, "(A,2(1X,I0))") "# matrix size on master:",&
          & size(matrix, dim=1), size(matrix, dim=2)
    end if

    allocate(iobuffer(nn))

    call collector%init(mygrid, matrixdesc, "c")

    do icol = 1, matrixdesc(N_)
      if (mygrid%master) then
        call collector%getline_master(mygrid, icol, matrix, iobuffer)
        write(*, *) iobuffer(:)
      else
        call collector%getline_slave(mygrid, icol, matrix)
      end if
    end do

  end subroutine main


  !> Process command line arguments.
  subroutine getarguments(nprow, npcol, bsize, matrixfile)
    integer, intent(out) :: npcol, nprow, bsize
    character(*), intent(out) :: matrixfile

    integer :: narg
    character(128) :: buffer

    narg = command_argument_count()
    if (narg /= 4) then
      write(stdout, *) "Incorrect number of arguments"
      write(stdout, *) "Need: nprow npcol bsize matrix"
      stop
    end if
    call get_command_argument(1, buffer)
    read(buffer, *) nprow
    call get_command_argument(2, buffer)
    read(buffer, *) npcol
    call get_command_argument(3, buffer)
    read(buffer, *) bsize
    call get_command_argument(4, matrixfile)

  end subroutine getarguments


end program test_linecomm
