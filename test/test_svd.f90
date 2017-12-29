!> Testing SVD
!! \details
!! Invoke it with the following parameters:
!!
!!   # number of processor rows
!!   # number of processor columns
!!   # block size
!!   # file with the matrix to decompose
!!
!! The file with the matrix should contain two integer numbers (m,n) in its first line, indicating
!! the dimension of the problem. Then n lines with m entries should follow.
program test_diag
  use, intrinsic :: iso_fortran_env, stdout => output_unit
  use test_common_module
  use libscalapackfx_module
  implicit none

  call main()
  

contains

  subroutine main()
    type(blacsgrid) :: mygrid
    real(dp), allocatable :: matrix(:,:), U(:,:), sigma(:), vt(:,:)
    integer :: matrixdesc(DLEN_), Udesc(DLEN_), vtdesc(DLEN_)
    integer :: mm, nn
    integer :: tsys1, tsys2, tcount
    real(dp) :: tcpu1, tcpu2
    character(100) :: matfile
    integer :: nprow, npcol, bsize
    integer :: iproc, nproc

    call getarguments(nprow, npcol, bsize, matfile)

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

    call readfromfile(mygrid, matfile, bsize, bsize, matrix, matrixdesc)

    mm = matrixdesc(M_)
    nn = matrixdesc(N_)
    
    if (mygrid%master) then
      write(stdout, "(A,2(1X,I0))") "# matrix size:", mm, nn
      write(stdout, "(A,1X,A)") "# singular value decomposition:", "standard"
      write(stdout, "(A,2(1X,I0))") "# matrix size on master:",&
          & size(matrix, dim=1), size(matrix, dim=2)      
    end if

    call system_clock(tsys1, tcount)
    call cpu_time(tcpu1)

    ! Allocate matrices    
    call scalafx_creatematrix(mygrid, mm, mm, bsize, bsize, U, Udesc)
    call scalafx_creatematrix(mygrid, nn, nn, bsize, bsize, Vt, Vtdesc)
    allocate(sigma(min(mm,nn)))

    call scalafx_pgesvd(matrix, matrixdesc, U, Udesc, sigma, vt, vtdesc)

    call cpu_time(tcpu2)
    call system_clock(tsys2)

    if (mygrid%master) then
      write(stdout, "(A,F8.1)") "# cpu time on master:", tcpu2 - tcpu1
      write(stdout, "(A,F8.1)") "# system time on master:",&
          & real(tsys2 - tsys1, dp) / real(tcount, dp)
    end if

    ! Write singular values and vectors
    if (mygrid%master) then
      open(12, file="singularvals.dat", form="formatted", status="replace")
      write(12, "(ES23.15)") sigma
      close(12)
      write(stdout, "(A,A,A)") "# Singular values written to '", "singularvals.dat", "'"
    end if
    call writetofile(mygrid, "singularU.dat", U, Udesc)
    call writetofile(mygrid, "singularVt.dat", vT, Vtdesc)
    if (mygrid%master) then
      write(stdout, "(A,A,A)") "# Singular vectors written to '", "singularU/singularVt.dat", "'"
    end if

    ! Destroy blacs communication layer.
    call blacsfx_exit()

  end subroutine main


  !> Process command line arguments.
  subroutine getarguments(nprow, npcol, bsize, matfile)
    integer, intent(out) :: npcol, nprow, bsize
    character(*), intent(out) :: matfile

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
    call get_command_argument(4, matfile)

  end subroutine getarguments


end program test_diag
