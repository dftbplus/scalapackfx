!> Testing transposition
program test_ptran
  use, intrinsic :: iso_fortran_env, stdout => output_unit
  use test_common_module
  use libscalapackfx_module
  implicit none


  ! Block size (using an extremely small value for test purposes)
  integer, parameter :: bsize = 2

  call main()

contains

  subroutine main()
    type(blacsgrid) :: mygrid
    real(dp), allocatable :: xx(:,:), res(:,:)
    complex(dp), allocatable :: xxc(:,:), resc(:,:)
    integer :: descx(DLEN_), descres(DLEN_)
    integer :: nprow, npcol, mm, nn, iproc, nproc

    ! Initialize blas and create a square processor grid
    call blacsfx_pinfo(iproc, nproc)
    do nprow = int(sqrt(real(nproc, dp))), nproc
      if (mod(nproc, nprow) == 0) then
        exit
      end if
    end do
    npcol = nproc / nprow
    call mygrid%initgrid(nprow, npcol)
    if (mygrid%lead) then
      write(stdout, "(A,2(1X,I0))") "# processor grid:", nprow, npcol
      write(stdout, "(A,1X,I0)") "# block size:", bsize
    end if

    ! Set up matrix
    if (mygrid%lead) then
      write(stdout, "(A)") "Matrix read from file 'hamsqr1.dat'."
    end if
    call readfromfile(mygrid, "hamsqr1.dat", bsize, bsize, xx, descx)
    mm = descx(M_)
    nn = descx(N_)
    if (mygrid%lead) then
      write(stdout, "(A,2(1X,I0))") "# global matrix size:", mm, nn
      write(stdout, "(A,2(1X,I0))") "# local matrix size on leader:",&
          & size(xx, dim=1), size(xx, dim=2)
    end if

    call scalafx_creatematrix(mygrid, nn, mm, bsize, bsize, res, descres)
    res = 0.0_dp

    call pblasfx_ptran(xx,descx,res,descres)

    ! Write results to disc.
    call writetofile(mygrid, "ptran_realresult.dat", res, descres)
    if (mygrid%lead) then
      write(stdout, "(A)") "Result written to file 'ptran_realresult.dat'."
    end if

    allocate(xxc(size(xx,dim=1),size(xx,dim=2)))
    allocate(resc(size(res,dim=1),size(res,dim=2)))

    xxc = cmplx(0,1,dp) * xx + xx

    call pblasfx_ptranu(xxc,descx,resc,descres)

    ! Write results to disc.
    call writetofile(mygrid, "ptran_cmplxresult.dat", resc, descres)
    if (mygrid%lead) then
      write(stdout, "(A)") "Result written to file 'ptran_cmplxresult.dat'."
    end if

    call pblasfx_ptranc(xxc,descx,resc,descres)

    ! Write results to disc.
    call writetofile(mygrid, "ptran_cmplxHresult.dat", resc, descres)
    if (mygrid%lead) then
      write(stdout, "(A)") "Result written to file 'ptran_cmplxHresult.dat'."
    end if

    deallocate(resc)
    deallocate(xxc)

    ! Finish blacs.
    call blacsfx_exit()

  end subroutine main

end program test_ptran
