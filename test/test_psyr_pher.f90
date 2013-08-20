!> Testing rank one updates.
program test_psyr_pher
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
    integer :: descx(DLEN_), descres(DLEN_)
    real(dp) :: alpha
    integer :: nprow, npcol, mm, nn, iproc, nproc
    integer :: ii

    ! Initialize blas and create a square processor grid
    call blacsfx_pinfo(iproc, nproc)
    do nprow = int(sqrt(real(nproc, dp))), nproc
      if (mod(nproc, nprow) == 0) then
        exit
      end if
    end do
    npcol = nproc / nprow
    call mygrid%initgrid(nprow, npcol)
    if (mygrid%master) then
      write(stdout, "(A,2(1X,I0))") "# processor grid:", nprow, npcol
      write(stdout, "(A,1X,I0)") "# block size:", bsize
    end if

    ! Read in matrix from disc.
    call readfromfile(mygrid, "hamsqr1.dat", bsize, bsize, xx, descx)
    mm = descx(M_)
    nn = descx(N_)
    if (mygrid%master) then
      write(stdout, "(A,2(1X,I0))") "# global matrix size:", mm, nn
      write(stdout, "(A,2(1X,I0))") "# local matrix size on master:",&
          & size(xx, dim=1), size(xx, dim=2)
    end if

    ! Do rank one update with all column vectors of the matrix.
    call scalafx_creatematrix(mygrid, mm, mm, bsize, bsize, res, descres)
    res(:,:) = 0.0_dp
    do ii = 1, mm
      alpha = real(ii, dp)
      call pblasfx_psyr(xx, descx, res, descres, alpha=alpha, ix=1, jx=ii,&
          & incx=1)
    end do

    ! Write results to disc.
    call writetofile(mygrid, "psyr_result.dat", res, descres)
    if (mygrid%master) then
      write(stdout, "(A)") "Result written to file 'psyr_result.dat'."
    end if

    ! Finish blacs.
    call blacsfx_exit()

  end subroutine main


end program test_psyr_pher

