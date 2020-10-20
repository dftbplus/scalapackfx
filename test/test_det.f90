!> Testing determinant evaluation via getrf
program test_pdet
  use, intrinsic :: iso_fortran_env, stdout => output_unit
  use test_common_module
  use libscalapackfx_module
  implicit none

  ! Block size (using an extremely small value for test purposes)
  integer, parameter :: bsize = 1

  call main()

contains

  subroutine main()
    type(blacsgrid) :: mygrid
    real(dp), allocatable :: xx(:,:)
    complex(dp), allocatable :: xxc(:,:)
    real(dp) :: det
    integer :: exponent
    complex(dp) :: detc
    integer :: descx(DLEN_)
    integer :: nprow, npcol, mm, nn, iproc, nproc, ii, jj, iLoc, jLoc
    integer, allocatable :: ipiv(:)
    logical :: tLocal

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

    allocate(ipiv(min(mm,nn)))
    ipiv = 0

    if (mygrid%lead) then
      write(stdout, "(A,2(1X,I0))") "# global matrix size:", mm, nn
      write(stdout, "(A,2(1X,I0))") "# local matrix size on leader:",&
          & size(xx, dim=1), size(xx, dim=2)
    end if

    call scalafx_pgetrf(xx,descx,ipiv)

    ! note, this includes under-/over-flow protection similar to LINPACK routine dgedi.f
    det = 1.0_dp
    exponent = 0
    do ii = 1, size(xx,dim=2)
      jj = scalafx_indxl2g(ii, descx(NB_), mygrid%mycol, descx(CSRC_), mygrid%ncol)
      call scalafx_islocal(mygrid, descx, jj, jj, tLocal, iLoc, jLoc)
      if (tLocal) then
        if (jj /= ipiv(ii)) then
          det = -det * xx(iLoc,jLoc)
        else
          det = det * xx(iLoc,jLoc)
        end if
        do while (abs(det) > 2)
          det = det / 2.0_dp
          exponent = exponent + 1
        end do
        do while (abs(det) < 0.5_dp)
          det = det * 2.0_dp
          exponent = exponent - 1
        end do
      end if
    end do
    det = det * 2.0_dp ** exponent

    ! Would normally accumulate product via mpi calls, but messy to do via blacs operations only
    write(stdout,*)'Det part from process', mygrid%iproc, ' : ', det

    ! Finish blacs.
    call blacsfx_exit()

  end subroutine main

end program test_pdet
