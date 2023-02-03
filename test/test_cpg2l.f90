!> Testing rank one updates.
program test_cpg2l
  use, intrinsic :: iso_fortran_env, stdout => output_unit
  use test_common_module
  use libscalapackfx_module
  implicit none


  ! Block size (using an extremely small value for test purposes)
  integer, parameter :: bsize = 2

  call main()

contains

  subroutine main()
    type(blacsgrid) :: grid1, grid2

    integer :: nprow, npcol, iproc, nproc

    ! Initialize blas and create a square processor grid
    call blacsfx_pinfo(iproc, nproc)
    do nprow = int(sqrt(real(nproc, dp))), nproc
      if (mod(nproc, nprow) == 0) then
        exit
      end if
    end do
    npcol = nproc / nprow

    call grid1%initgrid(nprow, npcol)
    if (grid1%lead) then
      write(stdout, "(A,2(1X,I0))") "# processor grid:", nprow, npcol
    end if

    call grid2%initgrid(1, nproc)
    if (grid2%lead) then
      write(stdout, "(A,2(1X,I0))") "# processor grid:", 1, nproc
    end if

    if (.not. readMatrixAndTest(grid1, 2, 2, 1, 2)) then
      write(stdout, "(A)") "Test 1 failed"
    end if
    if (.not. readMatrixAndTest(grid1, 5, 5, 1, 1)) then
      write(stdout, "(A)") "Test 2 failed"
    end if
    if (.not. readMatrixAndTest(grid2, 2, 2, 1, 2)) then
      write(stdout, "(A)") "Test 3 failed"
    end if
    if (.not. readMatrixAndTest(grid2, 3, 5, 3, 1)) then
      write(stdout, "(A)") "Test 4 failed"
    end if

    call grid1%destruct()
    call grid1%initgrid(nproc, 1)
    if (.not. readMatrixAndTest(grid1, 2, 2, 1, 2)) then
      write(stdout, "(A)") "Test 5 failed"
    end if

    ! Finish blacs.
    call blacsfx_exit()

  end subroutine main

  function readMatrixAndTest(mygrid, iSize, jSize, i0, j0) result(success)
    type(blacsgrid), intent(inout) :: mygrid
    integer, intent(in) :: iSize, jSize, i0, j0
    logical :: success

    real(dp), allocatable :: glob(:,:), localTest(:,:), localRef(:,:)
    integer :: desc(DLEN_)
    integer :: mm, nn, i, j
    integer :: iloc, jloc, prow, pcol

    ! Read in matrix from disc.
    call readfromfile(mygrid, "hamsqr1.dat", bsize, bsize, glob, desc)
    mm = desc(M_)
    nn = desc(N_)
    if (mygrid%lead) then
      write(stdout, "(A,2(1X,I0))") "# global matrix size:", mm, nn
      write(stdout, "(A,2(1X,I0))") "# local matrix size on leader:",&
          & size(glob, dim=1), size(glob, dim=2)
    end if

    allocate(localRef(iSize,jSize), localTest(iSize,jSize))

    localRef(:,:) = 0.0_dp
    do j = 1, jSize
      do i = 1, iSize
        call scalafx_infog2l(mygrid, desc, i + i0 - 1, j + j0 - 1, iloc, jloc,&
            & prow, pcol)
        if (prow == mygrid%myrow .and. pcol == mygrid%mycol) then
          localRef(i, j) = glob(iloc, jloc)
        end if
      end do
    end do

    localTest(:,:) = 0.0_dp
    call scalafx_cpg2l(mygrid, desc, i0, j0, glob, localTest)

    success = all(localTest == localRef)

  end function readMatrixAndTest


end program test_cpg2l

