!> Test app driving Fortuno unit tests.
module test_scalapackfx
  use libscalapackfx_module
  use fortuno_mpi, only : global_comm, is_equal, test => mpi_case_item, check => mpi_check,&
      & test_list, this_rank
  implicit none

contains

  function tests()
    type(test_list) :: tests

    tests = test_list([&
        test("desc", test_desc)&
    &])

  end function tests


  subroutine test_desc()
    type(blacsgrid) :: myGrid
    type(blocklist) :: block_List
    integer :: desc(DLEN_)
    integer :: ii, iGlob, iLoc, nElem
    integer :: iProc, nProc, nGrid

    call blacsfx_pinfo(iProc, nProc)

    nGrid = floor(sqrt(real(nProc)))

    call myGrid%initgrid(2, 2)

    call check(.false.)

    if (myGrid%mycol == -1) then

      print *, "IDLE Node:", iProc, " Exiting..."

    else

      call scalafx_getdescriptor(myGrid, 5, 5, 2, 2, desc)
      call block_List%init(myGrid, desc, "c")
      do ii = 1, size(block_List)
        call block_List%getblock(ii, iGlob, iLoc, nElem)
        print "(6(A,I3,2X))", "PROW:", myGrid%myrow, "PCOL:", myGrid%mycol,&
            & "II:", ii, "GLOB:", iGlob, "LOC:", iLoc, "NEL:", nElem
      end do

    end if

    ! THEN each rank must contain source rank's value
    !call check(is_equal(buffer, sourceval))
    call check(.false.)

    call blacsfx_exit()

  end subroutine test_desc

end module test_scalapackfx


program testapp
  use fortuno_mpi, only : execute_mpi_cmd_app
  use test_scalapackfx, only : tests
  implicit none

  call execute_mpi_cmd_app(tests())

end program testapp
