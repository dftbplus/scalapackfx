!> Test app driving Fortuno unit tests.
module test_scalapackfx
  use libscalapackfx_module
  use fortuno_mpi, only : check => mpi_check, failed => mpi_failed, skip => mpi_skip, test_list
  use blacstestutils, only : blacs_grid_env, get_grid_or_fail, test => blacs_test,&
      & num_procs
  implicit none

contains


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        test("desc", test_desc)&
    &])

  end function tests


  subroutine test_desc()
    type(blacs_grid_env) :: env
    type(blocklist) :: blist
    integer :: desc(DLEN_)
    integer :: ii, iglob, iloc, nelem

    if (num_procs() /= 4) then
      call skip()
      return
    end if

    call get_grid_or_fail(env, nrow=2, ncol=2)
    if (failed()) return

    if (env%mycol == -1) then
      print *, "IDLE Node:", env%iproc, " Exiting..."
    else
      call scalafx_getdescriptor(env%blacsgrid, 5, 5, 2, 2, desc)
      call blist%init(env%blacsgrid, desc, "c")
      do ii = 1, size(blist)
        call blist%getblock(ii, iglob, iloc, nelem)
        print "(6(A,I3,2X))", "PROW:", env%myrow, "PCOL:", env%mycol,&
            & "II:", ii, "GLOB:", iglob, "LOC:", iloc, "NEL:", nelem
      end do
    end if

    call check(.false.)

  end subroutine test_desc

end module test_scalapackfx


program testapp
  use fortuno_mpi, only : execute_mpi_cmd_app
  use test_scalapackfx, only : tests
  implicit none

  call execute_mpi_cmd_app(tests())

end program testapp
