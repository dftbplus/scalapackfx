!> Test BLACS subgrids
module test_scalapackfx
  use libscalapackfx_module
  use fortuno_mpi, only : fortuno_check => mpi_check, fortuno_failed => mpi_failed,&
      & fortuno_skip => mpi_skip, test_list
  use blacstestutils, only : blacs_grid_env, get_grid_or_fail, test => blacs_test,&
      & num_procs
  implicit none

contains

  !> Tests to run
  function tests()
    type(test_list) :: tests

    tests = test_list([&
        test("desc_subgrids", test_subgrids)&
    &])

  end function tests


  !> Test splitting up a BLACS grid
  subroutine test_subgrids()

    type(blacs_grid_env) :: allProc, grpLeaders, grpProc
    integer :: iProc, nProc, sqrtProcs

    sqrtProcs = gridSide(num_procs())

    call get_grid_or_fail(allProc, nrow=sqrtProcs, ncol=sqrtProcs, includeall=.false.)
    if (fortuno_failed()) return

    call blacsfx_pinfo(iProc, nProc)

    call fortuno_check(iProc == allProc%iProc)
    call fortuno_check(nProc == allProc%nProc)
    print *, "iProc=", iProc, " nProc=", nProc

    call grpproc%initsplitgrids(sqrtProcs, sqrtProcs, sqrtProcs, leadgrid=grpLeaders)

    print "(A,4(A,1X,I4,5X))", "ROW-MAJOR| ", "GLOBAL ID:", iProc,&
        & "GROUP ID:", grpProc%iProc, "ROW:", grpProc%myrow,&
        & "COL:", grpProc%mycol

    if (grpproc%iproc == -1) then
      print *, "IDLE:", iproc
    else
      if (grpLeaders%iproc /= -1) then
        print *, "GRIDLEAD:", iproc, grpproc%iproc, grpLeaders%iproc
        !call blacsfx_barrier(grpLeaders)
        call grpLeaders%destruct()
      else
        print *, "NORMAL:", iproc, grpproc%iproc
      end if

      !call blacsfx_barrier(grpproc)
      call grpproc%destruct()

    end if

    call grpproc%initsplitgrids(sqrtProcs, sqrtProcs, sqrtProcs, colmajor=.true.)
    print "(A,4(A,1X,I4,5X))", "COL-MAJOR| ", "GLOBALID:", iproc, &
        & "GROUPID:", grpproc%iproc, "ROW:", grpproc%myrow,&
        & "COL:", grpproc%mycol
    if (grpproc%iproc == -1) then
      print *, "IDLE:", iproc
    else
      !call blacsfx_barrier(grpproc)
    end if

    ! call blacsfx_exit()

    call fortuno_check(.false.)

  end subroutine test_subgrids


  !> Calculate side of a square grid containing nProc processors
  function gridSide(nProc)

    !> Number of processors
    integer, intent(in) :: nProc

    !> sqrt(nProc)
    integer :: gridSide

    gridSide = floor(sqrt(real(nProc) + epsilon(0.0)))

  end function gridSide

end module test_scalapackfx


program testapp
  use fortuno_mpi, only : execute_mpi_cmd_app
  use test_scalapackfx, only : tests
  implicit none

  call execute_mpi_cmd_app(tests())

end program testapp
