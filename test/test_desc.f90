!> Test app driving Fortuno unit tests.
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
        test("desc_square", test_desc_square),&
        test("desc_nonsquare", test_desc_nonsquare)&
    &])

  end function tests


  !> Test cases where all processors fit onto a square BLACS grid
  subroutine test_desc_square()
    type(blacs_grid_env) :: env
    type(blocklist) :: blist
    integer :: desc(DLEN_)
    integer :: ii, iGlob, iLoc, nElem, sqrtProcs, info

    sqrtProcs = gridSide(num_procs())
    if (num_procs() /= sqrtProcs**2) then
      call fortuno_skip()
      return
    end if

    print *, "Testing square grid of processors"

    call get_grid_or_fail(env, nrow=sqrtProcs, ncol=sqrtProcs)
    if (fortuno_failed()) return

    if (env%mycol == -1) then
      print *, "IDLE Node:", env%iproc, " Should be a square grid, so should not get here..."
      call fortuno_check(.false.)
    else
      print *, "WORKING Node:", env%iproc
      call scalafx_getdescriptor(env%blacsgrid, num_procs()+1, num_procs()+1, sqrtProcs, sqrtProcs,&
          & desc, info=info)
      if (info == 0) then
        call fortuno_check(.true.)
      else
        print *, "Node:", env%iproc, " ScaLAPACK descriptor returned info=", info
        call fortuno_check(.false.)
      end if
      call blist%init(env%blacsgrid, desc, "c")
      do ii = 1, size(blist)
        ! Loop over local blocks of the matrix of this proc.
        call blist%getblock(ii, iGlob, iLoc, nElem)
        print "(6(A,1X,I0,2X))", "PROW: ", env%myrow, "PCOL:", env%mycol, "II:", ii,&
            & "GLOB:", iGlob, "LOC:", iLoc, "NEL:", nElem
      end do

    end if

  end subroutine test_desc_square


  !> Test cases where not all processors neccessarily fit onto a square BLACS grid
  subroutine test_desc_nonsquare()
    type(blacs_grid_env) :: env
    type(blocklist) :: blist
    integer :: desc(DLEN_)
    integer :: ii, iGlob, iLoc, nElem, sqrtProcs, info

    print *, "Testing non-square grid of processors"

    sqrtProcs = gridSide(num_procs())

    call get_grid_or_fail(env, nrow=sqrtProcs, ncol=sqrtProcs, includeall=.false.)
    if (fortuno_failed()) return

    if (env%mycol == -1) then
      print *, "IDLE Node:", env%iproc, " Idle so exiting successfully..."
      call fortuno_check(.true.)
    else
      print *, "WORKING Node:", env%iproc
      call scalafx_getdescriptor(env%blacsgrid, num_procs()+1, num_procs()+1, sqrtProcs, sqrtProcs,&
          & desc, info=info)
      if (info == 0) then
        call fortuno_check(.true.)
      else
        print *, "Node:", env%iproc, " ScaLAPACK descriptor returned info=", info
        call fortuno_check(.false.)
      end if
      call blist%init(env%blacsgrid, desc, "c")
      do ii = 1, size(blist)
        ! Loop over local blocks of the matrix of this proc.
        call blist%getblock(ii, iGlob, iLoc, nElem)
        print "(6(A,1X,I0,2X))", "PROW: ", env%myrow, "PCOL:", env%mycol, "II:", ii,&
            & "GLOB:", iGlob, "LOC:", iLoc, "NEL:", nElem
      end do

    end if

  end subroutine test_desc_nonsquare


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
