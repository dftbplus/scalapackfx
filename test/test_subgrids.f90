program test_subgrid
  use libscalapackfx_module
  implicit none

  type(blacsgrid) :: grpproc, allproc, grpmasters
  integer :: iproc, nproc

  call blacsfx_pinfo(iproc, nproc)
  call allproc%initgrid(2, 4)
  call grpproc%initsplitgrids(2, 2, 2, mastergrid=grpmasters)
  print "(A,4(A,1X,I4,5X))", "ROW-MAJOR| ", "GLOBALID:", iproc, &
      & "GROUPID:", grpproc%iproc, "ROW:", grpproc%myrow,&
      &  "COL:", grpproc%mycol
  if (grpproc%iproc == -1) then
    print *, "IDLE:", iproc
  else
    if (grpmasters%iproc /= -1) then
      print *, "GRIDMASTER:", iproc, grpproc%iproc, grpmasters%iproc
      call blacsfx_barrier(grpmasters)
    else
      print *, "NORMAL:", iproc, grpproc%iproc
    end if
    call blacsfx_barrier(grpproc)
    call grpproc%destruct()
    if (grpmasters%iproc /= -1) then
      call grpmasters%destruct()
    end if
  end if
  call grpproc%initsplitgrids(2, 2, 2, colmajor=.true.)
  print "(A,4(A,1X,I4,5X))", "COL-MAJOR| ", "GLOBALID:", iproc, &
      & "GROUPID:", grpproc%iproc, "ROW:", grpproc%myrow,&
      &  "COL:", grpproc%mycol
  if (grpproc%iproc == -1) then
    print *, "IDLE:", iproc
  else
    call blacsfx_barrier(grpproc)
  end if
  call blacsfx_exit()

end program test_subgrid
