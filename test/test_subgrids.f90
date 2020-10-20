program test_subgrid
  use libscalapackfx_module
  implicit none

  type(blacsgrid) :: grpproc, allproc, grpleaders
  integer :: iproc, nproc

  call blacsfx_pinfo(iproc, nproc)
  call allproc%initgrid(2, 4)
  call grpproc%initsplitgrids(2, 2, 2, leadgrid=grpleaders)
  print "(A,4(A,1X,I4,5X))", "ROW-MAJOR| ", "GLOBALID:", iproc, &
      & "GROUPID:", grpproc%iproc, "ROW:", grpproc%myrow,&
      &  "COL:", grpproc%mycol
  if (grpproc%iproc == -1) then
    print *, "IDLE:", iproc
  else
    if (grpleaders%iproc /= -1) then
      print *, "GRIDLEAD:", iproc, grpproc%iproc, grpleaders%iproc
      call blacsfx_barrier(grpleaders)
    else
      print *, "NORMAL:", iproc, grpproc%iproc
    end if
    call blacsfx_barrier(grpproc)
    call grpproc%destruct()
    if (grpleaders%iproc /= -1) then
      call grpleaders%destruct()
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
