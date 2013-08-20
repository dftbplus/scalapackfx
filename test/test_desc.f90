program test_descriptors
  use libscalapackfx_module
  implicit none

  type(blacsgrid) :: mygrid
  type(blocklist) :: blist
  integer :: desc(DLEN_)
  integer :: ii, iglob, iloc, nelem
  integer :: iproc, nproc

  call blacsfx_pinfo(iproc, nproc)
  call mygrid%initgrid(2, 2)
  if (mygrid%mycol == -1) then
    print *, "IDLE Node:", iproc, " Exiting..."
    call blacsfx_exit()
    stop
  end if
  call scalafx_getdescriptor(mygrid, 5, 5, 2, 2, desc)
  call blist%init(mygrid, desc, "c")
  do ii = 1, size(blist)
    call blist%getblock(ii, iglob, iloc, nelem)
    print "(6(A,I3,2X))", "PROW:", mygrid%myrow, "PCOL:", mygrid%mycol,&
        & "II:", ii, "GLOB:", iglob, "LOC:", iloc, "NEL:", nelem
  end do
  call blacsfx_exit()
  
  
end program test_descriptors
