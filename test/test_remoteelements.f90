!> Testing copy between matrix patterns
program test_remoteelements
  use, intrinsic :: iso_fortran_env, stdout => output_unit, stdin => input_unit
  use test_common_module
  use libscalapackfx_module
  use blacsfx_module
  implicit none

  call main()

contains

  subroutine main()
    type(blacsgrid) :: myGrid
    integer, allocatable :: AA(:,:), AAreceive(:,:), All(:,:)
    integer :: desc(DLEN_), nprow, npcol, mm, nn, mb, nb, iproc, nproc, ii, jj, kk, iGlob, jGlob
    integer :: iprow, ipcol, iLoc, jLoc, locrow, loccol
    logical :: isLocal

    ! grid information
    call blacsfx_pinfo(iproc, nproc)

    ! approximately square grid for group
    do npcol = floor(sqrt(real(nproc, dp))), nproc
      if (mod(nproc, npcol) == 0) then
        exit
      end if
    end do
    nprow = nproc / npcol
    call myGrid%initgrid(nprow, npcol)

    if (myGrid%lead) then
      write(stdout, "(A,2(1X,I0))") "# source processor grid:", nprow, npcol
      write(stdOut, "(A)")'# Matrix size to re-distrbute?'
      read(stdin, *) mm, nn
      write(stdout,"(A)")'# Block sizes of matrix?'
      read(stdin,*) mb, nb
      call blacsfx_gebs(myGrid, mm)
      call blacsfx_gebs(myGrid, nn)
      call blacsfx_gebs(myGrid, mb)
      call blacsfx_gebs(myGrid, nb)
    else
      call blacsfx_gebr(myGrid, mm)
      call blacsfx_gebr(myGrid, nn)
      call blacsfx_gebr(myGrid, mb)
      call blacsfx_gebr(myGrid, nb)
    end if

    call scalafx_creatematrix(myGrid, mm, nn, mb, nb, AA, desc)

    if (myGrid%lead) then
      do ii = 0, nproc-1
        call scalafx_getremoteshape(myGrid, desc, ii, locrow, loccol)
        write(stdout, "(A,I0,A,1X,I0,1X,I0)") "# A block (data:",ii,"):", locrow, loccol
        call scalafx_getremoteshape(myGrid, desc, ii, locrow, loccol, .true.)
        write(stdout, "(T8,A,I0,A,1X,I0,1X,I0)") "(storage:",ii,"):", locrow, loccol
      end do

      do ii = 1, nn
        do jj = 1, mm
          call scalafx_infog2p(mygrid, desc, jj, ii, iprow, ipcol)
          write(stdout, "(A,I0,1X,I0,A,I0,1X,I0,A,I0)") "# Matrix elements (",jj, ii,&
              & ') stored on BLACS grid at ', iprow, ipcol, " proc. ",&
              & blacsfx_pnum(myGrid, iprow, ipcol)
        end do
      end do

    end if

    call blacsfx_barrier(myGrid)

    AA(:,:) = -1
    call scalafx_getremoteshape(myGrid, desc, iProc, locrow, loccol, .true.)
    call blacsfx_pcoord(myGrid, iProc, iprow, ipcol)
    do ii = 1, loccol
      iGlob = scalafx_indxl2g(ii, desc(NB_), ipcol, desc(CSRC_), myGrid%ncol)
      do jj = 1, locrow
        jGlob = scalafx_indxl2g(jj, desc(MB_), iprow, desc(RSRC_), myGrid%nrow)
        AA(jj,ii) = jGlob + (iGlob-1) * mm
      end do
    end do
    write(stdout,*)iproc,':', AA(:locrow,:loccol)

    call blacsfx_barrier(myGrid)

    if (myGrid%lead) then

      allocate(All(mm,nn), source = 0)

      do kk = 0, nproc-1

        if (kk /=  myGrid%leadproc) then
          write(stdout,*)'Follow', kk

          call scalafx_getremoteshape(myGrid, desc, kk, locrow, loccol, .true.)
          allocate(AAreceive(locrow, loccol), source=-1)

          call blacsfx_pcoord(myGrid, kk, iprow, ipcol)
          write(*,*)kk, iprow, ipcol,'Expecting', shape(AAreceive)
          call blacsfx_gerv(myGrid, AAreceive, iprow, ipcol)

        else

          write(stdout,*)'Lead', kk
          call scalafx_getlocalshape(mygrid, desc, locrow, loccol)
          AAreceive = AA(:locrow,:loccol)

        end if

        call scalafx_getremoteshape(myGrid, desc, kk, locrow, loccol, .true.)
        call blacsfx_pcoord(myGrid, kk, iprow, ipcol)
        do ii = 1, loccol
          iGlob = scalafx_indxl2g(ii, desc(NB_), ipcol, desc(CSRC_), myGrid%ncol)
          do jj = 1, locrow
            jGlob = scalafx_indxl2g(jj, desc(MB_), iprow, desc(RSRC_), myGrid%nrow)
            write(stdout,*)'Proc',kk,',', iprow, ipcol,': (', jGlob,iGlob,')(',jj,ii,')'
            All(jGlob,iGlob) = AAreceive(jj,ii)
          end do
        end do

        deallocate(AAreceive)

      end do

    else

      write(*,*)iproc, myGrid%myrow, myGrid%mycol,'sending', shape(AA)
      call blacsfx_gesd(myGrid, AA, myGrid%leadrow, myGrid%leadcol)

    end if

    call blacsfx_barrier(myGrid)

    if (myGrid%lead) then
      do jj = 1, nn
        write(stdout,*)All(:,jj)
      end do
    end if

    ! Finish blacs.
    call blacsfx_exit()

  end subroutine main

end program test_remoteelements
