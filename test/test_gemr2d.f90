!> Testing copy between matrix patterns
program test_gemr2d
  use, intrinsic :: iso_fortran_env, stdout => output_unit, stdin => input_unit
  use test_common_module
  use libscalapackfx_module
  use blacsfx_module
  implicit none

  ! Block sizes for matrix A (using an extremely small value for test purposes)
  integer, parameter :: mbA = 2
  integer, parameter :: nbA = 2

  call main()

contains

  ! number up the elements in the global matrix
  pure function numelem(iGlob, jGlob, mm)
    integer, intent(in) :: iGlob, jGlob, mm
    integer :: numelem
    numelem = jGlob + (iGlob-1) * mm
  end function numelem

  subroutine main()
    type(blacsgrid) :: myGridA, myGridB
    integer, allocatable :: AA(:,:), BB(:,:)
    integer :: descA(DLEN_), descB(DLEN_)
    integer :: nprow, npcol, mm, nn, iproc, nproc, ii, jj, iGlob, jGlob, iprow, ipcol
    integer, allocatable :: iErr(:)

    ! Block size for matrix B
    integer :: mbB
    integer :: nbB

    ! grid information
    call blacsfx_pinfo(iproc, nproc)

    ! approximately square grid for group A
    do npcol = int(sqrt(real(nproc, dp))), nproc
      if (mod(nproc, npcol) == 0) then
        exit
      end if
    end do
    nprow = nproc / npcol
    call myGridA%initgrid(nprow, npcol)
    call myGridB%initgrid(1, nproc)
    if (myGridA%lead) then
      write(stdout, "(A,2(1X,I0))") "# source processor grid:", nprow, npcol
    end if
    if (myGridB%lead) then
      write(stdout, "(A,2(1X,I0))") "# destination processor grid:", 1, nproc
    end if

    call blacsfx_barrier(myGridA)
    call blacsfx_barrier(myGridB)

    if (myGridA%lead) then
      write(stdOut, *)'# Matrix size to re-distrbute?'
      read(stdin, *) mm, nn
      call blacsfx_gebs(myGridA, mm)
      call blacsfx_gebs(myGridA, nn)
    else
      call blacsfx_gebr(myGridA, mm)
      call blacsfx_gebr(myGridA, nn)
    end if

    ! Block size for matrix B
    mbB = mm
    nbB = 1

    if (myGridA%lead) then
      write(stdout, "(A,2(1X,I0))") "# A processor grid:", nprow, npcol
      write(stdout, "(A,1X,I0,1X,I0)") "# A block size:", mbA, nbA
    end if
    if (myGridB%lead) then
      write(stdout, "(A,2(1X,I0))") "# B processor grid:", nprow, npcol
      write(stdout, "(A,1X,I0,1X,I0)") "# B block size:", mbB, nbB
    end if

    call scalafx_creatematrix(myGridA, mm, nn, mbA, nbA, AA, descA)
    call scalafx_creatematrix(myGridB, mm, nn, mbB, nbB, BB, descB)

    if (myGridA%lead) then
      do ii = 0, nproc-1
        call blacsfx_pcoord(myGridA, ii, iprow, ipcol)
        write(stdout, "(A,I0,A,1X,I0,1X,I0)") "# A block (proc:",ii,"):",&
            & max(1, scalafx_numroc(descA(M_), descA(MB_), iprow, descA(RSRC_), myGridA%nrow)),&
            & scalafx_numroc(descA(N_), descA(NB_), ipcol, descA(CSRC_), myGridA%ncol)
      end do

      write(stdout,*)'Processor locations in the A BLACS grid'
      do ii = 0, nproc - 1
        call blacsfx_pcoord(myGridA, ii, iprow, ipcol)
        write(stdout, "(X,A,I0,A,I0,X,I0)")'Proc ', ii, ' is at BLACS grid loc. ', iprow, ipcol
      end do

      write(stdout,*)'Matrix element locations in the A BLACS grid'
      do ii = 1, nn
        do jj = 1, mm
          iprow = scalafx_indxg2p(jj, descA(MB_), descA(RSRC_), myGridA%nrow)
          ipcol = scalafx_indxg2p(ii, descA(NB_), descA(CSRC_), myGridA%ncol)
          write(stdout, "(I0,X,I0,A,I0,X,I0,A,I0)")jj, ii, ' is at BLACS location ', ipcol, iprow,&
              & ' proc:', blacsfx_pnum(myGridA, iprow, ipcol)
        end do
      end do
    end if

    call blacsfx_barrier(myGridA)
    call blacsfx_barrier(myGridB)

    if (myGridB%lead) then
      do ii = 0, nproc - 1
        call blacsfx_pcoord(myGridB, ii, iprow, ipcol)
        write(stdout, "(A,I0,A,1X,I0,1X,I0)") "# B block (proc:",ii,"):",&
            & max(1, scalafx_numroc(descB(M_), descB(MB_), iprow, descB(RSRC_), myGridB%nrow)),&
            & scalafx_numroc(descB(N_), descB(NB_), ipcol, descB(CSRC_), myGridB%ncol)
      end do

      write(stdout,*)'Processor locations in the B BLACS grid'
      do ii = 0, nproc - 1
        call blacsfx_pcoord(myGridB, ii, iprow, ipcol)
        write(stdout, "(X,A,I0,A,I0,X,I0)")'Proc ', ii, ' is at BLACS grid loc. ', iprow, ipcol
      end do

      write(stdout,*)'Matrix element locations in the B BLACS grid'
      do ii = 1, nn
        do jj = 1, mm
          iprow = scalafx_indxg2p(jj, descB(MB_), descB(RSRC_), myGridB%nrow)
          ipcol = scalafx_indxg2p(ii, descB(NB_), descB(CSRC_), myGridB%ncol)
          write(stdout, "(I0,X,I0,A,I0,X,I0,A,I0)")jj, ii, ' is at BLACS location ', ipcol, iprow,&
              & ' proc:', blacsfx_pnum(myGridB, iprow, ipcol)
        end do
      end do

    end if

    call blacsfx_barrier(myGridA)
    call blacsfx_barrier(myGridB)

    AA(:,:) = 0
    BB(:,:) = 0

    write(stdOut,"(A,I0,A,I0,',',I0,A,I0,X,I0)")'Processor ', myGridA%iproc,' (grid A) holds A(',&
        & shape(AA), ') elements vs ',&
        & max(1, scalafx_numroc(descA(M_), descA(MB_), myGridA%myrow, descA(RSRC_), myGridA%nrow)),&
        & scalafx_numroc(descA(N_), descA(NB_), myGridA%mycol, descA(CSRC_), myGridA%ncol)

    write(stdOut,"(A,I0,A,I0,',',I0,A,I0,X,I0)")'Processor ', myGridB%iproc,' (grid B) holds B(',&
        & shape(BB), ') elements vs ',&
        & max(1, scalafx_numroc(descB(M_), descB(MB_), myGridB%myrow, descB(RSRC_), myGridB%nrow)),&
        & scalafx_numroc(descB(N_), descB(NB_), myGridB%mycol, descB(CSRC_), myGridB%ncol)

    ! number off the elements in A
    do ii = 1, size(AA, dim=2)
      iGlob = scalafx_indxl2g(ii, descA(NB_), myGridA%mycol, descA(CSRC_), myGridA%ncol)
      do jj = 1, size(AA, dim=1)
        jGlob = scalafx_indxl2g(jj, descA(MB_), myGridA%myrow, descA(RSRC_), myGridA%nrow)
        AA(jj,ii) = numelem(iGlob, jGlob, mm)
      end do
    end do

    call blacsfx_barrier(myGridA)
    call blacsfx_barrier(myGridB)

    ! A -> B
    call blacsfx_gemr2d(mm, nn, aa, 1, 1, descA, bb, 1, 1, descB, myGridA%ctxt)

    aa(:,:) = 0

    ! B -> A
    call blacsfx_gemr2d(mm, nn, bb, 1, 1, descB, aa, 1, 1, descA, myGridA%ctxt)

    allocate(iErr(0:nProc-1))
    iErr(:) = 0
    ! check it's what is expected
    if (all(shape(AA) > [0,0])) then
      lpCol: do ii = 1, size(AA, dim=2)
        iGlob = scalafx_indxl2g(ii, descA(NB_), myGridA%mycol, descA(CSRC_), myGridA%ncol)
        if (iGlob > nn) then
          cycle lpCol
        end if
        lpRow: do jj = 1, size(AA, dim=1)
          jGlob = scalafx_indxl2g(jj, descA(MB_), myGridA%myrow, descA(RSRC_), myGridA%nrow)
          if (jGlob > mm) then
            cycle lpRow
          end if
          if (AA(jj,ii) /= numelem(iGlob, jGlob, mm)) then
            write(stdOut,*) "Error on processor ", myGridA%iproc, " Mismatch in element ", jGlob,&
                & iGlob
            iErr(myGridA%iproc) = numelem(iGlob, jGlob, mm)
            exit lpcol
          end if
        end do lpRow
      end do lpCol
    end if

    call blacsfx_gsum(myGridA, iErr)

    if (myGridB%lead) then
      if (any(iErr /= 0)) then
        write(stdOut,*)"Errors for matrix elements"
        write(stdOut,*) iErr
      end if
    end if

    ! Finish blacs.
    call blacsfx_exit()

  end subroutine main

end program test_gemr2d
