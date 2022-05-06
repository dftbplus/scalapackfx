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

  ! number elements in global matrix
  pure function numelem(iGlob, jGlob, mm)
    integer, intent(in) :: iGlob, jGlob, mm
    integer :: numelem
    numelem = jGlob + (iGlob-1) * mm
  end function numelem

  subroutine main()
    type(blacsgrid) :: mygridA, mygridB
    integer, allocatable :: AA(:,:), BB(:,:)
    integer :: descA(DLEN_), descB(DLEN_)
    integer :: nprow, npcol, mm, nn, iproc, nproc, ii, jj, iGlob, jGlob
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
    call mygridA%initgrid(nprow, npcol)
    call mygridB%initgrid(1, nproc)
    if (mygridA%lead) then
      write(stdout, "(A,2(1X,I0))") "# source processor grid:", nprow, npcol
    end if
    if (mygridB%lead) then
      write(stdout, "(A,2(1X,I0))") "# destination processor grid:", 1, nproc
    end if

    if (mygridA%lead) then
      write(stdOut, *)'# Matrix size to re-distrbute?'
      read(stdin, *) mm, nn
      call blacsfx_gebs(mygridA, mm)
      call blacsfx_gebs(mygridA, nn)
    else
      call blacsfx_gebr(mygridA, mm)
      call blacsfx_gebr(mygridA, nn)
    end if

    ! Block size for matrix B
    mbB = mm
    nbB = 1

    if (mygridA%lead) then
      write(stdout, "(A,2(1X,I0))") "# A processor grid:", nprow, npcol
      write(stdout, "(A,1X,I0,1X,I0)") "# A block size:", mbA, nbA
    end if
    if (mygridB%lead) then
      write(stdout, "(A,2(1X,I0))") "# B processor grid:", nprow, npcol
      write(stdout, "(A,1X,I0,1X,I0)") "# B block size:", mbB, nbB
    end if

    call scalafx_creatematrix(mygridA, mm, nn, mbA, nbA, AA, descA)
    call scalafx_creatematrix(mygridB, mm, nn, mbB, nbB, BB, descB)

    AA(:,:) = 0
    BB(:,:) = 0

    write(stdOut,"(A,I0,A,I0,',',I0,A)")'Processor ', mygridA%iproc,' holds A(', shape(AA),&
        & ') elements'

    write(stdOut,"(A,I0,A,I0,',',I0,A)")'Processor ', mygridB%iproc,' holds B(', shape(BB),&
        & ') elements'

    ! number off elements in A
    do ii = 1, size(AA, dim=2)
      iGlob = scalafx_indxl2g(ii, descA(NB_), mygridA%mycol, descA(CSRC_), mygridA%ncol)
      do jj = 1, size(AA, dim=1)
        jGlob = scalafx_indxl2g(jj, descA(MB_), mygridA%myrow, descA(RSRC_), mygridA%nrow)
        AA(jj,ii) = numelem(iGlob, jGlob, mm)
      end do
    end do

    ! A -> B
    call blacsfx_gemr2d(mm, nn, aa, 1, 1, descA, bb, 1, 1, descB, mygridA%ctxt)

    aa(:,:) = 0

    ! B -> A
    call blacsfx_gemr2d(mm, nn, bb, 1, 1, descB, aa, 1, 1, descA, mygridA%ctxt)

    allocate(iErr(0:nProc-1))
    iErr(:) = 0
    ! check it's what is expected
    if (all(shape(AA) > [0,0])) then
      lpCol: do ii = 1, size(AA, dim=2)
        iGlob = scalafx_indxl2g(ii, descA(NB_), mygridA%mycol, descA(CSRC_), mygridA%ncol)
        if (iGlob > nn) then
          cycle lpCol
        end if
        lpRow: do jj = 1, size(AA, dim=1)
          jGlob = scalafx_indxl2g(jj, descA(MB_), mygridA%myrow, descA(RSRC_), mygridA%nrow)
          if (jGlob > mm) then
            cycle lpRow
          end if
          if (AA(jj,ii) /= numelem(iGlob, jGlob, mm)) then
            write(stdOut,*) "Error on processor ", mygridA%iproc, " Mismatch in element ", jGlob,&
                & iGlob
            iErr(mygridA%iproc) = numelem(iGlob, jGlob, mm)
            exit lpcol
          end if
        end do lpRow
      end do lpCol
    end if

    call blacsfx_gsum(mygridA, iErr)

    if (mygridB%lead) then
      if (any(iErr /= 0)) then
        write(stdOut,*)"Errors for matrix elements"
        write(stdOut,*) iErr
      end if
    end if

    ! Finish blacs.
    call blacsfx_exit()

  end subroutine main

end program test_gemr2d
