!> Testing matrix inversion via getrf and getri
program test_pmatinv
  use, intrinsic :: iso_fortran_env, stdout => output_unit
  use test_common_module
  use libscalapackfx_module
  implicit none

  ! Block size (using an extremely small value for test purposes)
  integer, parameter :: bsize = 1

  call main()

contains

  subroutine main()
    type(blacsgrid) :: mygrid
    real(dp), allocatable :: r_xx(:,:), r_xxOld(:,:), r_matI(:,:)
    complex(dp), allocatable :: c_xx(:,:), c_xxOld(:,:), c_matI(:,:)

    integer :: descx(DLEN_)
    integer :: nprow, npcol, mm, nn, iproc, nproc, ii, jj, iLoc, jLoc, kLoc, jGlob
    logical :: isLocal, failed

    ! Initialize blas and create a square processor grid
    call blacsfx_pinfo(iproc, nproc)
    do nprow = int(sqrt(real(nproc, dp))), nproc
      if (mod(nproc, nprow) == 0) then
        exit
      end if
    end do
    npcol = nproc / nprow
    call mygrid%initgrid(nprow, npcol)
    if (mygrid%lead) then
      write(stdout, "(A,2(1X,I0))") "# processor grid:", nprow, npcol
      write(stdout, "(A,1X,I0)") "# block size:", bsize
    end if

    ! Set up matrix
    if (mygrid%lead) then
      write(stdout, "(A)") "Matrix read from file 'oversqr.dat'."
    end if
    call readfromfile(mygrid, "oversqr.dat", bsize, bsize, r_xx, descx)
    mm = descx(M_)
    nn = descx(N_)

    if (mygrid%lead) then
      write(stdout, "(A,2(1X,I0))") "# global matrix size:", mm, nn
      write(stdout, "(A,2(1X,I0))") "# local matrix size on leader:",&
          & size(r_xx, dim=1), size(r_xx, dim=2)
    end if

    !!!!!!!!!!!!!
    ! Real case !
    !!!!!!!!!!!!!

    r_xxOld = r_xx
    r_matI = r_xx ! just to get sizing

    ! Set r_xx to an eye matrix
    r_xx(:,:) = 0.0_dp
    do ii = 1, size(r_xx,dim=2)
      jGlob = scalafx_indxl2g(ii, descx(NB_), mygrid%mycol, descx(CSRC_), mygrid%ncol)
      ! where is the global?
      call scalafx_islocal(mygrid, descx, jGlob, jGlob, isLocal, iLoc, jLoc)
      if (isLocal) then
        ! a global diagonal element is stored here
        r_xx(iLoc,jLoc) = 1.0_dp
      end if
    end do

    ! Solve for xx (*) xx^-1 = 1
    call scalafx_pposv(r_matI, descx, r_xx, descx)


    ! xx (*) xx^-1
    call pblasfx_pgemm(r_xx, descx, r_xxOld, descx, r_matI, descx)

    failed = .false.
    do ii = 1, size(r_matI,dim=2)
      jGlob = scalafx_indxl2g(ii, descx(NB_), mygrid%mycol, descx(CSRC_), mygrid%ncol)
      ! where is the global
      call scalafx_islocal(mygrid, descx, jGlob, jGlob, isLocal, iLoc, jLoc)
      if (isLocal) then
        ! a global diagonal element is stored here
        if (abs(r_matI(iLoc,jLoc) - 1.0_dp) > epsilon(0.0)) failed = .true.
      else
        do kLoc = 1, size(r_matI,dim=1)
          if (abs(r_matI(iLoc,jLoc)) > epsilon(0.0)) failed = .true.
        end do
      end if
    end do

    deallocate(r_xx)
    deallocate(r_matI)

    ! Would normally accumulate product via mpi calls, but messy to do via blacs operations only
    if (failed) then
      write(stdout,*)'Real matrix element(s) on processor', mygrid%iproc, ' non-identity matrix'
    else
      write(stdout,*)'Real matrix elements on processor', mygrid%iproc, ' are OK'
    end if

    !!!!!!!!!!!!!!!!
    ! Complex case !
    !!!!!!!!!!!!!!!!

    c_xx = r_xxOld
    deallocate(r_xxOld)
    c_xxOld = c_xx
    c_matI = c_xx

    ! Set c_xx to an eye matrix
    c_xx(:,:) = cmplx(0,0,dp)
    do ii = 1, size(c_matI,dim=2)
      jGlob = scalafx_indxl2g(ii, descx(NB_), mygrid%mycol, descx(CSRC_), mygrid%ncol)
      ! where is the global?
      call scalafx_islocal(mygrid, descx, jGlob, jGlob, isLocal, iLoc, jLoc)
      if (isLocal) then
        ! a global diagonal element is stored here
        c_xx(iLoc,jLoc) = cmplx(1,0,dp)
      end if
    end do

    ! Solve for xx (*) xx^-1 = 1
    call scalafx_pposv(c_matI, descx, c_xx, descx)

    ! xx (*) xx^-1
    call pblasfx_pgemm(c_xx, descx, c_xxOld, descx, c_matI, descx)

    failed = .false.
    do ii = 1, size(c_matI,dim=2)
      jGlob = scalafx_indxl2g(ii, descx(NB_), mygrid%mycol, descx(CSRC_), mygrid%ncol)
      ! where is the global
      call scalafx_islocal(mygrid, descx, jGlob, jGlob, isLocal, iLoc, jLoc)
      if (isLocal) then
        ! a global diagonal element is stored here
        if (abs(c_matI(iLoc,jLoc) - cmplx(1,0,dp)) > epsilon(0.0)) failed = .true.
      else
        do kLoc = 1, size(c_matI,dim=1)
          if (abs(c_matI(iLoc,jLoc)) > epsilon(0.0)) failed = .true.
        end do
      end if
    end do

    ! Would normally accumulate product via mpi calls, but messy to do via blacs operations only
    if (failed) then
      write(stdout,*)'Complex matrix element(s) on processor', mygrid%iproc, ' non-identity matrix'
    else
      write(stdout,*)'Complex matrix elements on processor', mygrid%iproc, ' are OK'
    end if

    deallocate(c_xx)
    deallocate(c_matI)
    deallocate(c_xxOld)

    ! Finish blacs.
    call blacsfx_exit()

  end subroutine main

end program test_pmatinv
