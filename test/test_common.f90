module test_common_module
  use blacsfx_module
  use scalapackfx_module
  use scalapackfx_tools_module
  implicit none
  private

  public :: dp
  public :: readfromfile, writetofile

  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.0d0)

  interface writetofile
    module procedure writetofile_real
    module procedure writetofile_cmplx
  end interface writetofile

contains


  !> Read distributed matrix from file.
  !! \param mygrid  BLACS descriptor
  !! \param fname  Name of the file to read the matrix from.
  !! \param mb  Row block size of the matrix.
  !! \param nb  Column block size of the matrix.
  !! \param mtxloc  Allocated local matrix on exit.
  !! \param desc  Matrix descriptor
  subroutine readfromfile(mygrid, fname, mb, nb, mtxloc, desc)
    type(blacsgrid), intent(inout) :: mygrid
    character(*), intent(in) :: fname
    integer, intent(in) :: mb, nb
    real(dp), allocatable, intent(out) :: mtxloc(:,:)
    integer, intent(out) :: desc(DLEN_)

    integer :: mm, nn

    if (mygrid%lead) then
      open(12, file=fname, action="read", status="old")
      read(12, *) mm, nn
      call blacsfx_gebs(mygrid, mm)
      call blacsfx_gebs(mygrid, nn)
    else
      call blacsfx_gebr(mygrid, mm)
      call blacsfx_gebr(mygrid, nn)
    end if
    call scalafx_creatematrix(mygrid, mm, nn, mb, nb, mtxloc, desc)
    if (mygrid%lead) then
      call readarray_lead(mygrid, 12, desc, mtxloc, formatted=.true.)
      close(12)
    else
      call readarray_follow(mygrid, desc, mtxloc)
    end if

  end subroutine readfromfile


  !> Write distributed matrix to file.
  !! \param mygrid  BlACS descriptor
  !! \param fname  Name of the file to write the matrix to.
  !! \param mtxloc  Local part of the matrix.
  !! \param desc  Matrix descriptor
  subroutine writetofile_real(mygrid, fname, mtxloc, desc)
    type(blacsgrid), intent(in) :: mygrid
    character(*), intent(in) :: fname
    real(dp), intent(in) :: mtxloc(:,:)
    integer, intent(in) :: desc(DLEN_)

    if (mygrid%lead) then
      open(12, file=fname, form="formatted", status="replace")
      write(12, "(I0,1X,I0)") desc(M_), desc(N_)
      call writearray_lead(mygrid, 12, desc, mtxloc, elemformat="(ES23.15)")
      close(12)
    else
      call writearray_follow(mygrid, desc, mtxloc)
    end if

  end subroutine writetofile_real

  !> Write distributed matrix to file.
  !! \param mygrid  BlACS descriptor
  !! \param fname  Name of the file to write the matrix to.
  !! \param mtxloc  Local part of the matrix.
  !! \param desc  Matrix descriptor
  subroutine writetofile_cmplx(mygrid, fname, mtxloc, desc)
    type(blacsgrid), intent(in) :: mygrid
    character(*), intent(in) :: fname
    complex(dp), intent(in) :: mtxloc(:,:)
    integer, intent(in) :: desc(DLEN_)

    if (mygrid%lead) then
      open(12, file=fname, form="formatted", status="replace")
      write(12, "(I0,1X,I0)") desc(M_), desc(N_)
      call writearray_lead(mygrid, 12, desc, mtxloc, elemformat="(ES23.15)")
      close(12)
    else
      call writearray_follow(mygrid, desc, mtxloc)
    end if

  end subroutine writetofile_cmplx


end module test_common_module

