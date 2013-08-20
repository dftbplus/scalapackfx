include(blacsfx.m4)

!> High level Fortran wrappers for the BLACS library.
module blacsfx_module
  use scalapackfx_common_module
  use blacsgrid_module
  use blacs_module
  implicit none
  private

  ! Public names.
  public :: DLEN_, DT_, CTXT_, M_, N_, MB_, NB_, RSRC_, CSRC_, LLD_
  public :: blacsgrid
  public :: blacsfx_gebs, blacsfx_gebr
  public :: blacsfx_gesd, blacsfx_gerv
  public :: blacsfx_gsum
  public :: blacsfx_barrier
  public :: blacsfx_pinfo, blacsfx_pcoord, blacsfx_pnum, blacsfx_exit


  !> Wrapper around ?gebs2d for data of rank 0, 1, 2.
  interface blacsfx_gebs
    module procedure blacsfx_gebs_i0, blacsfx_gebs_i1, blacsfx_gebs_i2
    module procedure blacsfx_gebs_s0, blacsfx_gebs_s1, blacsfx_gebs_s2
    module procedure blacsfx_gebs_d0, blacsfx_gebs_d1, blacsfx_gebs_d2
    module procedure blacsfx_gebs_c0, blacsfx_gebs_c1, blacsfx_gebs_c2
    module procedure blacsfx_gebs_z0, blacsfx_gebs_z1, blacsfx_gebs_z2
  end interface blacsfx_gebs

  !> Wrapper around ?gebr2d for data of rank 0, 1, 2.
  interface blacsfx_gebr
    module procedure blacsfx_gebr_i0, blacsfx_gebr_i1, blacsfx_gebr_i2
    module procedure blacsfx_gebr_s0, blacsfx_gebr_s1, blacsfx_gebr_s2
    module procedure blacsfx_gebr_d0, blacsfx_gebr_d1, blacsfx_gebr_d2
    module procedure blacsfx_gebr_c0, blacsfx_gebr_c1, blacsfx_gebr_c2
    module procedure blacsfx_gebr_z0, blacsfx_gebr_z1, blacsfx_gebr_z2
  end interface blacsfx_gebr

  !> Wrapper around ?gesd2d for data of rank 0, 1, 2.
  interface blacsfx_gesd
    module procedure blacsfx_gesd_i0, blacsfx_gesd_i1, blacsfx_gesd_i2
    module procedure blacsfx_gesd_s0, blacsfx_gesd_s1, blacsfx_gesd_s2
    module procedure blacsfx_gesd_d0, blacsfx_gesd_d1, blacsfx_gesd_d2
    module procedure blacsfx_gesd_c0, blacsfx_gesd_c1, blacsfx_gesd_c2
    module procedure blacsfx_gesd_z0, blacsfx_gesd_z1, blacsfx_gesd_z2
  end interface blacsfx_gesd

  !> Wrapper around ?gerv2d for data of rank 0, 1, 2.
  interface blacsfx_gerv
    module procedure blacsfx_gerv_i0, blacsfx_gerv_i1, blacsfx_gerv_i2
    module procedure blacsfx_gerv_s0, blacsfx_gerv_s1, blacsfx_gerv_s2
    module procedure blacsfx_gerv_d0, blacsfx_gerv_d1, blacsfx_gerv_d2
    module procedure blacsfx_gerv_c0, blacsfx_gerv_c1, blacsfx_gerv_c2
    module procedure blacsfx_gerv_z0, blacsfx_gerv_z1, blacsfx_gerv_z2
  end interface blacsfx_gerv

  !> Wrapper around ?gsum2d for data of rank 0, 1, 2.
  interface blacsfx_gsum
    module procedure blacsfx_gsum_i0, blacsfx_gsum_i1, blacsfx_gsum_i2
    module procedure blacsfx_gsum_s0, blacsfx_gsum_s1, blacsfx_gsum_s2
    module procedure blacsfx_gsum_d0, blacsfx_gsum_d1, blacsfx_gsum_d2
    module procedure blacsfx_gsum_c0, blacsfx_gsum_c1, blacsfx_gsum_c2
    module procedure blacsfx_gsum_z0, blacsfx_gsum_z1, blacsfx_gsum_z2
  end interface blacsfx_gsum
  

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! GEBS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  _subroutine_blacsfx_gebs_r2(integer, i, integer)
  _subroutine_blacsfx_gebs_r2(real, s, real(sp))
  _subroutine_blacsfx_gebs_r2(dreal, d, real(dp))
  _subroutine_blacsfx_gebs_r2(complex, c, complex(sp))
  _subroutine_blacsfx_gebs_r2(dcomplex, z, complex(dp))

  _subroutine_blacsfx_gebs_r1(integer, i, integer)
  _subroutine_blacsfx_gebs_r1(real, s, real(sp))
  _subroutine_blacsfx_gebs_r1(dreal, d, real(dp))
  _subroutine_blacsfx_gebs_r1(complex, c, complex(sp))
  _subroutine_blacsfx_gebs_r1(dcomplex, z, complex(dp))

  _subroutine_blacsfx_gebs_r0(integer, i, integer)
  _subroutine_blacsfx_gebs_r0(real, s, real(sp))
  _subroutine_blacsfx_gebs_r0(dreal, d, real(dp))
  _subroutine_blacsfx_gebs_r0(complex, c, complex(sp))
  _subroutine_blacsfx_gebs_r0(dcomplex, z, complex(dp))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! GEBR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  _subroutine_blacsfx_gebr_r2(integer, i, integer)
  _subroutine_blacsfx_gebr_r2(real, s, real(sp))
  _subroutine_blacsfx_gebr_r2(dreal, d, real(dp))
  _subroutine_blacsfx_gebr_r2(complex, c, complex(sp))
  _subroutine_blacsfx_gebr_r2(dcomplex, z, complex(dp))
  
  _subroutine_blacsfx_gebr_r1(integer, i, integer)
  _subroutine_blacsfx_gebr_r1(real, s, real(sp))
  _subroutine_blacsfx_gebr_r1(dreal, d, real(dp))
  _subroutine_blacsfx_gebr_r1(complex, c, complex(sp))
  _subroutine_blacsfx_gebr_r1(dcomplex, z, complex(dp))

  _subroutine_blacsfx_gebr_r0(integer, i, integer)
  _subroutine_blacsfx_gebr_r0(real, s, real(sp))
  _subroutine_blacsfx_gebr_r0(dreal, d, real(dp))
  _subroutine_blacsfx_gebr_r0(complex, c, complex(sp))
  _subroutine_blacsfx_gebr_r0(dcomplex, z, complex(dp))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! GESD2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  _subroutine_blacsfx_gesd_r2(integer, i, integer)
  _subroutine_blacsfx_gesd_r2(real, s, real(sp))
  _subroutine_blacsfx_gesd_r2(dreal, d, real(dp))
  _subroutine_blacsfx_gesd_r2(complex, c, complex(sp))
  _subroutine_blacsfx_gesd_r2(dcomplex, z, complex(dp))
  
  _subroutine_blacsfx_gesd_r1(integer, i, integer)
  _subroutine_blacsfx_gesd_r1(real, s, real(sp))
  _subroutine_blacsfx_gesd_r1(dreal, d, real(dp))
  _subroutine_blacsfx_gesd_r1(complex, c, complex(sp))
  _subroutine_blacsfx_gesd_r1(dcomplex, z, complex(dp))

  _subroutine_blacsfx_gesd_r0(integer, i, integer)
  _subroutine_blacsfx_gesd_r0(real, s, real(sp))
  _subroutine_blacsfx_gesd_r0(dreal, d, real(dp))
  _subroutine_blacsfx_gesd_r0(complex, c, complex(sp))
  _subroutine_blacsfx_gesd_r0(dcomplex, z, complex(dp))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! GERV2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  _subroutine_blacsfx_gerv_r2(integer, i, integer)
  _subroutine_blacsfx_gerv_r2(real, s, real(sp))
  _subroutine_blacsfx_gerv_r2(dreal, d, real(dp))
  _subroutine_blacsfx_gerv_r2(complex, c, complex(sp))
  _subroutine_blacsfx_gerv_r2(dcomplex, z, complex(dp))

  _subroutine_blacsfx_gerv_r1(integer, i, integer)
  _subroutine_blacsfx_gerv_r1(real, s, real(sp))
  _subroutine_blacsfx_gerv_r1(dreal, d, real(dp))
  _subroutine_blacsfx_gerv_r1(complex, c, complex(sp))
  _subroutine_blacsfx_gerv_r1(dcomplex, z, complex(dp))

  _subroutine_blacsfx_gerv_r0(integer, i, integer)
  _subroutine_blacsfx_gerv_r0(real, s, real(sp))
  _subroutine_blacsfx_gerv_r0(dreal, d, real(dp))
  _subroutine_blacsfx_gerv_r0(complex, c, complex(sp))
  _subroutine_blacsfx_gerv_r0(dcomplex, z, complex(dp))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! GSUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  _subroutine_blacsfx_gsum_r2(integer, i, integer)
  _subroutine_blacsfx_gsum_r2(real, s, real(sp))
  _subroutine_blacsfx_gsum_r2(dreal, d, real(dp))
  _subroutine_blacsfx_gsum_r2(complex, c, complex(sp))
  _subroutine_blacsfx_gsum_r2(dcomplex, z, complex(dp))

  _subroutine_blacsfx_gsum_r1(integer, i, integer)
  _subroutine_blacsfx_gsum_r1(real, s, real(sp))
  _subroutine_blacsfx_gsum_r1(dreal, d, real(dp))
  _subroutine_blacsfx_gsum_r1(complex, c, complex(sp))
  _subroutine_blacsfx_gsum_r1(dcomplex, z, complex(dp))

  _subroutine_blacsfx_gsum_r0(integer, i, integer)
  _subroutine_blacsfx_gsum_r0(real, s, real(sp))
  _subroutine_blacsfx_gsum_r0(dreal, d, real(dp))
  _subroutine_blacsfx_gsum_r0(complex, c, complex(sp))
  _subroutine_blacsfx_gsum_r0(dcomplex, z, complex(dp))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Barrier
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Holds up execution of all processes within given scope.
  !!
  !! \param self  BLACS group descriptor
  !! \param scope  Scope of the barrier (default: "A")
  !!
  subroutine blacsfx_barrier(mygrid, scope)
    type(blacsgrid), intent(in) :: mygrid
    character, intent(in), optional :: scope

    character :: scope0

    _handle_inoptflag(scope0, scope, "A")
    call blacs_barrier(mygrid%ctxt, scope0)

  end subroutine blacsfx_barrier

  
  !> Delivers process information.
  !!
  !! \param iproc  Id of the process (0 <= iproc < nproc)
  !! \param nproc  Nr. of processes.
  !!
  subroutine blacsfx_pinfo(iproc, nproc)
    integer, intent(out) :: iproc, nproc

    call blacs_pinfo(iproc, nproc)

  end subroutine blacsfx_pinfo


  !> Delivers row and column of a given process in a grid.
  !!
  !! \param mygrid  BLACS grid.
  !! \param iproc  Process of which position should be determined.
  !! \param prow  Row of the process.
  !! \param pcol  Column of the process.
  !!
  subroutine blacsfx_pcoord(mygrid, iproc, prow, pcol)
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: iproc
    integer, intent(out) :: prow, pcol

    call blacs_pcoord(mygrid%ctxt, iproc, prow, pcol)
    
  end subroutine blacsfx_pcoord


  !> Delivers process number for a given process in the grid.
  !!
  !! \param mygrid BLACS grid.
  !! \param prow  Row of the process.
  !! \param pcol  Column of the process.
  !! \return  Process number (id) of the process with the given coordinates.
  !!
  function blacsfx_pnum(mygrid, prow, pcol) result(pnum)
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: prow, pcol
    integer :: pnum

    pnum = blacs_pnum(mygrid%ctxt, prow, pcol)

  end function blacsfx_pnum

  
  !> Stops BLACS communication.
  !!
  !! \param keepmpi If set to yes, the MPI framework will kept alive after
  !!     BLACS is switched off (default: .false.)
  !!
  subroutine blacsfx_exit(keepmpi)
    logical, intent(in), optional :: keepmpi

    logical :: keepmpi0

    _handle_inoptflag(keepmpi0, keepmpi, .false.)
    if (keepmpi0) then
      call blacs_exit(1)
    else
      call blacs_exit(0)
    end if

  end subroutine blacsfx_exit


end module blacsfx_module
