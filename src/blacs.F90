include(blacs.m4)

!> Interfaces to BLACS routines.
module blacs_module
  use scalapackfx_common_module
  implicit none
  private

  public :: blacs_pinfo, blacs_get, blacs_gridinfo, blacs_gridinit
  public :: blacs_barrier, blacs_exit, blacs_abort, blacs_pnum
  public :: gebs2d, gebr2d, gesd2d, gerv2d, gsum2d
  
  interface
    !> Returns the number of processes available for use.
    !! \see BLACS documentation for details.
    subroutine blacs_pinfo(id, nproc)
      integer, intent(out) :: id, nproc
    end subroutine blacs_pinfo
    
    !> Gets values that BLACS uses for internal defaults.
    !! \see BLACS documentation for details.
    subroutine blacs_get(ictxt, what, val)
      integer, intent(in) :: ictxt, what
      integer, intent(out) :: val
    end subroutine blacs_get

    !> Delivers information about the grid.
    !! \see BLACS documentation for details.
    subroutine blacs_gridinfo(ictxt, nprow, npcol, myrow, mycol)
      integer, intent(in) :: ictxt,nprow, npcol
      integer, intent(out) :: myrow, mycol
    end subroutine blacs_gridinfo

    !> Assigns available processes into BLACS process grid. 
    !! \see BLACS documentation for details.
    subroutine blacs_gridinit(ictxt, order, nprow, npcol)
      integer, intent(inout) :: ictxt
      character, intent(in) :: order
      integer, intent(in) :: nprow, npcol
    end subroutine blacs_gridinit

    !> Creates a customized process grid.
    !! \see BLACS documentation for details.
    subroutine blacs_gridmap(ictxt, umap, ldumap, nprow, npcol)
      integer, intent(inout) :: ictxt
      integer, intent(in) :: ldumap
      integer, intent(in) :: umap(ldumap, *), nprow, npcol
    end subroutine blacs_gridmap

    !> Calls a barrier.
    !! \see BLACS documentation for details.
    subroutine blacs_barrier(ictxt, scope)
      integer, intent(in) :: ictxt
      character, intent(in) :: scope
    end subroutine blacs_barrier

    !> Exits blacs communicator.
    !! \see BLACS documentation for details.
    subroutine blacs_exit(cont)
      integer, intent(in) :: cont
    end subroutine blacs_exit
    
    !> Frees all BLACS contexts and releases all allocated memory.
    !! \see BLACS documentation for details.
    subroutine blacs_abort(ictxt, errornum)
      integer, intent(in) :: ictxt, errornum
    end subroutine blacs_abort

    !> Returns the system process number of the process in the process grid.
    !! \see BLACS documentation for details.
    function blacs_pnum(ictxt, prow, pcol) result(res)
      integer, intent(in) :: ictxt, prow, pcol
      integer :: res
    end function blacs_pnum

    !> Returns row and column of a process in the grid
    !! \see BLACS documentation for details.
    subroutine blacs_pcoord(ictxt, pnum, prow, pcol)
      integer, intent(in) :: ictxt, pnum
      integer, intent(out) :: prow, pcol
    end subroutine blacs_pcoord
  end interface

  !> Wrapper around type specific ?gebs2d subroutines.
  !! \see BLACS documentation for details.
  interface gebs2d
    _subroutine_interface_gebs2d(integer, i, integer)
    _subroutine_interface_gebs2d(real, s, real(sp))
    _subroutine_interface_gebs2d(dreal, d, real(dp))
    _subroutine_interface_gebs2d(complex, c, complex(sp))
    _subroutine_interface_gebs2d(dcomplex, z, complex(dp))
  end interface gebs2d

  !> Wrapper around type specific ?gebr2d subroutines.
  !! \see BLACS documentation for details.
  interface gebr2d
    _subroutine_interface_gebr2d(integer, i, integer)
    _subroutine_interface_gebr2d(real, s, real(sp))
    _subroutine_interface_gebr2d(dreal, d, real(dp))
    _subroutine_interface_gebr2d(complex, c, complex(sp))
    _subroutine_interface_gebr2d(dcomplex, z, complex(dp))
  end interface gebr2d

  !> Wrapper around type specific ?gesd2d subroutines.
  !! \see BLACS documentation for details.
  interface gesd2d
    _subroutine_interface_gesd2d(integer, i, integer)
    _subroutine_interface_gesd2d(real, s, real(sp))
    _subroutine_interface_gesd2d(dreal, d, real(dp))
    _subroutine_interface_gesd2d(complex, c, complex(sp))
    _subroutine_interface_gesd2d(dcomplex, z, complex(dp))
  end interface gesd2d

  !> Wrapper around type specific ?gerv2d subroutines.
  !! \see BLACS documentation for details.
  interface gerv2d
    _subroutine_interface_gerv2d(integer, i, integer)
    _subroutine_interface_gerv2d(real, s, real(sp))
    _subroutine_interface_gerv2d(dreal, d, real(dp))
    _subroutine_interface_gerv2d(complex, c, complex(sp))
    _subroutine_interface_gerv2d(dcomplex, z, complex(dp))
  end interface gerv2d

  !> Wrapper around type specific ?gsum2d subroutines.
  !! \see BLACS documentation for details.
  interface gsum2d
    _subroutine_interface_gsum2d(integer, i, integer)
    _subroutine_interface_gsum2d(real, s, real(sp))
    _subroutine_interface_gsum2d(dreal, d, real(dp))
    _subroutine_interface_gsum2d(complex, c, complex(sp))
    _subroutine_interface_gsum2d(dcomplex, z, complex(dp))
  end interface gsum2d

end module blacs_module
