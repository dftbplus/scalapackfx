include(scalapackfx_tools.m4)
  
!> Some extension routines to the scalapack library making it more usable.
!!
!! \details Some of the types exported by this module are defined in other
!! modules. Therefore, see also the documentation of the following modules:
!!
!!   * \ref linecomm_module "linecomm_module"
!!
module scalapackfx_tools_module
  use scalapackfx_common_module
  use blacsfx_module
  use scalapackfx_module
  use linecomm_module
  implicit none
  private

  public :: scalafx_cpl2g, scalafx_cpg2l, scalafx_addl2g, scalafx_addg2l
  public :: writearray_master, writearray_slave
  public :: readarray_master, readarray_slave
  public :: blocklist, size
  public :: linecomm

  interface scalafx_cpl2g
    module procedure cpl2g_real, cpl2g_dreal, cpl2g_complex, cpl2g_dcomplex
    module procedure cpl2g_int
  end interface scalafx_cpl2g

  interface scalafx_cpg2l
    module procedure cpg2l_real, cpg2l_dreal, cpg2l_complex, cpg2l_dcomplex
    module procedure cpg2l_int
  end interface scalafx_cpg2l

  interface scalafx_addl2g
    module procedure addl2g_real, addl2g_dreal, addl2g_complex, addl2g_dcomplex
    module procedure addl2g_int
  end interface scalafx_addl2g

  interface scalafx_addg2l
    module procedure addg2l_real, addg2l_dreal, addg2l_complex, addg2l_dcomplex
    module procedure addg2l_int
  end interface scalafx_addg2l

  interface writearray_master
    module procedure writearray_master_int
    module procedure writearray_master_real, writearray_master_dreal
    module procedure writearray_master_complex, writearray_master_dcomplex
  end interface writearray_master

  interface writearray_slave
    module procedure writearray_slave_int
    module procedure writearray_slave_real, writearray_slave_dreal
    module procedure writearray_slave_complex, writearray_slave_dcomplex
  end interface writearray_slave

  interface readarray_master
    module procedure readarray_master_int
    module procedure readarray_master_real, readarray_master_dreal
    module procedure readarray_master_complex, readarray_master_dcomplex
  end interface readarray_master

  interface readarray_slave
    module procedure readarray_slave_int
    module procedure readarray_slave_real, readarray_slave_dreal
    module procedure readarray_slave_complex, readarray_slave_dcomplex
  end interface readarray_slave
  

  !> List of the local blocks of a distributed matrix.
  !!
  !! \details This structure can be helpful when modifying a distributed
  !! matrix directly on the local nodes. For example, in order to multiply
  !! every column of a distributed matrix by a column dependent factor, 
  !! you could use the blocklist the following way:
  !!
  !!     type(blocklist) :: blocks
  !!     integer :: ii, jj, jglob, jloc, bsize
  !!     :
  !!     call blocks%init(mygrid, descaa, "c")
  !!     do ii = 1, size(blocks)
  !!       call blocks%getblock(ii, jglob, jloc, bsize)
  !!       do jj = 0, bsize  - 1
  !!         aa(:,jloc + jj) =  aa(:,jloc + jj) * ff(jglob + jj)
  !!       end do
  !!     end do
  !!
  type :: blocklist
    private
    integer :: nn, nb, nproc, myproc, srcproc, nblock
  contains
    !> Initializes the instance.
    procedure :: init => blocklist_init

    !> Returns the size of the blocklist.
    procedure :: getsize => blocklist_getsize

    !> Returns the indices (local and global) of a local block.
    procedure :: getblock => blocklist_getblock
  end type blocklist

  interface size
    module procedure blocklist_getsize
  end interface size

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Blocklist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initializes a blocklist instance.
  !! \param self  Initialized instance on exit.
  !! \param mygrid  BLACS descriptor
  !! \param desc  Descriptor of the distributed matrix.
  !! \param rowcol  "C" for column, "R" for row blocks.
  subroutine blocklist_init(self, mygrid, desc, rowcol)
    class(blocklist), intent(inout) :: self
    type(blacsgrid), intent(in) :: mygrid
    integer, intent(in) :: desc(DLEN_)
    character, intent(in) :: rowcol

    integer :: nblockall, nextrablock, mydist

    if (rowcol == "c" .or. rowcol == "C") then
      self%nn = desc(N_)
      self%nb = desc(NB_)
      self%nproc = mygrid%ncol
      self%myproc = mygrid%mycol
      self%srcproc = desc(CSRC_)
    else
      self%nn = desc(M_)
      self%nb = desc(MB_)
      self%nproc = mygrid%nrow
      self%myproc = mygrid%myrow
      self%srcproc = desc(RSRC_)
    end if
    nblockall = self%nn / self%nb
    self%nblock = nblockall / self%nproc
    nextrablock = mod(nblockall, self%nproc)
    mydist = mod(self%nproc + self%myproc - self%srcproc, self%nproc)
    if (mydist < nextrablock) then
      self%nblock = self%nblock + 1
    elseif (mydist == nextrablock .and. mod(self%nn, self%nb) > 0) then
      self%nblock = self%nblock +1
    end if

  end subroutine blocklist_init


  !> Returns the size of the blocklist.
  !! \param self  Instance.
  !! \returns Number of local blocks of the distributed matrix.
  function blocklist_getsize(self) result(res)
    class(blocklist), intent(in) :: self
    integer :: res

    res = self%nblock

  end function blocklist_getsize


  !> Returns the indices (local and global) of a local block.
  !! \param self  Blocklist instance.
  !! \param iblock  Number of local block.
  !! \param iglob  Index of the first element of the block in the global matrix.
  !! \param iloc  Index of the first element of the block in the local matirx.
  !! \param bsize  Size of the block (number of elements in the block).
  subroutine blocklist_getblock(self, iblock, iglob, iloc, bsize)
    class(blocklist), intent(in) :: self
    integer, intent(in) :: iblock
    integer, intent(out) :: iglob, iloc, bsize

    integer :: mydist

    if (iblock >= 1 .and. iblock <= self%nblock) then
      mydist = mod(self%nproc + self%myproc - self%srcproc, self%nproc)
      iglob = ((iblock - 1) * self%nproc + mydist) * self%nb + 1
      iloc = (iblock - 1) * self%nb + mod(iglob - 1, self%nb) + 1
      bsize = min(self%nb, self%nn - iglob + 1)
    else
      iglob = 0
      iloc = 0
      bsize = 0
    end if

  end subroutine blocklist_getblock


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Block copy/adding routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  _subroutine_cpl2g(real, real(sp))
  _subroutine_cpl2g(dreal, real(dp))
  _subroutine_cpl2g(complex, complex(sp))
  _subroutine_cpl2g(dcomplex, complex(dp))
  _subroutine_cpl2g(int, integer)

  _subroutine_cpg2l(real, real(sp))
  _subroutine_cpg2l(dreal, real(dp))
  _subroutine_cpg2l(complex, complex(sp))
  _subroutine_cpg2l(dcomplex, complex(dp))
  _subroutine_cpg2l(int, integer)

  _subroutine_addl2g(real, real(sp))
  _subroutine_addl2g(dreal, real(dp))
  _subroutine_addl2g(complex, complex(sp))
  _subroutine_addl2g(dcomplex, complex(dp))
  _subroutine_addl2g(int, integer)

  _subroutine_addg2l(real, real(sp))
  _subroutine_addg2l(dreal, real(dp))
  _subroutine_addg2l(complex, complex(sp))
  _subroutine_addg2l(dcomplex, complex(dp))
  _subroutine_addg2l(int, integer)
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! writearray/readarray
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  _subroutine_writearray_master(int, integer)
  _subroutine_writearray_master(real, real(sp))
  _subroutine_writearray_master(dreal, real(dp))
  _subroutine_writearray_master(complex, complex(sp))
  _subroutine_writearray_master(dcomplex, complex(dp))

  _subroutine_writearray_slave(int, integer)
  _subroutine_writearray_slave(real, real(sp))
  _subroutine_writearray_slave(dreal, real(dp))
  _subroutine_writearray_slave(complex, complex(sp))
  _subroutine_writearray_slave(dcomplex, complex(dp))

  _subroutine_readarray_master(int, integer)
  _subroutine_readarray_master(real, real(sp))
  _subroutine_readarray_master(dreal, real(dp))
  _subroutine_readarray_master(complex, complex(sp))
  _subroutine_readarray_master(dcomplex, complex(dp))

  _subroutine_readarray_slave(int, integer)
  _subroutine_readarray_slave(real, real(sp))
  _subroutine_readarray_slave(dreal, real(dp))
  _subroutine_readarray_slave(complex, complex(sp))
  _subroutine_readarray_slave(dcomplex, complex(dp))
  

end module scalapackfx_tools_module
