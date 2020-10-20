Using SCALAPACKFX
=================

Before you can use the SCALAPACKFX routines, you need basically the following
two steps.

#. Import the module `libscalapackfx_module` in your routines.

#. Initialize a grid for the BLACS-communcation using `type(blacsgrid)`.

Below you find an example draft for reading in a matrix from a file,
distributing it across the processes and finally diagonalizing it::

    program test_diag
      use libscalapackfx_module
    
      integer, parameter :: nprow = 4, npcol = 4   ! process rows/cols
      integer, parameter :: bsize = 64             ! block size
      integer, parameter :: nn = 1000              ! matrix size
    
      type(blacsgrid) :: mygrid  ! BLACS grid descriptor
      real(dp), allocatable :: aa(:,:), bb(:,:), eigvecs(:,:), eigvals(:)
      integer :: desc(DLEN_)    ! matrix descriptor
      integer :: mloc, nloc     ! nr. of local rows/columns of the matrices
      type(linecomm) :: distr               ! distributes matrix when read
      real(dp), allocatable :: iobuffer(:)  ! buffer during read
    
      ! Initialize your BLACS grid
      call mygrid%gridinit(nprow, npcol)
    
      ! Allocate the local part of the distributed matrices
      call scalafx_getdescriptor(mygrid, nn, nn, bsize, bsize, desc)
      call scalafx_getlocalshape(mygrid, desc, mloc, nloc)
      allocate(aa(mloc, nloc))
      allocate(bb(mloc, nloc))
      allocate(eigvecs(mloc, nloc))
      allocate(eigvals(nn))
      ...
      ! Here comes the code which distributes your matrix
      ! You can use the linecomm type to distribute it if you read it from file
      if (mygrid%lead) then
          allocate(iobuffer(nn))
      end if
      call distr%init(desc, "c")
      do icol = 1, nn
        if (mygrid%lead) then
          read(12, *) iobuffer(:)
          call distr%setline_lead(mygrid, icol, iobuffer, aa)
        else
          call distr%setline_follow(mygrid, icol, mtxloc)
        end if
      end do
      ...
      ! Get eigenvalues (on all nodes) and eigenvectors (distributed)
      call psygvd(aa, desc, bb, desc, eigvals, eigvecs, desc, jobz="V", uplo="L")
      ...
    end program test_diag


Have a look at the test folder in the source tree for further examples.
