include(common.m4)

dnl ****************************************************************************
dnl *** blacsfx_gebs
dnl ****************************************************************************

define(`_subroutine_blacsfx_gebs_r2',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Starts broadcast ($1, rank 2).
!! \param mygrid  BLACS descriptor.
!! \param aa  Matrix to broadcast.
!! \param scope  Scope of the broadcast (default: "A").
!! \param top  Topology of the broadcast (default " ").
!! \see BLACS documentation (routine ?gebs2d).
subroutine blacsfx_gebs_$2`2'(mygrid, aa, scope, top)
  class(blacsgrid), intent(in) :: mygrid
  $3, intent(in) :: aa(:,:)
  character, intent(in), optional :: scope, top

  character :: scope0, top0

  _handle_inoptflag(scope0, scope, "A")
  _handle_inoptflag(top0, top, " ")  
  call gebs2d(mygrid%ctxt, scope0, top0, size(aa, dim=1), size(aa, dim=2),&
      & aa, size(aa, dim=1))
  
end subroutine blacsfx_gebs_$2`2'
')

define(`_subroutine_blacsfx_gebs_r1',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Starts broadcast ($1, rank 1).
!! \param mygrid  BLACS descriptor.
!! \param aa  Vector to broadcast.
!! \param scope  Scope of the broadcast (default: "A").
!! \param top  Topology of the broadcast (default " ").
!! \see BLACS documentation (routine ?gebs2d).
subroutine blacsfx_gebs_$2`1'(mygrid, aa, scope, top)
  class(blacsgrid), intent(in) :: mygrid
  $3, intent(in), target :: aa(:)
  character, intent(in), optional :: scope, top

  $3, pointer :: buffer(:,:)

  buffer(1:size(aa), 1:1) => aa(1:size(aa))
  call blacsfx_gebs(mygrid, buffer, scope, top)

end subroutine blacsfx_gebs_$2`1'
')

define(`_subroutine_blacsfx_gebs_r0',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Starts broadcast ($1, rank 0).
!! \param mygrid  BLACS descriptor.
!! \param aa  Scalar to broadcast.
!! \param scope  Scope of the broadcast (default: "A").
!! \param top  Topology of the broadcast (default " ").
!! \see BLACS documentation (routine ?gebs2d).
subroutine blacsfx_gebs_$2`0'(mygrid, aa, scope, top)
  class(blacsgrid), intent(in) :: mygrid
  $3, intent(in) :: aa
  character, intent(in), optional :: scope, top

  $3 :: buffer(1,1)

  buffer(1,1) = aa
  call blacsfx_gebs(mygrid, buffer, scope, top)
  
end subroutine blacsfx_gebs_$2`0'
')

dnl ****************************************************************************
dnl *** blacsfx_gebr
dnl ****************************************************************************

define(`_subroutine_blacsfx_gebr_r2',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Receives broadcast ($1, rank 2).
!! \param mygrid  BLACS descriptor.
!! \param aa  Matrix to receive.
!! \param scope  Scope of the broadcast (default: "A").
!! \param top  Topology of the broadcast (default: " ").
!! \param rsrc  Row of the source (default: row of master process).
!! \param csrc  Column of the source (default: column of master process).
!! \see BLACS documentation (routine ?gebr2d).
subroutine blacsfx_gebr_$2`2'(mygrid, aa, scope, top, rsrc, csrc)
  class(blacsgrid), intent(in) :: mygrid
  $3, intent(out) :: aa(:,:)
  character, intent(in), optional :: scope, top
  integer, intent(in), optional :: rsrc, csrc

  character :: scope0, top0
  integer :: rsrc0, csrc0

  _handle_inoptflag(scope0,scope, "A")
  _handle_inoptflag(top0, top, " ")
  _handle_inoptflag(rsrc0, rsrc, mygrid%masterrow)
  _handle_inoptflag(csrc0, csrc, mygrid%mastercol)
  call gebr2d(mygrid%ctxt, scope0, top0, size(aa, dim=1), size(aa, dim=2),&
      & aa, size(aa, dim=1), rsrc0, csrc0)

end subroutine blacsfx_gebr_$2`2'
')

define(`_subroutine_blacsfx_gebr_r1',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Receives broadcast ($1, rank 1).
!! \param mygrid  BLACS descriptor.
!! \param aa  Vector to receive.
!! \param scope  Scope of the broadcast (default: "A").
!! \param top  Topology of the broadcast (default: " ").
!! \param rsrc  Row of the source (default: row of master process).
!! \param csrc  Column of the source (default: column of master process).
!! \see BLACS documentation (routine ?gebr2d).
subroutine blacsfx_gebr_$2`1'(mygrid, aa, scope, top, rsrc, csrc)
  class(blacsgrid), intent(in) :: mygrid
  $3, intent(out), target :: aa(:)
  character, intent(in), optional :: scope, top
  integer, intent(in), optional :: rsrc, csrc

  $3, pointer :: buffer(:,:)

  buffer(1:size(aa), 1:1) => aa(1:size(aa))
  call blacsfx_gebr(mygrid, buffer, scope, top, rsrc, csrc)
  
end subroutine blacsfx_gebr_$2`1'
')

define(`_subroutine_blacsfx_gebr_r0',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Receives broadcast ($1, rank 0).
!! \param mygrid  BLACS descriptor.
!! \param aa  Scalar to receive.
!! \param scope  Scope of the broadcast (default: "A").
!! \param top  Topology of the broadcast (default: " ").
!! \param rsrc  Row of the source (default: row of master process).
!! \param csrc  Column of the source (default: column of master process).
!! \see BLACS documentation (routine ?gebr2d).
subroutine blacsfx_gebr_$2`0'(mygrid, aa, scope, top, rsrc, csrc)
  class(blacsgrid), intent(in) :: mygrid
  $3, intent(out) :: aa
  character, intent(in), optional :: scope, top
  integer, intent(in), optional :: rsrc, csrc

  $3 :: buffer(1,1)

  call blacsfx_gebr(mygrid, buffer, scope, top, rsrc, csrc)
  aa = buffer(1,1)

end subroutine blacsfx_gebr_$2`0'
')

dnl ****************************************************************************
dnl *** blacsfx_gesd
dnl ****************************************************************************

define(`_subroutine_blacsfx_gesd_r2',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Sends general rectangular matrix to destination process ($1, rank 2).
!! \param mygrid  BLACS descriptor.
!! \param aa  Object to send.
!! \param rdest  Row of the destination process.
!! \param cdest  Column of the destination proces.
!! \see BLACS documentation (routine ?gesd2d).
subroutine blacsfx_gesd_$2`2'(mygrid, aa, rdest, cdest)
  type(blacsgrid), intent(in) :: mygrid
  $3, intent(in) :: aa(:,:)
  integer, intent(in) :: rdest, cdest

  call gesd2d(mygrid%ctxt, size(aa, dim=1), size(aa, dim=2), aa,&
      & size(aa, dim=1), rdest, cdest)
  
end subroutine blacsfx_gesd_$2`2'
')

define(`_subroutine_blacsfx_gesd_r1',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Sends vector to destination process ($1, rank 1).
!! \param mygrid  BLACS descriptor
!! \param aa  Object to send.
!! \param rdest  Row of the destination process.
!! \param cdest  Column of the destination proces.
!! \see BLACS documentation (routine ?gesd2d).
subroutine blacsfx_gesd_$2`1'(mygrid, aa, rdest, cdest)
  type(blacsgrid), intent(in) :: mygrid
  $3, intent(in), target :: aa(:)
  integer, intent(in) :: rdest, cdest

  $3, pointer :: buffer(:,:)

  buffer(1:size(aa), 1:1) => aa(1:size(aa))
  call blacsfx_gesd(mygrid, buffer, rdest, cdest)
  
end subroutine blacsfx_gesd_$2`1'
')

define(`_subroutine_blacsfx_gesd_r0',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Sends scalar to destination process ($1, rank 0).
!! \param mygrid  BLACS descriptor.
!! \param aa  Object to send.
!! \param rdest  Row of the destination process.
!! \param cdest  Column of the destination proces.
!! \see BLACS documentation (routine ?gesd2d).
subroutine blacsfx_gesd_$2`0'(mygrid, aa, rdest, cdest)
  type(blacsgrid), intent(in) :: mygrid
  $3, intent(in) :: aa
  integer, intent(in) :: rdest, cdest

  $3 :: buffer(1,1)

  buffer(1,1) = aa
  call blacsfx_gesd(mygrid, buffer, rdest, cdest)

end subroutine blacsfx_gesd_$2`0'
')

dnl ****************************************************************************
dnl *** blacsfx_gerv
dnl ****************************************************************************

define(`_subroutine_blacsfx_gerv_r2',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Receives general rectangular matrix from source process ($1, rank 2).
!! \param mygrid  BLACS descriptor
!! \param aa  Object to receive.
!! \param rdest  Row of the destination process (default: master row).
!! \param cdest  Column of the destination proces (default: master col).
!! \see BLACS documentation (routine ?gerv2d).
subroutine blacsfx_gerv_$2`2'(mygrid, aa, rsrc, csrc)
  type(blacsgrid), intent(in) :: mygrid
  $3, intent(out) :: aa(:,:)
  integer, intent(in), optional :: rsrc, csrc

  integer :: rsrc0, csrc0

  _handle_inoptflag(rsrc0, rsrc, mygrid%masterrow)
  _handle_inoptflag(csrc0, csrc, mygrid%mastercol)
  call gerv2d(mygrid%ctxt, size(aa, dim=1), size(aa, dim=2), aa,&
      & size(aa, dim=1), rsrc0, csrc0)

end subroutine blacsfx_gerv_$2`2'
')

define(`_subroutine_blacsfx_gerv_r1',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Receives general vector from source process ($1, rank 1).
!! \param mygrid  BLACS descriptor
!! \param aa  Object to receive.
!! \param rdest  Row of the destination process (default: master row).
!! \param cdest  Column of the destination proces (default: master col).
!! \see BLACS documentation (routine ?gerv2d).
subroutine blacsfx_gerv_$2`1'(mygrid, aa, rsrc, csrc)
  type(blacsgrid), intent(in) :: mygrid
  $3, intent(out), target :: aa(:)
  integer, intent(in), optional :: rsrc, csrc

  $3, pointer :: buffer(:,:)

  buffer(1:size(aa), 1:1) => aa(1:size(aa))
  call blacsfx_gerv(mygrid, buffer, rsrc, csrc)

end subroutine blacsfx_gerv_$2`1'
')

define(`_subroutine_blacsfx_gerv_r0',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Receives scalar from source process ($1, rank 0).
!! \param mygrid  BLACS descriptor.
!! \param aa  Object to receive.
!! \param rdest  Row of the destination process (default: master row).
!! \param cdest  Column of the destination proces (default: master col).
!! \see BLACS documentation (routine ?gerv2d).
subroutine blacsfx_gerv_$2`0'(mygrid, aa, rsrc, csrc)
  type(blacsgrid), intent(in) :: mygrid
  $3, intent(out) :: aa
  integer, intent(in), optional :: rsrc, csrc

  $3 :: buffer(1,1)

  call blacsfx_gerv(mygrid, buffer, rsrc, csrc)
  aa = buffer(1,1)

end subroutine blacsfx_gerv_$2`0'
')

dnl ****************************************************************************
dnl *** blacsfx_gsum
dnl ****************************************************************************

define(`_subroutine_blacsfx_gsum_r2',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Performs element-wise summation($1, rank 2).
!! \param mygrid  BLACS descriptor.
!! \param aa  Matrix to sum up.
!! \param scope  Scope of the broadcast (default: "A").
!! \param top  Topology of the broadcast (default: " ").
!! \param rdest  Row of the destination (default: row of master process).
!! \param rcol  Column of the destination (default: column of master process).
!! \see BLACS documentation (routine ?gsum2d).
subroutine blacsfx_gsum_$2`2'(mygrid, aa, scope, top, rdest, cdest)
  class(blacsgrid), intent(in) :: mygrid
  $3, intent(inout) :: aa(:,:)
  character, intent(in), optional :: scope, top
  integer, intent(in), optional :: rdest, cdest

  character :: scope0, top0
  integer :: rdest0, cdest0

  _handle_inoptflag(scope0, scope, "A")
  _handle_inoptflag(top0, top, " ")  
  _handle_inoptflag(rdest0, rdest, mygrid%masterrow)
  _handle_inoptflag(cdest0, cdest, mygrid%mastercol)
  call gsum2d(mygrid%ctxt, scope0, top0, size(aa, dim=1), size(aa, dim=2), aa,&
      & size(aa, dim=1), rdest0, cdest0)
  
end subroutine blacsfx_gsum_$2`2'
')

define(`_subroutine_blacsfx_gsum_r1',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Performs element-wise summation($1, rank 1).
!! \param mygrid  BLACS descriptor.
!! \param aa  Vector to sum up.
!! \param scope  Scope of the broadcast (default: "A").
!! \param top  Topology of the broadcast (default: " ").
!! \param rdest  Row of the destination (default: row of master process).
!! \param rcol  Column of the destination (default: column of master process).
!! \see BLACS documentation (routine ?gsum2d).
subroutine blacsfx_gsum_$2`1'(mygrid, aa, scope, top, rdest, cdest)
  class(blacsgrid), intent(in) :: mygrid
  $3, intent(inout), target :: aa(:)
  character, intent(in), optional :: scope, top
  integer, intent(in), optional :: rdest, cdest

  $3, pointer :: buffer(:,:)

  buffer(1:size(aa), 1:1) => aa(1:size(aa))
  call blacsfx_gsum(mygrid, buffer, scope, top, rdest, cdest)

end subroutine blacsfx_gsum_$2`1'
')

define(`_subroutine_blacsfx_gsum_r0',`
dnl $1 comment
dnl $2 subroutine suffix
dnl $3 dummy argument type
!> Performs element-wise summation($1, rank 0).
!! \param mygrid  BLACS descriptor
!! \param aa  Scalar to sum up.
!! \param scope  Scope of the broadcast (default: "A").
!! \param top  Topology of the broadcast (default: " ").
!! \param rdest  Row of the destination (default: row of master process).
!! \param rcol  Column of the destination (default: column of master process).
!! \see BLACS documentation (routine ?gsum2d).
subroutine blacsfx_gsum_$2`0'(mygrid, aa, scope, top, rdest, cdest)
  class(blacsgrid), intent(in) :: mygrid
  $3, intent(inout) :: aa
  character, intent(in), optional :: scope, top
  integer, intent(in), optional :: rdest, cdest

  $3 :: buffer(1,1)

  call blacsfx_gsum(mygrid, buffer, scope, top, rdest, cdest)
  aa = buffer(1,1)

end subroutine blacsfx_gsum_$2`0'
')

