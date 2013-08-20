include(common.m4)

define(`_subroutine_interface_gebs2d',`
dnl $1: comment
dnl $2: subroutine prefix
dnl $3: dummy arguments type
!> Starts broadcast for general rectangular matrix ($1).
!! \see BLACS documentation for details.
subroutine $2gebs2d(ictxt, scope, top, mm, nn, aa, lda)
  import
  integer, intent(in) :: ictxt
  character, intent(in) :: scope, top
  integer, intent(in) :: mm, nn
  integer, intent(in) :: lda
  $3, intent(in) :: aa(lda,*)
end subroutine $2gebs2d
')

define(`_subroutine_interface_gebr2d',`
dnl $1: comment
dnl $2: subroutine prefix
dnl $3: dummy arguments type
!> Receives broadcast for general rectangular matrix ($1).
!! \see BLACS documentation for details.
subroutine $2gebr2d(ictxt, scope, top, mm, nn, aa, lda, rsrc, csrc)
  import
  integer, intent(in) :: ictxt
  character, intent(in) :: scope, top
  integer, intent(in) :: mm, nn
  integer, intent(in) :: lda
  $3, intent(out) :: aa(lda,*)
  integer, intent(in) :: rsrc, csrc
end subroutine $2gebr2d
')dnl

define(`_subroutine_interface_gesd2d',`
dnl $1: comment
dnl $2: subroutine prefix
dnl $3: dummy arguments type
!> Sends general rectangular matrix to destination ($1).
!! \see BLACS documentation for details.
subroutine $2gesd2d(ictxt, mm, nn, aa, lda, rdest, cdest)
  import
  integer, intent(in) :: ictxt, mm, nn
  integer, intent(in) :: lda, rdest, cdest
  $3, intent(in) :: aa(lda,*)
end subroutine $2gesd2d
')dnl

define(`_subroutine_interface_gerv2d',`
dnl $1: comment
dnl $2: subroutine prefix
dnl $3: dummy arguments type
!> Receives general rectangular matrix from process ($1).
!! \see BLACS documentation for details.
subroutine $2gerv2d(ictxt, mm, nn, aa, lda, rsrc, csrc)
  import
  integer, intent(in) :: ictxt, mm, nn
  integer, intent(in) :: lda, rsrc, csrc
  $3, intent(out) :: aa(lda,*)
end subroutine $2gerv2d
')

define(`_subroutine_interface_gsum2d',`
dnl $1: comment
dnl $2: subroutine prefix
dnl $3: dummy arguments type
!> Performs element-wise summation ($1).
!! \see BLACS documentation for details.
subroutine $2gsum2d(ictxt, scope, top, mm, nn, aa, lda, rdest, cdest)
  import
  integer, intent(in) :: ictxt
  character, intent(in) :: scope, top
  integer, intent(in) :: mm, nn
  integer, intent(in) :: lda
  $3, intent(inout) :: aa(lda,*)
  integer, intent(in) :: rdest, cdest
end subroutine $2gsum2d
')
