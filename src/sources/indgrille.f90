! Function:
!    indgrille
!
! Purpose:
!    none
!
! Input:
!    none
!
! Output:
!    none
!
! Restrictions:
!    none
!
! Comments:
!    none
!
! History
!    MIMOSA was written by A. Hauchecorne
!    201004   M.-A. DROUIN/LPC2E conversion into FORTRAN 95
!    201103   M.-A. DROUIN/LPC2E add comments
!
!======================================================================
integer(kind=1) function indgrille(i,j,gind0,ng)

    implicit none

!======================================================================
! Declarations & Initialisation
!======================================================================

    ! Input variables
    !------------------------------------------------------------------
    integer,                                   intent(in) :: i
    integer,                                   intent(in) :: j
    integer,                                   intent(in) :: ng
    integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind0

!======================================================================
! Program
!======================================================================

    indgrille = 1

    if (i < -ng .or. i >= ng) then
        indgrille = 0
    endif
    if (j < -ng .or. j >= ng) then
        indgrille = 0
    endif
    if (indgrille == 1) then
        if (gind0(i,j) == 0) then
            indgrille = 0
        endif
        if (gind0(i+1,j) == 0) then
            indgrille = 0
        endif
        if (gind0(i,j+1) == 0) then
            indgrille = 0
        endif
        if (gind0(i+1,j+1) == 0) then
            indgrille = 0
        endif
    endif

    return

end function indgrille
