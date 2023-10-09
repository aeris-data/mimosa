! Function:
!    poids
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
real function poids(x,i,ib)

    implicit none

    !======================================================================
    ! Declarations & Initialisation
    !======================================================================

    ! Input variables
    !------------------------------------------------------------------
    integer, intent(in) :: i
    integer, intent(in) :: ib
    real, intent(in) :: x

    !======================================================================
    ! Program
    !======================================================================

    if (ib == 1)then

        if (i == -1) then
            poids = -x*(x-1)*(x-2)/6.
        endif
        if (i == 0) then
            poids = (x-2)*(x-1)*(x+1)/2.
        endif
        if (i == 1) then
            poids = -x*(x-2)*(x+1)/2.
        endif
        if (i == 2) then
            poids = x*(x-1)*(x+1)/6.
        endif

    else if (ib == 0) then

        if (i == 0) then
            poids = 1-x
        endif
        if (i == 1) then
            poids = x
        endif

    endif

    return

end function poids
