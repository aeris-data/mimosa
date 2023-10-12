! Function:
!    njour
!
! Purpose:
!    calcule le numero du jour (le 1 1 1901 est le jour 1)
!   si l'annee est avant 51, on est apres 2000
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
!    201109   M.-A. DROUIN/LPC2E add comments
!
!======================================================================
integer function njour(j, m, ia)

    implicit none

    !==================================================================
    ! Declarations & Initialisation
    !==================================================================

    ! Input variables
    !------------------------------------------------------------------
    integer, intent(in) :: ia
    integer, intent(in) :: j
    integer, intent(in) :: m

    ! Local variables
    !------------------------------------------------------------------
    integer :: mm
    integer :: nbis
    integer, dimension(13) :: ical

    data ical/0,31,59,90,120,151,181,212,243,273,304,334,365/

    !==================================================================
    ! Program
    !==================================================================
    mm = m + ia*12
    nbis = (mm-3)/48

    if (mm < 3) then
        nbis = -1
    endif

    njour = j+ical(m)+365*(ia-1)+nbis
    if (ia < 51) then
        njour = njour+36525
    endif

    return

end function njour
