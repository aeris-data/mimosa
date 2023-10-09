! Subroutine:
!    njouri
!
! Purpose:
!    calcule la date a partir du numero du jour
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
subroutine njouri(j,m,ia,nj)

    implicit none

    !======================================================================
    ! Declarations & Initialisation
    !======================================================================

    ! Input variables
    !------------------------------------------------------------------
 	integer, intent(in) :: nj

    ! Output variables
    !------------------------------------------------------------------
    integer, intent(out) :: ia
    integer, intent(out) :: j
    integer, intent(out) :: m

    ! Local variables
    !------------------------------------------------------------------
    integer:: mm
    integer:: nbis
    integer:: nj1
    integer, dimension(13) :: ical
    integer, dimension(13) :: ical1

 	data ical/0,31,59,90,120,151,181,212,243,273,304,334,365/

 	!======================================================================
    ! Program
    !======================================================================

	nj1  = nj
	nbis = (nj1-1)/1461
	nj1  = nj1-nbis*1461
	ia   = (nj1-1)/365+1

	if (ia == 5) then
		ia = 4
	endif

	nj1 = nj1-365*(ia-1)
	ia = ia + nbis*4
	if (ia >= 100) then
		ia = ia-100
	endif

	do mm = 1,13
		ical1(mm) = ical(mm)
		if ( mod(ia,4) == 0 .and. mm >= 3) then
	 		ical1(mm) = ical1(mm)+1
		endif
	enddo

	m = 0
	m = m +1
	do while (nj1 > ical1(m+1))
		m = m+1
	enddo
	j = nj1-ical1(m)

	return

end subroutine njouri
