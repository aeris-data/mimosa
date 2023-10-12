! Subroutine:
!    xyinterp
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
subroutine xyinterp(inter,wxp,wyp,xxp,yyp,xxinterp,yyinterp)

	implicit none

  	!======================================================================
   	! Declarations & Initialisation
   	!======================================================================

   	! Input variables
   	!------------------------------------------------------------------
	real, intent(in) :: xxp
	real, intent(in) :: yyp
	real, intent(in) :: wxp
	real, intent(in) :: wyp

    ! Output variables
    !------------------------------------------------------------------
    integer, intent(out) :: inter
    real, intent(out) :: xxinterp
    real, intent(out) :: yyinterp

   	! Local variables
   	!------------------------------------------------------------------
	integer:: i
	real:: ax
	real:: ay
	real:: xx0
	real:: xx1
	real:: yy0
	real:: yy1
	real:: eps
	real:: xymin

    !======================================================================
    ! Program
    !======================================================================

	inter = 0
	eps = 1.e-6
 	xymin = .1

	if ((xxp > -eps) .and. &
		(yyp > -eps) .and. &
		((wxp-1.)*yyp-(xxp-1.)*wyp > -eps) .and. &
		(wxp*(yyp-1.)-xxp*(wyp-1.) < eps)) then

		inter = 1
		ax = wxp-1.
		ay = wyp-1.

		if (max(abs(ax),abs(ay)) > xymin) then
			if (abs(ay) > abs(ax)) then
				xxinterp = (-1.+ay*xxp-ax*yyp+sqrt((-1.+ay*xxp-ax*yyp)**2 &
						 + 4.*ay*xxp))/(2.*ay)
				yyinterp = yyp/(1.+ay*xxinterp)
			else
				yyinterp = (-1.+ax*yyp-ay*xxp+sqrt((-1.+ax*yyp-ay*xxp)**2 &
						 + 4.*ax*yyp))/(2.*ax)
				xxinterp = xxp/(1.+ax*yyinterp)
			endif
		else
			xx0 = xxp
			yy0 = yyp

			i = 1
			xx1 = xxp-ax*xx0*yy0
			yy1 = yyp-ay*xx1*yy0

			do while(abs(xx1-xx0) > eps .and. abs(yy1-yy0) > eps .and. i <= 10)
				xx0 = xx1
				yy0 = yy1
				xx1 = xxp-ax*xx0*yy0
				yy1 = yyp-ay*xx1*yy0
				i = i+1

			enddo

			xxinterp = xx1
			yyinterp = yy1

		endif
	endif

	return

end
