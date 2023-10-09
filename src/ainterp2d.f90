! Function:
!    ainterp2d
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
!    201111   M.-A. DROUIN/LPC2E Add condition to prevent index to exceed a array bounds
!
!======================================================================
real function ainterp2D(a,glong,glat,nx,ny,ig,jg,ng,paslong,latmin,paslat)

    implicit none

!======================================================================
! Declarations & Initialisation
!======================================================================

    ! Input variables
    !------------------------------------------------------------------
    integer,                        intent(in) :: ng
	real, dimension(-ng:ng,-ng:ng), intent(in) :: glong
	real, dimension(-ng:ng,-ng:ng), intent(in) :: glat
	real, dimension(nx,ny),         intent(in) :: a
    integer,                        intent(in) :: nx
    integer,                        intent(in) :: ny
    integer,                        intent(in) :: ig
    integer,                        intent(in) :: jg
	real,                           intent(in) :: paslong
	real,                           intent(in) :: paslat
	real,                           intent(in) :: latmin

    ! Local variables
    !------------------------------------------------------------------
	real    :: fracx
	real    :: fracy
	real    :: x
	real    :: y
    integer :: intx
    integer :: intx1
    integer :: inty
    integer :: inty1 ! 201111 MAD

!======================================================================
! Program
!======================================================================

    x = gclong(ig,jg,glong,ng,paslong)
    intx = int(x)
    fracx = x-intx
    y = gclat(ig,jg,glat,ng,latmin,paslat)
    inty = int(y)
    fracy = y-inty
    intx1 = intx+1
    inty1 = inty+1 ! 201111 MAD

    if (intx1 > nx) then
        intx1 = intx1-nx
    end if
    if (inty1 > ny) then !
        inty1 = inty-1   ! 201111 MAD
    end if               !

    ainterp2D = a(intx,inty)*(1.-fracx)*(1.-fracy) + &
                a(intx1,inty)*fracx*(1.-fracy) + &
                a(intx,inty1)*(1.-fracx)*fracy + &
                a(intx1,inty1)*fracx*fracy

!======================================================================
! Internal functions
!======================================================================
    contains

    real function gclat(ig,jg,glat,ng,latmin,paslat)

        ! Input variables
        !--------------------------------------------------------------
        integer,                        intent(in) :: ig
        integer,                        intent(in) :: jg
        integer,                        intent(in) :: ng
		real,                           intent(in) :: paslat
		real,                           intent(in) :: latmin
		real, dimension(-ng:ng,-ng:ng), intent(in) :: glat

        gclat=1.+(glat(ig,jg)-latmin)/paslat

        return

    end function gclat

    real function gclong(ig,jg,glong,ng,paslong)

        ! Input variables
        !--------------------------------------------------------------
        integer,                        intent(in) :: ig
        integer,                        intent(in) :: jg
        integer,                        intent(in) :: ng
		real,                           intent(in) :: paslong
		real, dimension(-ng:ng,-ng:ng), intent(in) :: glong

        gclong = 1.+glong(ig,jg)/paslong

        return

    end function gclong

end function ainterp2D
