! Subroutine:
!    initgrid
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
subroutine initgrid(gind,gind0,gx,gy,glat,glong,ng,hem,ndeg)

    implicit none
    
!======================================================================
! Declarations & Initialisation
!======================================================================

    ! Input variables
    !------------------------------------------------------------------
    integer, intent(in) :: ng
    integer, intent(in) :: ndeg
    real,    intent(in) :: hem

    ! Output variables
    !------------------------------------------------------------------
    integer(kind=1), dimension(-ng:ng,-ng:ng), intent(out) :: gind
    integer(kind=1), dimension(-ng:ng,-ng:ng), intent(out) :: gind0
    real,            dimension(-ng:ng,-ng:ng), intent(out) :: gx
    real,            dimension(-ng:ng,-ng:ng), intent(out) :: gy
    real,            dimension(-ng:ng,-ng:ng), intent(out) :: glat
    real,            dimension(-ng:ng,-ng:ng), intent(out) :: glong

    ! Local variables
    integer          :: ig
    integer          :: jg
    double precision :: radeg

!======================================================================
! Program
!======================================================================

    radeg = 45./atan(1.)

    do jg = -ng,ng
        do ig = -ng,ng
            gx(ig,jg) = 0.
            gy(ig,jg) = 0.
            if (ig*ig+jg*jg > (ng+1)*(ng+1)) then
                gind0(ig,jg) = 0
            else
                gind0(ig,jg) = 1
            endif

            gind(ig,jg) = gind0(ig,jg)

            if (gind0(ig,jg) == 1) then
                gx(ig,jg) = real(ig)
                gy(ig,jg) = real(jg)

                if (ig /= 0 .or. jg /= 0) then
                    glat(ig,jg) = hem*(90.-sqrt(abs(real(ig**2+jg**2)))/real(ndeg))
                    glong(ig,jg) = atan2(float(ig),-hem*real(jg))*radeg
                else
                    glong(ig,jg) = 0.
                    glat(ig,jg) = 90.*hem
                endif

                if (glong(ig,jg) < 0.) then
                    glong(ig,jg) = glong(ig,jg)+360.
                endif
            endif
        enddo
    enddo

    return

end subroutine initgrid
