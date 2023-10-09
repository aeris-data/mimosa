! Subroutine:
!    raccord
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
subroutine raccord(gpv,gind,gind0,gxhem,gyhem,ng,lathem,ndeg,nlis2D)

    implicit none

    !======================================================================
    ! Declarations & Initialisation
    !======================================================================

    ! Input variables
    !------------------------------------------------------------------
    integer, intent(in) :: lathem
    integer, intent(in) :: ndeg
    integer, intent(in) :: ng
    integer, intent(in) :: nlis2D
    integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind0
    real, dimension(-ng:ng,-ng:ng), intent(in) :: gxhem
    real, dimension(-ng:ng,-ng:ng), intent(in) :: gyhem

    ! Input/Output variables
    !------------------------------------------------------------------
    real, dimension(-ng:ng,-ng:ng,2), intent(inout) :: gpv
    integer(kind=1), dimension(-ng:ng,-ng:ng,2), intent(inout) :: gind

    ! Local variables
    !------------------------------------------------------------------
    integer:: i
    integer:: ig
    integer:: ihem
    integer:: j
    integer:: jg
    integer:: nvide
    integer:: igxh
    integer:: jgyh
    real:: fgxh
    real:: fgyh

    !======================================================================
    ! Program
    !======================================================================

    nvide = 0
    do ihem = 1,2
        do j = -ng,ng
            do i = -ng,ng
                if (gind0(i,j) == 1) then

                    if (gind(i,j,ihem) == 0 .or. i*i+j*j > (lathem*ndeg-nlis2D)**2) then
                        nvide = nvide+1
                        igxh = nint(gxhem(i,j)-.5)
                        fgxh = gxhem(i,j)-igxh
                        jgyh = nint(gyhem(i,j)-.5)
                        fgyh = gyhem(i,j)-jgyh

                        if ( gind(igxh,jgyh,3-ihem)+gind(igxh+1,jgyh,3-ihem) &
                        + gind(igxh,jgyh+1,3-ihem)+gind(igxh+1,jgyh+1,3-ihem) /= 4 ) then
                            print *, 'ATTENTION VALEUR NULLE 2 HEMISPHERES',ihem,i,j
                        endif

                        gpv(i,j,ihem) = gpv(igxh,jgyh,3-ihem)*(1.-fgxh)*(1.-fgyh) &
                        + gpv(igxh+1,jgyh,3-ihem)*fgxh*(1.-fgyh) &
                        + gpv(igxh,jgyh+1,3-ihem)*(1.-fgxh)*fgyh &
                        + gpv(igxh+1,jgyh+1,3-ihem)*fgxh*fgyh
                    endif
                endif
            enddo
        enddo
    enddo

    print *,'nbre pts vides = ', nvide

    do jg=-ng,ng
        do ig=-ng,ng
            if (gind0(ig,jg) == 1) then
                gind(ig,jg,1) = 1
                gind(ig,jg,2) = 1
            endif
        enddo
    enddo

    return

end subroutine raccord
