! Subroutine:
!    deplacement
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
!    202302   C. BOONNE/IPSL migration SPIRIT1 ajout
!             use interfaces_mod, except_this_one => deplacement
!
!======================================================================
subroutine deplacement(gind0,gind,gu,gv,gx,gy,ng,nh)

    use interfaces_mod, except_this_one => deplacement

    implicit none

    !======================================================================
    ! Declarations & Initialisation
    !======================================================================

    ! Input variables
    !------------------------------------------------------------------
    integer,                                   intent(in) :: ng
    integer,                                   intent(in) :: nh
    integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind0
    real,            dimension(-ng:ng,-ng:ng), intent(in) :: gu
    real,            dimension(-ng:ng,-ng:ng), intent(in) :: gv

    ! Input/Output variables
    !------------------------------------------------------------------
    real,            dimension(-ng:ng,-ng:ng), intent(inout) :: gx
    real,            dimension(-ng:ng,-ng:ng), intent(inout) :: gy
    integer(kind=1), dimension(-ng:ng,-ng:ng), intent(inout) :: gind

    ! Local variables
    !------------------------------------------------------------------
    integer(kind=1) :: gind1
    integer(kind=1) :: gind2
    integer         :: intx
    integer         :: inty
    integer         :: intx2
    integer         :: inty2
    integer         :: ig
    integer         :: jg
    real            :: dt
    real            :: fracx
    real            :: fracx2
    real            :: fracy
    real            :: fracy2
    real            :: gx1
    real            :: gx2
    real            :: gy1
    real            :: gy2
    real            :: u1
    real            :: u2
    real            :: v1
    real            :: v2

    !======================================================================
    ! Program
    !======================================================================

    dt = 1./real(nh)

    do jg = -ng,ng
        do ig = -ng,ng

            ! Wind at starting point
            if (gind(ig,jg) == 1) then

                intx = nint(gx(ig,jg)-.5)
                fracx = gx(ig,jg)-intx
                inty = nint(gy(ig,jg)-.5)
                fracy = gy(ig,jg)-inty
                gind1 = indgrille(intx,inty,gind0,ng)

                if (gind1 == 1) then

                    u1 = gu(intx,inty)*(1.-fracx)*(1.-fracy) + &
                    gu(intx+1,inty)*fracx*(1.-fracy) + &
                    gu(intx,inty+1)*(1.-fracx)*fracy + &
                    gu(intx+1,inty+1)*fracx*fracy
                    v1 = gv(intx,inty)*(1.-fracx)*(1.-fracy) + &
                    gv(intx+1,inty)*fracx*(1.-fracy) + &
                    gv(intx,inty+1)*(1.-fracx)*fracy + &
                    gv(intx+1,inty+1)*fracx*fracy

                    gx1 = gx(ig,jg)+dt*u1
                    gy1 = gy(ig,jg)+dt*v1
                    intx = nint(gx1-.5)
                    fracx = gx1-intx
                    inty = nint(gy1-.5)
                    fracy = gy1-inty
                    gind1 = indgrille(intx,inty,gind0,ng)
                endif

                ! Wind at arrival point
                if (gind1 == 1) then

                    u2 = gu(intx,inty)*(1.-fracx)*(1.-fracy) &
                    + gu(intx+1,inty)*fracx*(1.-fracy)    &
                    + gu(intx,inty+1)*(1.-fracx)*fracy    &
                    + gu(intx+1,inty+1)*fracx*fracy
                    v2 = gv(intx,inty)*(1.-fracx)*(1.-fracy) &
                    + gv(intx+1,inty)*fracx*(1.-fracy)    &
                    + gv(intx,inty+1)*(1.-fracx)*fracy    &
                    + gv(intx+1,inty+1)*fracx*fracy

                    gx2 = gx(ig,jg)+dt*(u1+u2)*.5
                    gy2 = gy(ig,jg)+dt*(v1+v2)*.5
                    intx2 = nint(gx2-.5)
                    fracx2 = gx2-intx2
                    inty2 = nint(gy2-.5)
                    fracy2 = gy2-inty2
                    gind2 = indgrille(intx2,inty2,gind0,ng)
                else
                    gind2 = 0
                endif
            endif

            if (gind(ig,jg) == 1) then
                gind(ig,jg) = gind2

                if(gind2 == 1)then
                    gx(ig,jg) = gx2
                    gy(ig,jg) = gy2
                endif
            endif

        enddo
    enddo

    return

end subroutine deplacement
