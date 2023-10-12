! Subroutine:
!    regrid2
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
!    201109   M.-A. DROUIN/LPC2E add comments
!    201306   C. BOONNE/IPSL add nhmod
!    202302   C. BOONNE/IPSL migration SPIRIT1 ajout 
!             use interfaces_mod, except_this_one => regrid2 
!
!======================================================================
subroutine regrid2( gpv, gpv0, gpvn, gpoids, gx, gy, gind, gind0 &
                  , ng, nhrelax, nhmod, indrelax, dgpv, dgpvlis, nlis2D)

    use interfaces_mod, except_this_one => regrid2

    implicit none

    !======================================================================
    ! Declarations & Initialisation
    !======================================================================

    ! Input variables
    !------------------------------------------------------------------
    integer, intent(in) :: indrelax
    integer, intent(in) :: ng
    integer, intent(in) :: nhrelax
    integer, intent(in) :: nhmod
    integer, intent(in) :: nlis2D
    real, dimension(-ng:ng,-ng:ng), intent(in) :: gpv0
    integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind0

    ! Output variables
    !------------------------------------------------------------------
    integer(kind=2), dimension(-ng:ng,-ng:ng), intent(out) :: gpoids
    real, dimension(-ng:ng,-ng:ng), intent(out) :: dgpv
    real, dimension(-ng:ng,-ng:ng), intent(out) :: dgpvlis

    ! Input/Output variables
    !------------------------------------------------------------------
    real, dimension(-ng:ng,-ng:ng), intent(inout) :: gpv
    real, dimension(-ng:ng,-ng:ng), intent(inout) :: gpvn
    real, dimension(-ng:ng,-ng:ng), intent(inout) :: gx
    real, dimension(-ng:ng,-ng:ng), intent(inout) :: gy
    integer(kind=1), dimension(-ng:ng,-ng:ng) :: gind

    ! Local variables
    !------------------------------------------------------------------
    integer :: gxmin
    integer :: gxmax
    integer :: gymin
    integer :: gymax
    integer :: ibord
    integer :: inddif
    integer :: inter
    integer :: nmoy
    integer :: ii
    integer :: ig
    integer :: igg
    integer :: jj
    integer :: jg
    integer :: jgg

    real :: eps
    real :: poidsexp
    real :: poidsmin
    real :: poidsmod
    real :: axmax
    real :: axmoy
    real :: aymax
    real :: aymoy
    real :: ux
    real :: uy
    real :: vx
    real :: vy
    real :: wx
    real :: wxp
    real :: wy
    real :: wyp
    real :: alx
    real :: alx2
    real :: ax
    real :: aly
    real :: aly2
    real :: ay
    real :: xx
    real :: xxp
    real :: xxinterp
    real :: yy
    real :: yyp
    real :: yyinterp

    !======================================================================
    ! Program
    !======================================================================

    eps = 1.e-6
    poidsmin = 1.e-6

    !	inddif = 0 pas de diffusion, inddif = 1 diffusion
    inddif = 0

    if (indrelax == 1) then
!       poidsmod = 24./real(nhrelax)
!   Changement A. Hauchecorne au 31/01/2008
        poidsmod=real(nhmod)/real(nhrelax)
    else
        poidsmod = 1.e-6
    endif

    do jg = -ng,ng
        do ig = -ng,ng
            gpoids(ig,jg) = 0
        enddo
    enddo

    axmax = 0.; aymax = 0.; axmoy = 0.; aymoy = 0.
    nmoy = 0

    do jg = -ng,ng-1
        do ig=-ng,ng-1
            if (gind(ig,jg)+gind(ig+1,jg)+gind(ig,jg+1)+gind(ig+1,jg+1) ==4 ) then

                if ( inddif == 0) then

                    ibord = 0
                    if ((ig /= -ng) .and. (ig /= ng-1) .and. &
                    (jg /= -ng) .and. (jg /= ng-1)) then

                        if ( gind(ig-1,jg-1)+gind(ig,jg-1)+gind(ig+1,jg-1) &
                        + gind(ig+2,jg-1)+gind(ig-1,jg)+gind(ig+2,jg) &
                        + gind(ig-1,jg+1)+gind(ig+2,jg+1)+gind(ig-1,jg+2) &
                        + gind(ig,jg+2)+gind(ig+1,jg+2)+gind(ig+2,jg+2) == 12) then

                            ibord = 1
                        endif
                    endif
                endif

                ux = gx(ig+1,jg)-gx(ig,jg)
                uy = gy(ig+1,jg)-gy(ig,jg)
                vx = gx(ig,jg+1)-gx(ig,jg)
                vy = gy(ig,jg+1)-gy(ig,jg)
                wx = gx(ig+1,jg+1)-gx(ig,jg)
                wy = gy(ig+1,jg+1)-gy(ig,jg)
                wxp = (wx*vy-wy*vx)/(ux*vy-uy*vx)
                wyp = (wy*ux-wx*uy)/(ux*vy-uy*vx)
                alx2 = ux**2+uy**2
                aly2 = vx**2+vy**2
                alx = sqrt(alx2)
                aly = sqrt(aly2)
                axmoy = axmoy+(alx-1)**2
                aymoy = aymoy+(aly-1)**2
                nmoy = nmoy+1
                ax = wxp-1
                ay = wyp-1

                if (alx > axmax) then
                    axmax = alx
                endif
                if (aly > aymax) then
                    aymax = aly
                endif

                gxmin = min(nint(gx(ig,jg)+.5-eps),nint(gx(ig,jg+1)+.5-eps))
                gxmin = max(gxmin,-ng)
                gxmax = max(nint(gx(ig+1,jg)-.5+eps),nint(gx(ig+1,jg+1)-.5+eps))
                gxmax = min(gxmax,ng)
                gymin = min(nint(gy(ig,jg)+.5-eps),nint(gy(ig+1,jg)+.5-eps))
                gymin = max(gymin,-ng)
                gymax = max(nint(gy(ig,jg+1)-.5+eps),nint(gy(ig+1,jg+1)-.5+eps))
                gymax = min(gymax,ng)

                do jj = gymin,gymax
                    do ii = gxmin,gxmax
                        if (gpoids(ii,jj) /= 1) then
                            xx = real(ii)-gx(ig,jg)
                            yy = real(jj)-gy(ig,jg)
                            xxp = (xx*vy-yy*vx)/(ux*vy-uy*vx)
                            yyp = (yy*ux-xx*uy)/(ux*vy-uy*vx)

                            call xyinterp(inter,wxp,wyp,xxp,yyp,xxinterp,yyinterp)
                            if (inter == 1) then
                                gpoids(ii,jj) = 1
                                if (inddif == 1) then
                                    gpvn(ii,jj) = gpv(ig,jg)*(1.-xxinterp)*(1.-yyinterp) &
                                    + gpv(ig+1,jg)*xxinterp*(1.-yyinterp) &
                                    + gpv(ig,jg+1)*(1.-xxinterp)*yyinterp &
                                    + gpv(ig+1,jg+1)*xxinterp*yyinterp
                                else
                                    gpvn(ii,jj) = 0.
                                    do jgg = jg-ibord,jg+ibord+1
                                        do igg = ig-ibord,ig+ibord+1
                                            gpvn(ii,jj) = gpvn(ii,jj)+poids(xxinterp,igg-ig,ibord) &
                                            * poids(yyinterp,jgg-jg,ibord)*gpv(igg,jgg)
                                        enddo
                                    enddo
                                endif
                            endif
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo

    print '("nmoy  = ",I7)', nmoy
    print '("axmoy = ",E12.7," aymoy = ",E12.7)', sqrt(axmoy/nmoy), sqrt(aymoy/nmoy)
    print '("axmax = ",E12.7," aymax = ",E12.7)', axmax, aymax
    print *, ""

    do jg = -ng,ng
        do ig = -ng,ng
            if (gind0(ig,jg) == 1) then
                if (gpoids(ig,jg) == 1) then
                    gpvn(ig,jg) = (gpvn(ig,jg)+gpv0(ig,jg) &
                    *  poidsmin)/(1.+poidsmin)
                else
                    gpvn(ig,jg) = gpv0(ig,jg)
                endif
                if (indrelax /= 2) then
                    if (indrelax == 0) then
                        gpv(ig,jg) = gpvn(ig,jg)
                    else if (indrelax == 1) then
                        dgpv(ig,jg) = gpv0(ig,jg)-gpvn(ig,jg)
                        ! if (ig == jg) then
                            ! print *, ig, jg, dgpv(ig,jg), gpv0(ig,jg), gpvn(ig,jg)
                        ! endif
                    endif
                    gx(ig,jg) = real(ig)
                    gy(ig,jg) = real(jg)
                    gind(ig,jg) = 1
                endif
            else
                gpv(ig,jg) = -99.9
            endif
        enddo
    enddo

    if (indrelax == 1) then

        call lisse2D(ng,dgpv,dgpvlis,nlis2D,gind0)
        poidsexp = 1.-exp(-poidsmod)

        do jg = -ng,ng
            do ig = -ng,ng
                if (gind0(ig,jg) == 1) then
                    gpv(ig,jg) = gpvn(ig,jg)+dgpvlis(ig,jg)*poidsexp
                    gpvn(ig,jg) = gpv(ig,jg)
                endif
            enddo
        enddo
    endif

    return

end subroutine regrid2
