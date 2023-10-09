! Subroutine:
!    pvp
!
! Purpose:
!    Routine pour le calcule de la vorticite potentielle
!   sur un  surface isentrope a partir de champs 3D
!   a pas constants en lat, long de U,V et T en niveaux de pression
!
! Input:
!    theta0  niveau isentrope en K  Input
!    pres    tableau des pressions en hPa les pressions sont en ordre croissant (altitudes croissantes)
!    nx      dimension tableaux en x
!    ny      dimension tableaux en y
!    np      nbre niveaux de pression dimesion tableaux
!    np0     nbre niveaux pression calcul en principe np=np0
!    ul      champ 3D U en m/s
!    vl      champ 3D V  en m/s
!    tl      champ 3D T  en K
!    alatmin latitude mini des champs
!    alatmax latitude maxi des champs
!
! Output:
!    thetal  champ 3D Theta en K
!    u       champ U sur surface theta
!    v       champ V sur surface theta
!    t       champ T sur surface theta
!    p       champ P sur surface theta
!    pv      champ PV sur surface theta
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
!    2012 --> 2022 C. BOONNE Maintenance code
!    202212   C. BOONNE Migration sur le cluter SPIRIT IPSL
!             5.4.0-148-generic #165-Ubuntu
!
!======================================================================
subroutine pvp( theta0, nx, ny, np, np0 &
			  , pl, ul, vl, tl, thetal, u, v, t, p, dtdp, pv &
			  , alatmin, alatmax)

	implicit none

	!======================================================================
	! Declarations & Initialisation
	!======================================================================

	! Input variables
	!------------------------------------------------------------------
	integer, intent(in) :: nx
	integer, intent(in) :: ny
	integer, intent(in) :: np
	integer, intent(in) :: np0
	real, intent(in) :: theta0
	real, intent(in) :: alatmin
	real, intent(in) :: alatmax
	real, dimension(nx,ny,np), intent(in) :: ul
	real, dimension(nx,ny,np), intent(in) :: vl
	real, dimension(nx,ny,np), intent(in) :: tl
	real, dimension(nx,ny,np), intent(in) :: pl

	! Output variables
	!------------------------------------------------------------------
	real, dimension(nx,ny), intent(out) :: u
	real, dimension(nx,ny), intent(out) :: v
	real, dimension(nx,ny), intent(out) :: t
	real, dimension(nx,ny), intent(out) :: p
	real, dimension(nx,ny), intent(out) :: pv
	real, dimension(nx,ny), intent(out) :: dtdp
	real, dimension(nx,ny,np), intent(out) :: thetal

    ! Local variables
    !------------------------------------------------------------------
	integer :: i
	integer :: j
	integer :: k
	integer :: im1
	integer :: ip1
	integer :: j1
	real :: g
	real :: gamma
	real :: p0
	real :: rt
	real :: dlambda
	real :: dphi
	real :: dx
	real :: dy
	real :: phi
	real :: utphi
	real :: alfa
	real :: beta
	real :: dltdlp
	real :: dudy
	real :: dvdx
	real :: f
	real :: pvmoy
	real :: pi
	real :: pio180

    !======================================================================
    ! Program
    !======================================================================

    ! Constantes physiques
	pi = 4.*atan(1.)
	pio180 = pi / 180.
	gamma = 2./7.
	rt = 6.37e6
	g = 9.81
	p0 = 1013.25
	dphi = (alatmax-alatmin)/real(ny-1)
	dlambda = 360./real(nx)
	dy = pi/180.*rt*dphi
	dx = pi/180.*rt*dlambda

    ! Calcul de theta sur la grille pression
	do k = 1,np0
		do j = 1,ny
			do i = 1,nx

				thetal(i,j,k) = tl(i,j,k)*(p0/pl(i,j,k))**gamma

				if (k /= 1) then
					if (thetal(i,j,k) < thetal(i,j,k-1)+1.e-3) then
						thetal(i,j,k) = thetal(i,j,k-1)+1.e-3
					endif
				endif
			enddo
		enddo
	enddo

    ! interpolation sur surface theta0
	do j = 1,ny
		do i = 1,nx
			k = 0
			k = k+1
			do while (k < np0-1 .and. thetal(i,j,k+1) < theta0)
				k = k+1
			enddo

			alfa = log(theta0/thetal(i,j,k))/ &
				   log(thetal(i,j,k+1)/thetal(i,j,k))

! Pb points manquants chgt A. Hauchecorne le 20/06/2013
                    if(alfa.lt.0.)then
!       print *,'premier niveau au dessus de theta0'
                      alfa=0.
                    else
                      if (alfa.gt.1.) then
!       print *,'dernier niveau en dessous de theta0'
                      alfa=1.
                      endif
                    endif

			u(i,j) = (1.-alfa)*ul(i,j,k)+alfa*ul(i,j,k+1)
			v(i,j) = (1.-alfa)*vl(i,j,k)+alfa*vl(i,j,k+1)
			t(i,j) = (1.-alfa)*tl(i,j,k)+alfa*tl(i,j,k+1)
			p(i,j) = pl(i,j,k)**(1.-alfa)*pl(i,j,k+1)**alfa
			t(i,j) = theta0*(p(i,j)/p0)**gamma

            ! Calcul de dtheta/dp
			if (alfa <= .5) then
				if (k == 1)then

	    			dltdlp = log(thetal(i,j,k+1)/thetal(i,j,k))/ &
							 log(pl(i,j,k+1)/pl(i,j,k))

				else
					beta = 1.-2.*(alfa-.5)**2
					dltdlp = beta*log(thetal(i,j,k+1)/thetal(i,j,k))/ &
							 log(pl(i,j,k+1)/pl(i,j,k)) &
						   + (1.-beta)*alog(thetal(i,j,k)/thetal(i,j,k-1))/ &
							 log(pl(i,j,k)/pl(i,j,k-1))
				endif
			else
				if (k == np0-1) then
					dltdlp = log(thetal(i,j,k+1)/thetal(i,j,k))/ &
							 log(pl(i,j,k+1)/pl(i,j,k))
				else
					beta = 1.-2.*(alfa-.5)**2
					dltdlp = beta*log(thetal(i,j,k+1)/thetal(i,j,k))/ &
							 log(pl(i,j,k+1)/pl(i,j,k)) &
						   + (1.-beta)*log(thetal(i,j,k+2)/thetal(i,j,k+1))/ &
							 log(pl(i,j,k+2)/pl(i,j,k+1))
				endif
			endif

			dtdp(i,j) = -dltdlp*theta0/(100.*p(i,j))
		enddo
	enddo

   ! Calcul de la PV
	do j = 1,ny
		phi=alatmin+(j-1)*dphi
		if (abs(abs(phi)-90.) > 1.e-3) then

			do i = 1,nx
				im1 = i-1
				if (im1 == 0) then
					im1 = nx
				endif
				ip1=i+1
				if (ip1 == nx+1) then
					ip1 = 1
				endif
				dvdx = (v(ip1,j)-v(im1,j))/(2.*dx*cos(phi*pio180))
				if (j == ny) then
					dudy = (u(i,j)-u(i,j-1))/dy
				else if (j == 1) then
					dudy = (u(i,j+1)-u(i,j))/dy
				else
					dudy = (u(i,j+1)-u(i,j-1))/(2.*dy)
				endif

				utphi = u(i,j)*tan(phi*pio180)/rt
				f = 4.*pi/86164.*sin(phi*pio180)
				pv(i,j) = g*(dvdx-dudy+utphi+f)*dtdp(i,j)*1.e6
			enddo
		endif
	enddo

    ! Si on est a un pole la PV ne depend pas de la longitude
	do j = 1,ny,ny-1
		phi = alatmin+(j-1)*dphi
		if (abs(abs(phi)-90.) <= 1.e-3) then
			pvmoy = 0.

			if (j == 1) then
				j1 = 2
			endif
			if (j == ny) then
				j1 = ny-1
			endif

			do i = 1,nx
				pvmoy = pvmoy+pv(i,j1)
			enddo

			pvmoy = pvmoy/real(nx)
			do i = 1,nx
				pv(i,j) = pvmoy
			enddo
		endif
	enddo

	return

end subroutine pvp

