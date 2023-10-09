! Program:
!    filament
!
! Purpose:
!    This is the main program of MIMOSA.
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
!    201003   M.-A. DROUIN/LPC2E conversion into FORTRAN 95
!    201104   M.-A. DROUIN/LPC2E add comments
!
!======================================================================
program filament

    use constantes
    use stations
    use interfaces_mod

    implicit none

    !==================================================================
    ! Declarations & Initialisation
    !==================================================================

    ! Config of the simulation
    !------------------------------------------------------------------

    ! Type of run
    integer:: zone
    integer:: ihemd
    integer:: ihemf
    integer:: runtype
    integer:: intype

    ! Type of output Files
    integer:: wind_out
    integer:: t_out
    integer:: stations_out

    ! ECMWF grid
    integer:: lathem
    integer:: latraccord
    integer:: latminmim
    integer:: latmaxmim
    integer:: nbdegmim
    integer:: nx
    integer:: ny
    integer:: np
    real:: latminecmr
    real:: latmaxecmr
    real:: paslat
    real:: paslong

    ! MIMOSA grid
    integer:: ndeg
    integer:: ndeg2
    integer:: ng

    ! ECMWF read fields
    real, dimension(1:50) :: pres = 0.
    real, dimension(:,:,:), allocatable :: pl
    real, dimension(:,:,:), allocatable :: thetal
    real, dimension(:,:,:), allocatable :: tl
    real, dimension(:,:,:), allocatable :: ul
    real, dimension(:,:,:), allocatable :: vl

    ! ECMWF fields on chosen isentropic level
    real, dimension(:,:), allocatable :: p
    real, dimension(:,:), allocatable :: pv
    real, dimension(:,:), allocatable :: t
    real, dimension(:,:), allocatable :: u
    real, dimension(:,:), allocatable :: v
    real, dimension(:,:), allocatable :: dtdp

    ! Array for pvp routine (allocatable)
    integer:: alloc_stat   ! var to know if allocation is a sucess
    integer(kind=1), dimension(:,:)  , allocatable :: gind0  ! gind0 = 1 grid point in ECMWF field
    integer(kind=1), dimension(:,:,:), allocatable :: gind   ! gind  = 1 advected grid point in ECMWF field
    real, dimension(:,:,:), allocatable :: gu                ! ECMWF wind component on grid (step grid/hour)
    real, dimension(:,:,:), allocatable :: gv                ! ECMWF wind component on grid (step grid/hour)
    real, dimension(:,:,:), allocatable :: gpv0              ! ECMWF PV on grid points
    real, dimension(:,:,:), allocatable :: gpv               ! model PV on grid points
    real, dimension(:,:)  , allocatable :: gpvn
    real, dimension(:,:)  , allocatable :: dgpv
    real, dimension(:,:)  , allocatable :: dgpvlis           ! model PV smoothed on grid points
    real, dimension(:,:,:), allocatable :: gpvll
    integer(kind=2), dimension(:,:)  , allocatable :: gpoids
    real, dimension(:,:,:), allocatable :: gx                ! coordinate of advected grid points
    real, dimension(:,:,:), allocatable :: gy                ! coordinate of advected grid points
    real, dimension(:,:,:), allocatable :: glong             ! lat and long of grid points
    real, dimension(:,:,:), allocatable :: glat              ! lat and long of grid points
    real, dimension(:,:)  , allocatable :: gxhem
    real, dimension(:,:)  , allocatable :: gyhem
    real, dimension(:,:,:), allocatable :: gt0
    real, dimension(:,:,:), allocatable :: gtll	             ! T ECMWF aux points de grille

    !	Config
    integer:: nlis2D
    integer:: nhgrid 	! nbre heure remaillage
    integer:: nwrite 	! nbre heure fichier pvf
    integer:: nhmod    	! nbre heures entre 2 ch. ECMWF
    integer:: nhrelax	! nbre heures relaxation
    integer:: nwtemp	! nbre heures fichier temperature
    integer:: nprhmod	! heure premier fichier ECMWF
    integer:: nprwrite 	! heure premier fichier pvf
    integer:: nrun     	! numero de dossier de sortie
    integer:: indifexpl ! 0 pas de diffusion, 1 diffusion
    real:: diff			! valeur de la diffusion

    !	parametre de simulation
    integer:: iand		! date initialisation simulation
    integer:: moisd		! date initialisation simulation
    integer:: jourd		! date initialisation simulation
    integer:: iheured	! date initialisation simulation
    integer:: ianf		! date fin simulation
    integer:: moisf		! date fin simulation
    integer:: jourf		! date fin simulation
    integer:: iheuref	! date fin simulation
    integer:: initpv	! 1 initialisation, 0 reprise de fichier ph*
    real:: teta			! niveau isentrope

    !
    !	Entete fichier de sortie
    !
    !	  iparam     :parametres ecrits dans fichier resultat (majuscules)
    !	  0 1 2 3    : Annee Mois Jour Heure
    !	  4 5 6 7    : Annee Mois Jour Heure initialisation
    !	  8          : niveau Teta
    !	  9 10 11 12 : coins W S E N ecmwf
    !	  13 14      : nombre points Longit Latit ecmwf
    !	  15 16      : Pas en heure Premiere heure ecmwf
    !     17 18	     : Points par degre calcul, N echantillonage ecriture
    !	  19 20      : nombre points grille X Y
    !	  21 22      : Pas par heure, Nombre heures avant  regrillage
    !	  23         : Relaxation (heures)
    !	  24         : Simulation filaments=1  analyse ecmwf=2  temperature=3
    !	  25	     : Echeance
    !	  26-29      : Reserve
    integer, dimension(0:29) :: iparam
    integer, dimension(0:29) :: iparam1

    !	chemins et nom de fichiers
    character(len=20):: NomSAap
    character(len= 7):: Sisobaric
    character(len=22):: nompvf
    character(len=21):: nomu
    character(len=21):: nomv
    character(len=21):: nomt
    character(len=23):: nompvh
	
    character(len= 2), dimension(3) :: nomhem
    data nomhem/'s', 'n', 'g'/

    ! Local Variables
    !------------------------------------------------------------------
    integer:: nfull
    integer:: i
    integer:: ig
    integer:: ixx
    integer:: j
    integer:: jg
    integer:: jyy
    integer:: ihem
    integer:: ihemwr
    integer:: njourd
    integer:: njourf
    integer:: njourc
    integer:: njourmod
    integer:: njour0
    integer:: njouran
    integer:: icent
    integer:: ian
    integer:: ianmod
    integer:: jour
    integer:: jourmod
    integer:: mois
    integer:: moismod
    integer:: ntimed
    integer:: ntimef
    integer:: it
    integer:: itmod
    integer:: itmod_old
    integer:: iheure
    integer:: iheuremod
    integer:: iheurenom
    integer:: indrelax
    integer:: nhrelax1

    real, dimension(2) 	:: hemtab
    data hemtab/-1.,1./
    real:: hem
    real:: alfadiff
    real:: latmin
    real:: latmax
    real:: alfax
    real:: alfay
    real:: alongx
    real:: alaty
    real:: xx
    real:: yy
    real:: guint
    real:: gvint
    real:: coef
    real:: rayon
    real:: rapcos
    real:: sum
    real:: sumpv2
    real:: sumt
    real:: sumt2

    ! ERROR
    logical:: error

    ! Declaration of the namelists
    !------------------------------------------------------------------
    logical:: nml_exists
    namelist /grid/   nx, ny, np, paslat, paslong, pres, latminecmr, latmaxecmr
    namelist /config/ nlis2D, ndeg, nhgrid, nwrite, nhmod, nhrelax,nprhmod, nprwrite, indifexpl, diff
    namelist /run/	  iand, moisd, jourd, iheured, ianf, moisf, jourf, iheuref, initpv, teta, zone, runtype, intype
    namelist /output/ nrun, nwtemp, wind_out, t_out, stations_out

    ! verifier si toujour utile
    integer:: iech
    integer:: indat
    integer:: nh = 1

    !======================================================================
    ! Program
    !======================================================================

    pi	   = 4.*atan(1.)
    pio180 = 4.*atan(1.)/180.
    radeg  = 45./atan(1.)
    Iech   = 0
    indat  = 0

    ! Reading of the parameters in input.namelist
    inquire(file = "input.namelist", exist=nml_exists)
    if (.not. nml_exists) then
        print *, "ERROR : The file input.namelist doesn't exist"
        error = .true. ; stop
    end if
    open(unit = 20, file = "input.namelist")
    read(unit = 20, nml = run)
    read(unit = 20, nml = grid)
    read(unit = 20, nml = config)
    read(unit = 20, nml = output)
    close(unit = 20)

    ! Definition en fonction de la zone de calcul
    select case (zone)
        case(1)
            ihemd = 2
            ihemf = 2
            nfull = 1
        case(2)
            ihemd = 1
            ihemf = 1
            nfull = 1
        case(3)
            ihemd = 1
            ihemf = 2
            nfull = 0
        case default
            print *, "ERROR ", zone, " is not a value for zone"
            error = .true. ; stop
    end select

    ! Initialization of variables based on value in input.namelist
    lathem     = 100
    latraccord = 5
    latminmim  = -90+(180-lathem)*(ihemd-1)
    latmaxmim  =  90-(180-lathem)*(2-ihemf)
    nbdegmim   = latmaxmim - latminmim
    ndeg2      = ndeg
    ng         = ndeg*lathem
    alfadiff   = diff*ndeg**2*3600.*nhgrid/111111.**2

    ! Defining allocatable array
    !------------------------------------------------------------------
    allocate(gind0(-ng:ng,-ng:ng), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array gind0 "
        error = .true. ; stop
    endif
    allocate(gind(-ng:ng,-ng:ng, ihemd:ihemf), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array gind "
        error = .true. ; stop
    endif
    allocate(gu(-ng:ng,-ng:ng, ihemd:ihemf), gv(-ng:ng,-ng:ng, ihemd:ihemf), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array gu or gv "
        error = .true. ; stop
    endif
    allocate(gpv(-ng:ng,-ng:ng, ihemd:ihemf), gpv0(-ng:ng,-ng:ng, ihemd:ihemf), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array gpv or gpv0 "
        error = .true. ; stop
    endif
    allocate(gpvn(-ng:ng,-ng:ng), dgpv(-ng:ng,-ng:ng), dgpvlis(-ng:ng,-ng:ng), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array gpvn, dgpv or dgpvlis "
        error = .true. ; stop
    endif
    allocate(gpvll(-180*ndeg2:180*ndeg2,latminmim*ndeg2:latmaxmim*ndeg2,ihemd:ihemf), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array gpvll "
        error = .true. ; stop
    endif
    allocate(gpoids(-ng:ng,-ng:ng), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array gpoids "
        error = .true. ; stop
    endif
    allocate(gx(-ng:ng,-ng:ng,ihemd:ihemf), gy(-ng:ng,-ng:ng,ihemd:ihemf), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array gx or gy "
        error = .true. ; stop
    endif
    allocate(glong(-ng:ng,-ng:ng,ihemd:ihemf), glat(-ng:ng,-ng:ng,ihemd:ihemf), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array glong or glat "
        error = .true. ; stop
    endif
    allocate(gxhem(-ng:ng,-ng:ng), gyhem(-ng:ng,-ng:ng), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array gxhem or gyhem "
        error = .true. ; stop
    endif
    allocate(pl(nx,ny,np), ul(nx,ny,np), vl(nx,ny,np), tl(nx,ny,np), thetal(nx,ny,np), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array pl, ul, vl, tl or thetal "
        error = .true. ; stop
    endif
    allocate(u(nx,ny), v(nx,ny), pv(nx,ny), p(nx,ny), t(nx,ny), dtdp(nx,ny))
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array u, v, pv, p, t or detdp "
        error = .true. ; stop
    endif
    !	test sortie vents et temp
    allocate(gt0(-ng:ng,-ng:ng,ihemd:ihemf), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array gt0 "
        error = .true. ; stop
    endif
    allocate(gtll(-180*ndeg2:180*ndeg2,latminmim*ndeg2:latmaxmim*ndeg2,ihemd:ihemf), stat = alloc_stat)
    if (alloc_stat /= 0) then
        print *, "ERROR allocating array gtll "
        error = .true. ; stop
    endif

    ! Config based to ECMWF file
    select case (intype)
        case(1)
            do i=1,nx
                do j=1,ny
                    pl(i,j,:) = pres(1:np)
                enddo
            enddo
            print *, "Fichier ECMR"
        case(2)
            print *, "Fichier GRIB"
        case default
            print *, "ERROR ", intype, " is not a valid value for intype"
            error = .true. ; stop
    end select

    ! Header of outfiles
    if (initpv == 1) then
        iparam(4) = iand
        iparam(5) = moisd
        iparam(6) = jourd
        iparam(7) = iheured
    endif
    iparam(9)  = -180
    iparam(10) = latminmim
    iparam(11) = 180
    iparam(12) = latmaxmim
    iparam(13) = nx
    iparam(14) = ny
    iparam(15) = 24
    iparam(16) = 12
    iparam(17) = ndeg
    iparam(21) = nh
    iparam(22) = nhgrid
    iparam(23) = nhrelax

    iparam(8) = nint(teta)

    ! If stations output activated
    ! we open units for the file and write headers
    if ( stations_out == 1 .and. zone == 3 ) then
        do ista	= 1,nbsta
            ifilesta = ista+50
            write(unit=filesta,FMT='("RUN",I2.2,"/pv",i4.4,"_",i2.2,".",a3)')nrun,iparam(8),iand,nomsta(ista)
            open(unit=ifilesta,file=filesta,status='unknown')
            write(unit=ifilesta,FMT='(1x,a3,2f8.1,f6.0)')nomsta(ista),alatsta(ista),alongsta(ista),teta
        enddo
    endif

    ! Start initialization of the grid
    do ihem = ihemd,ihemf
        hem = hemtab(ihem)
        call initgrid(gind(-ng,-ng,ihem),gind0, &
        gx(-ng,-ng,ihem),gy(-ng,-ng,ihem),glat(-ng,-ng,ihem), &
        glong(-ng,-ng,ihem),ng,hem,ndeg)
    enddo

    ! Connection between hemispheres
    do jg = -ng,ng
        do ig = -ng,ng
            if (gind0(ig,jg) == 1) then
                rayon=sqrt(real(ig*ig+jg*jg))
                if (rayon > real((180-lathem)*ndeg-1)) then
                    gxhem(ig,jg)=real(ig)*(180.*real(ndeg)-rayon)/rayon
                    gyhem(ig,jg)=-real(jg)*(180.*real(ndeg)-rayon)/rayon
                endif
            endif
        enddo
    enddo
    ! End initialization of the grid

    ! Loop over time
    !------------------------------------------------------------------
    njourd = njour(jourd,moisd,iand)
    njourf = njour(jourf,moisf,ianf)
    ntimed = njourd*24+iheured+1
    ntimef = njourf*24+iheuref
    itmod_old = 0

    do it = ntimed,ntimef
        itmod = (it-nprhmod+(nhmod/2)-1)/nhmod
        itmod = itmod*nhmod+nprhmod
        njourc = it/24
        njourmod = itmod/24
        call njouri(jour,mois,ian,njourc)
        call njouri(jourmod,moismod,ianmod,njourmod)
        iheure = it-njourc*24
        iheuremod = itmod-njourmod*24
        print '("jour ",I2.2,1X,I2.2,1X,I2.2,1X,I2.2," heures | mod  ",I2.2,1X,I2.2,1X,I2.2,1X,I2.2," heures")', &
        jour,mois,ian,iheure,jourmod,moismod,ianmod,iheuremod
        print *, ""

        if (itmod /= itmod_old) then
            itmod_old = itmod

            ! Reading isobaric files
            if (intype == 1) then

                if (ian >= 51) then
                    Sisobaric="ECMR/19"
                else
                    Sisobaric="ECMR/20"
                endif

                write(UNIT=nomSAap, FMT='(I2.2,"/",I2.2,"/D",4I2.2,".ECMR")') ianmod,moismod,ianmod,moismod,jourmod,iheuremod
                open(UNIT=21, FILE=Sisobaric//nomSAap, status="old")
                print '(1x,a22,a20)', Sisobaric, nomsaap

                call readecmr(ianmod,moismod,jourmod,iech,ul,vl,tl,nx,ny,np)
                close(21)

            else

                if (ian >= 51) then
                    Sisobaric="GRIB/19"
                else
                    Sisobaric="GRIB/20"
                endif

                write(UNIT=nomSAap, FMT='(I2.2,"/",I2.2,"/D",4I2.2,".grib")') ianmod,moismod,ianmod,moismod,jourmod,iheuremod
                print '(1x,a22,a20)', Sisobaric, nomsaap
                call readgrib(Sisobaric,nomSAap,pl,tl,ul,vl,nx,ny,np)

            endif

            call pvp(teta,nx,ny,np,np,pl,ul,vl,tl,thetal,u,v,t,p,dtdp,pv,-90.,90.)

            ! interpolation grid
            print '(1X,A20)', "interpolation grille"
            print *, ""

            do ihem = ihemd,ihemf
                hem = hemtab(ihem)
                do ig = -ng,ng
                    do jg = -ng,ng
                        if (gind0(ig,jg) /= 1) then
                            if (it == ntimed) then
                                gpv0(ig,jg,ihem) = -99.9
                            endif
                        else
                            guint = ainterp2D(u,glong(-ng,-ng,ihem),glat(-ng,-ng,ihem) &
                            ,nx,ny,ig,jg,ng,paslong,latminecmr,paslat)
                            gvint = ainterp2D(v,glong(-ng,-ng,ihem),glat(-ng,-ng,ihem) &
                            ,nx,ny,ig,jg,ng,paslong,latminecmr,paslat)
                            if (ig /= 0 .or. jg /= 0) then
                                rapcos=(90.-hem*glat(ig,jg,ihem))/cos(glat(ig,jg,ihem)*pio180)
                            else
                                rapcos=radeg
                            endif

                            gu(ig,jg,ihem)   = ndeg/rt*(-hem*radeg*gvint &
                            * sin(glong(ig,jg,ihem)*pio180) &
                            + rapcos*guint*cos(glong(ig,jg,ihem)*pio180))*3600.
                            gv(ig,jg,ihem)   = ndeg/rt*(radeg*gvint*cos(glong(ig,jg,ihem)*pio180) &
                            + hem*rapcos*guint*sin(glong(ig,jg,ihem)*pio180))*3600.
                            gpv0(ig,jg,ihem) = ainterp2D(pv,glong(-ng,-ng,ihem) &
                            , glat(-ng,-ng,ihem),nx,ny,ig,jg,ng,paslong,latminecmr,paslat)
                            gt0(ig,jg,ihem)  = ainterp2D(t,glong(-ng,-ng,ihem) &
                            , glat(-ng,-ng,ihem),nx,ny,ig,jg,ng,paslong,latminecmr,paslat)

                            if (it == ntimed .and. initpv == 1) then
                                gpv(ig,jg,ihem) = gpv0(ig,jg,ihem)
                            endif
                        endif
                    enddo
                enddo
            enddo

            !	Restart PV file
            if (it == ntimed .and. initpv /= 1) then
                do ihem = ihemd,ihemf

                    write(UNIT=nompvh, FMT='("RUN",I2.2,"/ph",A1,3I2.2,I3.3,".",I4.4)') nrun,nomhem(ihem) &
                                                                                      ,ian,mois,jour, iheure-1,int(teta)
                    print '(1x,A20)', nompvh

                    open(22,file=nompvh,form='unformatted',status='unknown')
                    read(22)iparam1,((gpv(i,j,ihem),i=-ng,ng),j=-ng,ng)
                    print *, "iparam1 ", iparam1

                    iparam(4)=iparam1(4)
                    iparam(5)=iparam1(5)
                    iparam(6)=iparam1(6)
                    iparam(7)=iparam1(7)
                    close(22)
                enddo
            endif

            print *,'      Fin lecture fichier'

        end if

        ! Little loop over time
        do ihem = ihemd,ihemf
            hem = hemtab(ihem)
            call deplacement(gind0,gind(-ng,-ng,ihem) &
            ,gu(-ng,-ng,ihem),gv(-ng,-ng,ihem),gx(-ng,-ng,ihem) &
            ,gy(-ng,-ng,ihem),ng,nh)

            ! indrelax=0      remaillage sans relaxation
            ! indrelax=1      remaillage avec relaxation
            ! indrelax=2      stockage fichier pvf reduit
            if (mod(it,nhgrid) == 0) then
                if (mod(it,nhmod) == nprhmod) then
                    indrelax = 1
                else if (mod(it,nhgrid) == 0) then
                    indrelax = 0
                else
                    indrelax = 2
                endif

                print '("appel a regrid ",I2.2,1X,I2.2,1X,I2.2,1X,I2.2," indrelax = ",I2.2)',ian,mois,jour,iheure,indrelax
                print *, ""

                nhrelax1 = nhrelax
                if (it-ntimed < nhgrid .and. initpv == 1) then
                    nhrelax1 = 1
                endif
                call regrid2( gpv(-ng,-ng,ihem),gpv0(-ng,-ng,ihem) &
                , gpvn,gpoids,gx(-ng,-ng,ihem),gy(-ng,-ng,ihem) &
                , gind(-ng,-ng,ihem),gind0,ng &
                , nhrelax,nhmod,indrelax,dgpv,dgpvlis,nlis2D )
                if(indifexpl == 1) then
                    print *,'appel diffusion'
                    call diffusion( gpv(-ng,-ng,ihem),gpvn,gind(-ng,-ng,ihem) &
                    ,	gind0, gx(-ng,-ng,ihem),gy(-ng,-ng,ihem),ng,alfadiff)
                endif
            endif
        enddo

        if (nfull == 0 .and. mod(it,nhgrid) == 0) then
            call raccord(gpv,gind,gind0,gxhem,gyhem,ng,lathem,ndeg,nlis2D)
        endif
        !
        if (mod(it,nwrite) == nprwrite)then

            !	Interpolation sur une grille latitude -longitude
            do ihem = ihemd,ihemf
                hem = hemtab(ihem)
                latmin = -90+(180-lathem)*(ihem-1)
                latmax = latmin+lathem

                do j = latmin*ndeg2,latmax*ndeg2
                    alaty = real(j)/real(ndeg2)

                    do i = -180*ndeg2,180*ndeg2

                        alongx = real(i)/real(ndeg2)
                        xx = ndeg*sin(alongx*pio180)*(90.-hem*alaty)*.999999
                        yy = -hem*ndeg*cos(alongx*pio180)*(90.-hem*alaty)*.999999
                        ixx = nint(xx-.5)
                        jyy = nint(yy-.5)
                        alfax = xx-ixx
                        alfay = yy-jyy
                        gpvll(i,j,ihem) = gpv(ixx,jyy,ihem)*(1.-alfax)*(1.-alfay) &
                        + gpv(ixx+1,jyy,ihem)*alfax*(1.-alfay) &
                        + gpv(ixx,jyy+1,ihem)*(1.-alfax)*alfay &
                        + gpv(ixx+1,jyy+1,ihem)*alfax*alfay

                        gtll(i,j,ihem) = gt0(ixx,jyy,ihem)*(1.-alfax)*(1.-alfay) &
                        + gt0(ixx+1,jyy,ihem)*alfax*(1.-alfay) &
                        + gt0(ixx,jyy+1,ihem)*(1.-alfax)*alfay &
                        + gt0(ixx+1,jyy+1,ihem)*alfax*alfay

                    enddo
                enddo
            enddo

            ! Connecting the 2 hemispheres
            ! All data in ihem = 1

            ! Connecting PV
            if	(ihemd /= ihemf) then
                Sumpv2 = 0.
                Sum = 0.
                if (t_out == 1) then
                    sumt2 = 0.
                    sumt  = 0.
                end if
                do j = -latraccord*ndeg2+1,latraccord*ndeg2-1
                    coef = .5*(1.+sin(real(j)/real(ndeg2*latraccord)*pi*.5))

                    do i = -180*ndeg2,180*ndeg2

                        sumpv2 = sumpv2+(gpvll(i,j,1)-gpvll(i,j,2))**2
                        sum = sum+1.
                        gpvll(i,j,1) = (1.-coef)*gpvll(i,j,1)+coef*gpvll(i,j,2)

                        if (t_out == 1) then
                            sumt2 = sumt2+(gtll(i,j,1)-gtll(i,j,2))**2
                            sumt  = sumt+1.
                            gtll(i,j,1) = (1.-coef)*gtll(i,j,1)+coef*gtll(i,j,2)
                        end if

                    enddo
                enddo

                print '("ecart hemisphres pv = ",E12.7)', sqrt(sumpv2/sum)
                print *, ""

                if (t_out == 1) then
                    print '("ecart hemisphres t  = ",E12.7)', sqrt(sumt2/sumt)
                    print *, ""
                end if

                do j = latraccord*ndeg2,latmaxmim*ndeg2
                    do i = -180*ndeg2,180*ndeg2

                        gpvll(i,j,1) = gpvll(i,j,2)
                        if (t_out == 1) then
                            gtll(i,j,1) = gtll(i,j,2)
                        end if

                    enddo
                enddo
            endif

            ! Writing sampled file
            iheurenom = iheure

            !	Filename of North/south Hemisphere
            if (ihemd /= ihemf) then
                ihemwr = 3
            else
                ihemwr = ihemd
            endif
            write(UNIT=nompvf, FMT='("RUN",i2.2,"/pv",a1,4i2.2,".",i4.4)') nrun,nomhem(ihemwr),ian,mois,jour,iheure,int(teta)

            open(23,file=nompvf,form='unformatted',status='unknown')
            iparam(0) = ian
            iparam(1) = mois
            iparam(2) = jour
            iparam(3) = iheure
            iparam(24) = 1
            iparam(18) = ndeg2
            iparam(19) = 360*ndeg2
            iparam(20) = nbdegmim*ndeg2+1

            write(23) iparam, ((gpvll(i,j,ihemd),i=0,180*ndeg2) &
            ,(gpvll(i,j,ihemd) &
            ,i=-180*ndeg2+1,-1),j=latmaxmim*ndeg2,latminmim*ndeg2,-1)
            close(23)
            if (t_out == 1 .and. mod(it,nwtemp) == 0) then

                write(UNIT=nomt, FMT='("RUN",i2.2,"/t",a1,4i2.2,".",i4.4)') nrun,nomhem(ihemwr),ian,mois,jour,iheure,int(teta)
                open(26,file=nomt,form='unformatted',status='unknown')
                write(26) iparam, ((gtll(i,j,ihemd),i=0,180*ndeg2) &
                ,(gtll(i,j,ihemd) &
                ,i=-180*ndeg2+1,-1),j=latmaxmim*ndeg2,latminmim*ndeg2,-1)
                close(26)

            end if
            if (wind_out == 1 .and. mod(it,nwtemp) == 0) then

                iparam(19) = nx
                iparam(20) = ny
                write(UNIT=nomu, FMT='("RUN",i2.2,"/u",a1,4i2.2,".",i4.4)') nrun,nomhem(ihemwr),ian,mois,jour,iheure,int(teta)
                open(24,file=nomu,form='unformatted',status='unknown')
                write(24) iparam, ((u(i,j),i=1,nx),j=ny,1,-1)
                close(24)
                write(UNIT=nomv, FMT='("RUN",i2.2,"/v",a1,4i2.2,".",i4.4)') nrun,nomhem(ihemwr),ian,mois,jour,iheure,int(teta)
                open(25,file=nomv,form='unformatted',status='unknown')
                write(25) iparam, ((v(i,j),i=1,nx),j=ny,1,-1)
                close(25)

            end if

            ! Writing full resolution files for restart
            if (it == ntimef .or. (jour == 1 .and. iheure == 0)) then
                do ihem = ihemd,ihemf
                    hem = hemtab(ihem)
                    write(UNIT=nompvh, FMT = '("RUN",i2.2,"/ph",a1,3i2.2,i3.3,".",i4.4)') nrun,nomhem(ihem) &
                                                                                              ,ian,mois,jour,iheure, int(teta)
                    open(27,file=nompvh,form='unformatted',status='unknown')

                    iparam(0)=ian
                    iparam(1)=mois
                    iparam(2)=jour
                    iparam(3)=iheure
                    iparam(18)=ndeg
                    iparam(19)=2*ng+1
                    iparam(20)=2*ng+1
                    iparam(24)=1
                    write(27)iparam,((gpv(i,j,ihem),i=-ng,ng),j=-ng,ng)
                    close(27)
                enddo
            endif

            ! Samples at stations
            if (mod(it,24) == 12 .and. zone == 3 .and. stations_out == 1) then
                do ista = 1,nbsta
                    ifilesta = ista+50
                    if (alatsta(ista) > 0.) then
                        xsta =  ndeg*sin(alongsta(ista)*pio180)*(90.-alatsta(ista))
                        ysta = -ndeg*cos(alongsta(ista)*pio180)*(90.-alatsta(ista))
                    else
                        xsta =  ndeg*sin(alongsta(ista)*pio180)*(90.+alatsta(ista))
                        ysta =  ndeg*cos(alongsta(ista)*pio180)*(90.+alatsta(ista))
                    endif
                    ixsta = nint(xsta)
                    jysta = nint(ysta)
                    icent = 1900
                    if (ian < 51) then
                        icent = 2000
                    endif
                    njoursta = njour(jour,mois,ian)
                    njour0   = njour(1,1,ian)-1
                    njouran  = njoursta-njour0
                    if (alatsta(ista) > 0.) then
                        write(ifilesta,'(1x,2i4,3i3,f8.2,f8.3)')ian+icent,njouran,mois,jour,iheure &
                        ,gpv(ixsta,jysta,2),gtll(ixsta,jysta,2)
                    else
                        write(ifilesta,'(1x,2i4,3i3,f8.2,f8.3)')ian+icent,njouran,mois,jour,iheure &
                        ,gpv(ixsta,jysta,1),gtll(ixsta,jysta,1)
                    endif
                enddo
            endif

        endif
    enddo

    ! If stations outputs was chosen; we close the files
    if ( stations_out == 1 .and. zone == 3 ) then
        do ista	= 1,nbsta
            ifilesta = ista+50
            close(ifilesta)
        end do
    end if

    ! Deallocation of array
    deallocate(gind0, gind, gu, gv, gpv, gpv0, gpvn, dgpv, dgpvlis, gpvll)
    deallocate(gpoids, gx, gy, glong, glat, gxhem, gyhem)
    deallocate(pl, ul, vl, tl, thetal)
    deallocate(u, v, pv, p, t, dtdp)
    deallocate(gt0, gtll)

end program filament
