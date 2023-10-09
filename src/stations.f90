! Module:
!    stations
!
! Purpose:
!    This module allow to saved all the caracteristics of stations
!
!   ABE Aberystwyth          Wales        52.484   -4.067  138
!   AIR Aire sur Adour       France       43.700   -0.250   75
!   AND Andoya               Norway       69.280   16.180
!   BAU Bauru                Brazil      -22.347  -49.027  640
!   BOR Bordeaux             France       44.500   -0.534
!   DJO Djougou              Benin         9.647    1.741  439
!   DOM Dome C               Antartica   -75.100  123.350 3250
!   DDU Dumont d'Urville     France      -66.666  140.017   45
!   EUR Eureka               Canada       80.053  -86.416  623
!   GAP Gap                  France       44.450    6.030  606
!   HAR Harestua             Norway       60.200   10.750
!   JUN Jungfraujoch         Switzerland  46.548    7.984 3600
!   KER Kerguelen            France      -49.352   70.256   36
!   KIR Kiruna               Sweden       67.881   21.067  330
!   LEO Leon                 Spain        42.580   -5.650  900
!   MAR Marambio             Antartica   -64.200  -56.700
!   NYA Ny Alesund           Spitzberg    78.908   11.883
!   OHP Obsv. Haute Provence France       43.935    5.712  683
!   PAR Paris                France       48.846    2.346   63
!   REU Reunion              France      -20.901   55.484  110
!   ROT Rothera              Antartica   -67.570  -68.13
!   SAL Salekhard            Russia       66.667   66.533
!   SAN Sanae                Antartica   -71.670   -2.840  850
!   SCO Scoresby Sund        Greenland    70.485  -21.952   67
!   SOD Sodankyla            Finland      67.368   26.633  170
!   SOU South Pole           Antartica   -90.000    0.000
!   TAR Tarawa               Kiribati      1.355  172.923   10
!   THU Thule                Greenland    76.540  -68.780
!   TOR Toronto              Canada       43.660  -79.390  170
!   VER Verrieres le Buisson France       48.755    2.244  202
!   YAK Yakutsk              Russia       62.030  129.020
!   ZIG Zhigansk             Russia       66.793  123.351  200
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
module stations

    implicit none

	integer, parameter                 :: nbsta = 32        ! number of stations
	character(len=3), dimension(nbsta) :: nomsta
	real, dimension(nbsta)             :: alatsta, alongsta ! location of each station

	character(len=19)	               :: filesta
	integer							   :: ista				! compteur boucle stations
	integer							   :: ifilesta			! numero unite fichier station
	integer							   :: njoursta
	integer							   :: ixsta, jysta
	real							   :: xsta, ysta

	data nomsta/ 'ABE', 'AIR', 'AND', 'BAU', 'BOR' &
			   , 'DJO', 'DOM', 'DDU', 'EUR', 'GAP' &
			   , 'HAR', 'JUN', 'KER', 'KIR', 'LEO' &
			   , 'MAR', 'NYA', 'OHP', 'PAR', 'REU' &
			   , 'ROT', 'SAL', 'SAN', 'SCO', 'SOD' &
			   , 'SOU', 'TAR', 'THU', 'TOR', 'VER' &
			   , 'YAK', 'ZIG'/

	data alatsta/  52.484,  43.700,  69.280, -22.347,  44.500 &
				,   9.647, -75.100, -66.666,  80.053,  44.450 &
				,  60.200,  46.548, -49.352,  67.881,  42.580 &
				, -64.200,   8.908,  43.935,  48.846, -20.901 &
				, -67.570,  66.667, -71.670,  70.485,  67.368 &
				, -90.000,   1.355,  76.540,  43.660,  48.755 &
				,  62.030,  66.793/

	data alongsta/   -4.067,   -0.250,   16.180,  -49.027,   -0.534 &
				 ,    1.741,  123.350,  140.017,  -86.416,    6.030 &
				 ,   10.750,    7.984,   70.256,   21.067,   -5.650 &
				 ,  -56.700,   11.883,    5.712,    2.346,   55.484 &
				 ,  -68.130,   66.533,   -2.840,  -21.952,   26.633 &
				 ,    0.000,  172.923,  -68.780,  -79.390,    2.244 &
				 ,  129.020,  123.351/

end module stations
