#!/bin/bash

# ----------------------------------------------------------------------------------------------
# Description:
#   This script allows to do MIMOSA simulations for multiple isentropic levels 
# between 2 dates. This script will also create the list of MIMOSA outputed files.
#
# Categories:
#   MIMOSA tools
#
# Params:
#   none
# 
# Example:
#   To use this script you need to :
#      - copy this script in the directory containing MIMOSA executable (mimosa.x)
#      - make this script executable (chmod 755 mimosa_namelist.sh
#      - compile MIMOSA using the MakeFile
#      - configure the simulations in the parameters part below.
#      - execute the script using ./this_script.sh
# 
# History:
#   200903 MAD/LPC2E Creation
#   201109 MAD/LPC2E Group all parameters
#                    Change creation of input.namelist due to change in MIMOSA model
#                    Add creation and checking of directory needed for the simulation
#   2011 -> 2023 CB  Maintenance
#
# Author:
#   M.-A. Drouin & Cathy Boonne
# ----------------------------------------------------------------------------------------------

function info_msg(){
	txt=$1
	echo "$(date +'%d/%m/%Y %H:%M:%S')   [INFO]   ${txt}"
}
function err_msg(){
	txt=$1
	echo "$(date +'%d/%m/%Y %H:%M:%S')   [ERROR]   ${txt}"
}
function warn_msg(){
	txt=$1
	echo "$(date +'%d/%m/%Y %H:%M:%S')   [WARNING]   ${txt}"
}

############################################
# PARAMETERS TO DEFINE BY THE USER         #
############################################

SIMUDIR="/home/damali/Work/SEDOO/MIMOSA_wdir/MIMOSA_my_ecmwf"
RUN_ID_NUMBER="1" # number of the run for the directory (from 1 to anything)

# Start date
SYEAR=23
SMONTH=10
SDAY=03
SHOUR=12

# End date
EYEAR=23
EMONTH=10
EDAY=04
EHOUR=12

## RUN
ZONE=3
INTYPE=2
INITPV=1
THETA=( 380 475 550 675 )

## GRID
NX=288
NY=145
# NP=91
NP=137
PRES="1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,700.,650.,600.,550.,500.,450.,400.,350.,300.,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1."
PASLAT=1.25
PASLONG=1.25
LATMIN=-90
LATMAX=90

## CONFIG
NDEG=6
NLIS2D=15
NHGRID=6
NWRITE=6
NHMOD=6
NHRELAX=240
NPRHMOD=0
NPRWRITE=0
INDIFEXPL=0
DIFF=4050

#OUTPUT
NRUN=${RUN_ID_NUMBER}
NWTEMP=6
WINDOUT=1
TOUT=1
STATIONS=0

############################################
# Folders/file organization		           #
############################################

cd ${SIMUDIR}

if (( RUN_ID_NUMBER < 10 )) 
then
	RUNDIR="RUN0${RUN_ID_NUMBER}"
else
	RUNDIR="RUN${RUN_ID_NUMBER}"
fi

if [ ! -d "${SIMUDIR}/${RUNDIR}" ]; then mkdir ${SIMUDIR}/${RUNDIR}; fi
if [ ! -d "${SIMUDIR}/${RUNDIR}/DATA" ]; then mkdir ${SIMUDIR}/${RUNDIR}/DATA; fi
if [ ! -d "${SIMUDIR}/${RUNDIR}/reprise" ]; then mkdir ${SIMUDIR}/${RUNDIR}/reprise; fi

if (( INTYPE == 1 )); then
	if [ ! -d "${SIMUDIR}/ECMR" ]; then
		err_msg "Can't find ECMR directory"
		err_msg "End of program"
		exit 1
	fi
elif (( INTYPE == 2  )); then
	if [ ! -d "${SIMUDIR}/GRIB" ]; then
		err_msg "Can't find GRIB directory"
		err_msg "End of program"
		exit 1
	fi
fi

SRCDIR="/usr/local/MIMOSA"
cp ${SRCDIR}/src/mimosa.x ${SIMUDIR}/

info_msg "Start of MIMOSA simulation"

############################################
# 3. Simulation launch			           #
############################################

for i in ${THETA[*]}; do

	info_msg "run for theta = ${i}K"

	# Generating the input.namelist file for each isentropic level requested
	cat > ${SIMUDIR}/input.namelist <<EOF
!======================================================================!
!               __  __ _____ __  __  ____   _____                      !
!              |  \/  |_   _|  \/  |/ __ \ / ____|  /\                 !
!              | \  / | | | | \  / | |  | | (___   /  \                !
!              | |\/| | | | | |\/| | |  | |\___ \ / /\ \               !
!              | |  | |_| |_| |  | | |__| |____) / ____ \              !
!              |_|  |_|_____|_|  |_|\____/|_____/_/    \_\             !
!======================================================================!
!
! CARACTERISCS OF THE RUN
!
&run
!-----
! ZONE defines the geographical area :
!   -> 1 for the Northern Hemisphere [-10N, 90N]
!   -> 2 for the Southern hemisphere [-90N, 10N]
!   -> 3 for the both Hemisphere     [-90N, 90N]
zone = $ZONE
!-----
! INTYPE defines the type of input files :
!   -> 1 for ASCII isobaric files (*.ECMR) 
!   -> 2 for GRIB encoded model levels files (*.grib)
intype  = $INTYPE
!-----
! IAND, MOISD, JOURD, IHEURED define the starting date of the simulation
!  -> Year  (YY)
iand    = $SYEAR
!  -> Month (MM)
moisd   = $SMONTH
!  -> Day   (DD)
jourd   = $SDAY
!  -> Hour  (HH)
iheured = $SHOUR
!-----
! IANF, MOISF, JOURF, IHEUREF define the final date of the simulation
!  -> Year  (YY) 
ianf    = $EYEAR
!  -> Month (MM)
moisf   = $EMONTH
!  -> Day   (DD)
jourf   = $EDAY
!  -> Hour  (HH)
iheuref = $EHOUR
!-----
! INITPV defines if the simulation is new or a restart
!  -> 1 if initialization is need 
!  -> 0 if a ph* file from a previous run should be read
initpv  = $INITPV
!-----
! TETA defines the isentropic surface (K)
!   -> shouldn't be greater than 950K for isobaric input files (intype = 1)
teta    = $i
/
!
! CARACTERISCS OF ECMWF GRID
!
&grid
!-----
! NX, NY and NP define the number of points of the ECMWF grid
!   -> Number of grid points along longitudes 
nx = $NX
!   -> Number of grid points along latitudes
ny = $NY
!   -> Number of pressure levels (intype = 1) or number or model levels (inType = 2) 
np = $NP
!-----
! PRES allows to define the pressure levels of isobaric file (intype = 1)
!   -> this variable is not needed for GRIB files (intype = 2)
!   -> If there is more than 50 levels, declaration of pres variable in mimosa.f95 
pres(1) = $PRES
!-----
! PASLAT and PASLONG define the horizontal resolution of the ECMWF grid
!   -> resolution along latitude
paslat = $PASLAT
!   -> resolution along longitude
paslong = $PASLONG
!-----
! LATMINECMR and LATMAXECMR define the minumum and maximum latitude of ECMWF grid
!   -> minimum latitude 
latminecmr = $LATMIN
!   -> maximum latitude
latmaxecmr = $LATMAX
/
!
! CONFIGURATION OF THE SIMULATION
!
&config
!-----
!  NDEG defines the number of MIMOSA grid points per degree of latitude and longitude
!    -> value should be 3 or 6
ndeg   = $NDEG
!-----
!  NLIS2D defines the number of points use for the smooth
nlis2d = $NLIS2D
!-----
! NHGRID defines the number of hours between two call to regrid
nhgrid = $NHGRID
!-----
! NWRITE defines the number of hour between two outputs of PV files
nwrite = $NWRITE
!-----
! NHMOD defines the number of hours between two ECMWF files
nhmod = $NHMOD
!-----
! NHRELAX defines the number of hours of relaxation time
nhrelax = $NHRELAX
!-----
! NPRHMOD defines the hour of the first ECMWF file
nprhmod = $NPRHMOD
!-----
! NPRWRITE defines the first hour of PVF file
nprwrite = $NPRWRITE
!-----
! Explicit diffusion
! INDIFEXPL defines if explicit diffusion is activated
!   -> 0 no explicit diffusion
!   -> 1 explicit diffusion
indifexpl = $INDIFEXPL
! DIFF defines the value of the explicit diffusion if it is activated
diff = $DIFF
/
!
! OUTPUT OF THE SIMULATION
!
&output
!-----
! NRUN defines the name of the directory where MIMOSA outputs will be saved
!   -> if nrun =  5, files will be saved in RUN05 directory
!   -> if nrun = 16, files will be saved in RUN16 directory
nrun = $NRUN
!-----
! NWTEMP defines the time between two output of temperature or wind     
nwtemp = $NWTEMP
!-----
! WIND_OUT defines if wind horizontal components files will be saved
!   -> 0 no output of wind files
!   -> 1 output of wind files
wind_out = $WINDOUT
!-----
! T_OUT defines if temperature files will be saved
!   -> 0 no output of temperature files
!   -> 1 output of temperature files
t_out = $TOUT
!-----
! STATIONS_OUT defines if PV and temperature and PV profiles at stations files will be saved
!   -> 0 no output of stations files
!   -> 1 output of stations files
stations_out = $STATIONS
/ 
EOF

	${SIMUDIR}/mimosa.x
done

info_msg "Simulations for all thetas done"

# \rm ${SIMUDIR}/input.namelist
# \rm fdat zdat

info_msg "Moving data into DATA and reprise folders..."
mv ${SIMUDIR}/${RUNDIR}/phn* ${SIMUDIR}/${RUNDIR}/reprise
mv ${SIMUDIR}/${RUNDIR}/phs* ${SIMUDIR}/${RUNDIR}/reprise
mv ${SIMUDIR}/${RUNDIR}/pvg* ${SIMUDIR}/${RUNDIR}/DATA
mv ${SIMUDIR}/${RUNDIR}/tg* ${SIMUDIR}/${RUNDIR}/DATA
mv ${SIMUDIR}/${RUNDIR}/ug* ${SIMUDIR}/${RUNDIR}/DATA
mv ${SIMUDIR}/${RUNDIR}/vg* ${SIMUDIR}/${RUNDIR}/DATA
info_msg "DONE"

info_msg "Post processing simulation results"
mkdir -p ${SIMUDIR}/${RUNDIR}/IMAGES
python3 /usr/local/MIMOSA/post-process-mimosa.py --out-dir ${SIMUDIR}/${RUNDIR}/DATA --im-dir ${SIMUDIR}/${RUNDIR}/IMAGES

info_msg "End of MIMOSA simulation"

exit 0
