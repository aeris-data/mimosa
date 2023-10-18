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
# To launch the simulation you need to changes parameters below for your
# simulation, and launch it inside the MIMOSA Singularity container
# 
# History:
#   200903 MAD/LPC2E Creation
#   201109 MAD/LPC2E Group all parameters
#                    Change creation of input.namelist due to change in MIMOSA model
#                    Add creation and checking of directory needed for the simulation
#   2011 -> 2023 CB  Maintenance
#   2023   Updated by Daria MALIK (Magellium) to use only this one script
#          to launch simulation
#
# Author:
#   M.-A. Drouin & Cathy Boonne
#   (updated by D.Malik)
# ----------------------------------------------------------------------------------------------

set -e
PYTHON_SCRIPT="./post-process-mimosa.py"

function help() {
    bold=$(tput bold)
    normal=$(tput sgr0)
    echo "#       .-'';'-."
	echo "#     ,'   <_,-.\`.       MIMOSA tool for"
	echo "#    /)   ,--,_>\_\       high-resolution"
	echo "#   |'   (      \_ |       potential vorticity"
	echo "#   |_    \`-.    / |         advection model"
	echo "#    \\\`-.   ;  _(\`/"
	echo "#     \`.(    \/ ,' "
	echo "#       \`-....-'"
    echo "#"
    echo "# This script handles the MIMOSA simulation input"
    echo "# parameters and launches the simulation, as well as"
	echo "# the post-processing for the output additional"
	echo "# reformating and results visualization"
    echo "#"
    echo "# Usage: ${SCRIPT_NAME} [options] arguments"
    echo "# Options:"
    echo "#   ${bold}-h, --help${normal}     Show this help message and exit"
    echo "# Arguments:"
    echo "#   ${bold}--config conf_fielpath${normal}  This argument must correspond to the configuration"
    echo "# file where the user defines input parameters needed for the extraction"
    echo "#"
}

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

function main(){
	info_msg "!======================================================================!"
	info_msg "!               __  __ _____ __  __  ____   _____                      !"
	info_msg "!              |  \/  |_   _|  \/  |/ __ \ / ____|  /\                 !"
	info_msg "!              | \  / | | | | \  / | |  | | (___   /  \                !"
	info_msg "!              | |\/| | | | | |\/| | |  | |\___ \ / /\ \               !"
	info_msg "!              | |  | |_| |_| |  | | |__| |____) / ____ \              !"
	info_msg "!              |_|  |_|_____|_|  |_|\____/|_____/_/    \_\             !"
	info_msg "!======================================================================!"

	############################################
	# Folders/file organization		           #
	############################################

	cd ${SIMUDIR}

	if (( RUN_ID_NUMBER < 10 )) 
	then
		RUNDIR="RUN0${NRUN}"
	else
		RUNDIR="RUN${NRUN}"
	fi

	if [ ! -d "${SIMUDIR}/${RUNDIR}" ]; then mkdir ${SIMUDIR}/${RUNDIR}; fi
	if [ ! -d "${SIMUDIR}/${RUNDIR}/DATA" ]; then mkdir ${SIMUDIR}/${RUNDIR}/DATA; fi

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
	# Parameters configuration                 #
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

	\rm ${SIMUDIR}/input.namelist

	info_msg "Moving data into DATA folder..."
	if ! mv ${SIMUDIR}/${RUNDIR}/pvg* ${SIMUDIR}/${RUNDIR}/DATA; then
		info_msg "No PV output was found"
	fi
	if ! mv ${SIMUDIR}/${RUNDIR}/tg* ${SIMUDIR}/${RUNDIR}/DATA; then
		info_msg "No temperature output was found"
	fi
	if ! mv ${SIMUDIR}/${RUNDIR}/ug* ${SIMUDIR}/${RUNDIR}/DATA; then
		info_msg "No U-wind output was found"
	fi
	if ! mv ${SIMUDIR}/${RUNDIR}/vg* ${SIMUDIR}/${RUNDIR}/DATA; then
		info_msg "No V-wind output was found"
	fi
	info_msg "DONE"

	info_msg "Post processing simulation results"
	mkdir -p ${SIMUDIR}/${RUNDIR}/IMAGES
	python3 ${PYTHON_SCRIPT} --out-dir ${SIMUDIR}/${RUNDIR}/DATA --im-dir ${SIMUDIR}/${RUNDIR}/IMAGES

	info_msg "End of MIMOSA simulation"

	exit 0
}

# ----------------------------------------------------------------------------------------------
# BASH SCRIPT
# ----------------------------------------------------------------------------------------------

opts=$(getopt --longoptions "config:,help" --name "$(basename "$0")" --options "h" -- "$@")
eval set --$opts

while [[ $# -gt 0 ]]; do
	case "$1" in
		--config) shift; CONFIG_FILE=$1; shift;;
        -h|--help) help; exit 0; shift;;
		\?) shift; err_msg "Unrecognized options"; exit 1; shift;;
		--) break;;
	esac
done

if [[ -z ${CONFIG_FILE} ]]; then
    err_msg "No configuration file was passed. Exiting script."
    exit 1
fi

source ${CONFIG_FILE}
main
