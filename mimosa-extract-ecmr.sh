#!/bin/bash

########################################
#  SCRIPT FOR              __   _      #
#   THE MIMOSA TOOL      _(  )_( )_    #
#     ECMF DATA         (_   _    _)   #
#        EXTRACTION    / /(_) (__)     #
#                     / / / / / /      #
#                    / / / / / /       #
########################################
#
# This script handles the ECMWF data extraction needed for the MIMOSA simulations.
# The user have to define the start and end date of the data extraction, the spatial 
# resolution of the data and the paths to the working directory and data directory.
# The working directory will just contain the request file for the MARS API; the data
# directory will contain the final extracted data; the working and data directories can be the same.
# The data will be extracted on the pressire levels 1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500.
# The grid spacing needs to be an integer fraction of 90 degrees (latitude range from the Equator to
# the Pole), e.g. 0.225, 0.25, 0.28125, 0.3, 0.4 etc.
#
# PARAMETERS SYNTAXE:
#   START_DATE          = "YYYYMMDD"                    : start date of the data extraction
#   END_DATE            = "YYYYMMDD"                    : end date of the data extraction
#   SPATIAL_RESOLUTION  = "1.125"                       : spatial resolution of the data in degrees
#   DATA_DIR            = "/home/path/to/data/dir"      : path to the data directory
#   WORKING_DIR         = "/home/path/to/working/dir"   : path to the working directory
########################################
# DO NOT CHANGE ANYTHING BELOW

set -e
module load ecmwf-toolbox

function help() {
    bold=$(tput bold)
    normal=$(tput sgr0)
    echo "# ########################################"
    echo "# #                          __   _      #"
    echo "# #  SCRIPT FOR            _(  )_( )_    #"
    echo "# #   THE ECMWF           (_   _    _)   #"
    echo "# #    DATA EXTRACTION   / /(_) (__)     #"
    echo "# #                     / / / / / /      #"
    echo "# #                    / / / / / /       #"
    echo "# ########################################"
    echo "#"
    echo "# ${bold}${SCRIPT_NAME}${normal} script extracts and formats the ECMWF data into"
    echo "# files that are accepted by the simulation MIMOSA tool."
    echo "#"
    echo "# The script knows which meterological data and variables are needed; the user input"
    echo "# must only define the dates, the spatial resolution and additional logistical information,"
    echo "# like the data directory and the working directory."
    echo "# The extraction is performed on the MARS server, so the script should be launched on this server."
    echo "#"
    echo "# Usage: ${SCRIPT_NAME} [options] arguments"
    echo "# Options:"
    echo "#   ${bold}-h, --help${normal}     Show this help message and exit"
    echo "# Arguments:"
    echo "#   ${bold}--config conf_filepath${normal}  This argument must correspond to the configuration"
    echo "# file where the user defines input parameters needed for the extraction See an example of a"
    echo "# configuration file below."
    echo "#"
    echo "# +--------------------------------------------------------------------------------+"
    echo "# | Example of the content in the configuration file                               |"
    echo "# |                                                                                |"
    echo "# | Example filename : ${bold}my_parameters.conf${normal}                                          |"
    echo "# | Example content below :                                                        |"
    echo "# |                                                                                |"
    echo "# | START_DATE='20230101'                                                          |"
    echo "# | END_DATE='20230105'                                                            |"
    echo "# | SPATIAL_RESOLUTION='1.125'                                                     |"
    echo "# | DATA_DIR='/my/dir/for/data'                                                    |"
    echo "# | WORKING_DIR='/working/dir/for/aux/files'                                       |"
    echo "# +--------------------------------------------------------------------------------+"
}

function check_args() {
    if [[ -z ${START_DATE} ]]; then
        err_msg "No start date was defined. Exiting script."
        exit 1
    fi
    if [[ -z ${END_DATE} ]]; then
        err_msg "No end date was defined. Exiting script."
        exit 1
    fi
    if [[ -z ${SPATIAL_RESOLUTION} ]]; then
        err_msg "No spatial resolution was defined. Exiting script."
        exit 1
    fi
    if [[ -z ${DATA_DIR} ]]; then
        err_msg "No data directory was defined. Exiting script."
        exit 1
    fi
    if [[ -z ${WORKING_DIR} ]]; then
        err_msg "No working directory was defined. Exiting script."
        exit 1
    fi
    if [[ ${START_DATE} =~ ^[0-9]{8}$ ]]; then
        exit_status=0
    else
        err_msg "Starting date 'START_DATE' is not in the correct format YYYYMMDD"
        exit 1
    fi
    if [[ ${END_DATE} =~ ^[0-9]{8}$ ]]; then
        exit_status=0
    else
        err_msg "Ending date 'END_DATE' is not in the correct format YYYYMMDD"
        exit 1
    fi
    if [ $START_DATE -le $END_DATE ]; then
        exit_status=0
    else
        err_msg "Start date is greater than end date"
        exit 1
    fi
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

function check_the_date(){
    today_date_obj=$(date -d "$(date +'%Y%m%d')" +%s)
    end_fc_date=$(date -d "$(date +'%Y%m%d')+10 days" +%s)
    start_date_obj=$(date -d "${START_DATE}" +%s)
    end_date_obj=$(date -d "${END_DATE}" +%s)
    if [[ ${start_date_obj} -gt ${end_date_obj} ]]; then
        echo  "$(date +'%d/%m/%Y - %H:%M:%S') - ERROR : Your start date is greater that the end date"
        exit 1
    elif [[ ${start_date_obj} -gt ${today_date_obj} ]]; then
        echo  "$(date +'%d/%m/%Y - %H:%M:%S') - ERROR : Your start date is greater that today date"
        exit 1
    elif [[ ${end_date_obj} -gt ${end_fc_date} ]]; then
        echo  "$(date +'%d/%m/%Y - %H:%M:%S') - WARNING : Your end date exceeds TODAY+10 days, extracting data from your start date and up to 10 days from TODAY"
    fi
}

function write_fortran_code(){
    NX=$(bc <<< "scale=0; 360/${SPATIAL_RESOLUTION}")
    NY=$(bc <<< "scale=0; 180/${SPATIAL_RESOLUTION} + 1")
    NDATA=$((${NX}*${NY}))

    cat > ${WORKING_DIR}/grib_mimosa.f90 << EOF
PROGRAM grib_mimosa
     implicit none
!
     CHARACTER(LEN=*), PARAMETER             :: input_file='datafile'
     CHARACTER(LEN=*), PARAMETER             :: open_mode='r'

! Local type declaration statements
     integer, parameter :: np=17, nx=${NX}, ny=${NY}
     integer, parameter :: nparam=3, ndata=${NDATA}
     real*8, parameter :: alatdeb=-90., alongdeb=0.
     real*8, parameter :: pasd=${SPATIAL_RESOLUTION}

!    var: champs ecmwf lus sur fichiers isobaric
     real*8 :: tab(nparam,np,ndata)
!
     real*8 :: press(np),lat,lon
     data press /1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 50.0, &
     70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0/

     character(len=14) :: nomfic
     character(len=1) :: cparam(nparam)
     character(len=5) :: unit(nparam)
     character(len=26) :: entete
     integer :: idatmes, ian, mois, jour, &
                itimes,itime, i, j, l, k, err
     data cparam /'T', 'U', 'V'/
     data unit /'    K', '  m/s', '  m/s'/

     call get_command_argument(1, value=nomfic)

!
!     Open data file for reading.
!
      open(25,file=trim(input_file),iostat=err,status='old')
      if (err /= 0) print*, "Attention le fichier d entree n existe pas."
      read(25,*,iostat=err) idatmes,itimes
      do i=1,nparam
         do j=1,np
            read(25,*,iostat=err) entete
            do k=1,ndata
               read(25,*,iostat=err) lat,lon,tab(i,j,k)
            enddo
         enddo
      enddo
!
!
      close(25)
!
!     Creation nom du fichier resultat
      ian=(idatmes-20000000)/10000
      mois=((idatmes-20000000)-(ian*10000))/100
      jour=((idatmes-20000000)-(ian*10000))-(mois*100)
      itime=itimes/100
      !write(nomfic,2000) ian,mois,jour,itime
      open(20,file=nomfic)
!
!     icriture du fichier resultat
      do i=1,nparam
         do j=1,np
         write(20,3000)ian,mois,jour,itime,cparam(i),unit(i),press(j),&
           nx,alongdeb,pasd,ny,alatdeb,pasd
          do k=1,ny
            write(20,3001) (tab(i,j,l),l=((ny-k)*nx)+1,((ny-k)*nx)+nx)
          enddo
         enddo
      enddo

2000  format('D',4i2.2,'.ECMR')
3000  format(1x,'Echeance 20',i2,3i3,/,1x,'Parameter',1x, &
      a1,/,1x,'Unit ',a5,/,1x,'Level ',f6.1,' hPa',/,1x, &
      'grid longitude',i4,2f8.3,/,1x,'grid latitude ',i4,2f8.3)
3001  format(10f7.2)

END PROGRAM grib_mimosa
EOF

    gfortran -march=native -O3 -fno-sign-zero -o ${WORKING_DIR}/grib_mimosa ${WORKING_DIR}/grib_mimosa.f90
    if [[ "${WORKING_DIR}" != "${DATA_DIR}" ]]; then cp ${WORKING_DIR}/grib_mimosa ${DATA_DIR}/grib_mimosa; fi
}

function main(){
    ########################################
    check_the_date
    ########################################
    if [ ! -d ${DATA_DIR} ]; then
        mkdir -p ${DATA_DIR}
    fi
    if [ ! -d ${WORKING_DIR} ]; then
        mkdir -p ${WORKING_DIR}
    fi
    ########################################
    write_fortran_code
    ########################################
    DATE=$(date -d ${START_DATE} +%s)
    while [ ${DATE} -le $(date -d ${END_DATE} +%s) ] && [ ${DATE} -lt  $(date -d "$(date +'%Y-%m-%d') + 10 days" +%s) ]; do

        if [ ${DATE} -lt $(date -d "$(date +'%Y-%m-%d')" +%s) ]; then
            cat > ${WORKING_DIR}/data.req <<EOF
retrieve,
    class    = od,
    format   = packed,
    type     = an,
    levtype  = pl,
    level    = 1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500,
    time     = 00,
    param    = t/u/v,
    date     = $(date -d "@${DATE}" +'%Y%m%d'),
    grid     = ${SPATIAL_RESOLUTION}/${SPATIAL_RESOLUTION},
    area     = 90/0/-90/360,
    target   = "${DATA_DIR}/data00"
retrieve,
    time     = 12,
    target   = "${DATA_DIR}/data12"
EOF
        FILENAME_00="D$(date -d "@${DATE}" +'%y%m%d')00.ECMR"
        FILENAME_12="D$(date -d "@${DATE}" +'%y%m%d')12.ECMR"
        fi

        if [ ${DATE} -ge $(date -d "$(date +'%Y-%m-%d')" +%s) ]; then
            st=$(date -d "$(date +'%Y-%m-%d')" +%s)
            en=${DATE}
            fc_hour=$((((en-st)/86400)*24))
            cat > ${WORKING_DIR}/data.req <<EOF
retrieve,
    class    = od,
    format   = packed,
    type     = fc,
    levtype  = pl,
    level    = 1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500,
    time     = 00,
    step     = ${fc_hour},
    param    = t/u/v,
    date     = $(date +'%Y%m%d'),
    grid     = ${SPATIAL_RESOLUTION}/${SPATIAL_RESOLUTION},
    area     = 90/0/-90/360,
    target   = "${DATA_DIR}/data00"
retrieve,
    time     = 00,
    step     = $((${fc_hour}+12)),
    target   = "${DATA_DIR}/data12"
EOF
        FILENAME_00=D$(date -d "$(date +'%Y%m%d') + ${fc_hour} hours" +%y%m%d%H).ECMR
        FILENAME_12=D$(date -d "$(date +'%Y%m%d') + $((${fc_hour}+12)) hours" +%y%m%d%H).ECMR
        fi

        mars ${WORKING_DIR}/data.req

        grib_get -w level=1,shortName=t -p dataDate,dataTime ${DATA_DIR}/data00 > ${DATA_DIR}/datafile
        grib_get_data -w shortName=t ${DATA_DIR}/data00 >> ${DATA_DIR}/datafile
        grib_get_data -w shortName=u ${DATA_DIR}/data00 >> ${DATA_DIR}/datafile
        grib_get_data -w shortName=v ${DATA_DIR}/data00 >> ${DATA_DIR}/datafile
        ${DATA_DIR}/grib_mimosa ${FILENAME_00}

        grib_get -w level=1,shortName=t -p dataDate,dataTime ${DATA_DIR}/data12 > ${DATA_DIR}/datafile
        grib_get_data -w shortName=t ${DATA_DIR}/data12 >> ${DATA_DIR}/datafile
        grib_get_data -w shortName=u ${DATA_DIR}/data12 >> ${DATA_DIR}/datafile
        grib_get_data -w shortName=v ${DATA_DIR}/data12 >> ${DATA_DIR}/datafile
        ${DATA_DIR}/grib_mimosa ${FILENAME_12}

        DATE=`expr ${DATE} + 86400`
    done

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
check_args
main
