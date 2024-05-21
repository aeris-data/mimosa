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
# The user have to define the start and end date of the data extraction, the class of the ECMWF 
# reanalysis data (operational, ERA40, ERA5) and the paths to the working directory and data
# directory. The working directory will just contain the request file for the MARS API; the data
# directory will contain the final extracted data; the working and data directories can be the same.
# The data will be extracted on the model levels of the ECMWF, and elvel values will be
# adapted based on the requested data class. The data classes do not have the same temporal
# coverage of the data, so make sure to check that the requested class and start/end dates are
# coherent. Visit this link to navigate through MARS parameters and verify that your data
# class VS start/end dates are correct ---> https://apps.ecmwf.int/mars-catalogue/
#
# PARAMETERS SYNTAXE:
#   START_DATE          = "YYYYMMDD"                    : start date of the data extraction
#   END_DATE            = "YYYYMMDD"                    : end date of the data extraction
#   DATA_CLASS          = "class_name"                  : class name of the data ("od"=operational, "e4":ERA40, "ei":ERA-Interim, "ea":ERA5)
#   SPATIAL_RESOLUTION  = "1.125"                       : spatial resolution of the data in degrees
#   DATA_DIR            = "/home/path/to/data/dir"      : path to the data directory
#   WORKING_DIR         = "/home/path/to/working/dir"   : path to the working directory
########################################
# DO NOT CHANGE ANYTHING BELOW

set -e

SCRIPT_NAME=$(basename "$0")

function help() {
    bold=$(tput bold)
    normal=$(tput sgr0)
    echo ""
    echo "###   SCRIPT FOR            _(  )_( )_    "
    echo "###    THE ECMWF           (_   _    _)   "
    echo "###     DATA EXTRACTION   / /(_) (__)     "
    echo "###                      / / / / / /      "
    echo "###                     / / / / / /       "
    echo "###"
    echo "### ${bold}${SCRIPT_NAME}${normal} script extracts and formats the ECMWF data into"
    echo "### files that are accepted by the simulation MIMOSA tool."
    echo "###"
    echo "### The script knows which meterological data and variables are needed; the user input"
    echo "### must only define the dates, the spatial resolution and additional logistical information,"
    echo "### like the data directory and the working directory."
    echo "### The extraction is performed on the MARS server, hence the script must be launched on MARS server."
    echo "###"
    echo "### Usage: ${SCRIPT_NAME} [options] arguments"
    echo "### Options:"
    echo "###   ${bold}-h, --help${normal}              Show this help message and exit"
    echo "### Arguments:"
    echo "###   ${bold}--config conf_filepath${normal}  This argument must correspond to the configuration"
    echo "###                           file where the user defines input parameters needed"
    echo "###                           for the extraction See an example of a configuration file below."
    echo "###"
    echo "### +--------------------------------------------------------------------------------+"
    # echo "### | Example of the content in the configuration file                               |"
    # echo "### |                                                                                |"
    echo "### | Example filename : ${bold}my_parameters.conf${normal}                                          |"
    echo "### | Example content below :                                                        |"
    echo "### |                                                                                |"
    echo "### | START_DATE='20230101'                                                          |"
    echo "### | END_DATE='20230105'                                                            |"
    echo "### | DATA_CLASS='od'                                                                |"
    echo "### | SPATIAL_RESOLUTION='1.125'                                                     |"
    echo "### | DATA_DIR='/my/dir/for/data'                                                    |"
    echo "### | WORKING_DIR='/working/dir/for/aux/files'                                       |"
    echo "### +--------------------------------------------------------------------------------+"
    echo ""
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
    if [[ -z ${DATA_CLASS} ]]; then
        err_msg "No data class was defined. Exiting script."
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

function count_fc_hours(){
    local _N_DAYS=$1
    local _N_HOURS=$((${_N_DAYS}*24))
    local result=""
    if (( ${_N_HOURS} <= 144 )); then result=$(seq -s '/' 0 3 ${_N_HOURS}); fi
    if (( ${_N_HOURS} > 144 )); then result=$(seq -s '/' 0 6 ${_N_HOURS}); fi
    echo ${result}
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
        END_DATE=$(date -d "$(date +'%Y%m%d')+9 days" +%Y%m%d)
        end_date_obj=$(date -d "${END_DATE}" +%s)
    fi
    if [[ ${START_DATE} == $(date +%Y%m%d) ]]; then
        DATA_TYPE="fc"
        N_DAYS=$(( $((${end_date_obj} - ${start_date_obj})) / 86400 + 1 ))
    elif [[ ${start_date_obj} -lt ${today_date_obj} ]]; then
        if [[ ${end_date_obj} -lt ${today_date_obj} ]]; then
            DATA_TYPE="an"
        elif [[ ${end_date_obj} -ge ${today_date_obj} ]]; then
            DATA_TYPE="an+fc"
        fi
    fi
}

function main(){
    export MARS_MULTITARGET_STRICT_FORMAT=1
    module load ecmwf-toolbox
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
    if [[ ${DATA_TYPE} == "an" ]]; then
        FC_HOURS=""
        TODAY_DATE=""
        END_AN_DATE=""
        AN_N_DAYS=""
        FC_N_DAYS=""
    fi
    if [[ ${DATA_TYPE} == "fc" ]]; then
        FC_HOURS=$(count_fc_hours ${N_DAYS})
        IFS="/" read -ra NUMBERS <<< "${FC_HOURS}"
        FC_DELTA=$(( NUMBERS[1] - NUMBERS[0] ))
        TODAY_DATE=""
        END_AN_DATE=""
        AN_N_DAYS=""
        FC_N_DAYS=""
    fi
    if [[ ${DATA_TYPE} == "an+fc" ]]; then
        TODAY_DATE=$(date +'%Y%m%d')
        END_AN_DATE=$(date -d "${TODAY_DATE} - 1 day" "+%Y%m%d")
        AN_N_DAYS=$(( ($(date -d ${TODAY_DATE} +%s) - $(date -d ${START_DATE} +%s)) / 86400 ))
        FC_N_DAYS=$(( ($(date -d ${END_DATE} +%s) - $(date -d ${TODAY_DATE} +%s) + 86400) / 86400 ))
        FC_HOURS=$((${FC_N_DAYS}*24))
    fi
    ########################################
    REQ_FILEPATH="${WORKING_DIR}/${START_DATE}_${END_DATE}.req"
    if [[ ${DATA_TYPE} == "an" ]]; then
        cat > ${REQ_FILEPATH} <<EOF
retrieve,
    class    = ${DATA_CLASS},
    format   = packed,
    type     = an,
    levtype  = ml,
    levelist = 1,
    time     = 00/06/12/18,
    param    = 152,
    date     = ${START_DATE}/to/${END_DATE},
    grid     = ${SPATIAL_RESOLUTION}/${SPATIAL_RESOLUTION},
    area     = 90/0/-90/360,
    target   = "${DATA_DIR}/[date]_[time]_[step].grib"
retrieve,
    levelist = 1/to/137,
    param    = t/u/v,
    target = "${DATA_DIR}/[date]_[time]_[step].grib"
EOF
    fi

    if [[ ${DATA_TYPE} == "fc" ]]; then
        cat > ${REQ_FILEPATH} <<EOF
retrieve,
    class    = ${DATA_CLASS},
    format   = packed,
    type     = fc,
    levtype  = ml,
    levelist = 1,
    time     = 00,
    step     = ${FC_HOURS},
    param    = 152,
    date     = ${START_DATE},
    grid     = ${SPATIAL_RESOLUTION}/${SPATIAL_RESOLUTION},
    area     = 90/0/-90/360,
    target   = "${DATA_DIR}/[date]_[time]_[step].grib"
retrieve,
    levelist = 1/to/137,
    param    = t/u/v,
    target = "${DATA_DIR}/[date]_[time]_[step].grib"
EOF
    fi

    if [[ ${DATA_TYPE} == "an+fc" ]]; then
        cat > ${REQ_FILEPATH} <<EOF
retrieve,
    class    = ${DATA_CLASS},
    format   = packed,
    type     = an,
    levtype  = ml,
    levelist = 1,
    time     = 00/06/12/18,
    param    = 152,
    date     = ${START_DATE}/to/${END_AN_DATE},
    grid     = ${SPATIAL_RESOLUTION}/${SPATIAL_RESOLUTION},
    area     = 90/0/-90/360,
    target   = "${DATA_DIR}/[date]_[time]_[step].grib"
retrieve,
    levelist = 1/to/137,
    param    = t/u/v,
    target = "${DATA_DIR}/[date]_[time]_[step].grib"
retrieve,
    class    = ${DATA_CLASS},
    format   = packed,
    type     = fc,
    levtype  = ml,
    levelist = 1,
    time     = 00,
    step     = 0/to/${FC_HOURS}/by/6,
    param    = 152,
    date     = ${TODAY_DATE},
    grid     = ${SPATIAL_RESOLUTION}/${SPATIAL_RESOLUTION},
    area     = 90/0/-90/360,
    target   = "${DATA_DIR}/[date]_[time]_[step].grib"
retrieve,
    levelist = 1/to/137,
    param    = t/u/v,
    target = "${DATA_DIR}/[date]_[time]_[step].grib"
EOF
    fi

    mars ${REQ_FILEPATH}
    if [ $? != 0 ]; then
        echo " The MARS request failed"
        exit 1
    fi

    for FILE in ${DATA_DIR}/????????_????_*.grib; do
        FILENAME=$(basename ${FILE} ".grib")
        date_part=$(echo "${FILENAME}" | cut -d'_' -f1)
        time_part=$(echo "${FILENAME}" | cut -d'_' -f2)
        fcstep_part=$(echo "${FILENAME}" | cut -d'_' -f3)
        new_name=$(date --date "${date_part} ${time_part} + ${fcstep_part} hours" +"%y%m%d%H")
        mv ${FILE} ${DATA_DIR}/D${new_name}.grib
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
