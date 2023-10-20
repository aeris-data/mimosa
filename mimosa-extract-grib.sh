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
#   START_DATE  = "YYYYMMDD"                    : start date of the data extraction
#   END_DATE    = "YYYYMMDD"                    : end date of the data extraction
#   DATA_CLASS  = "class_name"                  : class name of the data ("od"=operational, "e4":ERA40, "ei":ERA-Interim, "ea":ERA5)
#   DATA_DIR    = "/home/path/to/data/dir"      : path to the data directory
#   WORKING_DIR = "/home/path/to/working/dir"   : path to the working directory
########################################
# PARAMETER DEFINITION ZONE FOR THE USER

START_DATE="20230101"
END_DATE="20230105"
DATA_CLASS="od"
DATA_DIR=$(pwd)
WORKING_DIR=$(pwd)

########################################
# DO NOT CHANGE ANYTHING BELOW

set -e
export MARS_MULTITARGET_STRICT_FORMAT=1
module load ecmwf-toolbox

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
    ########################################
    check_the_date
    ########################################
    DATA_DIR="${DATA_DIR}/${ID_NAME}/${START_DATE}_${END_DATE}"
    if [ ! -d ${DATA_DIR} ]; then
        mkdir -p ${DATA_DIR}
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
    grid     = 1.125/1.125,
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
    grid     = 1.125/1.125,
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
    grid     = 1.125/1.125,
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
    grid     = 1.125/1.125,
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

    for FILE in ${DATA_DIR}/D${DATE}*.grib; do
        FILENAME=$(basename ${FILE} ".grib")
        NEW_NAME="D${FILENAME:3:8}.grib"
        mv ${FILE} ${DATA_DIR}/${NEW_NAME} 
    done

    for FILE in ${DATA_DIR}/????????_????_*.grib; do
        FILENAME=$(basename ${FILE} ".grib")
        date_part=$(echo "${FILENAME}" | cut -d'_' -f1)
        time_part=$(echo "${FILENAME}" | cut -d'_' -f2)
        fcstep_part=$(echo "${FILENAME}" | cut -d'_' -f3)
        new_name=$(date --date "${date_part} ${time_part} + ${fcstep_part} hours" +"%y%m%d%H")
        mv ${FILE} ${DATA_DIR}/D${new_name}.grib
    done
}

main
