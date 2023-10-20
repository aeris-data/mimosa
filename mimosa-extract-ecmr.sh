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
#   SPATIAL_RESOLUTION  = "1.125"               : spatial resolution of the data in degrees
#   DATA_DIR            = "/home/path/to/data/dir"      : path to the data directory
#   WORKING_DIR         = "/home/path/to/working/dir"   : path to the working directory
########################################
# PARAMETER DEFINITION ZONE FOR THE USER

START_DATE="20230101"
END_DATE="20230105"
SPATIAL_RESOLUTION="1.125"
DATA_DIR=$(pwd)
WORKING_DIR=$(pwd)

########################################
# DO NOT CHANGE ANYTHING BELOW

set -e
module load ecmwf-toolbox

cd ${DATA_DIR}

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
      write(nomfic,2000) ian,mois,jour,itime
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

DATE=${START_DATE}
while [ ${DATE} -le ${END_DATE} ]; do

    cat > ${WORKING_DIR}/data.req <<EOF
retrieve,
    class    = od,
    format   = packed,
    type     = an,
    levtype  = pl,
    level    = 1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500,
    time     = 00,
    param    = t/u/v,
    date     = ${DATE},
    grid     = ${SPATIAL_RESOLUTION}/${SPATIAL_RESOLUTION},
    area     = 90/0/-90/360,
    target   = "${DATA_DIR}/data00"
retrieve,
    time     = 12,
    target   = "${DATA_DIR}/data12"
EOF

      mars ${WORKING_DIR}/data.req

      grib_get -w level=1,shortName=t -p dataDate,dataTime ${DATA_DIR}/data00 > ${DATA_DIR}/datafile
      grib_get_data -w shortName=t ${DATA_DIR}/data00 >> ${DATA_DIR}/datafile
      grib_get_data -w shortName=u ${DATA_DIR}/data00 >> ${DATA_DIR}/datafile
      grib_get_data -w shortName=v ${DATA_DIR}/data00 >> ${DATA_DIR}/datafile
      ${DATA_DIR}/grib_mimosa

      grib_get -w level=1,shortName=t -p dataDate,dataTime ${DATA_DIR}/data12 > ${DATA_DIR}/datafile
      grib_get_data -w shortName=t ${DATA_DIR}/data12 >> ${DATA_DIR}/datafile
      grib_get_data -w shortName=u ${DATA_DIR}/data12 >> ${DATA_DIR}/datafile
      grib_get_data -w shortName=v ${DATA_DIR}/data12 >> ${DATA_DIR}/datafile
      ${DATA_DIR}/grib_mimosa

      DATE=`expr ${DATE} + 1`
done

rm ${DATA_DIR}/data00 ${DATA_DIR}/data12 ${DATA_DIR}/datafile
