# MIMOSA

MIMOSA is a high-resolution potential vorticity advection model developed in Fortran by A. Hauchecorne (Hauchecorne et al., 2002). It is initialized at a time `t` from ECMWF data (horizontal wind fields U and V, temperature and pressure) on an orthogonal grid centered on the north pole. MIMOSA calculates then advects the potential vorticity on isentropic surfaces with a spacial resolution of 1/3 or 1/6 of a degree.

## Requirements
The MIMOSA tool is containerized into a Singularity container so one must have Singularity installed on the host system intended for simulations.

## Installation
1. `git clone https://github.com/aeris-data/mimosa.git`
2. `sudo singularity build ./mimosa.sif ./mimosa-container.def`
The `singularity build` command will build the container `mimosa.sif` from its definition file, using the source files got from the git repo; so for the build it is important to call the command from the git repo directory that one has made. ⚠️ ***The build requires sudo rights.*** Afterwards, the sif image can be placed anywhere (even on another system) independently of the source files. To run the image no sudo rights are required.

## Usage
The main script is `mimosa-user-script.sh` which needs the input configuration file user-config.conf (which can be renamed, the name is not important). This bash script handles user's input parameters, launch simulations and post-process simulation results. The main usage is `./mimosa-user-script.sh --config user-config.conf`. The script must be launched inside the Singularity container. In the simulation working directory the one must have a tree folder `GRIB/[year of the data in YYYY format]/[month of the data in MM format]` where the meteorological GRIB files with names `DYYMMDDHH.grib` must be stored. The outputs of the simulation are : estimated temperature and potential vorticity in binary and netCDF format + PNG map plots of the data. More details about input/output and folder structure are in the manual `XXX.pdf`.

## Input meteorological data extraction
The input data for the simulations is meteorological data : wind, temperature and logarithm of surface pressure, coming from the ECMWF database. To extract and prepare the data in the correct format, the script `mimosa-extract-grib.sh` or `mimosa-extract-ecmr.sh` should be used. These scripts extract the data either in the GRIB or ASCII (ECMR) format, respectively. The user can configure the start and end date of the data, as well as the spatial resolution, and the data class (only in the grib version). The configuration of two data extractions are as follows:
    - GRIB data
        - extracted on ECMWF 137 model levels
        - the timestep is 3 hours if the requested date range is up to J+6; if the end date exceeds the J+6 limit, the timestep is 6 hours
    - ECMR data
        - extracted on 17 pressure levels
        - the timestep is 12 hours
The script must be launched on the ECMWF MARS server (ecs, hpc or other). The data extraction was tested with a member-state user account. Other more public accounts might customize the script based on the MARS services or APIs available for their type of user. The data is extracted and stored in the directory requested in the input configuration; afterwards, the data can be used for the simulation.
