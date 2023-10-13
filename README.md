# MIMOSA

MIMOSA is a high-resolution potential vorticity advection model developed in Fortran by A. Hauchecorne (Hauchecorne et al., 2002). It is initialized at a time `t` from ECMWF data (horizontal wind fields U and V, temperature and pressure) on an orthogonal grid centered on the north pole. MIMOSA calculates then advects the potential vorticity on isentropic surfaces with a spacial resolution of 1/3 or 1/6 of a degree.

## Requirements
The MIMOSA tool is containerized into a Singularity container so one must have Singularity installed on the host system intended for simulations.

## Installation
1. `git clone https://github.com/aeris-data/mimosa.git`
2. `sudo singularity build ./mimosa.sif ./mimosa-container.def`
The `singularity build` command will build the container `mimosa.sif` from its definition file, using the source files got from the git repo; so for the build it is important to call the command from the git repo directory that one has made. Afterwards, the sif image can be placed anywhere independently of the source files.

## Usage
The main script is `mimosa-user-script.sh` which needs the input configuration file user-config.conf (which can be renamed, the name is not important). This bash script handles user's input parameters, launch simulations and post-process simulation results. The main usage is `./mimosa-user-script.sh --config user-config.conf`. The script must be launched inside the Singularity container. In the simulation working directory the one must have a tree folder `GRIB/[year of the data in YYYY format]/[month of the data in MM format]` where the meteorological GRIB files with names `DYYMMDDHH.grib` must be stored. More details in the manual `XXX.pdf`.

## Input meteorological data extraction
The input data for the simulations is meteorological data : wind, temperature and logarithm of surface pressure, coming from the ECMWF database. To extract and prepare the data in the correct format, the script `mimosa-extract-ecmwf.sh` should be used. The script must be launched on the ECMWF MARS server (ecs, hpc or other). The data extraction was tested with a member-state user account. Other more public accounts might customize the script based on the MARS services or APIs available for their type of user. The data is extracted and stored in the directory requested in the input configuration; afterwards, the data can be used for the simulation.
