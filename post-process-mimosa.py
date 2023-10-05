# #####################################################################
#                                                                     #
#                   MIMOSA Python post processing                     #
#                 =================================                   #
#                                                                     #
#   This Python script process the output files from the MIMOSA       #
#   simulation to create netCDF files with output results for         #
#   an easier use and create figures with results visualization.      #
#                                                                     #
#   Usage : python post-process-mimosa.py \                           #
#                 --out-dir dir_with_output_files \                   #
#                 --im-dir dir_where_to_save_images                   #
#                                                                     #
#   Created by Daria MALIK (Magellium, 2023)                          #
# #####################################################################

import numpy as np
from scipy.io import FortranFile
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import seaborn as sns
import scipy.signal
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
import copy
import matplotlib.path
import sys
import glob
import netCDF4 as nc
import datetime as dt
import os
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

sys.path.append('/usr/local/idl-colorbars-python')
import idl_colorbars
custom_cmap = idl_colorbars.getcmap(5)

VARIABLES_INFO = {"pvg":{"var_name":"pv",
                         "units":"-",
                         "standard_name":"potential_vorticity",
                         "long_name":"MIMOSA output potential vorticity"},
                  "tg":{"var_name":"t",
                        "units":"K",
                        "standard_name":"temperature",
                        "long_name":"MIMOSA output temperature"},
                  "ug":{"var_name":"u",
                        "units":"m.s-1",
                        "standard_name":"u_component_of_wind",
                        "long_name":"MIMOSA output U wind fields"},
                  "vg":{"var_name":"v",
                        "units":"m.s-1",
                        "standard_name":"v_component_of_wind",
                        "long_name":"MIMOSA output V wind fields"}}

PV_LEVELS = {380:np.arange(0,32,1.5),
             435:np.arange(0,53,2.5),
             475:np.arange(0,85,4),
             550:np.arange(0,190,9),
             675:np.arange(0,631,30),
             775:np.arange(0,946,45),
             950:np.arange(0,3361,160)}

T_LEVELS = np.arange(175,281,5)

# =============================================================================
# MY RESULT
# =============================================================================

def format_datetime(datetime_string, old_format, new_format):
    dt_object = dt.datetime.strptime(datetime_string, old_format)
    new_dt_string = dt.datetime.strftime(dt_object, new_format)
    return new_dt_string

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def create_netcdf_file(fortran_file_filepath: str, output_folder: str) -> str:
    filename = os.path.basename(fortran_file_filepath)
    print(f"Processing {filename}...")
    nc_filepath = f"{output_folder}/{filename}.nc"
    with FortranFile(fortran_file_filepath, "r") as f:
        all_data_bytes = f.read_record(dtype=np.uint8)
    header_size_bytes = 30*4
    data_bytes = all_data_bytes[header_size_bytes:]
    data_values = np.frombuffer(data_bytes, dtype=np.float32)
    header_bytes = all_data_bytes[:header_size_bytes]
    header_values = np.frombuffer(header_bytes, dtype="i4")
    Nx = header_values[19]
    Ny = header_values[20]
    res_x, res_y = 360/(Nx-1), 180/(Ny-1)
    lat = (90 - res_y/2) - np.arange(Ny)*res_y
    lon = (0 + res_x/2) + np.arange(Nx)*res_x
    var_name = filename.split(".")[0][:-8]
    data_values = data_values.reshape((Ny,Nx))
    with nc.Dataset(nc_filepath,"w",format="NETCDF4") as ncfile:
        # ------------------------------------------------------------------------------------
        ncfile.conventions              = "CF-1.0"
        ncfile.netcdf_version_id        = nc.__netcdf4libversion__
        ncfile.standard_name_vocabulary = "NetCDF Standard"
        # ------------------------------------------------------------------------------------
        ncfile.title                    = "MIMOSA (Isentropic Modeling of Mesoscale Transport of Stratospheric Ozone by Advection)"
        ncfile.summary                  = "This file contains output results of the MIMOSA simulation (temperature, potential vorticity or u/v wind components)"
        ncfile.institution              = "OMP / Magellium"
        # ------------------------------------------------------------------------------------
        ncfile.simulation_datetime = format_datetime(f"{header_values[2]}/{header_values[1]}/{header_values[0]} {header_values[3]}:00",
                                                        "%d/%m/%y %H:%M",
                                                        "%Y-%m-%d %H:%M")
        ncfile.initialization_datetime = format_datetime(f"{header_values[6]}/{header_values[5]}/{header_values[4]} {header_values[7]}:00",
                                                            "%d/%m/%y %H:%M",
                                                            "%Y-%m-%d %H:%M")
        # ------------------------------------------------------------------------------------
        ncfile.westernmost_longitude    = header_values[9]
        ncfile.easternmost_longitude    = header_values[11]
        ncfile.southernmost_latitude    = header_values[10]
        ncfile.northernmost_latitude    = header_values[12]
        # ------------------------------------------------------------------------------------
        ncfile.ecmwf_nx = header_values[13]
        ncfile.ecmwf_ny = header_values[14]
        ncfile.ecmwf_timestep_hours = header_values[15]
        ncfile.ecmwf_first_file_hour = header_values[16]
        # ------------------------------------------------------------------------------------
        ncfile.mimosa_nx = header_values[19]
        ncfile.mimosa_ny = header_values[20]
        ncfile.mimosa_calculations_timestep_hours = header_values[21]
        ncfile.output_timestep_hours = header_values[18]
        ncfile.number_of_points_in_1_degree_output = header_values[17]
        ncfile.number_of_hours_before_regridding = header_values[22]
        ncfile.relaxation_hours = header_values[23]
        # ------------------------------------------------------------------------------------
        ncfile.isentropic_level_Kelvin = header_values[8]
        # ------------------------------------------------------------------------------------
        ncfile.createDimension("time",1)
        ncfile.createDimension("lon",header_values[19])
        ncfile.createDimension("lat",header_values[20])
        ncfile.createDimension("theta",1)
        # ------------------------------------------------------------------------------------
        timedelta = dt.datetime.strptime(ncfile.simulation_datetime, "%Y-%m-%d %H:%M") - dt.datetime.strptime("1970-01-01 00:00", "%Y-%m-%d %H:%M")
        timedelta = timedelta.days + timedelta.seconds/(3600*24)
        var = ncfile.createVariable("time", np.float32, ("time",))
        var.units = "days since 1970-01-01 00:00:00"
        var.standard_name = 'time'
        var.long_name = 'Simulation date and time'
        var.calendar = "standard"
        var[:] = timedelta
        # ------------------------------------------------------------------------------------
        var = ncfile.createVariable("theta", np.float32, ("theta",))
        var.units = "K"
        var.standard_name = 'theta'
        var.long_name = 'Isentropic level'
        var[:] = float(header_values[8])
        # ------------------------------------------------------------------------------------
        var = ncfile.createVariable("lat", np.float32, ("lat",))
        var.units = "degrees_north"
        var.standard_name = 'latitude'
        var.long_name = 'Latitude of the output grid'
        var[:] = lat
        # ------------------------------------------------------------------------------------
        var = ncfile.createVariable("lon", np.float32, ("lon",))
        var.units = "degrees_east"
        var.standard_name = 'longitude'
        var.long_name = 'Longitude of the output grid'
        var[:] = lon
        # ------------------------------------------------------------------------------------
        var = ncfile.createVariable(VARIABLES_INFO[var_name]["var_name"], np.float32, ("time","theta","lat","lon"))
        var.units = VARIABLES_INFO[var_name]["units"]
        var.standard_name = VARIABLES_INFO[var_name]["standard_name"]
        var.long_name = VARIABLES_INFO[var_name]["long_name"]
        var[:] = data_values
        # ------------------------------------------------------------------------------------
        return nc_filepath


def create_figure_pv(netcdf_file: str, images_folder: str) -> None:
    gaussian_kernel_sigma = 4
    ds = xr.open_dataset(netcdf_file)
    ds.pv[0,0,:,:] = gaussian_filter(ds.pv[0,0,:,:], sigma=gaussian_kernel_sigma)
    theta_level = ds.isentropic_level_Kelvin
    if theta_level not in PV_LEVELS.keys():
        print(np.array(PV_LEVELS.keys()))
        theta_level = find_nearest(list(PV_LEVELS.keys()), theta_level)
    fig = plt.figure()
    p = (-ds.pv[0,0,:,:]).plot.contourf(levels=PV_LEVELS[theta_level],
                                    transform=ccrs.PlateCarree(),
                                    subplot_kws={"projection": ccrs.SouthPolarStereo()},
                                    cmap=custom_cmap)
    p1 = (-ds.pv[0,0,:,:]).plot.contour(levels=PV_LEVELS[theta_level],
                                    transform=ccrs.PlateCarree(),
                                    subplot_kws={"projection": ccrs.SouthPolarStereo()},
                                    colors="k", linestyles="-", linewidths=0.5)
    p.axes.coastlines(color="w", linewidth=3)
    p.axes.gridlines(linestyle="--", linewidth=0.5)
    p.axes.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
    old_title = p.axes.title.get_text()
    new_title = ds.pv.long_name + "\n" + old_title
    p.axes.set_title(new_title, fontsize=15)
    p.colorbar.ax.yaxis.label.set_fontsize(15)
    p.colorbar.set_ticks(PV_LEVELS[theta_level])
    p.colorbar.ax.tick_params(labelsize=10)
    fig.savefig(f"{images_folder}/{os.path.basename(netcdf_file)[:-3]}.png", dpi=200)
    ds.close()
    plt.close(fig)

def main(output_folder: str, images_folder: str) -> None:
    list_of_fortran_files = glob.glob(f"{output_folder}/*g*.0625")
    for file in list_of_fortran_files:
        nc_file = create_netcdf_file(file, output_folder)
        if "pv" in os.path.basename(nc_file):
            create_figure_pv(nc_file, images_folder)

if __name__=="__main__":
    """
    Main function

    Args:
        --out-dir  : Path to the directory with the Fortran binary output files
        --im-dir   : Path to the directory where to save visualization
    """

    import argparse
    
    parser = argparse.ArgumentParser(description="Post-processing of the MIMOSA Fortran output files for the creation of netCDF copies and results visualization",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--out-dir", type=str, help="Path to the directory with the Fortran binary output files")
    parser.add_argument("--im-dir",  type=str, help="Path to the directory where to save visualization")
    args = parser.parse_args()

    OUTPUT_DIR = args.out_dir
    IMAGES_DIR = args.im_dir

    main(OUTPUT_DIR, IMAGES_DIR)