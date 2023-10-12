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
import matplotlib.pyplot as plt
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
import logging
from scipy.io import FortranFile

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
# FUNCTIONS
# =============================================================================

def start_log() -> logging.Logger:
    log_handlers = [logging.StreamHandler()]
    logging.basicConfig(format="%(asctime)s   [%(levelname)s]   %(message)s",
                        datefmt="%d/%m/%Y %H:%M:%S",
                        handlers=log_handlers)
    logger = logging.getLogger('my_log')
    logger.setLevel(logging.INFO)
    return logger

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
    nc_filepath = f"{output_folder}/{filename}.nc"
    if os.path.exists(nc_filepath):
        return nc_filepath
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
    # ------------------------------------------------------------------------------------
    # GLOBAL PLOT
    # ------------------------------------------------------------------------------------
    global_levels = np.concatenate(((-PV_LEVELS[theta_level][::-1])[:-1], PV_LEVELS[theta_level]))
    image_filepath = f"{images_folder}/{os.path.basename(netcdf_file)[:-3]}.png"
    if not os.path.exists(image_filepath):
        fig = plt.figure()
        p = (ds.pv[0,0,:,:]).plot.contourf(levels=global_levels,
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.PlateCarree()},
                                        cmap=custom_cmap,
                                        add_colorbar=True,
                                        cbar_kwargs={"fraction": 0.023})
        p1 = (ds.pv[0,0,:,:]).plot.contour(levels=global_levels,
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.PlateCarree()},
                                        colors="k", linestyles="-", linewidths=0.5)
        p.axes.coastlines(color="w", linewidth=1.5)
        obj = p.axes.gridlines(linestyle="--", linewidth=0.5, draw_labels=True)
        obj.top_labels   = False
        obj.right_labels = False
        obj.xlabel_style = {"size":8}
        obj.ylabel_style = {"size":8}
        p.axes.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        title = f"MIMOSA Potential Vorticity\ntime = {format_datetime(str(ds.time.values[0]).split('.')[0],'%Y-%m-%dT%H:%M:%S','%d/%m/%Y %H:%M')}, theta = {theta_level} [K]"
        p.axes.set_title(title, fontsize=10)
        p.colorbar.ax.yaxis.label.set_fontsize(10)
        p.colorbar.set_ticks(global_levels[::2])
        p.colorbar.ax.tick_params(labelsize=8)
        fig.savefig(image_filepath, dpi=200, bbox_inches='tight')
        plt.close(fig)
    # ------------------------------------------------------------------------------------
    # NORTH POLE PLOT
    # ------------------------------------------------------------------------------------
    image_filepath = f"{images_folder}/{os.path.basename(netcdf_file)[:-3].replace('pvg','pvn')}.png"
    if not os.path.exists(image_filepath):
        fig = plt.figure()
        p = (ds.pv[0,0,:,:]).plot.contourf(levels=PV_LEVELS[theta_level],
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.NorthPolarStereo()},
                                        cmap=custom_cmap)
        p1 = (ds.pv[0,0,:,:]).plot.contour(levels=PV_LEVELS[theta_level],
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.NorthPolarStereo()},
                                        colors="k", linestyles="-", linewidths=0.5)
        p.axes.coastlines(color="w", linewidth=1.5)
        obj = p.axes.gridlines(linestyle="--", linewidth=0.5)
        p.axes.set_extent([-180, 180, 30, 90], ccrs.PlateCarree())
        title = f"MIMOSA Potential Vorticity\nNorth Pole (PV=PV)\ntime = {format_datetime(str(ds.time.values[0]).split('.')[0],'%Y-%m-%dT%H:%M:%S','%d/%m/%Y %H:%M')}, theta = {theta_level} [K]"
        p.axes.set_title(title, fontsize=15)
        p.colorbar.ax.yaxis.label.set_fontsize(15)
        p.colorbar.set_ticks(PV_LEVELS[theta_level])
        p.colorbar.ax.tick_params(labelsize=10)
        fig.savefig(image_filepath, dpi=200, bbox_inches='tight')
        ds.close()
        plt.close(fig)
    # ------------------------------------------------------------------------------------
    # SOUTH POLE PLOT
    # ------------------------------------------------------------------------------------
    image_filepath = f"{images_folder}/{os.path.basename(netcdf_file)[:-3].replace('pvg','pvs')}.png"
    if not os.path.exists(image_filepath):
        fig = plt.figure()
        p = (-ds.pv[0,0,:,:]).plot.contourf(levels=PV_LEVELS[theta_level],
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.SouthPolarStereo()},
                                        cmap=custom_cmap)
        p1 = (-ds.pv[0,0,:,:]).plot.contour(levels=PV_LEVELS[theta_level],
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.SouthPolarStereo()},
                                        colors="k", linestyles="-", linewidths=0.5)
        p.axes.coastlines(color="w", linewidth=1.5)
        obj = p.axes.gridlines(linestyle="--", linewidth=0.5)
        p.axes.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
        title = f"MIMOSA Potential Vorticity\nSouth Pole (PV=-PV)\ntime = {format_datetime(str(ds.time.values[0]).split('.')[0],'%Y-%m-%dT%H:%M:%S','%d/%m/%Y %H:%M')}, theta = {theta_level} [K]"
        p.axes.set_title(title, fontsize=15)
        p.colorbar.ax.yaxis.label.set_fontsize(15)
        p.colorbar.set_ticks(PV_LEVELS[theta_level])
        p.colorbar.ax.tick_params(labelsize=10)
        fig.savefig(image_filepath, dpi=200, bbox_inches='tight')
        ds.close()
        plt.close(fig)

def create_figure_t(netcdf_file: str, images_folder: str) -> None:
    gaussian_kernel_sigma = 4
    ds = xr.open_dataset(netcdf_file)
    ds.t[0,0,:,:] = gaussian_filter(ds.t[0,0,:,:], sigma=gaussian_kernel_sigma)
    theta_level = ds.isentropic_level_Kelvin
    # ------------------------------------------------------------------------------------
    # GLOBAL PLOT
    # ------------------------------------------------------------------------------------
    image_filepath = f"{images_folder}/{os.path.basename(netcdf_file)[:-3]}.png"
    if not os.path.exists(image_filepath):
        fig = plt.figure()
        p = (ds.t[0,0,:,:]).plot.contourf(levels=T_LEVELS,
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.PlateCarree()},
                                        cmap=custom_cmap,
                                        add_colorbar=True,
                                        cbar_kwargs={"fraction": 0.023})
        p1 = (ds.t[0,0,:,:]).plot.contour(levels=T_LEVELS,
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.PlateCarree()},
                                        colors="k", linestyles="-", linewidths=0.5)
        p.axes.coastlines(color="w", linewidth=1.5)
        obj = p.axes.gridlines(linestyle="--", linewidth=0.5, draw_labels=True)
        obj.top_labels   = False
        obj.right_labels = False
        obj.xlabel_style = {"size":8}
        obj.ylabel_style = {"size":8}
        p.axes.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        title = f"MIMOSA Temperature [K]\ntime = {format_datetime(str(ds.time.values[0]).split('.')[0],'%Y-%m-%dT%H:%M:%S','%d/%m/%Y %H:%M')}, theta = {theta_level} [K]"
        p.axes.set_title(title, fontsize=10)
        p.colorbar.ax.yaxis.label.set_fontsize(10)
        p.colorbar.set_ticks(T_LEVELS)
        p.colorbar.ax.tick_params(labelsize=8)
        fig.savefig(image_filepath, dpi=200, bbox_inches='tight')
        plt.close(fig)
    # ------------------------------------------------------------------------------------
    # NORTH POLE PLOT
    # ------------------------------------------------------------------------------------
    image_filepath = f"{images_folder}/{os.path.basename(netcdf_file)[:-3].replace('tg','tn')}.png"
    if not os.path.exists(image_filepath):
        fig = plt.figure()
        p = (ds.t[0,0,:,:]).plot.contourf(levels=T_LEVELS,
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.NorthPolarStereo()},
                                        cmap=custom_cmap)
        p1 = (ds.t[0,0,:,:]).plot.contour(levels=T_LEVELS,
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.NorthPolarStereo()},
                                        colors="k", linestyles="-", linewidths=0.5)
        p.axes.coastlines(color="w", linewidth=1.5)
        obj = p.axes.gridlines(linestyle="--", linewidth=0.5)
        p.axes.set_extent([-180, 180, 30, 90], ccrs.PlateCarree())
        title = f"MIMOSA Temperature [K]\nNorth Pole\ntime = {format_datetime(str(ds.time.values[0]).split('.')[0],'%Y-%m-%dT%H:%M:%S','%d/%m/%Y %H:%M')}, theta = {theta_level} [K]"
        p.axes.set_title(title, fontsize=15)
        p.colorbar.ax.yaxis.label.set_fontsize(15)
        p.colorbar.set_ticks(T_LEVELS)
        p.colorbar.ax.tick_params(labelsize=10)
        fig.savefig(image_filepath, dpi=200, bbox_inches='tight')
        plt.close(fig)
    # ------------------------------------------------------------------------------------
    # SOUTH POLE PLOT
    # ------------------------------------------------------------------------------------
    image_filepath = f"{images_folder}/{os.path.basename(netcdf_file)[:-3].replace('tg','ts')}.png"
    if not os.path.exists(image_filepath):
        fig = plt.figure()
        p = (ds.t[0,0,:,:]).plot.contourf(levels=T_LEVELS,
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.SouthPolarStereo()},
                                        cmap=custom_cmap)
        p1 = (ds.t[0,0,:,:]).plot.contour(levels=T_LEVELS,
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.SouthPolarStereo()},
                                        colors="k", linestyles="-", linewidths=0.5)
        p.axes.coastlines(color="w", linewidth=1.5)
        obj = p.axes.gridlines(linestyle="--", linewidth=0.5)
        p.axes.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
        title = f"MIMOSA Temperature [K]\nSouth Pole\ntime = {format_datetime(str(ds.time.values[0]).split('.')[0],'%Y-%m-%dT%H:%M:%S','%d/%m/%Y %H:%M')}, theta = {theta_level} [K]"
        p.axes.set_title(title, fontsize=15)
        p.colorbar.ax.yaxis.label.set_fontsize(15)
        p.colorbar.set_ticks(T_LEVELS)
        p.colorbar.ax.tick_params(labelsize=10)
        fig.savefig(image_filepath, dpi=200, bbox_inches='tight')
        plt.close(fig)
        ds.close()

def create_figure_wind(netcdf_file: str, images_folder: str) -> None:
    u_file = netcdf_file[0]
    v_file = netcdf_file[1]
    ds_u = xr.open_dataset(u_file)
    ds_v = xr.open_dataset(v_file)
    ds_wind = xr.merge([ds_u.sum(dim=["time","theta"]), ds_v.sum(dim=["time","theta"])])
    ds_wind["magnitude"] = np.sqrt(ds_wind.u*ds_wind.u + ds_wind.v*ds_wind.v).T
    theta_level = ds_u.isentropic_level_Kelvin
    if theta_level not in PV_LEVELS.keys():
        print(np.array(PV_LEVELS.keys()))
        theta_level = find_nearest(list(PV_LEVELS.keys()), theta_level)
    N=5
    subset_ds = ds_wind.sel(lon=slice(None, None, N), lat=slice(None, None, N))
    # ------------------------------------------------------------------------------------
    # GLOBAL PLOT
    # ------------------------------------------------------------------------------------
    image_filepath = f"{images_folder}/w{os.path.basename(u_file)[1:-3]}.png"
    if not os.path.exists(image_filepath):
        fig = plt.figure()
        p = ds_wind.magnitude.plot.imshow(x="lon",y="lat",
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.PlateCarree()},
                                        cmap="jet",
                                        add_colorbar=True,
                                        cbar_kwargs={"fraction": 0.023},
                                        vmin=0, vmax=100, extend="max")
        p1 = subset_ds.plot.quiver(x="lon",y="lat",u="u",v="v",
                                transform=ccrs.PlateCarree(),
                                subplot_kws={"projection": ccrs.PlateCarree()},
                                color="k",
                                scale=2500,
                                width=0.001,
                                headwidth=7)
        p.axes.coastlines(color="w", linewidth=1.5)
        obj = p.axes.gridlines(linestyle="--", linewidth=0.5, draw_labels=True)
        obj.top_labels   = False
        obj.right_labels = False
        obj.xlabel_style = {"size":8}
        obj.ylabel_style = {"size":8}
        p.axes.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        title = f"MIMOSA Wind [m.s-]\ntime = {format_datetime(str(ds_u.time.values[0]).split('.')[0],'%Y-%m-%dT%H:%M:%S','%d/%m/%Y %H:%M')}, theta = {ds_u.isentropic_level_Kelvin} [K]"
        p.axes.set_title(title, fontsize=10)
        p.colorbar.ax.yaxis.label.set_fontsize(10)
        p.colorbar.set_label("Wind [m.s-1]")
        p.colorbar.ax.tick_params(labelsize=8)
        fig.savefig(image_filepath, dpi=200, bbox_inches='tight')
        plt.close(fig)
    # ------------------------------------------------------------------------------------
    # NORTH POLE PLOT
    # ------------------------------------------------------------------------------------
    image_filepath = f"{images_folder}/w{os.path.basename(u_file)[1:-3].replace('g','n')}.png"
    if not os.path.exists(image_filepath):
        fig = plt.figure()
        p = ds_wind.magnitude.plot.imshow(x="lon",y="lat",
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.NorthPolarStereo()},
                                        cmap="jet",
                                        add_colorbar=True,
                                        cbar_kwargs={"fraction": 0.023},
                                        vmin=0, vmax=100, extend="max")
        p1 = subset_ds.plot.quiver(x="lon",y="lat",u="u",v="v",
                                transform=ccrs.PlateCarree(),
                                subplot_kws={"projection": ccrs.NorthPolarStereo()},
                                color="w",
                                scale=1500,
                                width=0.003,
                                headwidth=7,
                                headlength=7)
        p.axes.coastlines(color="w", linewidth=1.5)
        obj = p.axes.gridlines(linestyle="--", linewidth=0.5)
        p.axes.set_extent([-180, 180, 30, 90], ccrs.PlateCarree())
        title = f"MIMOSA Wind [m.s-]\nNorth Pole\ntime = {format_datetime(str(ds_u.time.values[0]).split('.')[0],'%Y-%m-%dT%H:%M:%S','%d/%m/%Y %H:%M')}, theta = {ds_u.isentropic_level_Kelvin} [K]"
        p.axes.set_title(title, fontsize=10)
        p.colorbar.ax.yaxis.label.set_fontsize(10)
        p.colorbar.set_label("Wind [m.s-1]")
        p.colorbar.ax.tick_params(labelsize=8)
        fig.savefig(image_filepath, dpi=200, bbox_inches='tight')
        plt.close(fig)
    # ------------------------------------------------------------------------------------
    # NORTH POLE PLOT
    # ------------------------------------------------------------------------------------
    image_filepath = f"{images_folder}/w{os.path.basename(u_file)[1:-3].replace('g','s')}.png"
    if not os.path.exists(image_filepath):
        fig = plt.figure()
        p = ds_wind.magnitude.plot.imshow(x="lon",y="lat",
                                        transform=ccrs.PlateCarree(),
                                        subplot_kws={"projection": ccrs.SouthPolarStereo()},
                                        cmap="jet",
                                        add_colorbar=True,
                                        cbar_kwargs={"fraction": 0.023},
                                        vmin=0, vmax=100, extend="max")
        p1 = subset_ds.plot.quiver(x="lon",y="lat",u="u",v="v",
                                transform=ccrs.PlateCarree(),
                                subplot_kws={"projection": ccrs.SouthPolarStereo()},
                                color="w",
                                scale=1500,
                                width=0.003,
                                headwidth=7,
                                headlength=7)
        p.axes.coastlines(color="w", linewidth=1.5)
        obj = p.axes.gridlines(linestyle="--", linewidth=0.5)
        p.axes.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
        title = f"MIMOSA Wind [m.s-]\nSouth Pole\ntime = {format_datetime(str(ds_u.time.values[0]).split('.')[0],'%Y-%m-%dT%H:%M:%S','%d/%m/%Y %H:%M')}, theta = {ds_u.isentropic_level_Kelvin} [K]"
        p.axes.set_title(title, fontsize=10)
        p.colorbar.ax.yaxis.label.set_fontsize(10)
        p.colorbar.set_label("Wind [m.s-1]")
        p.colorbar.ax.tick_params(labelsize=8)
        fig.savefig(image_filepath, dpi=200, bbox_inches='tight')
        plt.close(fig)

def combine_netcdfs(data_dir: str) -> None:
    for var in ["pvg", "tg", "ug", "vg"]:
        files = glob.glob(f"{data_dir}/{var}*.nc")
        files.sort()
        theta_levels_all = list(set([os.path.basename(elem).split(".")[1] for elem in files]))
        for theta in theta_levels_all:
            final_filepath = f"{data_dir}/{var}.{theta}.nc"
            if os.path.exists(final_filepath):
                if os.path.getsize(final_filepath)==0:
                    os.remove(f"{data_dir}/{var}.{theta}.nc")
                else:
                    try:
                        trash = xr.open_dataset(final_filepath)
                        if trash.size==0:
                            os.remove(f"{data_dir}/{var}.{theta}.nc")
                        else:
                            continue
                    except:
                        os.remove(f"{data_dir}/{var}.{theta}.nc")
            else:
                current_theta_files = [elem for elem in files if f"{theta}.nc" in elem]
                datasets = [xr.open_dataset(file) for file in current_theta_files]
                final_dataset = xr.concat(datasets, dim="time")
                final_dataset = final_dataset.sortby("time")
                try:
                    final_dataset.to_netcdf(f"{data_dir}/{var}.{theta}.nc")
                except:
                    os.remove(f"{data_dir}/{var}.{theta}.nc")
                    final_dataset.to_netcdf(f"{data_dir}/{var}.{theta}.nc")
                datasets = [elem.close() for elem in datasets]

def main(output_folder: str, images_folder: str) -> None:
    list_of_fortran_files = glob.glob(f"{output_folder}/*g*.????")
    list_of_fortran_files.sort()
    LOGGER.info("!======================================================================!")
    LOGGER.info("!               __  __ _____ __  __  ____   _____                      !")
    LOGGER.info("!              |  \/  |_   _|  \/  |/ __ \ / ____|  /\                 !")
    LOGGER.info("!              | \  / | | | | \  / | |  | | (___   /  \                !")
    LOGGER.info("!              | |\/| | | | | |\/| | |  | |\___ \ / /\ \               !")
    LOGGER.info("!              | |  | |_| |_| |  | | |__| |____) / ____ \              !")
    LOGGER.info("!              |_|  |_|_____|_|  |_|\____/|_____/_/    \_\             !")
    LOGGER.info("!======================================================================!")
    ii = 1
    for file in list_of_fortran_files:
        N = len(list_of_fortran_files)
        LOGGER.info(f"Processing file {os.path.basename(file)} ({ii}/{N})...")
        nc_file = create_netcdf_file(file, output_folder)
        if "pv" in os.path.basename(nc_file):
            create_figure_pv(nc_file, images_folder)
        elif "t" in os.path.basename(nc_file):
            create_figure_t(nc_file, images_folder)
        elif "u" in os.path.basename(nc_file):
            u_filepath = nc_file
            v_filepath = nc_file.replace("ug","vg")
            if not os.path.exists(v_filepath):
                LOGGER.warning(f"Cannot plot the wind for {os.path.basename(nc_file)}, the v component is missing")
                continue
            else:
                create_figure_wind([u_filepath, v_filepath], images_folder)
                other_file = file.replace("ug","vg")
                list_of_fortran_files.remove(file)
                list_of_fortran_files.remove(other_file)
        elif "v" in os.path.basename(nc_file):
            v_filepath = nc_file
            u_filepath = nc_file.replace("vg","ug")
            if not os.path.exists(u_filepath):
                LOGGER.warning(f"Cannot plot the wind for {os.path.basename(nc_file)}, the u component is missing")
                continue
            else:
                create_figure_wind([u_filepath, v_filepath], images_folder)
                other_file = file.replace("vg","ug")
                list_of_fortran_files.remove(file)
                list_of_fortran_files.remove(other_file)
        else:
            LOGGER.warning(f"Can't plot {os.path.basename(nc_file)}, unrecognized variable")
        ii+=1
    combine_netcdfs(output_folder)

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

    global LOGGER
    LOGGER = start_log()
    main(OUTPUT_DIR, IMAGES_DIR)
