#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 31 21:32:04 2020

@author: noelia
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from glob import glob
from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim)

wrf_file = sorted(glob("/data/noelia/modelo/wrf_chem_8/wrfout_d01_*"))  ## archivos wrfout
output = "/data/noelia/modelo/examples_plot/"

# Open the NetCDF file
ncfile = Dataset(wrf_file[10])

# Extract the pressure, geopotential height, and wind variables
p = getvar(ncfile, "pressure")
z = getvar(ncfile, "z", units="dm")
ua = getvar(ncfile, "ua", units="kt")
va = getvar(ncfile, "va", units="kt")
wspd = getvar(ncfile, "wspd_wdir", units="kts")[0,:]

# Interpolate geopotential height, u, and v winds to 500 hPa
ht_500 = interplevel(z, p, 500)
u_500 = interplevel(ua, p, 500)
v_500 = interplevel(va, p, 500)
wspd_500 = interplevel(wspd, p, 500)

# Get the lat/lon coordinates
lats, lons = latlon_coords(ht_500)

# Get the map projection information
cart_proj = get_cartopy(ht_500)

# Create the figure
fig = plt.figure(figsize=(12,9))
ax = plt.axes(projection=cart_proj)

# Download and add the states and coastlines
states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_shp")
ax.add_feature(states, linewidth=1.5, edgecolor="black")
ax.coastlines('50m', linewidth=0.8)

# Add the 500 hPa geopotential height contours
levels = np.arange(590., 592., 0.23)
contours = plt.contour(to_np(lons), to_np(lats), to_np(ht_500),
                       levels=levels, colors="red",
                       transform=crs.PlateCarree())
plt.clabel(contours, inline=1, fontsize=20, fmt="%i")

# Add the wind speed contours
levels1 = [0, 4, 8, 12, 16, 20, 24]
wspd_contours = plt.contourf(to_np(lons), to_np(lats), to_np(wspd_500),
                             levels=levels1,cmap=get_cmap("rainbow"),
                             transform=crs.PlateCarree())

plt.colorbar(wspd_contours, ax=ax, orientation="vertical", pad=.02)

# Add the 500 hPa wind barbs, only plotting every 125th data point.
plt.barbs(to_np(lons[::8,::8]), to_np(lats[::8,::8]),
          to_np(u_500[::8, ::8]), to_np(v_500[::8, ::8]),
          transform=crs.PlateCarree(), length=7)

# Set the map bounds
ax.set_xlim(cartopy_xlim(ht_500))
ax.set_ylim(cartopy_ylim(ht_500))

ax.gridlines(color="black", linestyle="dotted")
plt.title("500 MB Height (dm), Wind Speed (kt), Barbs (kt)",fontsize=18)
plt.savefig(output+"horizontal_interp_to_pressure.tif")
plt.show()