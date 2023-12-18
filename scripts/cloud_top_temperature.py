#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 08:48:30 2020

@author: noelia
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from netCDF4 import Dataset
from glob import glob
from wrf import (getvar, to_np, vertcross, smooth2d, CoordPair,
                 get_basemap, latlon_coords)

wrf_file = sorted(glob("/data/noelia/modelo/wrf_chem_8/wrfout_d01_*"))  ## archivos wrfout
output = "/data/noelia/modelo/examples_plot/"

# Open the NetCDF file
ncfile = Dataset(wrf_file[10])

# Get the WRF variables
slp = getvar(ncfile, "slp")
smooth_slp = smooth2d(slp, 3)
ctt = getvar(ncfile, "ctt")
z = getvar(ncfile, "z")
dbz = getvar(ncfile, "dbz")
Z = 10**(dbz/10.)
wspd =  getvar(ncfile, "wspd_wdir", units="kt")[0,:]

# Set the start point and end point for the cross section
start_point = CoordPair(lat=-24.5, lon=-48.5)
end_point = CoordPair(lat=-23., lon=-44.8)

# Compute the vertical cross-section interpolation.  Also, include the
# lat/lon points along the cross-section in the metadata by setting latlon
# to True.
z_cross = vertcross(Z, z, wrfin=ncfile, start_point=start_point,
                    end_point=end_point, latlon=True, meta=True)
wspd_cross = vertcross(wspd, z, wrfin=ncfile, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)
dbz_cross = 10.0 * np.log10(z_cross)

# Get the latitude and longitude points
lats, lons = latlon_coords(slp)

# Create the figure 
fig = plt.figure(figsize=(12,7))
ax = plt.axes()

# Get the basemap object
bm = get_basemap(slp)

# Convert the lat/lon points in to x/y points in the projection space
x, y = bm(to_np(lons), to_np(lats))

# Make the pressure contours
contour_levels = np.arange(to_np(smooth_slp).min(),to_np(smooth_slp).max(),0.7)

c1 = bm.contour(x, y, to_np(smooth_slp), levels=contour_levels,
                colors="white", zorder=3, linewidths=1.5, ax=ax)

# Create the filled cloud top temperature contours
contour_levels =np.arange(np.amin(ctt),np.amax(ctt)+1.5,1.5)
ctt_contours = bm.contourf(x, y, to_np(ctt), contour_levels,
                           cmap=get_cmap("binary"), zorder=2, ax=ax)

point_x, point_y = bm([start_point.lon, end_point.lon],
                      [start_point.lat, end_point.lat])
bm.plot([point_x[0], point_x[1]], [point_y[0], point_y[1]], color="red",
        marker="o", zorder=3, ax=ax)

# Create the color bar for cloud top temperature
cb_ctt = fig.colorbar(ctt_contours, ax=ax, shrink=.90)
cb_ctt.ax.tick_params(labelsize=10)

# Draw the oceans, land, and states
bm.drawcoastlines(linewidth=1.5, ax=ax)
bm.drawstates(linewidth=2.0, ax=ax)
bm.fillcontinents(color=np.array([ 0.9375 , 0.9375 , 0.859375]),ax=ax,
                                 lake_color=np.array([0.59375 ,0.71484375,0.8828125 ]))
bm.drawmapboundary(fill_color='None', ax=ax)

# Draw Parallels
parallels = np.arange(np.amin(lats), np.amax(lats), 0.9)
bm.drawparallels(parallels, ax=ax, color="blue", labels=[1, 0, 0, 0],fontsize=11)

merids = np.arange(np.amin(lons),np.amax(lons), 0.9)
bm.drawmeridians(merids, ax=ax, color="blue", labels=[0, 0, 0, 1],fontsize=11)
bm.readshapefile('/media/noelia/TOSHIBA EXT/doctorado/usp/shapefile/RM_Sao_Paulo/transformed','sp')

# Add titles
ax.set_title("Cloud Top Temperature (degC)", {"fontsize" : 15})

plt.savefig(output+"cloud_top_temperature.tif")
plt.show()