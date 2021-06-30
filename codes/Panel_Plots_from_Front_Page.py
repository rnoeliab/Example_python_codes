#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 18:12:21 2020

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

# Create the figure that will have 3 subplots
fig = plt.figure(figsize=(12,9))
ax_ctt = fig.add_subplot(1,2,1)
ax_wspd = fig.add_subplot(2,2,2)
ax_dbz = fig.add_subplot(2,2,4)

# Get the basemap object
bm = get_basemap(slp)

# Convert the lat/lon points in to x/y points in the projection space
x, y = bm(to_np(lons), to_np(lats))

# Make the pressure contours
contour_levels = [1016, 1017, 1018, 1019, 1020, 1021]
c1 = bm.contour(x, y, to_np(smooth_slp), levels=contour_levels,
                colors="white", zorder=3, linewidths=1.0, ax=ax_ctt)

# Create the filled cloud top temperature contours
contour_levels = [9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0]
ctt_contours = bm.contourf(x, y, to_np(ctt), contour_levels,
                           cmap=get_cmap("Greys"), zorder=2, ax=ax_ctt)

point_x, point_y = bm([start_point.lon, end_point.lon],
                      [start_point.lat, end_point.lat])
bm.plot([point_x[0], point_x[1]], [point_y[0], point_y[1]], color="yellow",
        marker="o", zorder=3, ax=ax_ctt)

# Create the color bar for cloud top temperature
cb_ctt = fig.colorbar(ctt_contours, ax=ax_ctt, shrink=.60)
cb_ctt.ax.tick_params(labelsize=5)

# Draw the oceans, land, and states
bm.drawcoastlines(linewidth=0.25, ax=ax_ctt)
bm.drawstates(linewidth=0.25, ax=ax_ctt)
bm.drawcountries(linewidth=0.25, ax=ax_ctt)
bm.fillcontinents(color=np.array([ 0.9375 , 0.9375 , 0.859375]),
                                 ax=ax_ctt,
                                 lake_color=np.array([0.59375 ,
                                                      0.71484375,
                                                      0.8828125 ]))
bm.drawmapboundary(fill_color='None', ax=ax_ctt)

# Draw Parallels
parallels = np.arange(np.amin(lats), -21., 1.0)
bm.drawparallels(parallels, ax=ax_ctt, color="white")

merids = np.arange(-49.0, -44.0, 1.0)
bm.drawmeridians(merids, ax=ax_ctt, color="white")

# Crop the image to the hurricane region
x_start, y_start = bm(-49.0, np.amin(lats))
x_end, y_end = bm(-44.0, -21.0)

ax_ctt.set_xlim([x_start, x_end])
ax_ctt.set_ylim([y_start, y_end])

# Make the contour plot for wspd
wspd_contours = ax_wspd.contourf(to_np(wspd_cross), cmap=get_cmap("jet"))
# Add the color bar
cb_wspd = fig.colorbar(wspd_contours, ax=ax_wspd)
cb_wspd.ax.tick_params(labelsize=5)

# Make the contour plot for dbz
levels = [-30 + 2.5*n for n in range(7)]
dbz_contours = ax_dbz.contourf(to_np(dbz_cross), levels=levels,
                               cmap=get_cmap("jet"))
cb_dbz = fig.colorbar(dbz_contours, ax=ax_dbz)
cb_dbz.ax.tick_params(labelsize=5)

# Set the x-ticks to use latitude and longitude labels.
coord_pairs = to_np(dbz_cross.coords["xy_loc"])
x_ticks = np.arange(coord_pairs.shape[0])
x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]
ax_wspd.set_xticks(x_ticks[::7])
ax_wspd.set_xticklabels([], rotation=45)
ax_dbz.set_xticks(x_ticks[::7])
ax_dbz.set_xticklabels(x_labels[::7], rotation=45, fontsize=4)

# Set the y-ticks to be height.
vert_vals = to_np(dbz_cross.coords["vertical"])
v_ticks = np.arange(vert_vals.shape[0])
ax_wspd.set_yticks(v_ticks[::7])
ax_wspd.set_yticklabels(vert_vals[::7], fontsize=4)
ax_dbz.set_yticks(v_ticks[::7])
ax_dbz.set_yticklabels(vert_vals[::7], fontsize=4)

# Set the x-axis and  y-axis labels
ax_dbz.set_xlabel("Latitude, Longitude", fontsize=5)
ax_wspd.set_ylabel("Height (m)", fontsize=5)
ax_dbz.set_ylabel("Height (m)", fontsize=5)

# Add titles
ax_ctt.set_title("Cloud Top Temperature (degC)", {"fontsize" : 7})
ax_wspd.set_title("Cross-Section of Wind Speed (kt)", {"fontsize" : 7})
ax_dbz.set_title("Cross-Section of Reflectivity (dBZ)", {"fontsize" : 7})

plt.savefig(output+"panel_plots_from_frontl_page.tif")
plt.show()
