#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 12:09:23 2019

@author: noelia
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from netCDF4 import Dataset
from glob import glob
from wrf import to_np, getvar, CoordPair, vertcross

wrf_file = sorted(glob("/data/noelia/modelo/wrf_chem_8/wrfout_d01_*"))  ## archivos wrfout
output = "/data/noelia/modelo/examples_plot/"

# Open the NetCDF file
ncfile = Dataset(wrf_file[10])

# Extract the model height and wind speed
z = getvar(ncfile, "z")
wspd =  getvar(ncfile, "uvmet_wspd_wdir", units="kt")[0,:]

# Create the start point and end point for the cross section
start_point = CoordPair(lat=-24.5, lon=-48.5)
end_point = CoordPair(lat=-23., lon=-44.8)


# Compute the vertical cross-section interpolation.  Also, include the
# lat/lon points along the cross-section.
wspd_cross = vertcross(wspd, z, wrfin=ncfile, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)

# Create the figure
fig = plt.figure(figsize=(12,10))
ax = plt.axes()

# Make the contour plot
wspd_contours = ax.contourf(to_np(wspd_cross[0:50,::]), cmap=get_cmap("jet"))

# Add the color bar
cbar=fig.colorbar(wspd_contours, ax=ax,shrink=1.0)
cbar.ax.tick_params(labelsize=15)

# Set the x-ticks to use latitude and longitude labels.
coord_pairs = to_np(wspd_cross.coords["xy_loc"])

x_ticks = np.arange(coord_pairs.shape[0])
x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}")
            for pair in to_np(coord_pairs)]

ax.set_xticks(x_ticks[::7])
ax.set_xticklabels(x_labels[::7], rotation=45, fontsize=10)

# Set the y-ticks to be height.
vert_vals = to_np(wspd_cross.coords["vertical"][0:50])
v_ticks = np.arange(vert_vals.shape[0])
ax.set_yticks(v_ticks[::7])
ax.set_yticklabels(vert_vals[::7], fontsize=12)

# Set the x-axis and  y-axis labels
ax.set_xlabel("Latitude, Longitude", fontsize=15)
ax.set_ylabel("Height (m)", fontsize=15)

plt.title("Vertical Cross Section of Wind Speed (kt)", fontsize=18)
plt.savefig(output+"vertical_cross_section.tif")
plt.show()




