#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 18:31:14 2021

@author: noelia
"""

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np

###################### posicion de las estaciones ###########################
sta_lat, sta_lon = ([-23.561,-22.413,-23.482,-22.689],[-46.735,-45.452,-46.500,-45.006])
station = ['Sao_Paulo','Itajuba','SP-EACH','Cachoeira_Paulista']

fig = plt.figure(figsize=(20, 15))
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3.0)
m = Basemap(projection='cyl', resolution='h', llcrnrlat=-30, urcrnrlat=-15,
            llcrnrlon=-53, urcrnrlon=-28)
m.drawmapboundary(fill_color='#7777ff')
m.fillcontinents(color='#ddaa66', lake_color='#7777ff', zorder=0)
m.drawcoastlines(linewidth=1.5)
m.drawstates(linewidth=1.5)    
m.drawparallels(np.arange(-90., 120., 2),labels=[1, 0, 0, 0],fontsize=22)
m.drawmeridians(np.arange(-180., 181., 2),labels=[0, 0, 0, 1],fontsize=22) 
m.readshapefile('../shapefile/RM_Sao_Paulo/transformed','sp')
x, y = m(sta_lon, sta_lat)
m.scatter(x, y, s=150, c='r', alpha=0.8)
axins = zoomed_inset_axes(ax, 2.35, loc=1)
axins.set_xlim(-49, -43)
axins.set_ylim(-25, -21)
# plt.xticks(visible=False)
# plt.yticks(visible=False)
m2 = Basemap(projection='cyl', resolution='h', llcrnrlat=-25, urcrnrlat=-21,
            llcrnrlon=-49, urcrnrlon=-43,ax=axins)
m2.drawmapboundary(fill_color='#7777ff')
m2.fillcontinents(color='#ddaa66', lake_color='#7777ff', zorder=0)
m2.drawcoastlines(linewidth=1.5)
m2.drawstates(linewidth=1.5)    
m2.drawparallels(np.arange(-90., 120., 1))
m2.drawmeridians(np.arange(-180., 181., 1)) 
m2.readshapefile('../shapefile/RM_Sao_Paulo/transformed','sp')
x, y = m(sta_lon, sta_lat)
m2.scatter(x, y, s=150, c='r', alpha=1.0)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", lw=3, ec="r")

#fig.savefig(output+'out_file.png',bbox_inches='tight')
plt.show()
