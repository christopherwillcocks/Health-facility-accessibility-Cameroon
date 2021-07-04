# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 10:54:22 2021

@author: Christopher Willcocks, christopher.willcocks@hotmail.com
Master Student in Environemental engineering in EPFL
"""

## MASTER THESIS ## Modelling health accessibility in urban Sub-Saharan Africa

# Aim of this code is to share the plots codes used in my final report


from matplotlib import cm
from matplotlib_scalebar.scalebar import ScaleBar
import pyproj
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt



mask=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\mask.shp')

Nearest_facility=gpd.read_file(r'C:\Users\Te Anau\Documents\Master thesis Christopher\Results\Nearest_facility.shp')
Nearest_pop_function=gpd.read_file(r'C:\Users\Te Anau\Documents\Master thesis Christopher\Results\Nearest_pop_function.shp')

Walking_scores=gpd.read_file(r'C:\Users\Te Anau\Documents\Master thesis Christopher\Results\Walking_scores.shp')
Accessibility_scores=gpd.read_file(r'C:\Users\Te Anau\Documents\Master thesis Christopher\Results\Accessibility_scores.shp')

# Frame of the coordinates
central_coordinate=[3.88, 11.51]
dlat=0.22
dlong=0.15



# %%

##### Plot of distance to nearest facility
# Need to adapt Nearest facility so that different healthcare facility are ploted


viridis_modified = cm.get_cmap('viridis', 256) # Colors sequence for colormap
newcmp = ListedColormap(viridis_modified(np.linspace(0.1, 0.92, 256)))

fig, axs = plt.subplots(2, 3, figsize=(16, 14))
Nearest_facility.plot(ax=axs[0,0], figsize=(6, 6), alpha=1, edgecolor='none', 
              column=Nearest_facility['Time to cl'].fillna(11000), cmap=newcmp, legend=False, markersize=2,  # time to cl is for set 1, need to be adapted 
             legend_kwds={'label': "Distance [km]", 'orientation': "vertical"},marker='o', vmax=10000)


mask.boundary.plot(ax=axs[0,0], color='k', alpha=0.4, linewidth = 0.5)

Nearest_facility.plot(ax=axs[0,1], figsize=(6, 6), alpha=1, edgecolor='none', 
              column=Nearest_facility['Time to cl'].fillna(11000), cmap=newcmp, legend=False, markersize=2, 
             legend_kwds={'label': "Distance [km]", 'orientation': "vertical"},marker='o', vmax=10000)
mask.boundary.plot(ax=axs[0,1], color='k', alpha=0.4, linewidth = 0.5)


axs[0,2].plot(Nearest_pop_function['abs popula'], color='blue', label='Basic structures')
axs[0,2].plot(Nearest_pop_function['prop popul'], color='red', label='Hospitals')
axs[0,2].plot(list([5,5]),list([0,100]),linewidth = 0.5, color='k' , label='5 km thresold')
leg = axs[0,2].legend(loc='upper right', markerscale=4)


Nearest_facility.plot(ax=axs[1,0], figsize=(6, 6), alpha=1, edgecolor='none', 
              column=Nearest_facility['Time to cl'].fillna(11000), cmap=newcmp, legend=False, markersize=2, 
             legend_kwds={'label': "Distance [km]", 'orientation': "vertical"},marker='o', vmax=10000)
mask.boundary.plot(ax=axs[1,0], color='k', alpha=0.4, linewidth = 0.5)

Nearest_facility.plot(ax=axs[1,1], figsize=(6, 6), alpha=1, edgecolor='none', 
              column=Nearest_facility['Time to cl'].fillna(11000), cmap=newcmp, legend=False, markersize=2, 
             legend_kwds={'label': "Distance [km]", 'orientation': "vertical"},marker='o', vmax=10000)
mask.boundary.plot(ax=axs[1,1], color='k', alpha=0.4, linewidth = 0.5)



axs[1,2].plot(Nearest_pop_function['abs popula'], color='blue', label='Basic structures')
axs[1,2].plot(Nearest_pop_function['prop popul'], color='red', label='Hospitals')
axs[1,2].plot(list([5,5]),list([0,100]),linewidth = 0.5, color='k' , label='5 km thresold')
leg = axs[1,2].legend(loc='upper right', markerscale=4)



cbar_ax = fig.add_axes([0.045, 0.2, 0.02, 0.6])
norm = cm.colors.Normalize(vmin=0, vmax=10)
fig.colorbar(cm.ScalarMappable(cmap=newcmp, norm=norm), cax=cbar_ax, label='Distance [m]', aspect=30, shrink=0.5, extend='max')


axs[0,0].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
axs[0,0].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)
axs[0,1].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
axs[0,1].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)


axs[1,0].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
axs[1,0].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)
axs[1,1].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
axs[1,1].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)

axs[0,2].set_xlim(0, 10)
axs[1,2].set_xlim(0, 10)
axs[0,2].set_ylim(0,50)
axs[1,2].set_ylim(0,50)

axs[0,0].set_ylabel('latitude')
axs[0,0].set_xlabel('longitude')
axs[0,1].set_xlabel('longitude')
#axs[0,2].set_xlabel('longitude')

axs[1,0].set_ylabel('latitude')
axs[1,0].set_xlabel('longitude')
axs[1,1].set_xlabel('longitude')
#axs[1,2].set_xlabel('longitude')

axs[0,2].set_ylabel('Population percentage %')
axs[1,2].set_ylabel('Population percentage %')
axs[0,2].set_xlabel('Distance [km]')
axs[1,2].set_xlabel('Distance [km]')



axs[0,0].set(title='(A) Distance to neartest basic structure \n mix Dataset')
axs[0,1].set(title='(B) Distance to neartest hospital \n mix Dataset')
axs[0,2].set(title='(C) Inaccessible population as a function of distance \n mix Dataset')
axs[1,0].set(title='(D) Distance to neartest basic structure \n gpd Dataset')
axs[1,1].set(title='(E) Distance to neartest hospital \n gpd Dataset')
axs[1,2].set(title='(F) Inaccessible population as a function of distance \n gpd Dataset')




g = pyproj.Geod(ellps='WGS84')
length=g.line_length([central_coordinate[0],central_coordinate[0]], [central_coordinate[1]-dlong, central_coordinate[1]+dlong])
scalebar = ScaleBar(length*4, "m", length_fraction=0.3, location='lower left',
                    border_pad=0.8, box_color='w', box_alpha=0.7)
axs[1,1].add_artist(scalebar)

# %%

# Plots  of accessibility scores

# Corlormap 
viridis_modified = cm.get_cmap('Spectral', 256)
newcmp = ListedColormap(viridis_modified(np.linspace(0.0, 0.49, 256)))

# Subplots
fig, ax = plt.subplots(1, 2, figsize=(8, 8), sharey=True, constrained_layout=True)
mask.plot(ax=ax[0],figsize=(8, 6), alpha=0.12, color='k')
Walking_scores.plot(ax=ax[0], alpha=0.85, edgecolor='none', 
              column=exp['ratio m1'], cmap=newcmp, legend=True, markersize=1,marker='s', scheme='quantiles', k=5, 
              legend_kwds={'fmt':'{:.5f}', 'loc':'upper left', 'title':'Accessibility score', 'title_fontsize':'medium',
                           'fontsize':'small','borderpad':0.2, 'borderaxespad':1})
mask.boundary.plot(ax=ax[0], linewidth = 0.5, alpha=0.4, color='k')

mask.plot(ax=ax[1],figsize=(8, 6), alpha=0.12, color='k')
Walking_scores.plot(ax=ax[1], alpha=0.85, edgecolor='none', 
              column=exp['ratio m1'], cmap=newcmp, legend=False, markersize=1,marker='s', scheme='quantiles', k=5, 
              legend_kwds={'fmt':'{:.5f}', 'loc':'upper left', 'title':'Accessibility score',
                           'fontsize':'small','borderpad':0.5, 'borderaxespad':1})
mask.boundary.plot(ax=ax[1], linewidth = 0.5, alpha=0.4, color='k')


# Axis labeling
ax[0].set(title='PeanutButter Accessibility \n dataset 1')
ax[1].set(title='PeanutButter Accessibility \n datset 2')
ax[0].set_ylabel('latitude')
ax[0].set_xlabel('longitude')
ax[1].set_xlabel('longitude')

# Including North arrow and comment
x, y, arrow_length = 0.44, 0.967, 0.06
ax[1].annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
            arrowprops=dict(facecolor='black', width=4, headwidth=10),
            ha='center', va='center', fontsize=10,
            xycoords=ax[1].transAxes)

x, y, arrow_length = 0.489, 0.837, 0.1
ax[0].annotate('Low', xy=(x, y), xytext=(x, y+arrow_length),
            arrowprops=dict(facecolor='black', width=0.013, headwidth=4),
            ha='center', va='center', fontsize=9,
            xycoords=ax[0].transAxes)
x, y = 0.489, 0.825
ax[0].annotate('High', xy=(x, y), xytext=(x, y),
            ha='center', va='center', fontsize=9,
            xycoords=ax[0].transAxes)


ax[0].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
ax[0].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)
ax[1].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
ax[1].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)

# calculate distance between points
g = pyproj.Geod(ellps='WGS84')
length=g.line_length([central_coordinate[0],central_coordinate[0]], [central_coordinate[1]-dlong, central_coordinate[1]+dlong])

scalebar = ScaleBar(length*4, "m", length_fraction=0.35, location='upper left',
                    border_pad=1.7, box_color='w', box_alpha=0.7)
ax[1].add_artist(scalebar)







