# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 10:20:28 2021

@author: Te Anau
"""
## MASTER THESIS ## Modelling health accessibility in urban Sub-Saharan Africa

# Aim of this code is share the codes used to to produced the plots in the methodological part of the project



import matplotlib.pyplot as plt
import georasters as gr
from geopandas import gpd
import contextily as ctx
from matplotlib_scalebar.scalebar import ScaleBar
import pyproj
import osmnx as ox



# Reading the data with geopandas
pop_centroid_UN_adj=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\pop_worldpop_UN_adj.shp')
pop_centroid_peanutButter=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\pop_peanutButter.shp')
raster_1_filepath = r'C:\Users\Te Anau\Documents\Data\Population_dataset_worldpop_UN_adj.tif'
grast=gr.from_file(raster_1_filepath)
gdf_worldpop_UN_adj=gr.to_geopandas(grast)
gdf_worldpop_UN_adj=gdf_worldpop_UN_adj.to_crs(epsg=4326)
gdf_worldpop_UN_adj['geometry']=gdf_worldpop_UN_adj.centroid
HStructures_clinic_geosm=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\HStructures_clinic_geosm.shp')
HStructures_hospital_geosm=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\HStructures_hospital_geosm.shp')
mask=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\mask.shp')
mask_transport=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\mask_transport.shp')

# For Yaounde, recangle extent used for the plots
central_coordinate=[3.88, 11.51]
dlat=0.22
dlong=0.15

# calculate distance between points for the scalebar
g = pyproj.Geod(ellps='WGS84')
length=g.line_length([central_coordinate[0],central_coordinate[0]], [central_coordinate[1]-dlong, central_coordinate[1]+dlong])


# %% 

# Plot presenting the urban footprint action

# Subbplots with 3 figures
fig, axs = plt.subplots(1, 3, figsize=(9, 4.2), sharey=True)
gdf_worldpop_UN_adj.plot(ax=axs[0], markersize=0.01, alpha=1, color='#fcc516')
mask.plot(ax=axs[1], markersize=4, alpha=0.3, color='w')
mask.boundary.plot(ax=axs[1], markersize=4, alpha=1, color='w')
pop_worldpop_UN_adj.plot(ax=axs[2], markersize=0.01, alpha=1, color='#fcc516')


axs[0].set_facecolor('#6a338d')
axs[1].set_facecolor('#6a338d')
axs[2].set_facecolor('#6a338d')
axs[0].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
axs[0].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)
axs[1].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
axs[1].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)
axs[2].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
axs[2].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)

axs[0].set_ylabel('latitude')

axs[0].set_xlabel('longitude')
axs[1].set_xlabel('longitude')
axs[2].set_xlabel('longitude')

axs[0].set_title('(A) Original population dataset', fontsize='medium')
axs[1].set_title('(B) Urban footprinter', fontsize='medium')
axs[2].set_title('(C) Urban area population dataset', fontsize='medium')

scalebar = ScaleBar(length*4, "m", length_fraction=0.3, location='lower left',
                    border_pad=0.8, box_color='w', box_alpha=0.7)
axs[2].add_artist(scalebar)

# %% 

# Plot presenting the road network in Yaounde

# Loading netowrk
G=ox.graph.graph_from_polygon(mask_transport['geometry'][0], network_type='walk', simplify=False, retain_all=False, truncate_by_edge=False, clean_periphery=False, custom_filter=None)




fig, ax = ox.plot_graph(G, figsize=(50,32), 
                        node_color='w', node_size=4, 
                        node_alpha=0.8, 
                        node_zorder=2, bbox=(central_coordinate[0]+dlat ,central_coordinate[0]-dlat, central_coordinate[1]+dlong ,central_coordinate[1]-dlong))

scalebar = ScaleBar(length*4, "m", length_fraction=0.4, location='lower left',
                    border_pad=1, box_color='w', box_alpha=0.9)
ax[2].add_artist(scalebar)


fig, ax = ox.plot_graph(G, figsize=(30,22), 
                        node_color=nc, node_size=ns, 
                        node_alpha=0.8, 
                        node_zorder=2)




# %%
 # plot presenting the population datasets
 
fig, axs = plt.subplots(1, 2, figsize=(9, 4.2), sharey=True)
pop_centroid_peanutButter.plot(ax=axs[0], markersize=1, alpha=1, column = 'value', cmap=plt.cm.BuPu_r)
im3=pop_centroid_UN_adj.plot(ax=axs[1], markersize=1, alpha=1, column = 'value', cmap=plt.cm.BuPu_r, legend=True)

#cbar1 = fig.colorbar(axs[0],axs[1], axs[2])

# Label x and y axis
axs[0].set_ylabel('latitude')
axs[0].set_xlabel('longitude')
axs[1].set_xlabel('longitude')

# Title
axs[0].set_title('(A) peanutButter database', fontsize='medium')
axs[1].set_title('(B) Worldpop UN adjusted', fontsize='medium')


axs[0].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
axs[0].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)
axs[1].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
axs[1].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)


scalebar = ScaleBar(length*4, "m", length_fraction=0.4, location='lower left',
                    border_pad=1, box_color='k', box_alpha=0.05)
axs[1].add_artist(scalebar)


# %%
 # plot presenting the health facilities datasets

# Subplots with two figures
fig, axs = plt.subplots(1, 2, figsize=(9, 8), sharey=True)
HStructures_clinic_geosm.plot(ax=axs[0], markersize=0.5, alpha=1, color='r', label='Geosm clinic')
HStructures_hospital_geosm.plot(ax=axs[1], markersize=0.5, alpha=1, color='r', label= 'Geosm hospitals')

mask.boundary.plot(ax=axs[0], markersize=1, alpha=1, color='w')
mask.boundary.plot(ax=axs[1], markersize=1, alpha=1, color='w')


axs[0].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
axs[0].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)
axs[1].set_xlim(central_coordinate[1]-dlong, central_coordinate[1]+dlong)
axs[1].set_ylim(central_coordinate[0]-dlat,central_coordinate[0]+dlat)
axs[0].set_ylabel('latitude')
axs[0].set_xlabel('longitude')
axs[1].set_xlabel('longitude')
axs[0].set_title('Dataset 1', fontsize='medium')
axs[1].set_title('Dataset 2', fontsize='medium')
leg = axs[0].legend(loc='upper left', markerscale=4)
leg = axs[1].legend(loc='upper left', markerscale=4)

# Adding OpenStreetMap basemap
ctx.add_basemap(axs[0], crs=mask.crs.to_string(),
               source=ctx.providers.OpenStreetMap.Mapnik,
              zoom=13, alpha=.6)
ctx.add_basemap(axs[1], crs=mask.crs.to_string(),
               source=ctx.providers.OpenStreetMap.Mapnik,
              zoom=13, alpha=.6)
ax[1].add_artist(scalebar)

scalebar = ScaleBar(length*4, "m", length_fraction=0.4, location='lower left',
                    border_pad=1, box_color='w', box_alpha=0.9)


