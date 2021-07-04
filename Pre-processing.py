# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 09:24:04 2021

@author: Christopher Willcocks
"""

## MASTER THESIS ## Modelling health accessibility in urban Sub-Saharan Africa

# Aim of this code is to delimit and convert population and healthfacilities datsets into shp




import rasterio as rio
import georasters as gr
from geopandas import gpd
import urban_footprinter as ufp
import os
print(os.getcwd())

# Reading population datasets
# ! To be changed according to files location
#Files imported from :https://www.worldpop.org/project/categories?id=3 and Accessible at :https://apps.worldpop.org/peanutButter/

raster_1_filepath = r'C:\Users\Te Anau\Documents\Data\Population_dataset_worldpop_UN_adj.tif'
raster_4_filepath = r'C:\Users\Te Anau\Documents\Data\Population_dataset_peanutButter.tif'

# Reading healthfacilities datasets
# ! To be changed according to files location
# Files imported from  http://geocameroun.cm/
HStructures_clinic_geosm = gpd.read_file(r'C:\Users\Te Anau\Documents\Data\Clinique_geosm.gpkg')
HStructures_hospital_geosm = gpd.read_file(r'C:\Users\Te Anau\Documents\Data\Hopital_geosm.gpkg')


# %%
# Opening and converting WorldPop population dataset to binary set (inhabited or not)
with rio.open(raster_1_filepath) as src:
    lulc_arr = src.read(1)
    lulc_arr[lulc_arr>0.1]=1 # if population is more that 0.1 person per hectare , we consider it inhabitated
    res = src.res[0]  # only square pixels are supported


# Writing binary set as reference
raster_base=r'raster_base.tif'
with rio.open(
    raster_base, 
    'w',
    driver='GTiff',
    height=lulc_arr.shape[0],
    width=lulc_arr.shape[1],
    count=1,
    dtype=lulc_arr.dtype,
    crs='epsg:4326', transform=src.transform,
) as dst:
    dst.write(lulc_arr, 1)


#%%
# urban footpinter function
# Bosch, M. (2020). Urban footprinter: a convolution-based approach to
# detect urban extents fromraster data. https://doi.org/10.5281/zenodo.3699310

# urban_mask is the urban area of interest of our study. A 500m radius was added

convert_ratio=100/res   # Our population data is in hectare (100 x 100), we prefer to work in meters
urban_mask=ufp.urban_footprint_mask_shp(raster_base, urban_classes=range(1,2), kernel_radius=500/convert_ratio, 
                                        urban_threshold=0.20, buffer_dist=500/convert_ratio)

# Creating a geodataframe (geopandas) file for mask
for geom in urban_mask:
    rt=geom.wkt
urban_mask_geometry={'geometry':[rt], 'Location': ['Yaounde']}
mask=gpd.GeoDataFrame(urban_mask_geometry)
mask['geometry'] = gpd.GeoSeries.from_wkt(mask['geometry'])
mask.crs = "EPSG:4326"



# Mask to delimit the health facilities
urban_mask_HStructures=ufp.urban_footprint_mask_shp(raster_base, urban_classes=range(1,2), kernel_radius=500/convert_ratio, 
                                        urban_threshold=0.20, buffer_dist=2000/convert_ratio)
for geom in urban_mask_HStructures:
    rt=geom.wkt
urban_HStructures_mask_geometry={'geometry':[rt], 'Location': ['Yaounde']}
mask_HStructures=gpd.GeoDataFrame(urban_HStructures_mask_geometry)
mask_HStructures['geometry'] = gpd.GeoSeries.from_wkt(mask_HStructures['geometry'])
mask_HStructures.crs = "EPSG:4326"



# Mask to delimit the transportation network
urban_mask_transportation=ufp.urban_footprint_mask_shp(raster_base, urban_classes=range(1,2), kernel_radius=500/convert_ratio, 
                                        urban_threshold=0.20, buffer_dist=8000/convert_ratio)
for geom in urban_mask_transportation:
    rt=geom.wkt
urban_tansportation_mask_geometry={'geometry':[rt], 'Location': ['Yaounde']}
mask_transport=gpd.GeoDataFrame(urban_tansportation_mask_geometry)
mask_transport['geometry'] = gpd.GeoSeries.from_wkt(mask_transport['geometry'])




# %%
# Conversion and delimitation of population datasets


grast=gr.from_file(raster_1_filepath)
gdf_worldpop_UN_adj=gr.to_geopandas(grast)
gdf_worldpop_UN_adj=gdf_worldpop_UN_adj.to_crs(epsg=4326)
gdf_worldpop_UN_adj['geometry']=gdf_worldpop_UN_adj.centroid # Converting polygon shapefiles to points shapfiles as the centroid is the point of interest
pop_worldpop_UN_adj=gpd.overlay(gdf_worldpop_UN_adj, mask, how='intersection') # Delimiting population dataset withthe urban footprint


grast=gr.from_file(raster_4_filepath)
gdf_peanutButter=gr.to_geopandas(grast)
gdf_peanutButter=gdf_peanutButter.to_crs(epsg=4326)
gdf_peanutButter['geometry']=gdf_peanutButter.centroid
pop_peanutButter=gpd.overlay(gdf_peanutButter, mask, how='intersection')



# Conversion and delimitation for health facilities datasets

HStructures_clinic_geosm=gpd.overlay(HStructures_clinic_geosm, mask_HStructures, how='intersection')
HStructures_hospital_geosm=gpd.overlay(HStructures_hospital_geosm, mask_HStructures, how='intersection')



# %%

#Save converted files

outfp = r'C:\Users\Te Anau\Documents\Data\pop_worldpop_UN_adj.shp'
pop_worldpop_UN_adj.to_file(outfp)


outfp = r'C:\Users\Te Anau\Documents\Data\HStructures_clinic_geosm.shp'
HStructures_clinic_geosm.to_file(outfp)


outfp = r'C:\Users\Te Anau\Documents\Data\HStructures_hospital_geosm.shp'
HStructures_hospital_geosm.to_file(outfp)

