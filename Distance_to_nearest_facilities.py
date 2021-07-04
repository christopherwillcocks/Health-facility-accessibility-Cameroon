# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 10:50:18 2021

@author: Christopher Willcocks, christopher.willcocks@hotmail.com
Master Student in Environemental engineering in EPFL
"""

## MASTER THESIS ## Modelling health accessibility in urban Sub-Saharan Africa

# Aim of this code is to construct Origine-destination matrix as a function of distance and then calculate the 
# distance to the nearest facility

import osmnx as ox
import networkx as nx
import pandas as pd
import geopandas as gpd




pop_peanutButter=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\pop_peanutButter.shp') # Population datasets converted with Pre-Processing steps (see Pre-processing code)
HStructures_clinic_geosm=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\HStructures_clinic_geosm.shp') # Healthcare datasets converted with Pre-Processing steps
HStructures_hospital_geosm=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\HStructures_hospital_geosm.shp') # Healthcare datasets converted with Pre-Processing steps
mask=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\mask.shp') # Mask created in the pre-processing 
mask_transport=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\mask_transport.shp') # Mask created in the pre-processing 


pop_centroid=pop_peanutButter # Only used to have hectaric results
HStructures = HStructures_clinic_geosm # Health facility dataeets used for the present code example

trip_times = list(range(0, 10000, 250)) # times in meters # goes from 0 to 10'000 meters with steps of 250 meters
# To increase speed, take bigger distance steps

# OSMnx walking network download
G=ox.graph.graph_from_polygon(mask_transport['geometry'][0], network_type='walk', simplify=False, retain_all=False, truncate_by_edge=False, clean_periphery=False, custom_filter=None)
print('Network downloaded')

# Finding node ID of the closest node to shapefile points
X = pop_centroid['geometry'].map(lambda pt: pt.coords[0][0])
Y = pop_centroid['geometry'].map(lambda pt: pt.coords[0][1])
pop_centroid["Closest node ID"] = ox.distance.nearest_nodes(G, X, Y)
del X,Y
X = HStructures['geometry'].map(lambda pt: pt.coords[0][0])
Y = HStructures['geometry'].map(lambda pt: pt.coords[0][1])
HStructures["Closest node ID"] = ox.distance.nearest_nodes(G, X, Y)
del X, Y

# %%
# Defining variables before entering the loop

HStructures['Name']=HStructures['name'] # In Name column, we will give thename of the facility from 0  len(HStructures)
HStructures_node = HStructures["Closest node ID"][0:len(HStructures)] # Number of nodes to compute and their ID
ind=-1 # When entering the lopp, starts at 0

print('Entering loop')

# %%

# Loop modifiied but inspired from https://github.com/gboeing/osmnx-examples/blob/main/notebooks/13-isolines-isochrones.ipynb, Boeing, G. (2016).Gboeing/osmnx.
# Retrieved Mars 19, 2021, from https://github.com/gboeing/osmnx


for n in HStructures_node: # Loop around all facility
    print(n, end=' ')
    node_distances={} # Node distance from this node to all others, within subgraph
    ind=ind+1
    colname = 'Infrastructure n°%d' % ind
    pop_ratio_name='Pop_ratio %d' % ind
    print(colname)
    HStructures.loc[ind,('Name')]=colname # Create a new column for each analyzed facilities
    
    for trip_time, trav_time in zip(sorted(trip_times, reverse=True), trip_times): # Loop over all distances
        subgraph = nx.ego_graph(G, n, radius=trip_time, distance='length')
        for node in subgraph.nodes():
            node_distances[node] = trip_time # Attribute a distance to all nodes within the subgraph
    pop_centroid[colname] =pop_centroid["Closest node ID"].apply(lambda x: node_distances.get(x)) # Present results in pop_centroid
# End of the loop  
    
    
pop_centroid['Time to closest inf']=pop_centroid.iloc[: , 8:].min(axis=1) # Calculation of the smallest distance value from all hectares = distance to the nearest facility

# Creation of dataframe with relation reachable population with distance to nearest facility
sum_pop=pd.DataFrame(trip_times, trip_times) 

for val in trip_times: # Loop around all distance (named time) to calculate population related
    print(val)
    sum_pop.loc[val,(0)]=sum(pop_centroid.loc[pop_centroid['Time to closest inf']<=val,'value'])
    sum_pop.loc[val,(1)]=100*(sum(pop_centroid['value'])-sum_pop.loc[val,(0)])/sum(pop_centroid['value']) # Calculation of population over this distance
    
# %%
# Saving results
# Location to be adapted=gpd.GeoDataFrame(pop_centroid['geometry']) # New GeoDataFrame
Nearest_facility=gpd.GeoDataFrame(pop_centroid['geometry'])
Nearest_facility['Time to closest inf']=pop_centroid['Time to closest inf']

Nearest_pop_function=gpd.GeoDataFrame(list(sum_pop[0]))
Nearest_pop_function['abs population - abs']=list(sum_pop[0])
Nearest_pop_function['prop population - unreachable']=list(sum_pop[1])
Nearest_pop_function['geometry']=pop_centroid['geometry'][0:4]
del Nearest_pop_function[0]

outfp = r'C:\Users\Te Anau\Documents\Master thesis Christopher\Results\Nearest_facility.shp' # Location to be adapted
Nearest_facility.to_file(outfp)

outfp = r'C:\Users\Te Anau\Documents\Master thesis Christopher\Results\Nearest_pop_function.shp' # Location to be adapted
Nearest_pop_function.to_file(outfp)

    