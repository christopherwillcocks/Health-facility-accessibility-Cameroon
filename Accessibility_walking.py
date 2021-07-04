# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 13:21:47 2021

@author: Christopher Willcocks, christopher.willcocks@hotmail.com
Master Student in Environemental engineering in EPFL
"""

## MASTER THESIS ## Modelling health accessibility in urban Sub-Saharan Africa

# Aim of this code is to construct the accessibility index for walking transportation mode
# more detailed comment presented in code Accessibility_index_multi-modal

import osmnx as ox
import networkx as nx
import pandas as pd
import geopandas as gpd



pop_peanutButter=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\pop_peanutButter.shp') # Population datasets converted with Pre-Processing steps (see Pre-processing code)
HStructures_clinic_geosm=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\HStructures_clinic_geosm.shp') # Healthcare datasets converted with Pre-Processing steps
HStructures_hospital_geosm=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\HStructures_hospital_geosm.shp') # Healthcare datasets converted with Pre-Processing steps
mask=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\mask.shp') # Mask created in the pre-processing 
mask_transport=gpd.read_file(r'C:\Users\Te Anau\Documents\Data\mask_transport.shp') # Mask created in the pre-processing 


pop_centroid=pop_peanutButter # Can be adapted
HStructures = HStructures_clinic_geosm # Depends ont the set wanted to be analyzed



Xpop = pop_centroid['geometry'].map(lambda pt: pt.coords[0][0]) 
Ypop = pop_centroid['geometry'].map(lambda pt: pt.coords[0][1])

Xstr = HStructures['geometry'].map(lambda pt: pt.coords[0][0])
Ystr = HStructures['geometry'].map(lambda pt: pt.coords[0][1])

r=0
travel_speed = 3.57 #walking speed in km/hour
trip_times = list(range(0,30,1))
meters_per_minute = travel_speed * 1000 / 60 #km per hour to meters per minute

# Download transportation network from osmnx 
G=ox.graph.graph_from_polygon(mask_transport['geometry'][0], network_type='walk', simplify=False, retain_all=False, truncate_by_edge=False, clean_periphery=False, custom_filter=None)
print('Network downloaded')

raster_paths = [r'C:\Users\Te Anau\Documents\Master thesis Christopher\Data\DEM\ALPSMLC30_N003E011_DSM.tif', r'C:\Users\Te Anau\Documents\Master thesis Christopher\Data\DEM\ALPSMLC30_N004E011_DSM.tif']
G = ox.elevation.add_node_elevations_raster(G, raster_paths) #DEM  # Brundson, C. (2018).Tobler’s hiking function. Retrieved May 25, 2021, from https://rpubs.com/chrisbrunsdon/hiking
assert not np.isnan(np.array(G.nodes(data="elevation"))[:, 1]).any()

G = ox.elevation.add_edge_grades(G, add_absolute=True)
print('Network slopes calculated')

for u, v, k, data in G.edges(data=True, keys=True):
    data['time'] = data['length'] / (meters_per_minute)
    data['time'] = data['length'] / (meters_per_minute*np.exp(-2.03*abs(data['grade']+0.133))) # Brundson, C. (2018).Tobler’s hiking function. Retrieved May 25, 2021, from https://rpubs.com/chrisbrunsdon/hiking
    
    
# Attribute nearest node to the network
pop_centroid["Closest node ID"] = ox.distance.nearest_nodes(G, Xpop, Ypop)
HStructures["Closest node ID"] = ox.distance.nearest_nodes(G, Xstr, Ystr)
HStructures['lon'] = HStructures['geometry'].x
 
    
HStructures['Name']=HStructures['lon']
HStructures_node = HStructures["Closest node ID"][0:len(HStructures)]
ind=-1

infra_radius={}


print('Entering loop')
for n in HStructures_node: # Loop over all healthcare facility
    print(n, end=' ')
    node_distances={}
    HStructures_u_30={}
    ind=ind+1
    colname = 'Infrastructure n°%d' % ind
    pop_ratio_name='Pop_ratio %d' % ind
    print(colname)

    HStructures.loc[ind,('Name')]=colname
    
    for trip_time, trav_time in zip(sorted(trip_times, reverse=True), trip_times):
        subgraph = nx.ego_graph(G, n, radius=trip_time, distance='time')
        for node in subgraph.nodes():
            node_distances[node] = trip_time
    pop_centroid[colname] =pop_centroid["Closest node ID"].apply(lambda x: node_distances.get(x)) # attributing distance to all hectares --> matrix Origin-destination filling

    HStructures_u_30=pop_centroid.loc[pop_centroid[colname]<=30,colname]
    infra_radius[colname]=HStructures_u_30.index.values.tolist() # End of the loop
    
    
HStructures_OD=pop_centroid.T
HStructures_OD = HStructures_OD.iloc[8:]



# gaussian distance decay weights 
# Other decay function proposed in code Accessibility_index_multi-modal
HStr_gaus=pd.DataFrame(HStructures_OD).copy().fillna(0)
HStr_gaus[HStr_gaus>30 ]=0 # weight=0 if more than 30 minutes away
HStr_gaus[(HStr_gaus>0) & (HStr_gaus<=30)]= ((2.71828**((-1/2)*((HStr_gaus/30)**2))-np.exp(-1/2))/(1 -np.exp(-1/2)))

OD_m1=pop_centroid.copy() # Oriigin_destination matrix for mode 1 which is walk

del OD_m1['row']# delete all relevant columns
del OD_m1['col']
del OD_m1['value']
del OD_m1['x']
del OD_m1['y']
del OD_m1['Location']
del OD_m1['geometry']
del OD_m1['Closest node ID']


Dec_m1=HStr_gaus.copy() # Distance decay weights
inf_rad_m1=infra_radius.copy()
    

    
#%% 2sfca part 1.1 calculation of population within all catchments



catch_pop_name= 'Catchment pop m1'
HStructures[catch_pop_name]=HStructures['lon']
HStructures.loc[:,(catch_pop_name)]=0 
catch_m1={}
cam1={}

f=0


while f<len(HStructures):
    colname = 'Infrastructure n°%d' % f
    catch_m1=OD_m1.loc[OD_m1[colname]<=30,'value']*Dec_m1.loc[colname] # Population considered within catchment
    print(f)
    catch_m1=catch_m1.fillna(0)
    cam1[colname]=sum(catch_m1)
    HStructures.loc[f, ('Catchment pop m1')]=cam1[colname] # population within one facility area
    f=f+1






# %% 2sfca part 1.2 and 2


pop_centroid['ratio m1']=pop_centroid['x']
pop_centroid['ratio m1'].values[:] = 0
value_hect=pop_centroid.index.values
Access={}
w=0
# Access relates to the HStructures within the time radius defined
for w in value_hect:
    list_of_keys = [key
                for key, list_of_values in inf_rad_m1.items()
                if value_hect[w] in list_of_values]
    Access[w]={w:list_of_keys} # Access list all facilities in the catchement constructed over one hectare

    if w%100==0:
        print(w/len(value_hect))
    yu=0
    ratio=0
    stop=len(Access[w][w])
    while yu < stop: # stop is for all hectares the number of facilities with an acceptable walking distance
        ratio=(Dec_m1[w][Access[w][w][yu]])*(1/(HStructures.loc[HStructures['Name']==Access[w][w][yu], 'Catchment pop m1'])) # Calulation of ratio facility per population
        ratio.fillna(0)
        yu=yu+1
        pop_centroid['ratio m1'][w]=pop_centroid['ratio m1'][w]+ratio # The addition of the ratio > 0 is the walking accessibility score
        yu=yu+1
        
        
# %%
#Saving results
Walking_scores=gpd.GeoDataFrame(pop_centroid['geometry'])
Walking_scores['ratio m1']=pop_centroid['ratio m1']


outfp = r'C:\Users\Te Anau\Documents\Master thesis Christopher\Results\Walking_scores.shp' # Location to be adapted
Walking_scores.to_file(outfp)

