# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 10:54:20 2021

@author: Christopher Willcocks, christopher.willcocks@hotmail.com
Master Student in Environemental engineering in EPFL
"""

## MASTER THESIS ## Modelling health accessibility in urban Sub-Saharan Africa

# Aim of this code is to construct the accessibility index as a function different modes of transport

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


r=0 # Number of transportation modes analyzed, initialized at 0
travel_speed = [1,2,3] #Number of transportation modes analyzed
# %%
#Entry in the multi-modal loop

for mode in travel_speed:
    r=r+1
    print(r)
    trip_times = list(range(0,30,1)) # 1 minutes precision, can be increased by changing 1 to 0.1 , 30 minutes = catchement threesold
    meters_per_minute = mode * 1000 / 60 #km per hour to meters per minute

    
    if mode==travel_speed[0]: # Walking loop
        G=ox.graph.graph_from_polygon(mask_transport['geometry'][0], network_type='walk', simplify=False, retain_all=False, truncate_by_edge=False, clean_periphery=False, custom_filter=None) # Walking network download
        print('Network downloaded')
        raster_paths = [r'C:\Users\Te Anau\Documents\Master thesis Christopher\Data\DEM\ALPSMLC30_N003E011_DSM.tif', r'C:\Users\Te Anau\Documents\Master thesis Christopher\Data\DEM\ALPSMLC30_N004E011_DSM.tif']
        G = ox.elevation.add_node_elevations_raster(G, raster_paths) # DEM files can be downloaded  https://www.eorc.jaxa.jp/ALOS/en/aw3d30/data/index.htm
        assert not np.isnan(np.array(G.nodes(data="elevation"))[:, 1]).any()
        
        G = ox.elevation.add_edge_grades(G, add_absolute=True)
        print('Network slopes calculated')
        
        # Brundson, C. (2018).Tobler’s hiking function. Retrieved May 25, 2021, from https://rpubs.com/chrisbrunsdon/hiking
        for u, v, k, data in G.edges(data=True, keys=True):
            data['time'] = data['length'] / (meters_per_minute*np.exp(-2.03*abs(data['grade']+0.133))) # Consider slope in the speed calculation
            
            
        # Nearest node from all shapefile points
        pop_centroid["Closest node ID"] = ox.distance.nearest_nodes(G, Xpop, Ypop)
        HStructures["Closest node ID"] = ox.distance.nearest_nodes(G, Xstr, Ystr)
        HStructures['lon'] = HStructures['geometry'].x
 
            
            
            
    
    if mode==travel_speed[1]: # Shared taxis loop
        G=G # Same as for walking as we consider walking if not on main axis

        # Attribute speed as function of the road considered, 30 km/h for motorway, 18 km/h for primary, secondary and tertiary, 3.57 for residential
        hwy_speeds = {'footway':3.57,
              'motorway': 30,
              'motorway_link': 30,
              'path':3.57,
              'pedestrian':3.57,
              'primary':18,
              'primary_link':18,
              'residential': 3.57,
              'road': 18,
              'secondary': 18,
              'secondary_link': 18,
              'service': 18,
              'steps': 3.57,
              'tertiary': 18,
              'tertiary_link': 18,
              'track': 18,
              'trunk': 18,
              'trunk_link': 18,
              'unclassified': 18}
        G=ox.speed.add_edge_speeds(G, hwy_speeds, fallback=None, precision=1) # Add chosen speeds to the network
        for u, v, k, data in G.edges(data=True, keys=True):
            data['time'] = data['length'] / (data['speed_kph']*1000/60)
        
            if data['speed_kph']==3.57:
                data['time'] = data['length'] / ((data['speed_kph']*1000/60)*np.exp(-2.03*abs(data['grade']+0.133)))
        
        # Closest node ID needs to be recalculated if G different than for Walking network
        pop_centroid["Closest node ID"] = ox.distance.nearest_nodes(G, Xpop, Ypop)
        HStructures["Closest node ID"] = ox.distance.nearest_nodes(G, Xstr, Ystr)
        HStructures['lon'] = HStructures['geometry'].x

        

    if mode==travel_speed[2]: # Motos-taxis network
        G=G
        #or G=ox.graph.graph_from_polygon(mask_transport['geometry'][0], network_type='drive', simplify=False, retain_all=False, truncate_by_edge=False, clean_periphery=False, custom_filter=None)
        print('Network downloaded')
        
        # Speeds of 23 km/h on main roads, 14 on others and 0.5 if walk is needed as it is considered particulliary unfavourable
        hwy_speeds = {'footway':3.5,
              'motorway': 23,
              'motorway_link': 23,
              'path':3.5,
              'pedestrian':3.5,
              'primary':23,
              'primary_link':23,
              'residential': 14,
              'road': 23,
              'secondary': 14,
              'secondary_link': 14,
              'service': 14,
              'steps': 0.1,
              'tertiary': 14,
              'tertiary_link': 14,
              'track': 14,
              'trunk': 14,
              'trunk_link': 14,
              'unclassified': 14}
        
        G=ox.speed.add_edge_speeds(G, hwy_speeds, fallback=None, precision=1)
        for u, v, k, data in G.edges(data=True, keys=True):
            data['time'] = data['length'] / (data['speed_kph']*1000/60)
            
            
        pop_centroid["Closest node ID"] = ox.distance.nearest_nodes(G, Xpop, Ypop)
        HStructures["Closest node ID"] = ox.distance.nearest_nodes(G, Xstr, Ystr)
        HStructures['lon'] = HStructures['geometry'].x
 

     #Variables initialization
    HStructures['Name']=HStructures['Source'] # Source column is not important, idea is to have a new column where infrastructure n can be transcribed
    HStructures_node = HStructures["Closest node ID"][0:len(HStructures)]
    ind=-1  # Number of loop in which we are, starts at 0 when coming in first facility
    infra_radius={}# List which will contain of all reachable facility in 30 minutes
    
    
    print('Entering loop')
    for n in HStructures_node: # Loop around all facilities
        print(n, end=' ')
        node_distances={} # In this variable, are stored the time needed to reach the facility for all locations (if they are within 30 minutes)
        HStructures_u_30={} # List which will contain of all reachable facility in 30 minutes, but is renewed on each loop
        ind=ind+1
        colname = 'Infrastructure n°%d' % ind
        pop_ratio_name='Pop_ratio %d' % ind
        print(colname)
        HStructures.loc[ind,('Name')]=colname
        
        # Loop modifiied but inspired from https://github.com/gboeing/osmnx-examples/blob/main/notebooks/13-isolines-isochrones.ipynb, Boeing, G. (2016).Gboeing/osmnx.
        # Retrieved Mars 19, 2021, from https://github.com/gboeing/osmnx
        for trip_time, trav_time in zip(sorted(trip_times, reverse=True), trip_times):
            subgraph = nx.ego_graph(G, n, radius=trip_time, distance='time')
            for node in subgraph.nodes():
                node_distances[node] = trip_time
                
        pop_centroid[colname] =pop_centroid["Closest node ID"].apply(lambda x: node_distances.get(x)) # Attribute the time needed to reach a facility from a closest node from the hectare centroid= Origin/destination matrix
        HStructures_u_30=pop_centroid.loc[pop_centroid[colname]<=30,colname]
        infra_radius[colname]=HStructures_u_30.index.values.tolist()
    
    # Origin destination matrix cleaning
    HStructures_OD=pop_centroid.T 
    HStructures_OD = HStructures_OD.iloc[8:]
    
    
    # Gaussian distance decay
    HStr_gen=pd.DataFrame(HStructures_OD).copy().fillna(0)
    HStr_gen[(HStr_gen>0) & (HStr_gen<=30)]=1
    HStr_gen[HStr_gen>30 ]=0
    
    # others distance decay functions proposed (not used in our example code)
    #HStr_gen can be replaced with one of the following distance decay functions
    
#Generalized 2SFCA

#HStructures_distance_decay_w_generalized=pd.DataFrame(HStructures_OD).copy().fillna(0)
#HStructures_distance_decay_w_generalized[(HStructures_distance_decay_w_generalized>0) & (HStructures_distance_decay_w_generalized<=30)]=1
#HStructures_distance_decay_w_generalized[HStructures_distance_decay_w_generalized>30 ]=0

#Enhanced 2SFCA
#HStructures_distance_decay_w_enhanced=pd.DataFrame(HStructures_OD).copy().fillna(0)
#HStructures_distance_decay_w_enhanced[(HStructures_distance_decay_w_enhanced>0) & (HStructures_distance_decay_w_enhanced<=5)]=1
#HStructures_distance_decay_w_enhanced[(HStructures_distance_decay_w_enhanced>5) & (HStructures_distance_decay_w_enhanced<=10)]=0.85
#HStructures_distance_decay_w_enhanced[(HStructures_distance_decay_w_enhanced>10) & (HStructures_distance_decay_w_enhanced<=15)]=0.70
#HStructures_distance_decay_w_enhanced[(HStructures_distance_decay_w_enhanced>15) & (HStructures_distance_decay_w_enhanced<=20)]=0.55
#HStructures_distance_decay_w_enhanced[(HStructures_distance_decay_w_enhanced>20) & (HStructures_distance_decay_w_enhanced<=25)]=0.40
#HStructures_distance_decay_w_enhanced[(HStructures_distance_decay_w_enhanced>25) & (HStructures_distance_decay_w_enhanced<=30)]=0.25
#HStructures_distance_decay_w_enhanced[HStructures_distance_decay_w_enhanced>30 ]=0

# G2SFCA
#HStructures_distance_decay_w_G2SFCA=pd.DataFrame(HStructures_OD).copy().fillna(0)
#HStructures_distance_decay_w_G2SFCA[HStructures_distance_decay_w_G2SFCA > 30 ]=0
#betha=2 # Need to be changed
#HStructures_distance_decay_w_G2SFCA[(HStructures_distance_decay_w_G2SFCA>0) & (HStructures_distance_decay_w_G2SFCA<=30)]= np.power(HStructures_distance_decay_w_G2SFCA[(HStructures_distance_decay_w_G2SFCA>0) & (HStructures_distance_decay_w_G2SFCA<=30)], -betha)

# Kernel density
#HStructures_distance_decay_w_kernel=pd.DataFrame(HStructures_OD).copy().fillna(0)
#HStructures_distance_decay_w_kernel[HStructures_distance_decay_w_kernel>30 ]=0
#HStructures_distance_decay_w_kernel[(HStructures_distance_decay_w_kernel>0) & (HStructures_distance_decay_w_kernel<=30)]= (3/4)*(1-(HStructures_distance_decay_w_kernel/30)**2)

# Gaussian density weigths
#HStructures_distance_decay_w_gaussian=pd.DataFrame(HStructures_OD).copy().fillna(0)
#HStructures_distance_decay_w_gaussian[HStructures_distance_decay_w_gaussian>30 ]=0
#HStructures_distance_decay_w_gaussian[(HStructures_distance_decay_w_gaussian>0) & (HStructures_distance_decay_w_gaussian<=30)]= ((2.71828**((-1/2)*((HStructures_distance_decay_w_gaussian/30)**2))-np.exp(-1/2))/(1 -np.exp(-1/2)))
    
    
    
    
    if mode==travel_speed[0]: # For walking scores
        OD_m1=pop_centroid.copy() # Proper Origin-destination matrix
        Dec_m1=HStr_gen.copy() # Weight distance_decay
        inf_rad_m1=infra_radius.copy() # List which show all included facilities
        
    if mode==travel_speed[1]:
        OD_m2=pop_centroid.copy()
        Dec_m2=HStr_gen.copy()
        inf_rad_m2=infra_radius.copy()
        
        
    if mode==travel_speed[2]:
        OD_m3=pop_centroid.copy()
        Dec_m3=HStr_gen.copy()
        inf_rad_m3=infra_radius.copy()        

    
    number_hectare_catchment_weighted={}
    catch_pop_name= 'Catchment pop m%d' %r
    HStructures[catch_pop_name]=HStructures['lon']
    HStructures.loc[:,(catch_pop_name)]=0 
  
# %%
# Step 1.1 of 2sfca Caclulation of catchement population


#   Variable initialization
HStructures['Catchment pop m']=HStructures['lon']
HStructures.loc[:,('Catchment pop m')]=0 
catch_m1={}
catch_m2={}
catch_m3={}
cam1={}
cam2={}
cam3={}
cam={}
f=0
#pop_centroid=OD_m1.fillna(0)
#pop_centroid=OD_m2.fillna(0)
#HStructures['Catchment pop m%d']=HStructures['lon'] %r
#HStructures['Catchment pop m%d'][:]=0 %r




while f<len(HStructures): # Loop again around all facilities
    colname = 'Infrastructure n°%d' % f
    catch_m1=OD_m1.loc[OD_m1[colname]<=30,'value']*Dec_m1.loc[colname] # Calculation of weighted population considered in catchement
    catch_m2=OD_m2.loc[OD_m2[colname]<=30,'value']*Dec_m2.loc[colname]
    catch_m3=OD_m3.loc[OD_m2[colname]<=30,'value']*Dec_m3.loc[colname]
    print(f)
    catch_m1=catch_m1.fillna(0) # replace nan value with 0
    catch_m2=catch_m2.fillna(0)
    catch_m3=catch_m3.fillna(0)
    cam1[colname]=sum(catch_m1)
    cam2[colname]=sum(catch_m2)
    cam3[colname]=sum(catch_m3)
    cam[colname]=sum(catch_m1)+sum(catch_m2)+sum(catch_m2)
    HStructures.loc[f, ('Catchment pop m1')]=cam1[colname]
    HStructures.loc[f, ('Catchment pop m2')]=cam2[colname]
    HStructures.loc[f, ('Catchment pop m3')]=cam3[colname]
    HStructures.loc[f, ('Catchment pop m')]=cam[colname]

    f=f+1
    print(f)



# %%

# Step 1.2 and 2 of the 2sfca

# initialization
pop_centroid['ratio m1']=pop_centroid['x']
pop_centroid['ratio m1'].values[:] = 0
pop_centroid['ratio m2']=pop_centroid['x']
pop_centroid['ratio m2'].values[:] = 0
pop_centroid['ratio m3']=pop_centroid['x']
pop_centroid['ratio m3'].values[:] = 0
pop_centroid['ratio m']=pop_centroid['x']
pop_centroid['ratio m'].values[:] = 0
value_hect=pop_centroid.index.values
Access={} # Access relates to the HStructures within the time radius defined. Access is the list of all facilities reachable from an hectare
w=0

###### For Walk
for w in value_hect: # Loop over all hectares
    list_of_keys = [key
                for key, list_of_values in inf_rad_m1.items()
                if value_hect[w] in list_of_values]
    Access[w]={w:list_of_keys}

    print(w/len(value_hect))
    yu=0
    ratio=0
    stop=len(Access[w][w])
    while yu < stop: # Calculate ratio for all facilities included in catchment , which means from 0 to stop 
        ratio=(Dec_m1[w][Access[w][w][yu]])*(1/(HStructures.loc[HStructures['Name']==Access[w][w][yu], 'Catchment pop m1'])) # calulation of facility to population/ratio
        ratio.fillna(0)
        pop_centroid['ratio m1'][w]=pop_centroid['ratio m1'][w]+ratio  # Addition of ratios = accessibility scores
        yu=yu+1

###### For shared taxis
Access={}
w=0
for w in value_hect:
    list_of_keys = [key
                for key, list_of_values in inf_rad_m2.items()
                if value_hect[w] in list_of_values]
    Access[w]={w:list_of_keys}

    if w%100==0:
        print(w/len(value_hect))
    yu=0
    ratio=0
    stop=len(Access[w][w])
    while yu < stop:
        ratio=(Dec_m2[w][Access[w][w][yu]])*(1/(HStructures.loc[HStructures['Name']==Access[w][w][yu], 'Catchment pop m2']))
        ratio.fillna(0)
        pop_centroid['ratio m2'][w]=pop_centroid['ratio m2'][w]+sum(ratio)
        yu=yu+1

##### For moto-taxis
Access={}
w=0
for w in value_hect:
    list_of_keys = [key
                for key, list_of_values in inf_rad_m3.items()
                if value_hect[w] in list_of_values]
    Access[w]={w:list_of_keys}


    if w%100==0:
        print(w/len(value_hect))
    yu=0
    ratio=0
    stop=len(Access[w][w])
    while yu < stop:
        ratio=(Dec_m3[w][Access[w][w][yu]])*(1/(HStructures.loc[HStructures['Name']==Access[w][w][yu], 'Catchment pop m3']))
        ratio.fillna(0)       
        pop_centroid['ratio m3'][w]=pop_centroid['ratio m3'][w]+sum(ratio)
        yu=yu+1
        
# Calculation of final score
# 0.4, 0.4 and 0.2 are the transportation splits of walking, shared taxis and motos-taxis in Yaounde        
pop_centroid['ratio m']=0.4*pop_centroid['ratio m1']+0.4*pop_centroid['ratio m2']+0.2*pop_centroid['ratio m3']


# Saving results in new geoDataFrame 

Accessibility_scores=gpd.GeoDataFrame(pop_centroid['geometry'])
Accessibility_scores['ratio m']=pop_centroid['ratio m']
Accessibility_scores['ratio m1']=pop_centroid['ratio m1']
Accessibility_scores['ratio m2']=pop_centroid['ratio m2']
Accessibility_scores['ratio m3']=pop_centroid['ratio m3']


outfp = r'C:\Users\Te Anau\Documents\Master thesis Christopher\Results\Accessibility_scores.shp' # Location to be adapted
Accessibility_scores.to_file(outfp)


