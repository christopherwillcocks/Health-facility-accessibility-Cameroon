# Health-facility-accessibility-Cameroon
This repository presents the python codes I have created to constructed my accessiiblity index in my 2021 Master project in EPFL (environmental engineering). The aim of the project was to model healthcare facility accessibility in the urban Sub-Saharan context. 


Four compuations codes can be found in this github repository and two plots examples codes. 

To use these codes , gridded population estimates and healthfacilities datasets are needed. DEM are optional They can be extracted:
gridded population estimates:  :https://apps.worldpop.org/peanutButter/, https://geopode.world/
healthfacilities: http://geocameroun.cm,  https://www.who.int/publications/m/item/who-cds-gmp-2019-01
DEM:  https://www.eorc.jaxa.jp/ALOS/en/aw3d30/data/index.htm


- Pre-processing : All pre-preprocessing steps performed to have cleaned shapefile data.
- Distance_to_nearest_facility: Computation of the distance between all hectares centroid and their nearest facility
- Accessibility_walking: Computation of accessibility scores for walking times 
- Accessibility_index_multi_modal: Computation of the mutli-modal scores using the gaussian 2sfca with 30 minutes. Three transportations modes are considered (walk, shared taxis and motostaxis)

Distance decay function can modified. Other functions are proposed in Accessibility_index_multi_modal code. The catchment size which we have fixed at 30 minutes can also be modified. 
To calculate the score of shared-taxis or motos-taxis individually, you'll need to to modify Accessibility_walking code with speeds proposed in Accessibility_index_multi_modal (line 79 to 103 or resp. 113 to 145).


I have also shared models of the plots I have done for my report. 
