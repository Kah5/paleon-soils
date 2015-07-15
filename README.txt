This repository contains soil parameters from STATSGO database that have been rasterized on the Paleon 8km albers grid. Code for reading them into R and extracting values for Kelly's Summer 2015 field sites are housed in this repository. 

R code in the R/ subfolder, while all data used is in the Data/ subfolder. 

R/ 
'surgo_soil_data_processing.R' processes the soil data into paleon grids
'field_sites.Rmd' reads in the kml/kmz files generated using google maps, then extracts point level environmental data at these sites for future use.


Citation for STATSGO data: Soil Survey Staff, Natural Resources Conservation Service, United States Department of Agriculture. Web Soil Survey. Available online at http://websoilsurvey.nrcs.usda.gov/. Accessed [3/20/2015].


Rasters were generated from the STATSGO2 database, projected to the paleon grid, and rasterized using the PalEON Albers grid in ArcGIS.

The following description of variables comes from the SSURGO_Metadata Column Descriptions:

ksat: The amount of water that would move vertically through a unit area of saturated soil in unit time under unit hydraulic gradient. The rasters in this directory include rasters on the 8km PalEON Albers grid over the Midwest (IL, IN, MN, MI).

paleon_slope: The difference in elevation between two points, expressed as a percentage of the distance between those points. (SSM)

paleon_awc: The amount of water that an increment of soil depth, inclusive of fragments, can store that is available to plants. AWC is expressed as a volume fraction, and is commonly estimated as the difference between the water contents at 1/10 or 1/3 bar (field capacity) and 15 bars (permanent wilting point) tension and adjusted for salinity, and fragments.

paleon_clay: Mineral particles less than 0.002mm in equivalent diameter as a weight percentage of the less than 2.0mm fraction.

paleon_sand: Mineral particles 0.05mm to 2.0mm in equivalent diameter as a weight percentage of the less than 2 mm fraction.

paleon_silt: Mineral particles 0.002 to 0.05mm in equivalent diameter as a weight percentage of the less than 2.0mm fraction.