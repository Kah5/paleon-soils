#This code reads the .asc soil rasters into R as rasters, plots them, and converts them to dataframes
library(raster)

#see README.txt for descriptions of the data and citation for STATSGO2 soil data


#simons code to create the 8km albers grid that I used to get the soil data on the 
#same grid scale
#input <- raster('paleon.unit.alb_01.tif')
#input <- setValues(input, 1:ncell(input))
#writeRaster(input, "paleon8kmalbers.tif")


#setwd('location of your soil directory')

#read the soil rasters into R
#example for ksat
ksat<-raster("paleon_ksat.asc")
#define Albers projcetion for this layer
proj4string(ksat)<-CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(ksat)
#convert this raster to a dataframe to work with the data
#adding xy=TRUE outputs the x-y albers coordinates with the data value
ksat.df<-as.data.frame(ksat, xy=TRUE)

#Repeat this process for the rest of the soil characters
#for AWC
awc<-raster("paleon_awc.asc")
proj4string(awc)<-CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(awc)
awc.df<-as.data.frame(awc, xy=TRUE)

#for sand
sand<-raster("paleon_sand.asc")
proj4string(sand)<-CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(sand)
sand.df<-as.data.frame(sand, xy=TRUE)

#for silt
silt<-raster("paleon_silt.asc")
proj4string(silt)<-CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(silt)
silt.df<-as.data.frame(silt, xy=TRUE)

#for clay
clay<-raster("paleon_clay.asc")
proj4string(clay)<-CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(clay)
clay.df<-as.data.frame(clay, xy=TRUE)

#for elev
elev<-raster("paleon_elev.asc")
proj4string(elev)<-CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(elev)
elev.df<-as.data.frame(elev, xy=TRUE)
