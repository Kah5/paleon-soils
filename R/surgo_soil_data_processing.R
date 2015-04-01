#This code reads the .asc soil rasters into R as rasters, plots them, and converts them to dataframes
library(raster)
library(sp)
#see README.txt for descriptions of the data and citation for STATSGO2 soil data

#simons code to create the 8km albers grid that I used to get the soil data on the 
#same grid scale
input <- raster('data/Abies_balsamea_diam_alb.tif') #this file came from biomass dropbox
input <- setValues(input, 1:ncell(input))
writeRaster(input, "paleon8km_midalbers1.tif")

input.df<-as.data.frame(input, xy=TRUE)
#setwd('location of your soil directory')


##read the soil rasters into R
##example for ksat
ksat<-raster("Data/paleon_ksat.asc")
#define Albers projcetion for this layer
proj4string(ksat)<-CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(ksat) #note that the extent of these data goes beyond statistical output
ksat <- resample(ksat, input) #resample ksat based on the psleon8km_midalbers1 grid created above

#need to convert this raster to a dataframe to work with the data

ksat.df <- as.data.frame(ksat, xy=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value



##Repeat this process for the rest of the soil characters
##for AWC
awc <- raster("Data/paleon_awc.asc")
proj4string(awc) <- CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(awc)
awc <- resample(awc, input)
awc.df <- as.data.frame(awc, xy=TRUE)

##for sand
sand <- raster("Data/paleon_sand.asc")
proj4string(sand) <- CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(sand)
sand <- resample(sand, input)
sand.df <- as.data.frame(sand, xy=TRUE)

##for silt
silt <- raster("Data/paleon_silt.asc")
proj4string(silt) <- CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(silt)
silt <- resample(silt, input)
silt.df <- as.data.frame(silt, xy=TRUE)

##for clay
clay <- raster("Data/paleon_clay.asc")
proj4string(clay) <- CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(clay)
clay <- resample(clay, input)
clay.df <- as.data.frame(clay, xy=TRUE)

#for elev from statsgo
elev<-raster("Data/paleon_elev.asc")
proj4string(elev) <- CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(elev)
elev <-resample(elev, input)
elev.cell <- extract(elev,extent(elev))
elev.df <- as.data.frame(elev, xy=TRUE)


##this biomass data does not have the same x and y values as the statistical model
#therefore, this biomass grid does not match up when I try to merge the biomass data with the soil data
biomass <- read.csv('Data/plss_biomass_alb_v0.9-3.csv')#this is gridded non-biomass model biomass
biomass.df<-data.frame(biomass)
coordinates(biomass) <- ~x+y
biomass <- raster(biomass)
values(biomass)<-biomass

proj4string(biomass) <- CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')


sample.merged.data<-merge(biomass.df, ksat.df, by =c("x", "y"))


