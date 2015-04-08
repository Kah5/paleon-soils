##This code reads the .asc soil rasters into R as rasters, plots them, and converts them to dataframes
library(raster)
library(sp)
#see README.txt for descriptions of the data and citation for STATSGO2 soil data

#simons code to create the 8km albers grid that I used to get the soil data on the 
#same grid scale
input <- raster('data/Abies_balsamea_diam_alb.tif') #this file came from biomass dropbox
input <- setValues(input, 1:ncell(input))
writeRaster(input, "paleon8km_midalbers1.tif", overwrite=TRUE)

input.df<-as.data.frame(input, xy=TRUE)
#setwd('location of your soil directory')


##read the soil rasters into R
##example for ksat
ksat<-raster("Data/paleon_ksat.asc")
#define Albers projcetion for this layer
proj4string(ksat)<-CRS('+init=epsg:3175')

#proj4string(ksat)<-CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
plot(ksat) #note that the extent of these data goes beyond statistical output
###can't use extract because ksat is a "RasterLayer", not a raster

ksat3 <- resample(ksat, input, method='ngb') #resample ksat based on the psleon8km_midalbers1 grid created above
ksat2<-stack(ksat3, input) #ksat has values of soil parameter, while input has values of "cell number of interest"
#ksat$cell <- setValues(ksat, 1:ncell(ksat)) #not sure if this is a valid thing to do, but this creates RasterBrick

#need to convert this raster to a dataframe to work with the data
ksat2.df <- as.data.frame(ksat2, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value

write.csv(ksat2.df, "ksat.soil.csv")




##this biomass data does not have the same x and y values as the statistical model
#therefore, this biomass grid does not match up when I try to merge the biomass data with the soil data
biomass <- read.csv('Data/plss_biomass_alb_v0.9-3.csv')#this is gridded non-biomass model biomass
biomass.df<-data.frame(biomass)
coordinates(biomass) <- ~x+y

biomass$soil <- extract(ksat3, biomass) # this extracts the ksat value for each biomass cell of the biomass dataframe
proj4string(biomass) <- CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')


