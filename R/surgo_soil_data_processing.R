##This code reads the .asc soil rasters into R as rasters, plots them, and converts them to dataframes
library(raster)
library(sp)
#see README.txt for descriptions of the data and citation for STATSGO2 soil data

#simons code to create the 8km albers grid that I used to get the soil data on the 
#same grid scale
input <- raster('data/Abies_balsamea_diam_alb.tif') #this file came from biomass dropbox
input <- setValues(input, 1:ncell(input))
writeRaster(input, "paleon8km_midalbers1.tif", overwrite=TRUE)

input.df <- as.data.frame(input, xy=TRUE)
#setwd('location of your soil directory')


##read the soil rasters into R
##example ksat
ksat <- raster("Data/paleon_ksat.asc")
#define Albers projcetion for this layer
proj4string(ksat) <- CRS('+init=epsg:3175')
plot(ksat) #note that the extent of these data goes beyond statistical output

##use resample to realign if the cells are off by a small amount
ksat3 <- resample(ksat, input, method='ngb') #resample ksat based on the psleon8km_midalbers1 grid created above
ksat2 <- stack(ksat3, input) #ksat has values of soil parameter, while input has values of "cell number of interest"
#ksat$cell <- setValues(ksat, 1:ncell(ksat)) #not sure if this is a valid thing to do, but this creates RasterBrick

##need to convert this raster to a dataframe to work with the data
ksat2.df <- as.data.frame(ksat2, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value
write.csv(ksat2.df, "ksat.soil.csv") #write to csv so you can haave it by itself

##########################
#Rest of the soil data
##########################
#clay
clay <- raster("Data/paleon_clay.asc")
#define Albers projcetion for this layer
proj4string(clay) <- CRS('+init=epsg:3175')
plot(clay) #note that the extent of these data goes beyond statistical output
clay3 <- resample(clay, input, method='ngb') #resample ksat based on the psleon8km_midalbers1 grid created above
clay2 <- stack(clay3, input) #ksat has values of soil parameter, while input has values of "cell number of interest"
#ksat$cell <- setValues(ksat, 1:ncell(ksat)) #not sure if this is a valid thing to do, but this creates RasterBrick

##need to convert this raster to a dataframe to work with the data
clay3.df <- as.data.frame(clay3, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value
write.csv(clay3.df, "clay.soil.csv") 

########################
#sand
#######################
sand <- raster("Data/paleon_sand.asc")
#define Albers projcetion for this layer
proj4string(sand) <- CRS('+init=epsg:3175')
plot(sand) #note that the extent of these data goes beyond statistical output
sand3 <- resample(sand, input, method='ngb') #resample ksat based on the psleon8km_midalbers1 grid created above
sand2 <- stack(sand3, input) #ksat has values of soil parameter, while input has values of "cell number of interest"
#ksat$cell <- setValues(ksat, 1:ncell(ksat)) #not sure if this is a valid thing to do, but this creates RasterBrick

##need to convert this raster to a dataframe to work with the data
sand3.df <- as.data.frame(sand3, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value
write.csv(sand3.df, "sand.soil.csv")

#####################
#silt
#####################
silt <- raster("Data/paleon_silt.asc")
#define Albers projcetion for this layer
proj4string(silt) <- CRS('+init=epsg:3175')
plot(silt) #note that the extent of these data goes beyond statistical output
silt3 <- resample(silt, input, method='ngb') #resample ksat based on the psleon8km_midalbers1 grid created above
silt2 <- stack(silt3, input) #ksat has values of soil parameter, while input has values of "cell number of interest"
#ksat$cell <- setValues(ksat, 1:ncell(ksat)) #not sure if this is a valid thing to do, but this creates RasterBrick

##need to convert this raster to a dataframe to work with the data
silt3.df <- as.data.frame(silt3, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value
write.csv(silt3.df, "silt.soil.csv")

#######################
#elev
#######################
elev <- raster("Data/paleon_elev.asc")
#define Albers projcetion for this layer
proj4string(elev) <- CRS('+init=epsg:3175')
plot(elev) #note that the extent of these data goes beyond statistical output
elev3 <- resample(elev, input, method='ngb') #resample ksat based on the psleon8km_midalbers1 grid created above
elev2 <- stack(elev3, input) #ksat has values of soil parameter, while input has values of "cell number of interest"
#ksat$cell <- setValues(ksat, 1:ncell(ksat)) #not sure if this is a valid thing to do, but this creates RasterBrick

##need to convert this raster to a dataframe to work with the data
elev3.df <- as.data.frame(elev3, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value
write.csv(elev3.df, "clay.soil.csv")
######################
#awc
#####################
awc <- raster("Data/paleon_awc.asc")
#define Albers projcetion for this layer
proj4string(awc) <- CRS('+init=epsg:3175')
plot(awc) #note that the extent of these data goes beyond statistical output
awc3 <- resample(awc, input, method='ngb') #resample ksat based on the psleon8km_midalbers1 grid created above
awc2 <- stack(awc3, input) #ksat has values of soil parameter, while input has values of "cell number of interest"
#ksat$cell <- setValues(ksat, 1:ncell(ksat)) #not sure if this is a valid thing to do, but this creates RasterBrick

##need to convert this raster to a dataframe to work with the data
awc3.df <- as.data.frame(awc3, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value
write.csv(awc3.df, "awc.soil.csv")


###################################
#Add soil to the biomass data grid#
###################################
biomass <- read.csv('biomass_v09_extract (1).csv')#this is gridded non-biomass model biomass
biomass$total <- rowSums(biomass[,5:33])
coordinates(biomass) <- ~x+y #assign coordinates to make biomass spatial points datframe
proj4string(biomass) <- CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
rast<-awc
Oak.biomass<- rasterize(biomass, rast, field = "Oak")
tot.biomass <- rasterize(biomass, rast, field = "total")
biomass$ksat <- extract(ksat3, biomass) # this extracts the ksat value for each biomass cell of the biomass dataframe
biomass$sand <- extract(sand3, biomass) 
biomass$silt <- extract(silt3, biomass) 
biomass$clay <- extract(clay3, biomass) 
biomass$awc <- extract(awc3, biomass) 
biomass$elev <- extract(elev3, biomass) 


spplot(biomass, "ksat")
spplot(biomass, "sand")
spplot(biomass, "silt")
spplot(biomass, "clay")
spplot(biomass, "awc")
spplot(biomass, "elev")
biomass.df <- data.frame(biomass) #export the data as a dataframe

#make a dataframe with only the soil data
soils.df <- biomass.df
soils<- c("x", "y", "sand", "silt", "clay", "awc", "ksat", "elev")
soils.df <- soils.df[,soils]

write.csv(soils.df, "StatsGO2.soil.covariates.csv")
