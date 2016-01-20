##This code reads the .asc soil rasters into R as rasters, plots them, and converts them to dataframes
library(raster)
library(sp)
library(rgdal)
#myShapeInR<-readOGR("C:/Users/Kelly/Documents/Eagleson/paleon-soils/merged_hydro_soils.shp", "merged_hydro_soils")
#see README.txt for descriptions of the data and citation for STATSGO2 soil data
#soils <- readOGR(dsn = "merged_hydro_soils/merged_hydro_soils.shp", layer = "merged_hydro_soils")
#simons code to create the 8km albers grid that I used to get the soil data on the 
#soils.il <- readOGR(dsn = "gsmsoilmu_a_il.shp", layer = "gsmsoilmu_a_il")


#same grid scale
input <- raster('Data/paleon.unit.alb_01.tif') #this file came from biomass dropbox
input <- setValues(input, 1:ncell(input))
writeRaster(input, "paleon8km_midalbers1.tif", overwrite=TRUE)

input.df <- as.data.frame(input, xy=TRUE)
#setwd('location of your soil directory')


##read the soil rasters into R
##example ksat
ksat <- raster("Data/ksat_raster1.tif")

#define Albers projcetion for this layer
ksat<-projectRaster(ksat, crs='+init=epsg:3175', method='bilinear')

plot(ksat) #note that the extent of these data goes beyond statistical output

##use resample to realign if the cells are off by a small amount
ksat3 <- resample(ksat, input, method='ngb') #resample ksat based on the psleon8km_midalbers1 grid created above
ksat2 <- stack(ksat3, input) #ksat has values of soil parameter, while input has values of "cell number of interest"
#ksat$cell <- setValues(ksat, 1:ncell(ksat)) #not sure if this is a valid thing to do, but this creates RasterBrick

##need to convert this raster to a dataframe to work with the data
ksat.df <- as.data.frame(ksat, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value
write.csv(ksat.df, "ksat.soil.csv") #write to csv so you can haave it by itself

##########################
#Rest of the soil data
##########################
#clay
clay <- raster("Data/clay_rast1.tif")
#project into Albers projcetion for this layer
clay<-projectRaster(clay, crs='+init=epsg:3175', method='bilinear')

plot(clay) #note that the extent of these data goes beyond statistical output
clay3 <- resample(clay, input, method='ngb') #resample ksat based on the psleon8km_midalbers1 grid created above
# to convert this raster to a dataframe to work with the data
clay.df <- as.data.frame(clay, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value
write.csv(clay.df, "clay.soil.csv") 

########################
#sand
#######################
sand <- raster("Data/sand_rast1.tif")
#define Albers projcetion for this layer
sand<-projectRaster(sand, crs='+init=epsg:3175', method='bilinear')

plot(sand) #note that the extent of these data goes beyond statistical output
sand3 <- resample(sand, input, method='bilinear') #resample ksat based on the psleon8km_midalbers1 grid created above
sand2 <- stack(sand3, input) #ksat has values of soil parameter, while input has values of "cell number of interest"
#ksat$cell <- setValues(ksat, 1:ncell(ksat)) #not sure if this is a valid thing to do, but this creates RasterBrick

##need to convert this raster to a dataframe to work with the data
sand.df <- as.data.frame(sand, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value
write.csv(sand.df, "sand.soil.csv")

#####################
#silt
#####################
silt <- raster("Data/silt_rast1.tif")
#define Albers projcetion for this layer
silt<-projectRaster(silt, crs='+init=epsg:3175', method='bilinear')

plot(silt) #note that the extent of these data goes beyond statistical output
silt3 <- resample(silt, input, method='bilinear') #resample ksat based on the psleon8km_midalbers1 grid created above
silt2 <- stack(silt, input) #ksat has values of soil parameter, while input has values of "cell number of interest"

##need to convert this raster to a dataframe to work with the data
silt.df <- as.data.frame(silt3, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value
write.csv(silt.df, "silt.soil.csv")

#######################
#elev
#######################
elev <- raster("Data/elev_rast1.tif")
#define Albers projcetion for this layer
elev<-projectRaster(elev, crs='+init=epsg:3175', method='bilinear')


plot(elev) #note that the extent of these data goes beyond statistical output
elev3 <- resample(elev, input, method='ngb') #resample ksat based on the psleon8km_midalbers1 grid created above
elev2 <- stack(elev3, input) #ksat has values of soil parameter, while input has values of "cell number of interest"
#ksat$cell <- setValues(ksat, 1:ncell(ksat)) #not sure if this is a valid thing to do, but this creates RasterBrick

##need to convert this raster to a dataframe to work with the data
elev.df <- as.data.frame(elev, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value
write.csv(elev.df, "elev.soil.csv")
######################
#awc
#####################
awc <- raster("Data/awc_rast1.tif")
#define Albers projcetion for this layer
awc<-projectRaster(awc, crs='+init=epsg:3175', method='bilinear')

plot(awc) #note that the extent of these data goes beyond statistical output
#awc <- resample(awc, input, method='ngb') #resample ksat based on the psleon8km_midalbers1 grid created above
#awc2= <- stack(awc3, input) #ksat has values of soil parameter, while input has values of "cell number of interest"

##need to convert this raster to a dataframe to work with the data
awc.df <- as.data.frame(awc, xy=TRUE, cellnumbers=TRUE) #adding xy=TRUE outputs the x-y albers coordinates with the data value
write.csv(awc.df, "awc.soil.csv")



###################################
#Modern Temp and precip covariates
###################################

#load rasters from IL_IN biomass directory
precip12 <- brick("C:/Users/Kelly/Documents/Indiana_Density_Biomass/precip12.grd")
precip13 <- brick("C:/Users/Kelly/Documents/Indiana_Density_Biomass/precip13.grd")

#tmean12 <- brick("C:/Users/Kelly/Documents/Indiana_Density_Biomass/tmean12.grd")
#tmean13 <- brick("C:/Users/Kelly/Documents/Indiana_Density_Biomass/tmean13.grd")


prec.12 <- extract(precip12, biomass)
total.pr12 <- as.vector(rowSums(prec.12[,1:12]))
pr12 <- total.pr12 #only works for half of data
prec.13 <- extract(precip13, biomass)
total.pr13 <- as.vector(rowSums(prec.13[,1:12]))
pr13 <- total.pr13
prs <- data.frame(pr13, pr12)

total <- rowSums(prs[,c("pr12", "pr13")], na.rm=TRUE)
jja <- c(6,7,8)
#precip13[[jja]]
djf <- c(12, 1, 2)
mam <- c(3, 4, 5)
son <- c(9, 10, 11)
jja.pr12 <- as.vector(rowSums(prec.12[,jja]))
djf.pr12 <- as.vector(rowSums(prec.12[,djf]))
mam.pr12 <- as.vector(rowSums(prec.12[,mam]))
son.pr12 <- as.vector(rowSums(prec.12[,son]))
total.pr12 <- as.vector(rowSums(prec.12[,1:12]))
#total.12 <- total.pr12

jja.pr13 <- as.vector(rowSums(prec.13[,jja]))
djf.pr13 <- as.vector(rowSums(prec.13[,djf]))
mam.pr13 <- as.vector(rowSums(prec.13[,mam]))
son.pr13 <- as.vector(rowSums(prec.13[,son]))
total.pr13 <- as.vector(rowSums(prec.13[,1:12]))

jja.prs <- rowSums(cbind(jja.pr12, jja.pr13), na.rm=TRUE)
djf.prs <-rowSums(cbind(djf.pr12, djf.pr13), na.rm=TRUE)
mam.prs <- rowSums(cbind(mam.pr12, mam.pr13), na.rm=TRUE)
son.prs <- rowSums(cbind(son.pr12, son.pr13), na.rm=TRUE)
total.prs <- rowSums(cbind(total.pr12, total.pr13), na.rm=TRUE)


biomass$jja_pr <- jja.prs
biomass$djf_pr <- djf.prs
biomass$mam_pr <- mam.prs
biomass$son_pr <- son.prs
biomass$total_pr <- total.prs

DEM <- raster("C:/Users/Kelly/Documents/Indiana_Density_Biomass/fullDEM.grd")


##########################
#Historical Precip & temp#
##########################

PT.dir<- "C:/Users/Kelly/Documents/biomodality"
precip.1900<- read.table("C:/Users/Kelly/Documents/biomodality/Global2011P.tar/Global2011P/precip.1900")
precip.1901<- read.table("C:/Users/Kelly/Documents/biomodality/Global2011P.tar/Global2011P/precip.1901")
precip.1902<- read.table("C:/Users/Kelly/Documents/biomodality/Global2011P.tar/Global2011P/precip.1902")
precip.1903<- read.table("C:/Users/Kelly/Documents/biomodality/Global2011P.tar/Global2011P/precip.1903")
precip.1904<- read.table("C:/Users/Kelly/Documents/biomodality/Global2011P.tar/Global2011P/precip.1904")
precip.1905<- read.table("C:/Users/Kelly/Documents/biomodality/Global2011P.tar/Global2011P/precip.1905")
precip.1906<- read.table("C:/Users/Kelly/Documents/biomodality/Global2011P.tar/Global2011P/precip.1906")
precip.1907<- read.table("C:/Users/Kelly/Documents/biomodality/Global2011P.tar/Global2011P/precip.1907")
precip.1908<- read.table("C:/Users/Kelly/Documents/biomodality/Global2011P.tar/Global2011P/precip.1908")
precip.1909<- read.table("C:/Users/Kelly/Documents/biomodality/Global2011P.tar/Global2011P/precip.1909")
precip.1910<- read.table("C:/Users/Kelly/Documents/biomodality/Global2011P.tar/Global2011P/precip.1910")
precip.2010<- read.table("C:/Users/Kelly/Documents/biomodality/Global2011P.tar/Global2011P/precip.2010")

Lat <- precip.1900[,2]
Long <- precip.1900[,1]

p.1900<- rowSums(precip.1900[,3:14], na.rm = TRUE)
p.1901<- rowSums(precip.1901[,3:14], na.rm = TRUE)
p.1902<- rowSums(precip.1902[,3:14], na.rm = TRUE)
p.1903<- rowSums(precip.1903[,3:14], na.rm = TRUE)
p.1904<- rowSums(precip.1904[,3:14], na.rm = TRUE)
p.1905<- rowSums(precip.1905[,3:14], na.rm = TRUE)
p.1906<- rowSums(precip.1906[,3:14], na.rm = TRUE)
p.1907<- rowSums(precip.1907[,3:14], na.rm = TRUE)
p.1908<- rowSums(precip.1908[,3:14], na.rm = TRUE)
p.1909<- rowSums(precip.1909[,3:14], na.rm = TRUE)
p.1910<- rowSums(precip.1910[,3:14], na.rm = TRUE)

avg.p<- rowMeans(data.frame(p.1900,
                            p.1901, 
                            p.1902,
                            p.1903, 
                            p.1904, 
                            p.1905, 
                            p.1906, 
                            p.1907,
                            p.1908,
                            p.1909, 
                            p.1910))/11

averages <- data.frame(Lat = Lat, 
                       Long = Long, 
                       avg = avg.p)

coordinates(averages) <- ~Long+Lat
gridded(averages) <- TRUE
avg.rast <- raster(averages)
projection(avg.rast) <- CRS("+init=epsg:4326")

avg.alb <- projectRaster(avg.rast, crs='+init=epsg:3175')
#extract total pr from AVG.alb
total.pr1910 <- extract(avg.alb, df.with.temp[,2:3] )

#for jja
p.1900<- rowSums(precip.1900[,8:10], na.rm = TRUE)
p.1901<- rowSums(precip.1901[,8:10], na.rm = TRUE)
p.1902<- rowSums(precip.1902[,8:10], na.rm = TRUE)
p.1903<- rowSums(precip.1903[,8:10], na.rm = TRUE)
p.1904<- rowSums(precip.1904[,8:10], na.rm = TRUE)
p.1905<- rowSums(precip.1905[,8:10], na.rm = TRUE)
p.1906<- rowSums(precip.1906[,8:10], na.rm = TRUE)
p.1907<- rowSums(precip.1907[,8:10], na.rm = TRUE)
p.1908<- rowSums(precip.1908[,8:10], na.rm = TRUE)
p.1909<- rowSums(precip.1909[,8:10], na.rm = TRUE)
p.1910<- rowSums(precip.1910[,8:10], na.rm = TRUE)

avg.pjja<- rowMeans(data.frame(p.1900,
                            p.1901, 
                            p.1902,
                            p.1903, 
                            p.1904, 
                            p.1905, 
                            p.1906, 
                            p.1907,
                            p.1908,
                            p.1909, 
                            p.1910))/11
averages <- data.frame(Lat = Lat, 
                       Long = Long, 
                       avg = avg.pjja)

coordinates(averages) <- ~Long+Lat
gridded(averages) <- TRUE
avg.rast <- raster(averages)
projection(avg.rast) <- CRS("+init=epsg:4326")

avg.alb <- projectRaster(avg.rast, crs='+init=epsg:3175')
#extract total pr from AVG.alb
jja.pr1910 <- extract(avg.alb, df.with.temp[,2:3] )



##not sure how useful these are, since they are on a fairly oarse grid
#have temperature and topography data here:
df.with.temp <- read.csv("C:/Users/Kelly/Documents/Indiana_Density_Biomass/biomass_full_extract_with_temp_pr_elev.csv")
df.temp<- df.with.temp[which(df.with.temp$cell %in% biom.export$cell),]
df.temp.cov <- data.frame(x= df.temp$x, 
                          y = df.temp$y, 
                          cell = df.temp$cell, 
                          df.temp[,54:62])

df.temp.cov<- df.temp.cov[order(df.temp.cov$cell),]
biom.t<- biom.t[order(biom.t$cell),]
#df.temp.cov1<- df.temp.cov[unique(df.temp.cov$cell),]
#biom.t1<- biom.t[unique(biom.t$cell),]


biom.export <- data.frame(biomass)
write.csv(biom.export, "biomass_and_covs.csv")
biom.t<- biom.export[which(biom.export$cell %in% df.temp.cov$cell),]
biom.merged<- merge( df.temp.cov, biom.t, 
                     by.x = "cell", by.y = "cell")

biom.merged1<- biom.merged[unique(biom.merged$cell),]

covars.export <- data.frame(x=biom.t$x, 
                            y=biom.t$y, 
                            cell = biom.t$cell, 
                            biom.t[,46:56],
                            df.temp.cov[,4:12])
ggplot(na.omit(covars.export), aes(x=x, y=y,colour=jja_tmean))+ 
  geom_point(shape=15, solid=TRUE,cex=1.5)
write.csv(covars.export, "C://Users/Kelly/Documents/BIOMASS/biomass_UW/IL_IN_covariates//data_v2/biomass_and_covs.csv")
###################################
#Add soil to the biomass data grid#
###################################
#read in the upper midwest data
#biomass <- read.csv('Data/biomass_prediction_v0.2-1.csv')#this is gridded non-biomass model biomass
#biomass$total <- rowSums(biomass[,4:25], na.rm = TRUE)
coordinates(biomass) <- ~x+y #assign coordinates to make biomass spatial points datframe
proj4string(biomass) <- CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')



#read in the indiana and illinois data
#biomass.inil<- read.csv('C:/Users/Kelly/Documents/Indiana_Density_Biomass/Data/outputsplss_biomass_indiana_illinois_v1.csv')
#biomass.combined<- rbind(biomass, biomass.inil)
#biomass.inil$total <- rowSums(biomass.inil[,5:36], na.rm = TRUE)
#coordinates(biomass.inil) <- ~x+y #assign coordinates to make biomass spatial points datframe
#proj4string(biomass.inil) <- CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')

#made a between way of loading in the full dataset
biomass.midwest <- read.csv('C:/Users/Kelly/Documents/Indiana_Density_Biomass/Data/outputsplss_biomass_full_midwest_v1.csv')
biomass.midwest$total <- rowSums(biomass.midwest[,4:40], na.rm = TRUE)
coordinates(biomass.midwest) <- ~x+y #assign coordinates to make biomass spatial points datframe
proj4string(biomass.midwest) <- CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')

biomass<- biomass.midwest

b <- read.csv('C:/Users/Kelly/Documents/BIOMASS/biomass_UW/IL_IN_covariates/partial_midwest_alb_biomass_v2.csv')
#b<- read.csv( "C:/Users/Kelly/Documents/BIOMASS/biomass_UW/IL_IN_covariates/inil_alb_biomass_v2.csv")

b$total <- rowSums(b[,-c(1,2,3,4,5,40,41,42)], na.rm = TRUE)
TBD_i = c("Alder", "Ash","Bald.cypress", "Basswood", "Beech", 
          "Birch", "Black.gum", "Black.gum.sweet.gum", "Buckeye",
          "Cherry","Chestnut", "Dogwood", "Elm", "Hackberry",
          "Other.hardwood","Hickory","Ironwood", "Locust","Maple",
          "Mulberry",
          "Oak","Poplar","Poplar.tulip.poplar", "Sweet.gum",
          "Sycamore","Tulip.poplar",
          "Walnut", "Willow")
TNE_i = c("Cedar.juniper","Tamarack","Bald.cypress", "Fir", "Hemlock", "Pine", "Spruce")

b$TBD = rowSums(b[,TBD_i])
b$TNE = rowSums(b[,TNE_i])
b$totwood = rowSums(b[,c(TNE_i, TBD_i)])
biomass.midwest <- b

coordinates(biomass.midwest) <- ~x+y #assign coordinates to make biomass spatial points datframe
proj4string(biomass.midwest) <- CRS('+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
library(corrgram)

#merge the biomass rasters
#test <- merge(biomass, biomass.inil, overlap=FALSE)
#rast<-awc
#Oak.biomass<- rasterize(biomass, rast, field = "Oak")
#tot.biomass <- rasterize(biomass, rast, field = "total")
#Oak.biomass.il<-rasterize(biomass.inil, rast, field="Oak")
#tot.biomass.il <- rasterize(biomass.inil, rast, field="total")
#biomass.tot<- biomass$total
biomass <- biomass.midwest

dens.all.test<- read.csv("full_density_alb.csv")

ksat.bio<-  extract(ksat, biomass.midwest, fun= mean, weights=TRUE, na.rm=TRUE) # this extracts the ksat value for each biomass cell of the bi
#ksat.bio[is.na(ksat.bio)]<- -9999
#ksat.dens <- extract(ksat3, density)
biomass$ksat <- ksat.bio
#biomass$ksat <- extract(ksat3, biomass, fun= mean, weights=TRUE, na.rm=TRUE) # this extracts the ksat value for each biomass cell of the biomass dataframe
biomass.sand <- extract(sand, biomass) 
#biomass.sand[is.na(biomass.sand)] <- -9999
biomass$sand <- biomass.sand

biomass.silt <- extract(silt, biomass) 
#biomass.silt[is.na(biomass.silt)] <- -9999
biomass$silt <- biomass.silt

biomass.clay <- extract(clay, biomass) 
#biomass.clay[is.na(biomass.clay)] <- -9999
biomass$clay <- biomass.clay

biomass.awc <- extract(awc, biomass) 
#biomass.awc[is.na(biomass.awc)] <- -9999
biomass$awc <- biomass.awc

biomass.elev <- extract(elev, biomass) 
#biomass.elev[is.na(biomass.elev)] <- -9999
biomass$elev <- biomass.elev

biomass$jja_pr <- jja.prs
biomass$djf_pr <- djf.prs
biomass$mam_pr <- mam.prs
biomass$son_pr <- son.prs
biomass$total_pr <- total.prs


test <- merge(dens.all.test, biomass.df, by="cell")

write.csv(test, "C://Users/Kelly/Documents/BIOMASS/biomass_UW/IL_IN_covariates//data_v2/density_biomass.csv")
test<- read.csv("C://Users/Kelly/Documents/BIOMASS/biomass_UW/IL_IN_covariates//data_v2/density_biomass.csv")

#basic total density plots
plot(test$jja_pr, test$total.x)
plot(test$total_pr, test$total.x)
plot(test$son_pr, test$total.x)
plot(test$djf_pr, test$total.x)
plot(test$mam_pr, test$total.x)
plot(test$sand, test$total.x)
plot(test$silt, test$total.x)
plot(test$clay, test$total.x)
plot(test$awc, test$total.x)
plot(test$ksat, test$total.x)
plot(test$elev, test$total.x)

#Oak density plots

plot(test$jja_pr, test$Oak.x)
plot(test$total_pr, test$Oak.x)
plot(test$son_pr, test$Oak.x)
plot(test$djf_pr, test$Oak.x)
plot(test$mam_pr, test$Oak.x)
plot(test$sand, test$Oak.x)
plot(test$silt, test$Oak.x)
plot(test$clay, test$Oak.x)
plot(test$awc, test$Oak.x)
plot(test$ksat, test$Oak.x)
plot(test$elev, test$Oak.x)

#BEech
plot(test$jja_pr, test$Beech.x)
plot(test$total_pr, test$Beech.x)
plot(test$son_pr, test$Beech.x)
plot(test$djf_pr, test$Beech.x)
plot(test$mam_pr, test$Beech.x)
plot(test$sand, test$Beech.x)
plot(test$silt, test$Beech.x)
plot(test$clay, test$Beech.x)
plot(test$awc, test$Beech.x)
plot(test$ksat, test$Beech.x)
plot(test$elev, test$Beech.x)

#Maple
plot(test$jja_pr, test$Maple.x)
plot(test$total_pr, test$Maple.x)
plot(test$son_pr, test$Maple.x)
plot(test$djf_pr, test$Maple.x)
plot(test$mam_pr, test$Maple.x)
plot(test$sand, test$Maple.x)
plot(test$silt, test$Maple.x)
plot(test$clay, test$Maple.x)
plot(test$awc, test$Maple.x)
plot(test$ksat, test$Maple.x)
plot(test$elev, test$Maple.x)

##subset by state and plot
indiana.f <- test[which(test$type == "Indiana forest"),]
indiana.p <- test[which(test$type == "Indiana prairie"),]
indiana <- test[which(test$states == "Indiana"),]


minnesota.f <- test[which(test$type == "Minnesota forest"),]
minnesota.p <- test[which(test$type == "Minnesota prairie"),]
minnesota <- test[which(test$states == "Minnesota"),]

wisconsin.f <- test[which(test$type == "Wisconsin forest"),]
wisconsin.p <- test[which(test$type == "Wisconsin prairie"),]
wisconsin <- test[which(test$states == "Wisconsin"),]

michigan.f <- test[which(test$type == "Michigan forest"),]
michigan.p <- test[which(test$type == "Michigan prairie"),]
michigan <- test[which(test$states == "Michigan"),]


illinois.f <- test[which(test$type == "Illinois forest"),]
illinois.p <- test[which(test$type == "Illinois prairie"),]
illinois <- test[which(test$states == "Illinois"),]

pdf("plots_nov_29.pdf")
par(mfrow=c(3,3))
#make state prairie vs. forest plots
plot(indiana.f$total_pr, indiana.f$total.x,
     xlim = c(0, 1600), ylim = c(0,600))
total_prinf <- lm(total.x ~ total_pr, data = indiana.f)
abline(total_prinf, col="red")
points(indiana.p$total_pr, indiana.p$total.x, col = "red")
total_prinp <- lm(total.x ~ total_pr, data = indiana.p)
abline(total_prinp, col="red")

plot(illinois.f$total_pr, illinois.f$total.x,
     xlim = c(700, 900), ylim = c(0,600))
total_prmif <- lm(total.x ~ total_pr, data = illinois.f)
abline(total_prmif, col="blue")
points(illinois.p$total_pr, illinois.p$total.x, col = "red")
total_prilp <- lm(total.x ~ total_pr, data = illinois.p)
abline(total_prilp, col="red")

plot(minnesota.f$total_pr, minnesota.f$total.x,
     xlim = c(700, 900), ylim = c(0,600))
total_prmnf <- lm(total.x ~ total_pr, data = minnesota.f)
abline(total_prmnf, col="blue")
points(minnesota.p$total_pr, minnesota.p$total.x, col = "red")
total_prmnp <- lm(total.x ~ total_pr, data = minnesota.p)
abline(total_prmnp, col="red")

plot(wisconsin.f$total_pr, wisconsin.f$total.x,
     xlim = c(700, 900), ylim = c(0,600))
total_prwif <- lm(total.x ~ total_pr, data = wisconsin.f)
abline(total_prwif, col="blue")
points(wisconsin.p$total_pr, wisconsin.p$total.x, col = "red")
total_prwip <- lm(total.x ~ total_pr, data = wisconsin.p)
abline(total_prwip, col="red")

plot(michigan.f$total_pr, michigan.f$total.x, 
     xlim = c(700, 900), ylim = c(0,600))
total_prmif <- lm(total.x ~ total_pr, data = michigan.f)
abline(total_prmif, col="blue")
points(michigan.p$total_pr, michigan.p$total.x, col = "red")
total_prmip <- lm(total.x ~ total_pr, data = michigan.p)
abline(total_prmip, col="red")

#just statewide plots
#for total_pr
plot(indiana$total_pr, indiana$total.x)
total_prlm1 <- lm(total.x ~ total_pr, data=indiana)
abline(total_prlm1, col="red")

plot(illinois$total_pr, illinois$total.x)
total_prlm2 <- lm(total.x ~ total_pr, data=illinois)
abline(total_prlm2, col="red")

plot(minnesota$total_pr, minnesota$total.x)
total_prlm3 <- lm(total.x ~ total_pr, data=minnesota)
abline(total_prlm3, col="red")

plot(wisconsin$total_pr, wisconsin$total.x)
total_prlm4 <- lm(total.x ~ total_pr, data=wisconsin)
abline(total_prlm4, col="red")

plot(michigan$total_pr, michigan$total.x)
total_prlm5 <- lm(total.x ~ total_pr, data=michigan)
abline(total_prlm5, col="red")

plot(test$total_pr, test$total.x)
total_prlmall <- lm(total.x ~ total_pr, data=test)
abline(total_prlmall, col="red")

##jja_pr
plot(indiana$jja_pr, indiana$total.x)
jja_lm1 <- lm(total.x ~jja_pr, data = indiana)
abline(jja_prlm1, col="red")

plot(illinois$jja_pr, illinois$total.x)
jja_prlm2 <- lm(total.x ~ jja_pr, data = illinois)
abline(jja_prlm2, col="red")

plot(minnesota$jja_pr, minnesota$total.x)
jja_prlm3 <- lm(total.x ~ jja_pr, data=minnesota)
abline(jja_prlm3, col="red")

plot(wisconsin$jja_pr, wisconsin$total.x)
jja_prlm4 <- lm(total.x ~ jja_pr, data=wisconsin)
abline(jja_prlm4, col="red")

plot(michigan$jja_pr, michigan$total.x)
jja_prlm5 <- lm(total.x ~ jja_pr, data=michigan)
abline(jja_prlm5, col="red")

plot(test$jja_pr, test$total.x)
jja_prlmall <- lm(total.x ~ jja_pr, data=test)
abline(jja_prlmall, col="red")

#son_pr
plot(indiana$son_pr, indiana$total.x)
son_prlm1 <- lm(total.x ~ son_pr, data=indiana)
abline(son_prlm1, col="red")

plot(illinois$son_pr, illinois$total.x)
son_prlm2 <- lm(total.x ~ son_pr, data=illinois)
abline(son_prlm2, col="red")

plot(minnesota$son_pr, minnesota$total.x)
son_prlm3 <- lm(total.x ~ son_pr, data=minnesota)
abline(son_prlm3, col="red")

plot(wisconsin$son_pr, wisconsin$total.x)
son_prlm4 <- lm(total.x ~ son_pr, data=wisconsin)
abline(son_prlm4, col="red")

plot(michigan$son_pr, michigan$total.x)
son_prlm5 <- lm(total.x ~ son_pr, data=michigan)
abline(son_prlm5, col="red")

plot(test$son_pr, test$total.x)
son_prlmall <- lm(total.x ~ son_pr, data=test)
abline(son_prlmall, col="red")

#mam pr
plot(indiana$mam_pr, indiana$total.x)
mam_prlm1 <- lm(total.x ~ mam_pr, data=indiana)
abline(mam_prlm1, col="red")

plot(illinois$mam_pr, illinois$total.x)
mam_prlm2 <- lm(total.x ~ mam_pr, data=illinois)
abline(mam_prlm2, col="red")

plot(minnesota$mam_pr, minnesota$total.x)
mam_prlm3 <- lm(total.x ~ mam_pr, data=minnesota)
abline(mam_prlm3, col="red")

plot(wisconsin$mam_pr, wisconsin$total.x)
mam_prlm4 <- lm(total.x ~ mam_pr, data=wisconsin)
abline(mam_prlm4, col="red")

plot(michigan$mam_pr, michigan$total.x)
mam_prlm5 <- lm(total.x ~ mam_pr, data=michigan)
abline(mam_prlm5, col="red")

plot(test$mam_pr, test$total.x)
mam_prlmall <- lm(total.x ~ mam_pr, data=test)
abline(mam_prlmall, col="red")

#elevation
plot(indiana$elev, indiana$total.x)
elevlm1 <- lm(total.x ~ elev, data=indiana)
abline(elevlm1, col="red")

plot(illinois$elev, illinois$total.x)
elevlm2 <- lm(total.x ~ elev, data=illinois)
abline(elevlm2, col="red")

plot(minnesota$elev, minnesota$total.x)
elevlm3 <- lm(total.x ~ elev, data=minnesota)
abline(elevlm3, col="red")

plot(wisconsin$elev, wisconsin$total.x)
elevlm4 <- lm(total.x ~ elev, data=wisconsin)
abline(elevlm4, col="red")

plot(michigan$elev, michigan$total.x)
elevlm5 <- lm(total.x ~ elev, data=michigan)
abline(elevlm5, col="red")

plot(test$elev, test$total.x)
elevlmall <- lm(total.x ~ elev, data=test)
abline(elevlmall, col="red")

#sand
plot(indiana$sand, indiana$total.x)
sandlm1 <- lm(total.x ~ elev, data=indiana)
abline(sandlm1, col="red")

plot(illinois$sand, illinois$total.x)
sandlm2 <- lm(total.x ~ elev, data=illinois)
abline(sandlm2, col="red")

plot(minnesota$sand, minnesota$total.x)
sandlm3 <- lm(total.x ~ elev, data=minnesota)
abline(sandlm3, col="red")

plot(wisconsin$sand, wisconsin$total.x)
sandlm4 <- lm(total.x ~ elev, data=wisconsin)
abline(sandlm4, col="red")

plot(michigan$sand, michigan$total.x)
sandlm5 <- lm(total.x ~ sand, data=michigan)
abline(sandlm5, col="red")

plot(test$sand, test$total.x)
sandlmall <- lm(total.x ~ sand, data=test)
abline(sandlmall, col="red")

##silt
plot(indiana$silt, indiana$total.x)
siltlm1 <- lm(total.x ~ silt, data=indiana)
abline(siltlm1, col="red")

plot(illinois$silt, illinois$total.x)
siltlm2 <- lm(total.x ~ silt, data=illinois)
abline(siltlm2, col="red")

plot(minnesota$silt, minnesota$total.x)
siltlm3 <- lm(total.x ~ silt, data=minnesota)
abline(siltlm3, col="red")

plot(wisconsin$silt, wisconsin$total.x)
siltlm4 <- lm(total.x ~ silt, data=wisconsin)
abline(siltlm4, col="red")

plot(michigan$silt, michigan$total.x)
siltlm5 <- lm(total.x ~ silt, data=michigan)
abline(siltlm5, col="red")

plot(test$silt, test$total.x)
siltlmall <- lm(total.x ~ silt, data=test)
abline(siltlmall, col="red")

##clay
plot(indiana$clay, indiana$total.x)
claylm1 <- lm(total.x ~ clay, data=indiana)
abline(claylm1, col="red")

plot(illinois$clay, illinois$total.x)
claylm2 <- lm(total.x ~ clay, data=illinois)
abline(claylm2, col="red")

plot(minnesota$clay, minnesota$total.x)
claylm3 <- lm(total.x ~ clay, data=minnesota)
abline(claylm3, col="red")

plot(wisconsin$clay, wisconsin$total.x)
claylm4 <- lm(total.x ~ clay, data=wisconsin)
abline(claylm4, col="red")

plot(michigan$clay, michigan$total.x)
claylm5 <- lm(total.x ~ clay, data=michigan)
abline(claylm5, col="red")

plot(test$clay, test$total.x)
claylmall <- lm(total.x ~ clay, data=test)
abline(claylmall, col="red")

#for awc
plot(indiana$awc, indiana$total.x)
awclm1 <- lm(total.x ~ awc, data=indiana)
abline(awclm1, col="red")

plot(illinois$awc, illinois$total.x)
awclm2 <- lm(total.x ~ awc, data=illinois)
abline(awclm2, col="red")

plot(minnesota$awc, minnesota$total.x)
awclm3 <- lm(total.x ~ awc, data=minnesota)
abline(awclm3, col="red")

plot(wisconsin$awc, wisconsin$total.x)
awclm4 <- lm(total.x ~ awc, data=wisconsin)
abline(awclm4, col="red")

plot(michigan$awc, michigan$total.x)
awclm5 <- lm(total.x ~ awc, data=michigan)
abline(awclm5, col="red")


plot(test$awc, test$total.x)
awclmall <- lm(total.x ~ awc, data=test)
abline(awclmall, col="red")

#now for ksat
plot(indiana$ksat, indiana$total.x)
ksatlm1 <- lm(total.x ~ ksat, data=indiana)
abline(ksatlm1, col="red")

plot(illinois$ksat, illinois$total.x)
ksatlm2 <- lm(total.x ~ ksat, data=illinois)
abline(ksatlm2, col="red")

plot(minnesota$ksat, minnesota$total.x)
ksatlm3 <- lm(total.x ~ ksat, data=minnesota)
abline(ksatlm3, col="red")

plot(wisconsin$ksat, wisconsin$total.x)
ksatlm4 <- lm(total.x ~ ksat, data=wisconsin)
abline(ksatlm4, col="red")

plot(michigan$ksat, michigan$total.x)
ksatlm5 <- lm(total.x ~ ksat, data=michigan)
abline(ksatlm5, col="red")

plot(test$ksat, test$total.x)
ksatlmall <- lm(total.x ~ ksat, data=test)
abline(ksatlmall, col="red")

dev.off()



#now for total biomass
pdf("biomass_plots_nov_28.pdf")
par(mfrow=c(3,3))
#make state prairie vs. forest plots
plot(indiana.f$total_pr, indiana.f$total.y,
     xlim = c(0, 1600), ylim = c(0,600))
total_prinf <- lm(total.y ~ total_pr, data = indiana.f)
abline(total_prinf, col="red")
points(indiana.p$total_pr, indiana.p$total.y, col = "red")
total_prinp <- lm(total.y ~ total_pr, data = indiana.p)
abline(total_prinp, col="red")

plot(illinois.f$total_pr, illinois.f$total.y,
     xlim = c(700, 900), ylim = c(0,600))
total_prmif <- lm(total.y ~ total_pr, data = illinois.f)
abline(total_prmif, col="blue")
points(illinois.p$total_pr, illinois.p$total.y, col = "red")
total_prilp <- lm(total.y ~ total_pr, data = illinois.p)
abline(total_prilp, col="red")

plot(minnesota.f$total_pr, minnesota.f$total.y,
     xlim = c(700, 900), ylim = c(0,600))
total_prmnf <- lm(total.y ~ total_pr, data = minnesota.f)
abline(total_prmnf, col="blue")
points(minnesota.p$total_pr, minnesota.p$total.y, col = "red")
total_prmnp <- lm(total.y ~ total_pr, data = minnesota.p)
abline(total_prmnp, col="red")

plot(wisconsin.f$total_pr, wisconsin.f$total.y,
     xlim = c(700, 900), ylim = c(0,600))
total_prwif <- lm(total.y ~ total_pr, data = wisconsin.f)
abline(total_prwif, col="blue")
points(wisconsin.p$total_pr, wisconsin.p$total.y, col = "red")
total_prwip <- lm(total.y ~ total_pr, data = wisconsin.p)
abline(total_prwip, col="red")

plot(michigan.f$total_pr, michigan.f$total.y, 
     xlim = c(700, 900), ylim = c(0,600))
total_prmif <- lm(total.y ~ total_pr, data = michigan.f)
abline(total_prmif, col="blue")
points(michigan.p$total_pr, michigan.p$total.y, col = "red")
total_prmip <- lm(total.y ~ total_pr, data = michigan.p)
abline(total_prmip, col="red")

#just statewide plots
#for total_pr
plot(indiana$total_pr, indiana$total.y)
total_prlm1 <- lm(total.y ~ total_pr, data=indiana)
abline(total_prlm1, col="red")

plot(illinois$total_pr, illinois$total.y)
total_prlm2 <- lm(total.y ~ total_pr, data=illinois)
abline(total_prlm2, col="red")

plot(minnesota$total_pr, minnesota$total.y)
total_prlm3 <- lm(total.y ~ total_pr, data=minnesota)
abline(total_prlm3, col="red")

plot(wisconsin$total_pr, wisconsin$total.y)
total_prlm4 <- lm(total.y ~ total_pr, data=wisconsin)
abline(total_prlm4, col="red")

plot(michigan$total_pr, michigan$total.y)
total_prlm5 <- lm(total.y ~ total_pr, data=michigan)
abline(total_prlm5, col="red")

plot(test$total_pr, test$total.y)
total_prlmall <- lm(total.y ~ total_pr, data=test)
abline(total_prlmall, col="red")

##jja_pr
plot(indiana$jja_pr, indiana$total.y)
jja_lm1 <- lm(total.y ~jja_pr, data = indiana)
abline(jja_prlm1, col="red")

plot(illinois$jja_pr, illinois$total.y)
jja_prlm2 <- lm(total.y ~ jja_pr, data = illinois)
abline(jja_prlm2, col="red")

plot(minnesota$jja_pr, minnesota$total.y)
jja_prlm3 <- lm(total.y ~ jja_pr, data=minnesota)
abline(jja_prlm3, col="red")

plot(wisconsin$jja_pr, wisconsin$total.y)
jja_prlm4 <- lm(total.y ~ jja_pr, data=wisconsin)
abline(jja_prlm4, col="red")

plot(michigan$jja_pr, michigan$total.y)
jja_prlm5 <- lm(total.y ~ jja_pr, data=michigan)
abline(jja_prlm5, col="red")

plot(test$jja_pr, test$total.y)
jja_prlmall <- lm(total.y ~ jja_pr, data=test)
abline(jja_prlmall, col="red")

#son_pr
plot(indiana$son_pr, indiana$total.y)
son_prlm1 <- lm(total.y ~ son_pr, data=indiana)
abline(son_prlm1, col="red")

plot(illinois$son_pr, illinois$total.y)
son_prlm2 <- lm(total.y ~ son_pr, data=illinois)
abline(son_prlm2, col="red")

plot(minnesota$son_pr, minnesota$total.y)
son_prlm3 <- lm(total.y ~ son_pr, data=minnesota)
abline(son_prlm3, col="red")

plot(wisconsin$son_pr, wisconsin$total.y)
son_prlm4 <- lm(total.y ~ son_pr, data=wisconsin)
abline(son_prlm4, col="red")

plot(michigan$son_pr, michigan$total.y)
son_prlm5 <- lm(total.y ~ son_pr, data=michigan)
abline(son_prlm5, col="red")

plot(test$son_pr, test$total.y)
son_prlmall <- lm(total.y ~ son_pr, data=test)
abline(son_prlmall, col="red")

#mam pr
plot(indiana$mam_pr, indiana$total.y)
mam_prlm1 <- lm(total.y ~ mam_pr, data=indiana)
abline(mam_prlm1, col="red")

plot(illinois$mam_pr, illinois$total.y)
mam_prlm2 <- lm(total.y ~ mam_pr, data=illinois)
abline(mam_prlm2, col="red")

plot(minnesota$mam_pr, minnesota$total.y)
mam_prlm3 <- lm(total.y ~ mam_pr, data=minnesota)
abline(mam_prlm3, col="red")

plot(wisconsin$mam_pr, wisconsin$total.y)
mam_prlm4 <- lm(total.y ~ mam_pr, data=wisconsin)
abline(mam_prlm4, col="red")

plot(michigan$mam_pr, michigan$total.y)
mam_prlm5 <- lm(total.y ~ mam_pr, data=michigan)
abline(mam_prlm5, col="red")

plot(test$mam_pr, test$total.y)
mam_prlmall <- lm(total.y ~ mam_pr, data=test)
abline(mam_prlmall, col="red")

#elevation
plot(indiana$elev, indiana$total.y)
elevlm1 <- lm(total.y ~ elev, data=indiana)
abline(elevlm1, col="red")

plot(illinois$elev, illinois$total.y)
elevlm2 <- lm(total.y ~ elev, data=illinois)
abline(elevlm2, col="red")

plot(minnesota$elev, minnesota$total.y)
elevlm3 <- lm(total.y ~ elev, data=minnesota)
abline(elevlm3, col="red")

plot(wisconsin$elev, wisconsin$total.y)
elevlm4 <- lm(total.y ~ elev, data=wisconsin)
abline(elevlm4, col="red")

plot(michigan$elev, michigan$total.y)
elevlm5 <- lm(total.y ~ elev, data=michigan)
abline(elevlm5, col="red")

plot(test$elev, test$total.y)
elevlmall <- lm(total.y ~ elev, data=test)
abline(elevlmall, col="red")

#sand
plot(indiana$sand, indiana$total.y)
sandlm1 <- lm(total.y ~ elev, data=indiana)
abline(sandlm1, col="red")

plot(illinois$sand, illinois$total.y)
sandlm2 <- lm(total.y ~ elev, data=illinois)
abline(sandlm2, col="red")

plot(minnesota$sand, minnesota$total.y)
sandlm3 <- lm(total.y ~ elev, data=minnesota)
abline(sandlm3, col="red")

plot(wisconsin$sand, wisconsin$total.y)
sandlm4 <- lm(total.y ~ elev, data=wisconsin)
abline(sandlm4, col="red")

plot(michigan$sand, michigan$total.y)
sandlm5 <- lm(total.y ~ sand, data=michigan)
abline(sandlm5, col="red")

plot(test$sand, test$total.y)
sandlmall <- lm(total.y ~ sand, data=test)
abline(sandlmall, col="red")

##silt
plot(indiana$silt, indiana$total.y)
siltlm1 <- lm(total.y ~ silt, data=indiana)
abline(siltlm1, col="red")

plot(illinois$silt, illinois$total.y)
siltlm2 <- lm(total.y ~ silt, data=illinois)
abline(siltlm2, col="red")

plot(minnesota$silt, minnesota$total.y)
siltlm3 <- lm(total.y ~ silt, data=minnesota)
abline(siltlm3, col="red")

plot(wisconsin$silt, wisconsin$total.y)
siltlm4 <- lm(total.y ~ silt, data=wisconsin)
abline(siltlm4, col="red")

plot(michigan$silt, michigan$total.y)
siltlm5 <- lm(total.y ~ silt, data=michigan)
abline(siltlm5, col="red")

plot(test$silt, test$total.y)
siltlmall <- lm(total.y ~ silt, data=test)
abline(siltlmall, col="red")

##clay
plot(indiana$clay, indiana$total.y)
claylm1 <- lm(total.y ~ clay, data=indiana)
abline(claylm1, col="red")

plot(illinois$clay, illinois$total.y)
claylm2 <- lm(total.y ~ clay, data=illinois)
abline(claylm2, col="red")

plot(minnesota$clay, minnesota$total.y)
claylm3 <- lm(total.y ~ clay, data=minnesota)
abline(claylm3, col="red")

plot(wisconsin$clay, wisconsin$total.y)
claylm4 <- lm(total.y ~ clay, data=wisconsin)
abline(claylm4, col="red")

plot(michigan$clay, michigan$total.y)
claylm5 <- lm(total.y ~ clay, data=michigan)
abline(claylm5, col="red")

plot(test$clay, test$total.y)
claylmall <- lm(total.y ~ clay, data=test)
abline(claylmall, col="red")

#for awc
plot(indiana$awc, indiana$total.y)
awclm1 <- lm(total.y ~ awc, data=indiana)
abline(awclm1, col="red")

plot(illinois$awc, illinois$total.y)
awclm2 <- lm(total.y ~ awc, data=illinois)
abline(awclm2, col="red")

plot(minnesota$awc, minnesota$total.y)
awclm3 <- lm(total.y ~ awc, data=minnesota)
abline(awclm3, col="red")

plot(wisconsin$awc, wisconsin$total.y)
awclm4 <- lm(total.y ~ awc, data=wisconsin)
abline(awclm4, col="red")

plot(michigan$awc, michigan$total.y)
awclm5 <- lm(total.y ~ awc, data=michigan)
abline(awclm5, col="red")


plot(test$awc, test$total.y)
awclmall <- lm(total.y ~ awc, data=test)
abline(awclmall, col="red")

#now for ksat
plot(indiana$ksat, indiana$total.y)
ksatlm1 <- lm(total.y ~ ksat, data=indiana)
abline(ksatlm1, col="red")

plot(illinois$ksat, illinois$total.y)
ksatlm2 <- lm(total.y ~ ksat, data=illinois)
abline(ksatlm2, col="red")

plot(minnesota$ksat, minnesota$total.y)
ksatlm3 <- lm(total.y ~ ksat, data=minnesota)
abline(ksatlm3, col="red")

plot(wisconsin$ksat, wisconsin$total.y)
ksatlm4 <- lm(total.y ~ ksat, data=wisconsin)
abline(ksatlm4, col="red")

plot(michigan$ksat, michigan$total.y)
ksatlm5 <- lm(total.y ~ ksat, data=michigan)
abline(ksatlm5, col="red")

plot(test$ksat, test$total.y)
ksatlmall <- lm(total.y ~ ksat, data=test)
abline(ksatlmall, col="red")

dev.off()




#lets compare the pls data across the transects










#corrgram(biomass)

library(maps)
library(rgdal)
#load us map data
library(maps)

#data from natural earth website
nat.earth <- readOGR(dsn='C:/Users/Kelly/Documents/BIOMASS/biomass_UW/IL_IN_covariates/data/ne_50m_admin_1_states_provinces_lakes/ne_50m_admin_1_states_provinces_lakes.shp', layer='ne_50m_admin_1_states_provinces_lakes')
nat.earth.alb <- spTransform(x = nat.earth, CRSobj = CRS('+init=epsg:3175'))
inil<- c("Minnesota", "Wisconsin", "Michigan", "Illinois", "Indiana")
#inil <- c("Illinois", "Indiana")
mw.alb = nat.earth.alb[match(toupper(inil), toupper(nat.earth.alb$name)),]
mw.alb <- spTransform(mw.alb,CRSobj = CRS('+init=epsg:3175'))

biomass <- spTransform(biomass,CRSobj = CRS('+init=epsg:3175'))






biomass$states <- as.character(over(biomass,mw.alb)$name)





biomass.df<-data.frame(biomass)

type <- rep(0, length(biomass.df$totwood) )
type[which(is.na(biomass.df[,"states"]))] = 'No data'
type[which(!is.na(biomass.df[,c("totwood")]) & biomass.df[,"states"]=="Minnesota" & biomass.df[, "totwood"]>= 5 )] = "Minnesota forest"
type[which(!is.na(biomass.df[,c("totwood")])& biomass.df[,"states"]=="Wisconsin" & biomass.df[, "totwood"]>= 50 )] = "Wisconsin forest"
type[which(!is.na(biomass.df[,c("totwood")]) & biomass.df[,"states"]=="Michigan" & biomass.df[, "totwood"]>= 50 )] = "Michigan forest"
type[which(!is.na(biomass.df[,c("totwood")]) & biomass.df[,"states"]=="Illinois" & biomass.df[, "totwood"]>= 50 )] = "Illinois forest"
type[which(!is.na(biomass.df[,c("totwood")])& biomass.df[,"states"]=="Indiana" & biomass.df[, "totwood"]>= 50 )] = "Indiana forest"

type[which(!is.na(biomass.df[,c("totwood")]) & biomass.df[,"states"]=="Minnesota"& biomass.df[, "totwood"]< 5 )] = "Minnesota prairie"
type[which(!is.na(biomass.df[,c("totwood")]) & biomass.df[,"states"]=="Wisconsin"& biomass.df[, "totwood"]< 50 )] = "Wisconsin prairie"
type[which(!is.na(biomass.df[,c("totwood")])& biomass.df[,"states"]=="Michigan"& biomass.df[, "totwood"]< 50 )] = "Michigan prairie"
type[which(!is.na(biomass.df[,c("totwood")]) & biomass.df[,"states"]=="Illinois"& biomass.df[, "totwood"]< 50 )] = "Illinois prairie"
type[which(!is.na(biomass.df[,c("totwood")]) & biomass.df[,"states"]=="Indiana"& biomass.df[, "totwood"]< 50 )] = "Indiana prairie"
type[which(type==0)]<- NA
#type[which(!is.na(biomass.df[,"states"]) & biomass.df[,"states"]=="Indiana" & biomass.df[, "totwood"]< 50 )] = "prairie"
df.with.temp <- read.csv("C:/Users/Kelly/Documents/Indiana_Density_Biomass/biomass_full_extract_with_temp_pr_elev.csv")
df.temp<- df.with.temp[which(df.with.temp$cell %in% biom.export$cell),]


                            



biomass.df <- data.frame(biomass.df,df.temp[,54:62])
                         
biomass.df$type <- type
bio<- biomass.df[,c("Oak", "Beech", "TBD", "TNE","totwood","jja_pr","djf_pr", "mam_pr","son_pr", "total_pr", "ksat", "sand",
                    "silt", "clay", "awc", "elev", "total_tmean",
                    "son_tmean", "djf_tmean", "mam_tmean", 
                    "jja_tmean", "TRI", "slope", "TPI")]
#write csv with soil covariates
write.csv(biomass, "biomass_v09_extract_with_soil.csv")

ggplot()+geom_point(data=biomass.df, aes(x=x,y=y, colour= type))

plot(log(bio$sand), bio$Beech)

pca.data<-biomass.df[,c("type","states","jja_pr","djf_pr", "mam_pr","son_pr", "total_pr", 
                        "jja_tmean", "djf_tmean", "mam_tmean", "son_tmean", 
                        "total_tmean","ksat", "sand", "silt", "clay", "awc", "elev", "TPI", "TRI", "slope")]
pca.data<- pca.data[which(!pca.data$type=="NA"),]
pca.data3 <-pca.data[complete.cases(pca.data),]
#index<- which(type=='NA')
#pca.data2<- pca.data[pca.data %in% index, ]

pc <- prcomp(na.omit(pca.data3[,3:ncol(pca.data3)]), scale = TRUE)
plot(pc)
p<- biplot(pc)

pdf("biplots.pdf")

g <- ggbiplot(pc, choices= c(1,2),obs.scale = 1, var.scale = 1, 
              groups = na.omit(pca.data3$type), ellipse = TRUE, 
              circle = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'vertical', 
               legend.position = 'right')
g

b <- ggbiplot(pc, choices= c(2,3),obs.scale = 1, var.scale = 1, 
              groups = na.omit(pca.data3$type), ellipse = TRUE, 
              circle = FALSE)
b <- b + scale_color_discrete(name = '')
b <- b + theme(legend.direction = 'vertical', 
               legend.position = 'right')
b

d <- ggbiplot(pc, choices= c(1,3),obs.scale = 1, var.scale = 1, 
              groups = na.omit(pca.data3$type), ellipse = TRUE, 
              circle = FALSE)
d <- d + scale_color_discrete(name = '')
d <- d + theme(legend.direction = 'vertical', 
               legend.position = 'right')
d
dev.off()




pdf("basic.plots.pdf")
plot(biomass.df$sand, biomass.df$totwood, ylim=c(0,500))
plot(biomass.df$silt, biomass.df$totwood, ylim=c(0,500))
plot(biomass.df$clay, biomass.df$totwood, ylim=c(0,500))
plot(biomass.df$awc, biomass.df$totwood, ylim=c(0,500))
plot(biomass.df$ksat, biomass.df$totwood, ylim=c(0,500))
plot(biomass.df$jja_pr, biomass.df$totwood,xlim=c(0,400), ylim=c(0,500))
plot(biomass.df$son_pr, biomass.df$totwood,xlim=c(0,400), ylim=c(0,500))
plot(biomass.df$djf_pr, biomass.df$totwood, xlim=c(0,400),ylim=c(0,500))
plot(biomass.df$mam_pr, biomass.df$totwood,xlim=c(0,400), ylim=c(0,500))
plot(biomass.df$total_pr, biomass.df$totwood,xlim=c(0,1100), ylim=c(0,500))


plot(log(biomass.df$sand), log(biomass.df$totwood))
plot(log(biomass.df$silt), log(biomass.df$totwood))
plot(log(biomass.df$clay), log(biomass.df$totwood))
plot(log(biomass.df$awc), log(biomass.df$totwood))
plot(log(biomass.df$ksat), log(biomass.df$totwood))
plot(log(biomass.df$jja_pr), log(biomass.df$totwood))
plot(log(biomass.df$son_pr), log(biomass.df$totwood))
plot(log(biomass.df$djf_pr), log(biomass.df$totwood))
plot(log(biomass.df$mam_pr), log(biomass.df$totwood))
plot(log(biomass.df$total_pr), log(biomass.df$totwood))

dev.off()

X11(width = 11)
install.packages('GGally')
library(GGally)

ggpairs(biomass.df)


#read indiana & illinois spec density
dens <- read.csv()

biomass.df$sand.10 <- biomass.df$sand/10

lm1 = lm(totwood~jja_pr,data=biomass.df, family=gaussian)



