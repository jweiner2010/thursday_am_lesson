# HSI LiDAR fusion lesson 

library(raster)
library(rhdf5)
library(rgdal)
library(neonAOP)

#import dsm
dsm <- raster("../NEONdata/D17-California/TEAK/2013/lidar/TEAK_lidarDSM.tif") 

#import dtm
dtm <- raster("../NEONdata/D17-California/TEAK/2013/lidar/TEAK_lidarDTM.tif")

#import chm
chm <- raster("../NEONdata/D17-California/TEAK/2013/lidar/TEAK_lidarCHM.tif")

#check out height distribution - seems reasonable?

plot(chm,
     main = "Canopy Height\n Lower Teakettle, CA")

#check out histogram
hist(chm,
     main = "Distribution of Canopy Height\n Teakettle, CA",
     xlab = "Tree Height (m)")

#look at overall image stats
cellStats(chm, max)
cellStats(chm, sd)

#stack rasters together
lidar.brick <- brick(dsm, dtm, chm)

#check it out
plot(lidar.brick)

#HYPERSPECTRAL FUN TIME

#identify file of interest
f <- "../NEONdata/D17-California/TEAK/2013/spectrometer/reflectance/Subset3NIS1_20130614_100459_atmcor.h5"

#id the projection
#define the CRS definition by EPSG
epsg <- 32611

#create a list of bands
bands <- c(60,83)

#look at all the wavelengths
wavelengths <- h5read(f, "wavelength")

#read in bands
ndvi.stack <- create_stack(f,
                           bands = bands,
                           epsg = epsg)

#plot to check
plot(ndvi.stack, 
     main = "Bands 60 and 83\n Teakettle, CA")

#calculate NDVI
ndvi <- (ndvi.stack[[2]]-ndvi.stack[[1]])/(ndvi.stack[[2]]+ndvi.stack[[1]])

names(ndvi) <- "TEAK_hsiNDVI"

#plot it!
plot(ndvi,
     main = "NDVI\n NEON Lower Teakettle Site")

#try stack
#all.data <- brick(ndvi, lidar.brick)

#check extents
extent(chm)
extent(ndvi)

#check the extents of the two rasters and crop to the smaller one

if (extent(chm) == extent(ndvi)) {
} else {
    print("Extents are different, cropping data")
    overlap <- raster::intersect(extent(ndvi),extent(lidar.brick))
    #now let's crop te lidar data to the HSI
    lidar.brick <- crop(lidar.brick, overlap)
    ndvi <- crop(ndvi, overlap)
}

#try stacking and see what happens
all.data <- brick(ndvi, lidar.brick)

#rename layers
names(all.data) <- c("NDVI", "DSM", "DTM", "CHM")

#read in the NEON NDVI product

ndvi.neon <- raster("../NEONdata/D17-California/TEAK/2013/spectrometer/veg_index/TEAK_NDVI.tif")

#compare neon to narrow band ndvi
ndvi.diff <- ndvi - ndvi.neon

plot(ndvi.diff,
     main = "NDVI Diff Map\n Lower Teakettle, CA")

#write function to compare extents and crop layers if they are different

same_extent <- function(raster1, raster2) {
  if (extent(raster1) ==  extent(raster2)) {
    print("Rasters have the same extent")
  } else {
    overlap <- raster::intersect(extent(raster1), extent(raster2))
    #crop both rasters
    #might be good to check which is larger and compare
    print("Extents are different, Cropping data")
    raster1 <- crop(raster1, overlap)
    raster2 <- crop(raster2, overlap)
  }
  #stack the rasters
  raster.stack <- stack(raster1, raster2)
  return(raster.stack)
}

#import NEON aspect
aspect <- raster("../NEONdata/D17-California/TEAK/2013/lidar/TEAK_lidarAspect.tif")

#crop aspect to other data
all.data <- same_extent(aspect, all.data)

#Error! How to fix
#extract stuff we want
all.data <- all.data[[1:5]]

names(all.data) <- c("Aspect", "NDVI", "DSM", "DTM", "CHM")

#plot aspect
plot(all.data$Aspect)

#plot how a mask works
plot(all.data$Aspect > 270)

#create a classified aspect intermediate output
#first create a matrix of values that represents the classification ranges
#north-facing = 1, south-facing = 2

class.m <- c(-0.1, 45, 1,
             45, 135, NA,
             135, 225, 2,
             225, 315, NA,
             315, 360, 1)

#reshape into matrix
rcl.m <- matrix(class.m,
                ncol = 3,
                byrow = TRUE)

#classify aspect using classification matrix
asp.ns <- reclassify(all.data$Aspect, rcl.m)

#get map extent
ns.extent <- extent(asp.ns)

plot(asp.ns,
     col = c("blue", "green"),
     axes = F,
     main = "North and South Facing Slopes\n NEON Lower Teakettle Field Site",
     bty = "n",
     legend = FALSE)

par(xpd = TRUE)
legend((par()$usr[2]+20), ns.extent@ymax-100,
       legend = c("North", "South"),
       fill = c("blue", "green"),
       bty = "n")
par(xpd = FALSE)

#create a northfacing slope mask
north.facing <- asp.ns == 1
north.facing[north.facing == 0] <- NA

#create a southfacing slope mask
south.facing <- asp.ns == 2
south.facing[south.facing == 0] <- NA

#histogram of tree height
hist(all.data$CHM,
     main = "Distribution of Canopy Height Model Values\n NEON Lower Teakettle")

#get mean, min, max values for all layers
all.data.stats <- data.frame(t(summary(values(all.data[[-1]]), na.rm = TRUE)))

ht.mean <- cellStats(all.data$CHM, mean, na.rm = T)

ht.sd <- cellStats(all.data$CHM, sd, na.rm = T)

#make a thresholds object
thresholds <- data.frame(id = 1)

#calling "tall trees" things that are mean + 1 sd tall
thresholds$height <- ht.mean + ht.sd

#TIME FOR NDVI
hist(all.data$NDVI,
     main = "Distribution of NDVI values\n Teakettle",
     col = "cornflowerblue")

#what counts as "green"?
greenRange <- cellStats(all.data$NDVI, range)
greenRange <- greenRange[2] - greenRange[1]

thresholds$green <- cellStats(all.data$NDVI, max) - greenRange/3

#frequency of north/south-facing slopes
north.count <- freq(asp.ns, value = 1)
south.count <- freq(asp.ns, value = 2)

#how many are tall and green?
north.tall.green <- asp.ns == 1 &
  all.data$NDVI >= thresholds$green &
  all.data$CHM >= thresholds$height

south.tall.green <- asp.ns == 2 &
  all.data$NDVI >= thresholds$green &
  all.data$CHM >= thresholds$height

#assign zeros to NAs
north.tall.green[north.tall.green == 0] <- NA

south.tall.green[south.tall.green == 0] <- NA

#count pixels
north.tall.green.count <- freq(north.tall.green, value = 1)

south.tall.green.count <- freq(south.tall.green, value = 1)

#fractional cover of tall, green trees
north.tall.green.frac <- north.tall.green.count/north.count

south.tall.green.frac <- south.tall.green.count/south.count

#NOW TO PLOT CIR-image with overlays
bands <- c(83, 60, 35)

cir.stack <- create_stack(f,
                          bands = bands,
                          epsg = epsg)

#plot everything in RGB world
plotRGB(cir.stack,
        scale = 1,
        stretch = "lin")

plot(north.tall.green,
     col = "cyan",
     add = T,
     legend = F)

plot(south.tall.green,
     col = "blue",
     add = T,
     legend = F)










