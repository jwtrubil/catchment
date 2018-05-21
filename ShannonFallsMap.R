# ShannonFallsMap.r
#
# Demonstration of working with spatial objects and mapping within R,
# including the following:
#  - creating Spatial* objects from data frames
#  - coordinate transformation
#  - use of RSAGA package to delineate catchment boundaries
#  - mapping using the plot() function
#
# 2016-Aug-11 RDM
#######################################################################

# empty work space
rm(list = ls())

########################################################################
#  You will need to download the following files to run this script:
#    - ShannonFallsDEM.tif (DEM file)
#    - canvec_151201_164206_shp (folder containing shape files)
#
#  After downloading, save the files in a working directory. 
########################################################################

# specify working directory; you will need to change this line
setwd("c:/Research/Shannon Falls/Map")

# reset plot window to a single frame 
par(mfrow = c(1, 1), mar = c(5, 5, 5, 5))

# load libraries of required functions
library(sp)
library(raster)
library(maptools)
library(rgdal)
library(RSAGA)
library(GISTools)

# set up environment for using RSAGA
myenv = rsaga.env(workspace = getwd(), path = 'C:/Program Files/saga_2.1.0_win32', version = '2.1.0')


########################################################################
# Read in DEM and hydrography shape files obtained from Geogratis
########################################################################

# read in dem as a geotiff and reproject to utm square grid with 25 m resolution
dem = raster("ShannonFallsDEM.tif")
dem.utm = projectRaster(dem, crs = "+proj=utm +zone=10 +ellps=WGS84", res = 25, method = "bilinear") 

# get all files with the .shp extension from directory
shpdir = "c:/Research/Shannon Falls/Map/canvec_151201_164206_shp"
shps = dir(shpdir, "*.shp")
nshps = length(shps)
used.shapes = c("hd_1440009_1.shp", "hd_1470009_1.shp")

nc = nchar(used.shapes)
layername = substr(used.shapes, 1, (nc - 4))
sf1 = readOGR(dsn = shpdir, layer = layername[1])
sf1.utm = spTransform(sf1, CRS("+proj=utm +zone=10 +ellps=WGS84"))

sf2 = readOGR(dsn = shpdir, layer = layername[2])
sf2.utm = spTransform(sf2, CRS("+proj=utm +zone=10 +ellps=WGS84"))

sf1l = as(sf1.utm, "SpatialLines")
sf2l = as(sf2.utm, "SpatialLines")

# plot to check
plot(dem.utm)
plot(sf1l, add = T)
plot(sf2l, add = T, col = "blue")


########################################################################
# Specify catchment outlet and sensor locations
#  - locations specified in long/lat, then transformed to utm using
#    spTransform()
#  - the coordinates() function converts a data frame with coordinate 
#    variables to a spatial points object
########################################################################

# Squamish airport
lat = 49 + 46/60 + 59.004/3600 
long = -(123 + 9/60 + 39.006/3600)
sapll = data.frame(long, lat)
coordinates(sapll) = ~long + lat
projection(sapll) = CRS("+proj=longlat +datum=WGS84")
sap = spTransform(sapll, CRSobj = CRS("+proj=utm +zone=10 +ellps=WGS84"))

# Upper falls weather station
lat = 49 + 40/60 + 1.15/3600
long = -(123 + 9/60 + 5.6/3600)
ufwxll = data.frame(long, lat)
coordinates(ufwxll) = ~long + lat
projection(ufwxll) = CRS("+proj=longlat +datum=WGS84")
ufwx = spTransform(ufwxll, CRSobj = CRS("+proj=utm +zone=10 +ellps=WGS84"))

# catchment outlet (point of interest, poi) 
poill = data.frame(x = -123.16139, y = 49.67105)  
coordinates(poill) = ~x + y
projection(poill) = CRS("+proj=longlat +datum=WGS84") 
poi = spTransform(poill, CRSobj = CRS("+proj=utm +zone=10 +ellps=WGS84"))

# tidbit locations (in UTM)
tbN = c(5501976, 5501645)
tbE = c(488425, 488904)
tb = data.frame(y = tbN, x = tbE)
coordinates(tb) = ~x + y
projection(tb) = CRS("+proj=utm +zone=10 +ellps=WGS84")

# read in dem as a geotiff and reproject to utm square grid with 25 m resolution
dem = raster("ShannonFallsDEM.tif")
dem.utm = projectRaster(dem, crs = "+proj=utm +zone=10 +ellps=WGS84", res = 25, method = "bilinear") 

# plot dem, tidbit locations and poi to check that everything looks okay
plot(dem.utm)
plot(sf1l, add = T)
plot(sf2l, add = T, col = "blue")
plot(poi, add = T, pch = 21, bg = "cyan")
plot(tb, add = T, pch = 21, bg = "yellow")
plot(sap, add = T, pch = 22, bg = "yellow")


########################################################################
# delineate catchment area polygon using SAGA
########################################################################

# write dem.utm to SAGA format
writeRaster(dem.utm, filename = "SFdemutm.sgrd", format = "SAGA", overwrite = T)

# fill sinks - note that the algorithm will create negative elevations in Howe Sound
rsaga.fill.sinks(
  in.dem = "SFdemutm.sgrd", out.dem = "SFdemutmfilled.sgrd", 
  method = "planchon.darboux.2001", minslope = 0.1,
  env = myenv
)

# read in filled dem, specify projection and convert negative values to NA
fill_dem = raster("SFdemutmfilled.sdat")
projection(fill_dem) <- CRS("+proj=utm +zone=10 +ellps=WGS84")
values(fill_dem) = ifelse(values(fill_dem) < 0, NA, values(fill_dem))

# determine upslope drainage areas and map
rsaga.parallel.processing("SFdemutmfilled.sgrd", out.carea = "SFda.sgrd", env = myenv)
SFda = raster("SFda.sdat")
projection(SFda) = CRS("+proj=utm +zone=10 +ellps=WGS84")
plot(log(SFda))
plot(poi, add = T)
contour(fill_dem, add = T)

# find maximum drainage area near the poi - assign as catchment outlet
# i'm using a buffer of 50 m, but a different size may be appropriate
# for different data sets
buff = raster::extract(x = SFda, y = poi, buffer = 50, cellnumbers = TRUE)
buff.df = as.data.frame(buff)
snap_cell <- buff.df$cell[which.max(buff.df$value)]
snap_xy <- xyFromCell(SFda, snap_cell)
snap_xy = data.frame(y = snap_xy[1, 2], x = snap_xy[1, 1])
coordinates(snap_xy) = ~ x + y

# determine pixels that drain into the poi - save as grid
rsaga.geoprocessor(
  lib = 'ta_hydrology', module = 4, env = myenv,
  param = list(TARGET_PT_X = snap_xy@coords[1],
               TARGET_PT_Y = snap_xy@coords[2],
               ELEVATION = 'SFdemutmfilled.sgrd',
               AREA = 'basin.sgrd',
               METHOD = 0
  )					 
)
					 
# convert grid version of catchment to a polygon
rsaga.geoprocessor(
  lib = 'shapes_grid', module = 6, env = myenv,
  param = list(GRID = 'basin.sgrd',
               POLYGONS = 'basin.shp',
               CLASS_ALL = 0,
               CLASS_ID = 100,
               SPLIT = 0
  )
)					 			 
					 
# read in the shapefile as a spatial polygons data frame
basin.spdf = readShapeSpatial('basin.shp', proj4string = CRS("+proj=utm +zone=10 +ellps=WGS84"))

# convert the spdf to a spatial lines object 
basin.sl = as(basin.spdf, "SpatialLines")

# create full map with all elements
bb = bbox(fill_dem)
xminmax = bb[1, ]
yminmax = bb[2, ]

xmin = min(xminmax)
xmax = max(xminmax)
ymin = min(yminmax)
ymax = max(yminmax)
xsb = xmax + 500
ysb = ymin + 700

# to generate map and store as a tiff file, uncomment the following line and the last line of the script
# tiff("Figure4.tiff")
par(mar = c(1, 1, 1, 1), cex = 1.1)
plot(sf1.utm, xlim = c(xmin, xmax + 700), col = "black", axes = F)
plot(fill_dem, add = T, 
  col = gray(seq(1, 0.2, by = -0.1), alpha = 1),
  legend.mar = 5,
  legend.width = 0.6, legend.shrink = 0.4,
  legend.args = list(text = '  Elevation (m)', side = 3, font = 2, line = 1, cex = 0.8)
)

plot(sf2.utm, add = T, col = "blue")
plot(tb, add = T, pch = 21, bg = "red")
plot(ufwx, add = T, pch = 21, bg = "yellow")
plot(sap, add = T, pch = 22, bg = "yellow")
plot(basin.sl, col = "yellow", add = T) 
plot(sf1.utm, xlim = c(xmin, xmax + 5000), col = "black", add = T)

scalebar(4000, xy = c(xsb, ysb), 
  divs = 4,
  label = c(0, 2, 4), below = "km",
  type = "bar"
)

north.arrow(xb = xmax + 2500, yb = ymin + 3500, len = 300, col = "lightgray")

legend("topright", bty = "n",
  pch = c(22, 21, 21), 
  legend = c("YOW", expression(italic(T[w])), "UF Wx"),
  pt.bg = c("yellow", "red", "yellow")
)

text(x = xmin + 2500, y = ymin + 5500, 
  labels = expression(italic("Howe Sound")),
  srt = 20, cex = 0.8
)

# dev.off()