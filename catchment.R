# code to illustrate application of the 'catchment' function
#
# 2018-May-24 Joel Trubilowicz and Dan Moore
#######################################################################################

library(RSAGA)
library(raster)
library(tidyverse)

#saga environment, in windows it seems to be able to find it automatically
saga.env <- rsaga.env(path = "/Applications/QGIS.app/Contents/MacOS/bin", 
                   modules = "/Applications/QGIS.app/Contents/MacOS/lib/saga")

#saga.env <- rsaga.env()
#rsaga.get.libraries(path = saga.env$modules)
#rsaga.get.modules('io_grid', env = saga.env)o


# input digital elevation model as a raster 
dem <- raster('./sourcedata/southcoastdem.sdat')

# coordinates of the catchment outlet
lat = 50.1194
long = -123.4361

# coordinate reference system - in this case, BC Albers
crs = 3005 

# buffer around catchment outlet to find location on digital stream network
pourpointsbuffer = 500 # units = m

# should sinks be filled
fillsinks = T

# name of output file
outname = 'elaho'

# 
#' Delineate a watershed
#'
#' @param dem Raster object of your dem.
#' @param lat A number.  Latitude of outlet in decimal degrees.
#' @param long A number. Longitude of outlet in decimal degrees.
#' @param buffsize A number, the buffer (m) around catchment outlet to find location on digital stream network.
#' @param crs A number.  The EPSG coordinate system number for the DEM and the output.
#' @param outname A character string.  The name for the output shapefile.
#' @param fillsinks Boolean.  Should sinks be filled?  Default is TRUE.
#' @param sinkmethod A character string. SAGA method for sink filling, options are "planchon.darboux.2001", "wang.liu.2006", or "xxl.wang.liu.2006" (default).
#' @param minslope A number.  Minimum slope angle preserved during sink filling, default is 0.01.
#' @param saga.env Saga environment object.  Default is to let saga find it on its own.
#' @return SpatialPolygonsDataFrame and a shapefile in working directory.
#' @export
catchment <- function(dem, 
                      lat, 
                      long, 
                      buffsize = 100, 
                      crs, 
                      outname, 
                      fillsinks = T,
                      sinkmethod = 'xxl.wang.liu.2006',
                      minslope = 0.01,
                      saga.env = rsaga.env()){
  library(RSAGA)
  library(raster)
  library(rgdal)
  
  #make crs string
  crs <- paste0("+init=epsg:", crs)
  
  #make a temporary directory 
  system('mkdir scratch')
  
  #put the dem object in there
  writeRaster(dem,"./scratch/dem.sdat",format="SAGA",NAflag=-9999, overwrite=T)
  
  #if you don't need to fill sinks, you can save a fair bit of processing time
  if (fillsinks == T) {     #fill sinks
    rsaga.fill.sinks("./scratch/dem.sgrd", './scratch/demfilled.sgrd', 
                     method = sinkmethod,
                     minslope = minslope,
                     env = saga.env)
    #calculate catchment area grid from filled dem
    rsaga.topdown.processing('./scratch/demfilled.sgrd', 
                             out.carea = './scratch/catchment_area.sgrd', 
                             env = saga.env)
  } else {
    #calculate catchment area grid direct from dem
    rsaga.topdown.processing("./scratch/dem.sgrd", 
                             out.carea = './scratch/catchment_area.sgrd', 
                             env = saga.env)
  }
  
  # make the base data frame, x is longitude and y is latitude
  gauge <- data.frame(y = lat, x = long)
  
  # turn into a spatial object
  coordinates(gauge) <- ~ x + y
  
  # assign the coordinate system (WGS84)
  projection(gauge) <- CRS("+init=epsg:4326")
  
  # reproject to BC Albers
  gauge <- spTransform(gauge, CRS(crs))
  
  # # plot it on the dem so I know it worked using the raster package
  # dem <- raster('./sourcedata/southcoastdem.sdat')
  # projection(dem) <- CRS("+init=epsg:3005")
  # plot(dem)
  # plot(gauge, add=T)
  
  # read in the catchment area grid
  catch_area <- raster('./scratch/catchment_area.sdat')
  
  # extract a window around around the gauge point, I am going to get the maximum value within 500 m of the gauge
  buffer <- raster::extract(catch_area, gauge, buffer = pourpointsbuffer, cellnumbers = T)[[1]] %>%
            as.data.frame
  
  # this is the location of the maximum catchment area on the grid, given as the id from the raster
  snap_loc <- buffer$cell[which.max(buffer$value)]
  
  # get the xy coordinates at that max location, which is now going to be the location of the gauge.
  snap_loc <- xyFromCell(catch_area, snap_loc)
  
  #make watershed as grid
  if (fillsinks == T){
    rsaga.geoprocessor(lib = 'ta_hydrology', 4,
    		param = list(TARGET_PT_X = snap_loc[1,1],
    			  				 TARGET_PT_Y = snap_loc[1,2],
    				  			 ELEVATION = './scratch/demfilled.sgrd',
    					  		 AREA = './scratch/bounds.sgrd',
    						  	 METHOD = 0),
    		env = saga.env)
  } else {
    rsaga.geoprocessor(lib = 'ta_hydrology', 4,
                       param = list(TARGET_PT_X = snap_loc[1,1],
                                    TARGET_PT_Y = snap_loc[1,2],
                                    ELEVATION = './scratch/dem.sgrd',
                                    AREA = './scratch/bounds.sgrd',
                                    METHOD = 0),
                       env = saga.env)
  }
  
  #convert shape to grid
  rsaga.geoprocessor(lib = 'shapes_grid', 6,
  		param = list(GRID = './scratch/bounds.sgrd',
  		             POLYGONS = outname,
  		             CLASS_ALL = 0,
  		             CLASS_ID = 100,
  		             SPLIT = 0),
  		env = saga.env)
  
  #return a spatialpolygonsdataframe and plot it onto the DEM
  basin <- readOGR('.', outname)
  projection(basin) <- CRS(crs)
  
  if (.Platform$OS.type == 'unix') {
    system('rm -r scratch/')
   } else {
     system('del /f /s /q scratch 1')
     system('rmdir /s /q scratch')
   } 
  return(basin)
}

basin <- catchment(dem, lat, long, pourpointsbuffer, crs, 
                   outname, fillsinks = T, sinkmethod = 'planchon.darboux.2001', 
                   minslope = 0.01, saga.env = saga.env)

plot(dem)
plot(basin, add = T) 
