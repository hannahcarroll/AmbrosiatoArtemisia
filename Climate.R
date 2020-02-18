### Temp and precip data ###
# Data from worldclim.org

library(raster)
library(rgeos)

####################################################
# Bioclim data in 2 degree buffer - uses buffer from Buffer.R

# Make sure it's loaded!
plot(buffer)

clim <- getData("worldclim",var="bio",res=2.5)
 # http://www.worldclim.org/bioclim

# Raster projection needs to match the buffer projection
clim.r <- projectRaster(clim, crs=(proj4string(buffer)))

# Clip climate data to buffer bounding box. Makes rgeos::over run much faster to start with less data.
clim.clip <- crop(clim.r, bbox(buffer))

# Convert to spatial pixels data frame
clim.spdf <- as(clim.clip, "SpatialPixelsDataFrame")

# This can be slow! 
clim.2degbuffer <- clim.spdf[!is.na(over(clim.spdf,as(buffer,"SpatialPolygons"))),]


