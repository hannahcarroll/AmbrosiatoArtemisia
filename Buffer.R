# Creating the 2 degree buffer
library(rgdal)
library(maptools)
library(rgeos)

# Download the shapefile
download.file("ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/na_cec_eco_l2.zip" , destfile="NA_CEC_Eco_Level2.zip")

# Unzip this file
unzip("NA_CEC_Eco_Level2.zip")

# Load the shapefile
ecor.shp <- readOGR(
  dsn = paste0(getwd()), 
  layer = "NA_CEC_Eco_Level2"
  )

# Transform to lat long
ecor.shp.2 <- spTransform(ecor.shp, CRS("+proj=longlat +datum=NAD83"))

# Subset to Great Plains

# Create points for ecoregions
ecor.shp.2@data$id <- rownames(ecor.shp.2@data)
ecor.shp.2.points <- as.data.frame(ecor.shp.2, region="id")

# Join
ecor.shp.2.df <- merge.data.frame(ecor.shp.2.points, ecor.shp.2@data, by="id", )

# Get just the Great Plains region
greatplains.sp <- ecor.shp.2[ecor.shp.2@data$NA_L1KEY == "9  GREAT PLAINS", ]

# Define the ID values
ID <- greatplains.sp@data$NA_L1KEY

# Dissolve the subregions within the great plains
greatplains.dissolved <- unionSpatialPolygons(greatplains.sp, ID)

# Create the buffer
buffer <- gBuffer(greatplains.dissolved, byid=FALSE, id=NULL, width=2.0, quadsegs=5, capStyle="ROUND",
 joinStyle="ROUND", mitreLimit=1.0)

# This feeds into Climate.R and ambtoartcode.R
