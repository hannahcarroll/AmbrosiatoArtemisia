packages <- c("grid", "Rcpp", "sp", "dtplyr", "ggmap", "rgdal", "rgeos", "maptools", 
              "plyr", "dplyr", "tidyr", "cowplot", "MASS", "caret",
       "RColorBrewer", "maps", "sqldf", "neotoma", "analogue", "data.table", "raster", 
       "spatialEco", "pca3d", "rgl", "plotly", "tictoc", "mda", "viridis",
       "htmlwidgets", "leaflet", "widgetframe", "MASS", "gstat", "mgcv")

lapply(packages, require, character.only = TRUE)

#######################################################################################################
# Some calculations are highly memory intensive and therefore all setup steps have been commented out.#
#######################################################################################################

# Original dataset from which presettlement set was split. Very large.
# pollen22k.buffer <- read.csv("pollen22k 2deg buffer whitmore full.csv", header=TRUE)

# Update raw radiocarbon ages
# correctionfactor <- read.table("intcal0914c.txt", sep=",")
#  names(correctionfactor) <- c("calBP", "age"," Error", "D14Cpermil", "sigma14C")
   # Reimer et al. (2009)
     # correction.curve <- lm(calBP ~ age, data=correctionfactor)

# pollen22k.buffer$calendar.ybp <- pollen22k.buffer$age

# splitpollen.calibrated <- subset(pollen22k.buffer, 
#              date.type=="Calibrated radiocarbon years BP" | date.type=="Calendar years BP" |
#                date.type=="Varve years BP")

# splpol.cal <- within(splitpollen.calibrated, 
#      {calendar.ybp[date.type=="Calendar years BP"] <- calendar.ybp[date.type=="Calendar years BP"]; 
#       calendar.ybp[date.type=="Calibrated radiocarbon years BP"] <-  calendar.ybp[date.type=="Calibrated radiocarbon years BP"];
#       calendar.ybp[date.type=="Varve years BP"] <-  calendar.ybp[date.type=="Varve years BP"]})

# splitpollen.ybp <- subset(pollen22k.buffer, date.type=="Radiocarbon years BP")
 
# splitpollen.ybp$calendar.ybp <- predict(correction.curve, splitpollen.ybp)

# The corrected dataset
# pollen22k.buffer <- rbind(splpol.cal, splitpollen.ybp)

# Subset by age
# presettlement.pollen <- pollen22k.buffer[pollen22k.buffer$calendar.ybp >=250 & 
#                                         pollen22k.buffer$calendar.ybp <=500, ]

############################# Clean dataset ###############################

# Drop spores

# pres.pollen.clean <- presettlement.pollen[ , -which(names(presettlement.pollen) %in% c("EQUISETU", "PTERIDIUM", "LYCOPODX", "SPHAGNUM", "x.counts.good_counts..."))]

#Reorder for the sake of readability
# pres.pollen.clean2 <- pres.pollen.clean[,c(1:8,118,9:117)]

# Drop columns that are all na
# pres.pollen.clean2 <- Filter(function(x)!all(is.na(x)), pres.pollen.clean2)

# Recode NA pollen observations to zero
# pres.pollen.clean2[, 13:112][is.na(pres.pollen.clean2[, 13:112])] <- 0

# Calculate the pollen sum
# pres.pollen.clean2$pollensum <- rowSums(pres.pollen.clean2[13:112])

# Drop cases where both Ambrosia and Artemisia are 0
# Keep a list of those that have been dropped

# pres.pollen.dropped <- pres.pollen.clean2 %>% 
#  filter(AMBROSIA == 0 & ARTEMISIA == 0)

# Tally of datasets reporting Artemisia but no Ambrosia
# count(pres.pollen.clean2 %>% 
#  filter(AMBROSIA == 0 & ARTEMISIA != 0)) # 77

# Tally of datasets reporting Ambrosia but no Artemisia
# count(pres.pollen.clean2 %>% 
#  filter(AMBROSIA != 0 & ARTEMISIA == 0)) # 13

# Tally of datasets reporting both taxa as nonzero
# count(pres.pollen.clean2 %>% 
#  filter(AMBROSIA != 0 & ARTEMISIA != 0)) # 680

# Tally of datasets reporting either taxon as nonzero
# count(pres.pollen.clean2 %>% 
#  filter(AMBROSIA != 0 | ARTEMISIA != 0)) # 770

# Keep datasets where either taxon is nonzero
# pres.pollen.clean2 <- pres.pollen.clean2 %>% 
#  filter(AMBROSIA != 0 | ARTEMISIA != 0)

# Calculate the proportion of each taxon
# pres.pollen.clean2 <- pres.pollen.clean2 %>% 
#  mutate(
#    ambrosia.proportion = AMBROSIA/pollensum, 
#    artemisia.proportion = ARTEMISIA/pollensum
#  )

# Recode zero ambrosia records
# pres.pollen.clean2$ambrosia.proportion[pres.pollen.clean2$ambrosia.proportion == 0] <- 0.00001
# pres.pollen.clean2$artemisia.proportion[pres.pollen.clean2$artemisia.proportion == 0] <- 0.00001

# Calculate ambrosia to artemesia ratio
# pres.pollen.clean2 <- pres.pollen.clean2 %>% 
#  mutate(
#    ambtoart = ambrosia.proportion/artemisia.proportion
#  )

# Double check that no weird ones snuck through   
# NaNcases2 <- sum(is.na(pres.pollen.clean2$ambtoart)) # 0 datasets
# Infcases2 <- sum(is.infinite(pres.pollen.clean2$ambtoart)) # 0 datasets


########################### Match datasets to their ecoregion ###################################

# Create a unique ID field to join on later
# pres.pollen.clean2$uniqueid <- 1:nrow(pres.pollen.clean2)

# Point in poly to assign ecoregions to datasets
# pres.pollen.sp <- pres.pollen.clean2[,c(10,11,116)]

# Load shapefile of ecoregions
# ecoregions <- readOGR(dsn=getwd(), layer="ecoregions_NAD83") 

# Assign coordinates
# coordinates(pres.pollen.sp) <- ~ long + lat

# Assign the ecoregions projection to the pollen datasets
# proj4string(pres.pollen.sp) <- CRS(proj4string(ecoregions))


# The datsets get the region they are on top of
# data.by.region <- over(pres.pollen.sp, ecoregions)


# Create a field to join on
# data.by.region$uniqueid <- 1:nrow(data.by.region)

# Merge dataset location with region
# pres.pollen.regions <- merge(pres.pollen.clean2, data.by.region, by="uniqueid")

######################### Join with climate data ################################
# Each presettlement set gets its modern-day climate
# Caution: very memory intensive

###################################################################################################################
# This code depends on the 2 decimal degree buffer created in Buffer.R and climate data from the Climate.R script #
###################################################################################################################

# Coordinates in long~lat order
# xy <- pres.pollen.regions[,c(12,11)]

# Create the spatial points data frame
#pres.pollen.regions.spdf <- SpatialPointsDataFrame(coords = xy, data = pres.pollen.regions,
#                               proj4string = CRS(proj4string(clim.2degbuffer)))

# Join with climate data

# Using pres.pollen.regions and clim.2degbuffer Climate script

#tic("compute distance matrix")
### compute the complete distance matrix between the two sets of points
#dist_mat <- pointDistance(pres.pollen.regions.spdf, clim.2degbuffer, lonlat = TRUE, allpairs = TRUE)
#toc()

#tic("find nearest point")
### identify nearest point in dataset B for every point in dataset A
#nearest <- apply(dist_mat, 1, which.min)
#toc()

#tic("bind data together")
### Bind the data
#pres.pollen.regions.spdf@data<- cbind(pres.pollen.regions.spdf@data, clim.2degbuffer@data[nearest,])
#toc()

# Flatten to data frame
#pres.pollen.regions.df <- as.data.frame(pres.pollen.regions.spdf)

# Set the names of the BioClim columns
#names(pres.pollen.regions.df)[names(pres.pollen.regions.df) %in%  c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8",                
#                                      "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",               
#                                      "bio17", "bio18", "bio19")] <- c("MeanT", 
#                    "Mean.Diurnal.RangeT", "Isothermality", "TSeasonality",
#                    "MaxTWarm", "MinTCold", "TAnnRange", "MeanTWetQ", "MeanTDryQ", "MeanTWarmQ",
#                    "MeanTColdQ", "AnnP", "PrecipWetMonth", "PrecipDryMonth", "PrecipSeasonality", 
#                    "PrecipWetQ", "PrecipDryQ", "PrecipWarmQ", "PrecipColdQ")

################################### Drop sites reflecting non-grassland taxa ############################################

# Remove Western Cordillera at reviewer request, and sites receiving pollen from non-grassland taxa
# pres.pollen.regions <- pres.pollen.regions.df[!grepl("UPPER GILA MOUNTAINS", pres.pollen.regions.df$NA_L2NAME),]    
# pres.pollen.regions2 <- pres.pollen.regions[!grepl("WESTERN CORDILLERA", pres.pollen.regions$NA_L2NAME),]
# pres.pollen.regions2 <-  pres.pollen.regions2[!grepl("COLD DESERTS", pres.pollen.regions2$NA_L2NAME),]

# write.csv("presettlement set.csv", pres.pollen.regions2, row.names = FALSE)

pres.pollen.regions2 <- read.csv("presettlement set.csv")



#######################################
# Set up surface samples for analysis #
#######################################

# Code used to download surface samples was generated in 8/2018
# This taxes *extremely* long to run. Recommend load from csv.

# bbox(buffer)

# Start with bbox of 2 deg buffer
# greatplains.surface <- get_dataset(loc = c(-116.5924, 21.1635, -88.14369, 56.11235), 
#                      datasettype = 'pollen surface sample')

# tic("download")
# gplains.surface.d <- get_download(greatplains.surface)
# tic("compile taxa")
# gplains.surface.taxa <- compile_taxa(gplains.surface.d, "WhitmoreFull")
# tic("compile downloads")
# gplains.surface.c <- compile_downloads(gplains.surface.taxa)
# tic("write to file")
# toc()

# Clip dataset to buffer.nad83
# Copy data to new object
# gplains.surface.sp <- gplains.surface.c

# Set coordinates and projection
# coordinates(gplains.surface.sp) <- ~ long + lat
# proj4string(gplains.surface.sp) <- CRS(proj4string(buffer.nad83))

# Clip to 2 degree buffer
# tic()
# gplains.surface.buffer <- gplains.surface.sp[!is.na(over(gplains.surface.sp,as(buffer.nad83,"SpatialPolygons"))),]
# write.csv(file = "gplains 2degbuffer surface whitmore full 10-23-2018.csv", x=gplains.surface.buffer, row.names=FALSE)
#toc()

# This is the raw surface set
# gplains.s.df <- read.csv("gplains 2degbuffer surface whitmore full 10-23-2018.csv")

# Drop spores

# surface.pollen.clean <- gplains.s.df[ , -which(names(gplains.s.df) %in% c(".id", "EQUISETU", "PTERIDIUM", "LYCOPODX", "SPHAGNUM", "x.counts.good_counts...", "optional"))]

# Drop columns that are all na
# surface.pollen.clean2 <- surface.pollen.clean[, unlist(lapply(surface.pollen.clean, function(x) !all(is.na(x))))]

# Recode NA pollen observations to zero
# surface.pollen.clean2[, 8:94][is.na(surface.pollen.clean2[, 8:94])] <- 0

# Calculate total pollen count
# surface.pollen.clean2$pollensum <- rowSums(surface.pollen.clean2[8:94])

# surface.pollen.clean2 <- surface.pollen.clean2 %>% 
#  mutate(
#    ambrosia.proportion = AMBROSIA/pollensum, 
#    artemisia.proportion = ARTEMISIA/pollensum
#  )

# Drop cases where both Ambrosia and Artemisia are 0
# surface.pollen.clean2 <- surface.pollen.clean2[(!(surface.pollen.clean2$AMBROSIA == 0) | 
#                                                !(surface.pollen.clean2$ARTEMISIA == 0)),]
#

# Recode zero ambrosia records
# surface.pollen.clean2$ambrosia.proportion[surface.pollen.clean2$ambrosia.proportion == 0] <- 0.00001
# surface.pollen.clean2$artemisia.proportion[surface.pollen.clean2$artemisia.proportion == 0] <- 0.00001

########################### Match datasets to their ecoregion ###################################

# Should already be environment from the presettlement set
# ecoregions <- readOGR(dsn=getwd(), layer="ecoregions_NAD83") 

# Point in poly to assign ecoregions to datasets
# surface.pollen.clean2$uniqueid <- 1:nrow(surface.pollen.clean2)
# surface.pollen.sp <- surface.pollen.clean2[,c(5,6,98)]
# coordinates(surface.pollen.sp) <- ~ long + lat
# proj4string(surface.pollen.sp) <- CRS(proj4string(ecoregions))

# Create the spatial polygons data frame

# The datsets get the region they are on top of
# surface.data.by.region <- over(surface.pollen.sp, ecoregions)

# Create a field to join on
# surface.data.by.region$uniqueid <- 1:nrow(surface.data.by.region)

# Merge dataset location with region
# surface.pollen.regions <- merge(surface.pollen.clean2, surface.data.by.region, by="uniqueid")

# Tally of observations by ecoregion
# s.obs.per.ecoregion <- table(droplevels(surface.pollen.regions)$NA_L2NAME)

################################### Drop sites reflecting non-grassland taxa ############################################

# Remove Western Cordillera at reviewer request, and sites receiving pollen from non-grassland taxa
# s.pollen.regions2 <- surface.pollen.regions[!grepl("WESTERN CORDILLERA", surface.pollen.regions$NA_L2NAME),]
# s.pollen.regions2 <- s.pollen.regions2[!grepl("COLD DESERTS", s.pollen.regions2$NA_L2NAME),]

# Create column to facet by
# s.pollen.regions2$category <- NULL
# s.pollen.regions2$category[s.pollen.regions2$AMBROSIA == 0] <- "NoAmbrosia"
# s.pollen.regions2$category[s.pollen.regions2$ARTEMISIA == 0] <- "NoArtemisia"
# s.pollen.regions2$category[s.pollen.regions2$AMBROSIA != 0 & s.pollen.regions2$ARTEMISIA !=0] <- "BothPresent"


######################### Join with climate data ################################
# Each presettlement set gets its modern-day climate
# Caution: very memory intensive

###################################################################################################################
# This code depends on the 2 decimal degree buffer created in Buffer.R and climate data from the Climate.R script #
###################################################################################################################

# Set up coordinates
# surface.xy <- s.pollen.regions2[,c(7,6)]

# Create theh spdf
# surface.pollen.regions.spdf <- SpatialPointsDataFrame(coords = surface.xy, data = s.pollen.regions2,
#                               proj4string = CRS(proj4string(clim.2degbuffer)))

# Join with climate data
# Should already be in memory from setting up presettlement pollen
# tic("compute distance matrix")
### compute the complete distance matrix between the two sets of points
# s.dist_mat <- pointDistance(surface.pollen.regions.spdf, clim.2degbuffer, lonlat = TRUE, allpairs = TRUE)
# toc()

# tic("find nearest point")
### identify nearest point in dataset B for every point in dataset A
# s.nearest <- apply(s.dist_mat, 1, which.min)
# toc()

# tic("bind data together")
### bind together 
# surface.pollen.regions.spdf@data<- cbind(surface.pollen.regions.spdf@data, clim.2degbuffer@data[s.nearest,])
# toc()

# Convert to data frame
# surface.pollen.regions.df <- as.data.frame(surface.pollen.regions.spdf)
# names(surface.pollen.regions.df)[names(surface.pollen.regions.df) %in%  c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8",                
#                                      "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",               
#                                      "bio17", "bio18", "bio19")] <- c("MeanT", 
#                    "Mean.Diurnal.RangeT", "Isothermality", "TSeasonality",
#                    "MaxTWarm", "MinTCold", "TAnnRange", "MeanTWetQ", "MeanTDryQ", "MeanTWarmQ",
#                    "MeanTColdQ", "AnnP", "PrecipWetMonth", "PrecipDryMonth", 
#                    "PrecipSeasonality", "PrecipWetQ", "PrecipDryQ", "PrecipWarmQ", "PrecipColdQ")

# write.csv(surface.pollen.regions.df, file="surface set.csv", row.names = FALSE)

surface.pollen.regions.df <- read.csv("surface set.csv")

# Apply the modern correction

summary(pres.pollen.regions2$ambrosia.proportion)
summary(surface.pollen.regions.df$ambrosia.proportion)

# Ambrosia has more than doubled since settlement - making modern climate appear artifically wet 
# without a correction (or, making past cliamte appear artificially dry without a correction)
median(surface.pollen.regions.df$ambrosia.proportion)/median(pres.pollen.regions2$ambrosia.proportion)
# Median is the correction: 2.201151

mean(surface.pollen.regions.df$ambrosia.proportion)/mean(pres.pollen.regions2$ambrosia.proportion)

# Artemisia
summary(pres.pollen.regions2$artemisia.proportion)
summary(surface.pollen.regions.df$artemisia.proportion)


# Calculate ambrosia to artemesia ratio
surface.pollen.corrected <- surface.pollen.regions.df %>% 
  mutate(
    ambtoart = (ambrosia.proportion/2.201151)/artemisia.proportion
  )

length(unique(droplevels(surface.pollen.corrected$site.name))) # 433

count(surface.pollen.corrected %>% 
  filter(AMBROSIA == 0 & ARTEMISIA != 0)) # 101

count(surface.pollen.corrected %>% 
  filter(AMBROSIA != 0 & ARTEMISIA == 0)) # 83

count(surface.pollen.corrected %>% 
  filter(AMBROSIA != 0 & ARTEMISIA != 0)) # 323

count(surface.pollen.corrected %>% 
  filter(AMBROSIA != 0 | ARTEMISIA != 0)) # 507















summary(lm(formula=log(ambtoart) ~ AnnP, data=surface.pollen.regions.df[grepl("NoAmbrosia", surface.pollen.regions.df$category),]))
summary(lm(formula=log(ambtoart) ~ AnnP, data=s.pollen.nowc[grepl("NoAmbrosia", s.pollen.nowc$category),]))

summary(lm(formula=log(ambtoart) ~ AnnP, data=surface.pollen.regions.df[grepl("NoArtemisia", surface.pollen.regions.df$category),]))
summary(lm(formula=log(ambtoart) ~ AnnP, data=s.pollen.nowc [grepl("NoArtemisia", s.pollen.nowc $category),]))


# AOV followed by Tukey's HSD
# Surface set

aov.surf <- (aov(log(ambtoart) ~ NA_L2NAME, data=s.pollen.regions2))
tukey.surf <- TukeyHSD(aov.surf)
t.surf <- as.data.frame(tukey.surf[["NA_L2NAME"]])
summary(aov.surf)
write.csv(t.surf, file="surfacetukey_nowc.csv")



