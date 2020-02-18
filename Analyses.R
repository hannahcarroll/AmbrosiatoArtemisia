##### Analyses #####
# To be run on pollen22k 2deg buffer.csv

packages <- c("grid", "Rcpp", "sp", "dtplyr", "ggmap", "rgdal", "rgeos", "maptools", 
              "plyr", "dplyr", "tidyr", "cowplot", "MASS", "plot3D", "caret",
       "RColorBrewer", "maps", "sqldf", "neotoma", "analogue", "data.table", "raster", 
       "spatialEco", "pca3d", "rgl", "plotly", "tictoc", "scatterplot3d", "mda", "viridis",
       "htmlwidgets", "leaflet", "widgetframe", "MASS", "gstat", "mgcv")

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, library, character.only = TRUE)
}

ipak(packages)

pollen22k.buffer <- read.csv("pollen22k 2deg buffer whitmore full.csv", header=TRUE)

####### Make data smaller #######

presettlement.pollen.raw <- pollen22k.buffer[pollen22k.buffer$calendar.ybp <=750, ]

# Convert radiocarbon ages to calendar YBP
correctionfactor <- read.table("intcal0914c.txt", sep=",")
names(correctionfactor) <- c("calBP", "age"," Error", "D14Cpermil", "sigma14C")
# Reimer et al. (2009)
correction.curve <- lm(calBP ~ age, data=correctionfactor)

correction.plot <- ggplot() + 
  geom_line(data=correctionfactor, aes(x=age, y=calBP)) +
  geom_point(data=correctionfactor, aes(x=age, y=calBP)) +
  xlab("Radiocarbon Years Before Present") + ylab("Calendar Years Before Present") +
  scale_y_continuous(breaks=seq(0, 25000, 5000), limits=c(5000,25000)) +
  scale_x_continuous(breaks=seq(0, 25000, 5000), limits=c(5000,25000))

# Apply correction factor to pollen22k.buffer

# Copy the age column to a new column
pollen22k.buffer$calendar.ybp <- pollen22k.buffer$age

splitpollen.calibrated <- subset(pollen22k.buffer, 
              date.type=="Calibrated radiocarbon years BP" | date.type=="Calendar years BP" |
                date.type=="Varve years BP")


splpol.cal <- within(splitpollen.calibrated, 
      {calendar.ybp[date.type=="Calendar years BP"] <- calendar.ybp[date.type=="Calendar years BP"]; 
       calendar.ybp[date.type=="Calibrated radiocarbon years BP"] <-  calendar.ybp[date.type=="Calibrated radiocarbon years BP"];
       calendar.ybp[date.type=="Varve years BP"] <-  calendar.ybp[date.type=="Varve years BP"]})

splitpollen.ybp <- subset(pollen22k.buffer, date.type=="Radiocarbon years BP")
 
splitpollen.ybp$calendar.ybp <- predict(correction.curve, splitpollen.ybp)

# The corrected dataset
pollen22k.buffer <- rbind(splpol.cal, splitpollen.ybp)

presettlement.pollen <- pollen22k.buffer[pollen22k.buffer$calendar.ybp >=250 & 
                                         pollen22k.buffer$calendar.ybp <=500, ]


presettlement.map <- allregionmap + 
  geom_polygon(data = gp.dissolve.f, aes(x=long, y = lat, group=group, color=id), fill=NA) + 
  geom_polygon(data=fortify(buffer), aes(x=long, y=lat, group=group, color = group), fill=NA) + 
   scale_color_manual(values=c("black", "blue"), 
                         labels=c("Great Plains Region", "2 Decimal Degree Buffer"), 
                         name="Outlines") +
  geom_point(data=presettlement.pollen, aes(x=long, y=lat), shape = 16, size=2, color="red3")

fossil.pollen <- pollen22k.buffer[pollen22k.buffer$calendar.ybp >500 & 
                                         pollen22k.buffer$calendar.ybp <=22000, ]

fossil.pollen.map <- allregionmap + geom_polygon(data = gp.dissolve.f, 
                        aes(x=long, y = lat, group=group), fill=NA, colour = "red") + 
                geom_point(data=fossil.pollen, aes(x=long, y=lat), shape = 16, size=1) +
                geom_polygon(data=fortify(buffer), aes(x=long, y=lat, group=group),
                           fill=NA, color = "blue") + xlab("Longitude") + ylab("Latitude") +
  ggtitle("Fossil Set - Late Quaternary (22,000-500 YBP)")
allregionmap + 
  geom_polygon(data = gp.dissolve.f, aes(x=long, y = lat, group=group, color=id), fill=NA) + 
  geom_polygon(data=fortify(buffer), aes(x=long, y=lat, group=group, color = group), fill=NA) + 
  geom_point(data=presettlement.pollen, aes(x=long, y=lat), shape = 16, size=2, color="red3") +
      scale_color_manual(values=c("black", "blue"), 
                         labels=c("Great Plains Region", "2 Decimal Degree Buffer"), 
                         name="Polygons")

plot_grid(presettlement.map, fossil.pollen.map, ncol=2)

#############################################################################
# Ambrosia to Artemisia Ratio

# Modern set: presettlement.pollen
# Drop spores

pres.pollen.clean2 <- presettlement.pollen[ , -which(names(presettlement.pollen) %in% c("EQUISETU", "PTERIDIUM", "LYCOPODX", "SPHAGNUM", "x.counts.good_counts..."))]
#reorder

pres.pollen.clean2 <- pres.pollen.clean2[,c(1:8,118,9:117)]

# Drop columns that are all na
pres.pollen.clean2 <- Filter(function(x)!all(is.na(x)), pres.pollen.clean2)

# Recode NA pollen observations to zero
pres.pollen.clean2[, 13:112][is.na(pres.pollen.clean2[, 13:112])] <- 0

pres.pollen.clean2$pollensum <- rowSums(pres.pollen.clean2[13:112])

# Drop cases where both Ambrosia and Artemisia are 0
pres.pollen.dropped <- pres.pollen.clean2 %>% 
  filter(AMBROSIA == 0 & ARTEMISIA == 0)

count(pres.pollen.clean2 %>% 
  filter(AMBROSIA == 0 & ARTEMISIA != 0)) # 75

count(pres.pollen.clean2 %>% 
  filter(AMBROSIA != 0 & ARTEMISIA == 0)) # 12

count(pres.pollen.clean2 %>% 
  filter(AMBROSIA != 0 & ARTEMISIA != 0)) # 706

count(pres.pollen.clean2 %>% 
  filter(AMBROSIA != 0 | ARTEMISIA != 0)) # 793

pres.pollen.clean2 <- pres.pollen.clean2 %>% 
  filter(AMBROSIA != 0 | ARTEMISIA != 0)

pres.pollen.clean2 <- pres.pollen.clean2 %>% 
  mutate(
    ambrosia.proportion = AMBROSIA/pollensum, 
    artemisia.proportion = ARTEMISIA/pollensum
  )

# Recode zero ambrosia records
pres.pollen.clean2$ambrosia.proportion[pres.pollen.clean2$ambrosia.proportion == 0] <- 0.00001
pres.pollen.clean2$artemisia.proportion[pres.pollen.clean2$artemisia.proportion == 0] <- 0.00001

# Calculate ambrosia to artemesia ratio
pres.pollen.clean2 <- pres.pollen.clean2 %>% 
  mutate(
    ambtoart = ambrosia.proportion/artemisia.proportion
  )

   
NaNcases2 <- sum(is.na(pres.pollen.clean2$ambtoart)) # 0 datasets
Infcases2 <- sum(is.infinite(pres.pollen.clean2$ambtoart)) # 0 datasets

# Remove infinite cases
# No longer needed
#pres.pollen.clean <- pres.pollen.clean[!is.infinite(pres.pollen.clean$ambtoart), ]

 ggplot() + geom_point(data=pres.pollen.clean2,
                       aes(x=long,y=lat,color=ambtoart)) +
  scale_color_gradient2(low="blue", mid="red", high=
                          "green", midpoint = mean(pres.pollen.clean2$ambtoart))

ggplot() + geom_point(data=pres.pollen.clean2,
                       aes(x=long,y=log(ambtoart),color=lat)) +
  scale_color_gradient(low="blue", high="green")

lm1 <- lm(ambrosia.proportion ~ artemisia.proportion, data=pres.pollen.clean2)
plot(lm1)

################### 
pres.pollen.clean2$uniqueid <- 1:nrow(pres.pollen.clean2)
# Point in poly to assign ecoregions to datasets
pres.pollen.sp <- pres.pollen.clean2[,c(10,11,116)]

ecoregions <- readOGR(dsn=getwd(), layer="ecoregions_NAD83") 

coordinates(pres.pollen.sp) <- ~ long + lat
proj4string(pres.pollen.sp) <- CRS(proj4string(ecoregions))

# Create the spatial polygons data frame


# The datsets get the region they are on top of

data.by.region <- over(pres.pollen.sp, ecoregions)

# Create a field to join on
data.by.region$uniqueid <- 1:nrow(data.by.region)
# Merge dataset location with region
pres.pollen.regions <- merge(pres.pollen.clean2, data.by.region, by="uniqueid")
obs.per.ecoregion <- table(droplevels(pres.pollen.regions$NA_L2NAME))

pres.pollen.regions$dummygr <- sample(1:13,793,replace=T) 

buffermap + geom_point(data=pres.pollen.regions, aes(x=long, y=lat), color="black") + 
     labs(title = "Presettlement pollen (250 - 500 CYA)", subtitle = "n = 793")
ggsave(file = "presettlementmap.png", dpi=900, scale=2)

xy <- pres.pollen.regions[,c(12,11)]
tic()
pres.pollen.regions.spdf <- SpatialPointsDataFrame(coords = xy, data = pres.pollen.regions,
                               proj4string = CRS(proj4string(clim.2degbuffer)))
toc()

# Join with climate data
# Using pres.pollen.regions and clim.2degbuffer from line 296 of Climate script
tic("compute distance matrix")
### compute the complete distance matrix between the two sets of points
dist_mat <- pointDistance(pres.pollen.regions.spdf, clim.2degbuffer, lonlat = TRUE, allpairs = TRUE)
toc()
tic("find nearest point")
### identify nearest point in dataset B for every point in dataset A
nearest <- apply(dist_mat, 1, which.min)
toc()
tic("bind data together")
### bind together the data from the dataset B (in your case the "red points")
### at the closest point to dataset A ("black points")
pres.pollen.regions.spdf@data<- cbind(pres.pollen.regions.spdf@data, clim.2degbuffer@data[nearest,])
toc()

pres.pollen.regions.df <- as.data.frame(pres.pollen.regions.spdf)
names(pres.pollen.regions.df)[names(pres.pollen.regions.df) %in%  c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8",                
                                      "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",               
                                      "bio17", "bio18", "bio19")] <- c("MeanT", 
                    "Mean.Diurnal.RangeT", "Isothermality", "TSeasonality",
                    "MaxTWarm", "MinTCold", "TAnnRange", "MeanTWetQ", "MeanTDryQ", "MeanTWarmQ",
                    "MeanTColdQ", "AnnP", "PrecipWetMonth", "PrecipDryMonth", 
                    "PrecipSeasonality", "PrecipWetQ", "PrecipDryQ", "PrecipWarmQ", "PrecipColdQ")


ratio <- as.matrix(pres.pollen.regions.df[,c(144,133,123)])

#############################
clim.2degbuffer.df <- as.data.frame(clim.2degbuffer)

templegend <- expression("Mean Annual Temperature " ( ~degree*C))
temp.plot <- ggplot() + theme_void() + geom_polygon(data=region2, 
                           aes(x=long, y=lat, group=group), colour="grey68", fill="grey95") + 
  geom_point(data=fortify(clim.2degbuffer.df), 
                        aes(x=x, y=y, colour=(bio1/10)), alpha=0.5) +
  scale_colour_gradient2(low = "#0571b0", mid = "#ffffbf", high = "#ca0020",
                         midpoint = mean((clim.2degbuffer.df$bio1/10)), guide = "colorbar", 
                         guide_legend("Mean Annual Temperature (°C)")) +
                          theme(legend.position="bottom") + 
  geom_polygon(data=surroundingregions, aes(x=long, y=lat, group=group), color="grey50", alpha=0.4, fill=NA) + 
  geom_polygon(data = gp.dissolve.f, 
                        aes(x=long, y = lat, group=group), fill=NA, colour = "black") + 
                geom_polygon(data=fortify(buffer), aes(x=long, y=lat, group=group),
                           fill=NA, color = "blue") +
  coord_map(xlim = c(-128, -68),ylim = c(25, 59)) +
  theme(legend.position = "bottom") + xlab("Longitude") + ylab("Latitude") +
            theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))

  ggsave(file = "tempplot.jpg", dpi = 400, scale=2)
          
  

ggplot() +
  geom_point(data=fortify(clim.2degbuffer.df), 
                        aes(x=x, y=y, colour=(bio12))) +
  scale_colour_gradient(low = "#ece7f2", high = "#023858", 
                         guide = "colorbar", 
                         guide_legend(title = "Mean Annual Precipitation (mm)")) + 
                        theme(legend.position="bottom") + 
  geom_polygon(data=region2, aes(x=long, y=lat, group=group), colour="grey60", fill=NA) + 

 coord_map(xlim = c(-128, -68),ylim = c(24, 58)) + theme_bw(base_size = 12) +
  xlab("Longitude") + ylab("Latitude") +
            theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_segment(aes(x=-69, y=26.1, xend=-69, yend=27.1), size = 2, arrow = arrow(length = unit(0.25, "cm"))) +
  annotate("text", x=-69, y=25.6, label="N")
 
  ggsave(file = "precipplot.jpg", dpi = 600, scale=2)

##################################################################
# Drop region with single case

#pres.pollen.regions2 <- pres.pollen.regions2[!grepl("CENTRAL USA PLAINS", pres.pollen.regions2$NA_L2NAME),]

# Remove Western Cordillera at reviewer request
pres.pollen.regions2 <- pres.pollen.regions.df[!grepl("UPPER GILA MOUNTAINS", pres.pollen.regions.df$NA_L2NAME),]    
pres.pollen.regions2 <- pres.pollen.regions2[!grepl("WESTERN CORDILLERA", pres.pollen.regions2$NA_L2NAME),]
pres.pollen.regions2 <-  pres.pollen.regions2[!grepl("COLD DESERTS", pres.pollen.regions2$NA_L2NAME),]

#pres.pollen.regions2 <- pres.pollen.regions[droplevels(pres.pollen.regions$NA_L2NAME)]
pca1 <- prcomp(pres.pollen.regions2[,c(11,12,117)], scale=T)
reduced.pollen <- (pres.pollen.regions[pres.pollen.regions$NA_L2NAME == "WEST-CENTRAL SEMIARID PRAIRIES" |
                                       pres.pollen.regions$NA_L2NAME == "TEMPERATE PRAIRIES", ])
              
pca.reduced <- prcomp(reduced.pollen[,c(11,12,117)], scale=F)
gr <- factor(pres.pollen.regions2[,119])
gr.reduced <- factor(reduced.pollen[,119])
gr <- factor(droplevels(pres.pollen.regions2$NA_L2NAME))
summary(gr)
pca3d(pca1, show.ellipses = TRUE, ellipse.ci = 0.95, group=gr, show.group.labels = TRUE)

pca2d(pca1, show.ellipses = TRUE, ellipse.ci = 0.95, group=gr)
pca3d(pca.reduced, show.ellipses = TRUE, ellipse.ci = 0.95, group=gr.reduced, show.group.labels = TRUE)
with(pres.pollen.clean, plot3d(x = long, y = lat, z = ambtoart))

##########################################################
# Repeat analyses for all surface samples in the region
greatplains.surface <- get_dataset(loc = c(-116.5924, 21.1635, -88.14369, 56.11235), 
                      datasettype = 'pollen surface sample')
tic("download")
gplains.surface.d <- get_download(greatplains.surface)
tic("compile taxa")
gplains.surface.taxa <- compile_taxa(gplains.surface.d, "WhitmoreFull")
tic("compile downloads")
gplains.surface.c <- compile_downloads(gplains.surface.taxa)
tic("write to file")
toc()

random <- get_dataset(loc = c(-100, 30, -90, 40), datasettype = "pollen surface sample")
random.d <- get_download(random)
r.comp <- compile_taxa(random.d, "WS64")
r.down <- compile_downloads(r.comp)

# Clip dataset to buffer.nad83
gplains.surface.sp <- gplains.surface.c
coordinates(gplains.surface.sp) <- ~ long + lat
proj4string(gplains.surface.sp) <- CRS(proj4string(buffer.nad83))

# Clip to 2 degree buffer
tic()
gplains.surface.buffer <- gplains.surface.sp[!is.na(over(gplains.surface.sp,as(buffer.nad83,"SpatialPolygons"))),]
write.csv(file = "gplains 2degbuffer surface whitmore full 10-23-2018.csv", x=gplains.surface.buffer, row.names=FALSE)
toc()
#################################33

# Ambrosia to Artemisia Ratio

# Surface set: gplains.surface.buffer
gplains.s.df <- read.csv("gplains 2degbuffer surface whitmore full 10-23-2018.csv")

# Drop spores

surface.pollen.clean <- gplains.s.df[ , -which(names(gplains.s.df) %in% c(".id", "EQUISETU", "PTERIDIUM", "LYCOPODX", "SPHAGNUM", "x.counts.good_counts...", "optional"))]

# Drop columns that are all na
surface.pollen.clean2 <- surface.pollen.clean[, unlist(lapply(surface.pollen.clean, function(x) !all(is.na(x))))]

# Recode NA pollen observations to zero
surface.pollen.clean2[, 8:94][is.na(surface.pollen.clean2[, 8:94])] <- 0

# Calculate total pollen count
surface.pollen.clean2$pollensum <- rowSums(surface.pollen.clean2[8:94])

surface.pollen.clean2 <- surface.pollen.clean2 %>% 
  mutate(
    ambrosia.proportion = AMBROSIA/pollensum, 
    artemisia.proportion = ARTEMISIA/pollensum
  )

# Drop cases where both Ambrosia and Artemisia are 0
surface.pollen.clean2 <- surface.pollen.clean2[(!(surface.pollen.clean2$AMBROSIA == 0) | 
                                                !(surface.pollen.clean2$ARTEMISIA == 0)),]

# Recode zero ambrosia records
surface.pollen.clean2$ambrosia.proportion[surface.pollen.clean2$ambrosia.proportion == 0] <- 0.00001
surface.pollen.clean2$artemisia.proportion[surface.pollen.clean2$artemisia.proportion == 0] <- 0.00001

# Do spatial join

ecoregions <- readOGR(dsn=getwd(), layer="ecoregions_NAD83") 

# Point in poly to assign ecoregions to datasets
surface.pollen.clean2$uniqueid <- 1:nrow(surface.pollen.clean2)
surface.pollen.sp <- surface.pollen.clean2[,c(5,6,98)]
coordinates(surface.pollen.sp) <- ~ long + lat
proj4string(surface.pollen.sp) <- CRS(proj4string(ecoregions))

# Create the spatial polygons data frame


# The datsets get the region they are on top of

surface.data.by.region <- over(surface.pollen.sp, ecoregions)

# Create a field to join on
surface.data.by.region$uniqueid <- 1:nrow(surface.data.by.region)
# Merge dataset location with region
surface.pollen.regions <- merge(surface.pollen.clean2, surface.data.by.region, by="uniqueid")
s.obs.per.ecoregion <- table(droplevels(surface.pollen.regions)$NA_L2NAME)

# Remove Western Cordillera at reviewer request
s.pollen.regions2 <- surface.pollen.regions[!grepl("WESTERN CORDILLERA", surface.pollen.regions$NA_L2NAME),]
# 538 to 512

# Drop regions with few cases

s.pollen.regions2 <- s.pollen.regions2[!grepl("COLD DESERTS", s.pollen.regions2$NA_L2NAME),]
# 512 to 507

# Create column to facet by
s.pollen.regions2$category <- NULL
s.pollen.regions2$category[s.pollen.regions2$AMBROSIA == 0] <- "NoAmbrosia"
s.pollen.regions2$category[s.pollen.regions2$ARTEMISIA == 0] <- "NoArtemisia"
s.pollen.regions2$category[s.pollen.regions2$AMBROSIA != 0 & s.pollen.regions2$ARTEMISIA !=0] <- "BothPresent"

surface.xy <- s.pollen.regions2[,c(7,6)]

surface.pollen.regions.spdf <- SpatialPointsDataFrame(coords = surface.xy, data = s.pollen.regions2,
                               proj4string = CRS(proj4string(clim.2degbuffer)))


# Join with climate data
# Using surface.pollen.regions and clim.2degbuffer from line 296 of Climate script
tic("compute distance matrix")
### compute the complete distance matrix between the two sets of points
s.dist_mat <- pointDistance(surface.pollen.regions.spdf, clim.2degbuffer, lonlat = TRUE, allpairs = TRUE)
toc()
tic("find nearest point")
### identify nearest point in dataset B for every point in dataset A
s.nearest <- apply(s.dist_mat, 1, which.min)
toc()
tic("bind data together")
### bind together 
surface.pollen.regions.spdf@data<- cbind(surface.pollen.regions.spdf@data, clim.2degbuffer@data[s.nearest,])
toc()

surface.pollen.regions.df <- as.data.frame(surface.pollen.regions.spdf)
names(surface.pollen.regions.df)[names(surface.pollen.regions.df) %in%  c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8",                
                                      "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",               
                                      "bio17", "bio18", "bio19")] <- c("MeanT", 
                    "Mean.Diurnal.RangeT", "Isothermality", "TSeasonality",
                    "MaxTWarm", "MinTCold", "TAnnRange", "MeanTWetQ", "MeanTDryQ", "MeanTWarmQ",
                    "MeanTColdQ", "AnnP", "PrecipWetMonth", "PrecipDryMonth", 
                    "PrecipSeasonality", "PrecipWetQ", "PrecipDryQ", "PrecipWarmQ", "PrecipColdQ")

write.csv(surface.pollen.regions.df, file="surpollenregions.csv", row.names = FALSE)

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


buffermap + geom_point(data=surface.pollen.corrected, aes(x=long, y=lat), color="black") + 
     labs(title = "Surface pollen", subtitle = "n = 507")
ggsave(file = "surfacemap.png", dpi=900, scale=2)
   
NaNcases2 <- sum(is.na(surface.pollen.corrected$ambtoart))
Infcases2 <- sum(is.infinite(surface.pollen.corrected$ambtoart))



 ggplot() + geom_point(data=surface.pollen.corrected,
                       aes(x=long,y=lat,color=ambtoart)) +
  scale_color_gradient2(low="blue", mid="red", high=
                          "green", midpoint = median(surface.pollen.corrected$ambtoart))

ggplot() + geom_point(data=surface.pollen.corrected,
                       aes(x=long,y=ambtoart,color=lat)) +
  scale_color_gradient(low="blue", high="green")







allregionmap.fade <- ggplot() + geom_polygon(data=region2, 
                           aes(x=long, y=lat, group=group), colour="grey68", fill="grey98") + 
  geom_polygon(data=surroundingregions, aes(x=long, y=lat, group=group, 
               fill=surroundingregions$NA_L2KEY), color="gray70", alpha=0.2) +
  scale_fill_manual(values = c("5.2  MIXED WOOD SHIELD" = "#aedee4",
                               "5.4  BOREAL PLAIN" = "#bed7c0",
                               "6.2  WESTERN CORDILLERA" = "#5ebc55",
                               "8.1  MIXED WOOD PLAINS" = "#a4cb4b",
                               "8.2  CENTRAL USA PLAINS" = "#e1f1e7",
                               "8.3  SOUTHEASTERN USA PLAINS" = "#cfe4ae",
                               "8.4  OZARK/OUACHITA-APPALACHIAN FORESTS" = "#51bd83",
                               "8.5  MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS" = "#b9d87e",
                               "9.2  TEMPERATE PRAIRIES" = "#fee5ca",
                               "9.3  WEST-CENTRAL SEMIARID PRAIRIES" = "#fadca1",
                               "9.4  SOUTH CENTRAL SEMIARID PRAIRIES" = "#f8ce9e",
                               "9.6  TAMAULIPAS-TEXAS SEMIARID PLAIN"  = "#fcb94b",
                               "10.1  COLD DESERTS" = "#ebed90",
                               "10.2  WARM DESERTS" = "#fadd69",
                               "13.1  UPPER GILA MOUNTAINS" = "#b9d87e"),
                    breaks = c("5.2  MIXED WOOD SHIELD",
                               "5.4  BOREAL PLAIN",
                               "6.2  WESTERN CORDILLERA",
                               "8.1  MIXED WOOD PLAINS",
                               "8.2  CENTRAL USA PLAINS",
                               "8.3  SOUTHEASTERN USA PLAINS",
                               "8.4  OZARK/OUACHITA-APPALACHIAN FORESTS",
                               "8.5  MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS",
                               "9.2  TEMPERATE PRAIRIES",
                               "9.3  WEST-CENTRAL SEMIARID PRAIRIES",
                               "9.4  SOUTH CENTRAL SEMIARID PRAIRIES",
                               "9.6  TAMAULIPAS-TEXAS SEMIARID PLAIN",
                               "10.1  COLD DESERTS",
                               "10.2  WARM DESERTS",
                               "13.1  UPPER GILA MOUNTAINS"),
                    labels = c("5.2  Mixed wood shield",
                               "5.4  Boreal plain",
                               "6.2  Western cordillera",
                               "8.1  Mixed wood plain",
                               "8.2  Central USA plains",
                               "8.3  Southeastern USA plains",
                               "8.4  Ozark/Ouachita-Appalachian forests",
                               "8.5  Mississippi alluvial and southeast USA coastal plains",
                               "9.2  Temperate prairies",
                               "9.3  West-central semiarid prairies",
                               "9.4  South central semiarid prairies",
                               "9.6  Tamaulipas-Texas semiarid plain",
                               "10.1  Cold deserts",
                               "10.2  Warm deserts",
                               "13.1  Upper Gila Mountains"), 
                    name="EPA Level II Ecoregions") + 
  geom_text(data=droplevels(centroids.df3), aes(x=long, y=lat, label=code)) +
  coord_quickmap(xlim=c(-117, -87), ylim=c(29,57)) + theme_bw(base_size = 12) +
  theme(legend.position = "right") + xlab("Longitude") + ylab("Latitude") +
            theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14)) + 
  geom_segment(aes(x=-86.1, y=28.2, xend=-86.1, yend=28.7), size = 2, arrow = arrow(length = unit(0.25, "cm"))) +
  annotate("text", x=-86.1, y=28, label="N")
allregionmap.fade

allregionmap.fade + geom_point(data=surface.pollen.corrected, aes(x=long, y=lat, color=log(ambtoart), shape=category), size=3.25) +
  scale_color_viridis(name=expression(paste("Log", italic(" Ambrosia"), " to", italic(" Artemisia"), " Ratio")), option="magma", end=0.9,
                      guide = guide_colourbar(direction = "horizontal",title.position = "top", barwidth=12)) + 
  scale_shape_manual(name="Category", values=c(19,15,17), labels=c("Both Present", expression(paste("No", italic(" Ambrosia"))),
                                                                   expression(paste("No", italic(" Artemisia")))))
ggsave(file="allregion.ratiomap.png", dpi=600, scale=2)

category.by.region <- pres.pollen.regions2 %>%
  group_by(category) %>%
  summarize(n())

noartcat <- pres.pollen.regions2 %>%
 # dplyr::filter(category == "NoArtemisia") %>%
  group_by(dataset) %>%
  summarize(n())

write.csv(category.by.region, file="catbyregion.csv", row.names = FALSE)
surface.xy <- surface.pollen.regions[,c(7,6)]
tic()
surface.pollen.regions.spdf <- SpatialPointsDataFrame(coords = surface.xy, data = surface.pollen.regions,
                               proj4string = CRS(proj4string(clim.2degbuffer)))
toc()

# Join with climate data
# Using surface.pollen.regions and clim.2degbuffer from line 296 of Climate script
tic("compute distance matrix")
### compute the complete distance matrix between the two sets of points
s.dist_mat <- pointDistance(surface.pollen.regions.spdf, clim.2degbuffer, lonlat = TRUE, allpairs = TRUE)
toc()
tic("find nearest point")
### identify nearest point in dataset B for every point in dataset A
s.nearest <- apply(s.dist_mat, 1, which.min)
toc()
tic("bind data together")
### bind together the data from the dataset B (in your case the "red points")
### at the closest point to dataset A ("black points")
surface.pollen.regions.spdf@data<- cbind(surface.pollen.regions.spdf@data, clim.2degbuffer@data[s.nearest,])
toc()

surface.pollen.regions.df <- as.data.frame(surface.pollen.regions.spdf)
names(surface.pollen.regions.df)[names(surface.pollen.regions.df) %in%  c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8",                
                                      "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",               
                                      "bio17", "bio18", "bio19")] <- c("MeanT", 
                    "Mean.Diurnal.RangeT", "Isothermality", "TSeasonality",
                    "MaxTWarm", "MinTCold", "TAnnRange", "MeanTWetQ", "MeanTDryQ", "MeanTWarmQ",
                    "MeanTColdQ", "AnnP", "PrecipWetMonth", "PrecipDryMonth", 
                    "PrecipSeasonality", "PrecipWetQ", "PrecipDryQ", "PrecipWarmQ", "PrecipColdQ")

write.csv(surface.pollen.regions.df, file="surpollenregions.csv", row.names = FALSE)



# Ratio in relationship to precip/temp
ggplot() + geom_point(data=surface.pollen.corrected,
                       aes(x=MeanT,y=ambtoart,color=AnnP)) +
  scale_color_gradient(low="blue", high="green") + ggtitle("Surface Pollen")
ggplot() + geom_point(data=surface.pollen.corrected,
                       aes(x=MeanT,y=ambtoart,color=AnnP)) +
   scale_color_gradient(low="blue", high="green") + ggtitle("Modern Presettlement Pollen")


# Drop regions with few observations to prepare for Discriminant Function Analysis

s.pollen.regions2 <- surface.pollen.corrected[!grepl("TEXAS-LOUISIANA COASTAL PLAIN", surface.pollen.corrected$NA_L2NAME),]

s.pollen.regions2 <- s.pollen.regions2[!grepl("SOFTWOOD SHIELD", s.pollen.regions2$NA_L2NAME),]

#s.pollen.regions2 <- s.pollen.regions2[!grepl("MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS", s.pollen.regions2$NA_L2NAME),]

colvar<- as.array(c(s.pollen.regions2$MeanT, s.pollen.regions2$AnnP, log(s.pollen.regions2$ambtoart)))
group.cols <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")



p.colors <- viridis(10, option = "D")
regions <- factor(unique(droplevels(s.pollen.regions2$NA_L2NAME)))

f <- list(
  size = 16,
  color = "black"
)
x <- list(
  title = "MAP (mm)",
  titlefont = f
)
y <- list(
  title = "MAT (\u00B0C)",
  titlefont = f
)
z <- list(
  title = "Log Ratio",
  titlefont = f
)
scene = list(xaxis=x,
             yaxis=y,
             zaxis=z)

colormatch <- data.frame(regions,p.colors, stringsAsFactors = FALSE)
colormatch2 <- colormatch[match(colormatch$regions, colormatch$p.colors), 'p.colors']
presettlement.set$plotcol <- colormatch[match(colormatch$regions, colormatch$p.colors), 'p.colors']

legend.layout <- list(size = 14)

p.presettlement <- plot_ly(presettlement.set, x=~(AnnP/10), y=~(MeanT/10), z=~ambtoart,
       #    width=900, height=600, 
        type="scatter3d", mode="markers", color=~droplevels(NA_L2NAME), colors=p.colors[c(1,3,4,6:8)],
        hoverinfo = 'text',
        text = ~paste('</br>', NA_L2NAME,
                      '</br> Precip: ', paste(round((AnnP),2)),
                      '</br> Temp: ', paste(round((MeanT),2)),
                      '</br> Ratio: ', paste(round(ambtoart,2)))) %>%
    layout(title="Presettlement Set (250-500 cal BP)", scene=scene, legend=list(font =list(size=18)), autosize=T)

p.surface <- plot_ly(s.pollen.regions2, x=~(AnnP), y=~(MeanT/10), z=~log(ambtoart),
     #   width=1200, height=600, 
        type="scatter3d", mode="markers", color=~droplevels(NA_L2NAME), colors=p.colors,
        hoverinfo = 'text',
        text = ~paste('</br>', NA_L2NAME,
                      '</br> Precip: ', paste(round((AnnP),2), '(mm/year)'),
                      '</br> Temp: ', paste(round((MeanT),2), '(\u00B0C)'),
                      '</br> Log Ratio: ', paste(round(ambtoart,2)))) %>%
    layout(title=NULL, scene=scene, legend=list(font =list(size=14)), autosize=T)

p.surface.annp <- plot_ly(surface.set2, x=~(AnnP), y=~log(ambtoart),
     #   width=1200, height=600, 
        type="scatter", mode="markers", color=~droplevels(NA_L2NAME), colors=p.colors,
        hoverinfo = 'text',
        text = ~paste('</br>', NA_L2NAME,
                      '</br> Precip: ', paste(round((AnnP),2), '(mm/year)'),
                      '</br> Temp: ', paste(round((MeanT),2), '(\u00B0C)'),
                      '</br> Log Ratio: ', paste(round(ambtoart,2)))) %>%
    layout(scene=scene, legend=list(font =list(size=14)), autosize=T)

htmlwidgets::saveWidget(p.surface, 'Appendix S4 Fig S1.html')

annp <- ggplot(s.pollen.regions2) + geom_point(aes(x=AnnP, y=log(ambtoart), color=droplevels(NA_L2NAME))) +
  scale_color_manual(values=p.colors, name="EPA Level II Ecoregion",
                     labels = c(
                               "Boreal plain",
                               "Central USA plains",
                               "MS alluvial and SE USA coastal plains",
                               "Mixed wood plains",
                               "Mixed wood shield",
                               "Ozark/Ouachita-Appalachian forests",
                               "South central semiarid prairies",
                               "Southeastern USA plains",
                               "Temperate prairies",
                               "West-central semiarid prairies"
                                )) + theme_bw() +
  labs(y= expression(paste("Log", italic(" Ambrosia")," to ", italic("Artemisia")," ratio")),
       x= "Mean Annual Precipitation (mm)") + ylim(c(-12,10.5))+ theme(plot.margin = margin(0, 0.5, 0, 0)) + xlim(c(200,1650))

meant <- ggplot(s.pollen.regions2) + geom_point(aes(x=MeanT/10, y=log(ambtoart), color=droplevels(NA_L2NAME))) +
  scale_color_manual(values=p.colors, name="EPA Level II Ecoregion",
                     labels = c(
                               "Boreal plain",
                               "Central USA plains",
                               "MS alluvial and SE USA coastal plains",
                               "Mixed wood plains",
                               "Mixed wood shield",
                               "Ozark/Ouachita-Appalachian forests",
                               "South central semiarid prairies",
                               "Southeastern USA plains",
                               "Temperate prairies",
                               "West-central semiarid prairies"
                               )) + theme_bw() +
  labs(y= NULL, x= expression("Mean Annual Temperature " ( degree*C))) + ylim(c(-12,10.5)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + theme(plot.margin = margin(0, 0, 0, -0.8))


library(cowplot)
prow <- plot_grid(annp + theme(legend.position="none"),
                  meant + theme(legend.position="none"),
  align = 'vh',
  #labels = c("b", "c"),
  hjust= -5, vjust= 2,
  nrow = 1)
legend <- get_legend(
  # create some space to the left of the legend
  meant + theme(legend.box.margin = margin(0, 0, 0, 2))
)
# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
comboplot <- plot_grid(prow, legend, rel_widths = c(7,3))
save_plot("comboplot.png", comboplot, base_width = 9)

precip.ratio <- plot_ly(surface.set, x=~AnnP, y=~MeanT, 
        type="scatter", mode="markers", color=~(MeanT/AnnP),
        hoverinfo = 'text',
        text = ~paste('</br> Precip: ', paste(round(AnnP,2)),
                      '</br> Temp: ', paste(round(MeanT,2)))) %>%
       layout(title="Climate space of the Great Plains region", xaxis = x, yaxis = y, titlefont=f,
              width=100, height=700,
              legend=list(title="Ratio",font =list(size=16)))

htmlwidgets::saveWidget(p.surface,"surface.html", selfcontained = TRUE)
htmlwidgets::saveWidget(p.presettlement,"presettlement.html", selfcontained = TRUE)

library(processx)
export(precip.ratio, "precip.ratio.png")

htmlwidgets::createWidget(p.surface,"surfacepollen.html",
  sizingPolicy = htmlwidgets::sizingPolicy(
    viewer.padding = 0,
    viewer.paneHeight = 500,
    browser.fill = TRUE
  ))
htmlwidgets::saveWidget(widget=p.presettlement, "presettlementpollen.html", selfcontained = TRUE)

### Rotating graph

args <- commandArgs(TRUE)
df <- read.csv(args[1], sep=",", check.names=FALSE, row.names = 1)
 args <- c("Documents/Warwick/PhD/Structural\ Work/Tail\ Fibres/CD\ Spectra/lumt.csv")
file <- strsplit(basename(args[1]), '[.]')[[1]]

for(i in seq(0,6.3,by=0.1)){
outfile <- paste(file[1],"surface3d",round(i,digits=2), sep = "_")
cam.zoom = 2
ver.angle = 0
cat("Now rendering iteration:", i,"\n")
plotly_IMAGE(p.surface,
         width = 1200,
         height = 1050,
         format = "gif",
         username="xxx",
         key="xxxx",
         scale = 1,
         out_file = paste(outfile,"png", sep="."))
}

## Subplots (equivalent to facet)
subplot(
        plot_ly(x = 1:3, y = 10*(1:3), name = "slope of 10") %>%
            add_trace(x = 2:4, y = 1:3, name = "slope of 1", yaxis = "y2") %>%
            layout(title = "Double Y Axis", yaxis2 = ay, xaxis = ax, yaxis = ay1), 
        plot_ly(x = 1:3, y = 10*(1:3), name = "slope of 10") %>%
            add_trace(x = 2:4, y = 1:3, name = "slope of 1", yaxis = "y2") %>%
            layout(title = "Double Y Axis", yaxis2 = ay, xaxis = ax, yaxis = ay1), nrows = 2)
##################################################################3
# Discriminant analysis

# Surface Pollen
# Use s.pollen.regions2

# Get only the columns we're worred about
surface.set <- s.pollen.regions2[,c(6,7,99,101,108,119)]
# Split surface samples into training set and test set
set.seed(123)
training.samples <- droplevels(surface.set$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- surface.set[training.samples, ]
test.data <- surface.set[-training.samples, ]

# Estimate preprocessing parameters
preproc.param <- train.data %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)

# Fit the model
model <- fda(droplevels(NA_L2NAME)~., data = train.transformed)
# Make predictions
predicted.classes <- model %>% predict(test.transformed)
# Model accuracy
mean(predicted.classes == droplevels(test.transformed$NA_L2NAME))


######################## Try this again without lat long data
# Get only the columns we're worred about
surface.set2 <- s.pollen.regions2[,c(99,101,108,119)]
# Split surface samples into training set and test set
set.seed(123)
training.samples2 <- droplevels(surface.set2$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data2 <- surface.set2[training.samples2, ]
test.data2 <- surface.set2[-training.samples2, ]
summary(test.data2)

# Estimate preprocessing parameters
preproc.param2 <- train.data2 %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed2 <- preproc.param2 %>% predict(train.data2)
test.transformed2 <- preproc.param2 %>% predict(test.data2)

# Fit the model
model2 <- fda(droplevels(NA_L2NAME)~., data = train.transformed2)
# Make predictions
predicted.classes2 <- model2 %>% predict(test.transformed2)
# Model accuracy
mean(predicted.classes2 == droplevels(test.transformed2$NA_L2NAME))
plot(model2, train.transformed2, group="pred")
plot(model2, train.transformed2, group="true")


######################## Try this again with only climate data
# Get only the columns we're worred about
surface.set3 <- s.pollen.regions2[,c(99,101,129)]
surface.set3$ambtoart <- log(surface.set3$ambtoart)

summary(droplevels(surface.set3$NA_L2NAME))
# Split surface samples into training set and test set
set.seed(123)
training.samples3 <- droplevels(surface.set3$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data3 <- surface.set3[training.samples3, ]
test.data3 <- surface.set3[-training.samples3, ]

# Estimate preprocessing parameters
preproc.param3 <- train.data3 %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed3 <- preproc.param3 %>% predict(train.data3)
test.transformed3 <- preproc.param3 %>% predict(test.data3)

# Fit the model
model3 <- fda(droplevels(NA_L2NAME)~., data = train.transformed3)
# Make predictions
predicted.classes3 <- model3 %>% predict(test.transformed3)
# Model accuracy
mean(predicted.classes3 == droplevels(test.transformed3$NA_L2NAME))

surfacetrainmap <- buffermap + geom_point(data=train.data, aes(x=long, y=lat), color="black") +
  ggtitle("Training Set (n=348)")

ggsave(file = "surfacetrainmap.jpg", dpi = 600, scale=2)

surfacetestmap <- buffermap + geom_point(data=test.data, aes(x=long, y=lat), color="black") +
  ggtitle("Test Set (n=85)")

ggsave(file = "surfacetestmap.jpg", dpi = 600, scale=2)

###################

aov.surf <- (aov(log(ambtoart) ~ NA_L2NAME, data=s.pollen.regions2))
tukey.surf <- TukeyHSD(aov.surf)
df.tukey <- as.data.frame(tukey.surf[["NA_L2NAME"]])
summary(aov.surf)
write.csv(df.tukey, file="surfacetukey.csv")


##################################################################3
# Discriminant analysis

# Presettlement Pollen
# Use pres.pollen.regions2, but drop Ozarks and SE Plains
presettlement.set <- pres.pollen.regions2
presettlement.set <- pres.pollen.regions2[!grepl("OZARK/OUACHITA-APPALACHIAN FORESTS", pres.pollen.regions2$NA_L2NAME),]
presettlement.set <- presettlement.set[!grepl("SOUTHEASTERN USA PLAINS", presettlement.set$NA_L2NAME),]
# Get only the columns we're worred about
presettlement.set <- presettlement.set[,c(117,119,126,137)]
# Split surface samples into training set and test set
set.seed(123)
p.training.samples <- droplevels(presettlement.set$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
p.train.data <- presettlement.set[p.training.samples, ]
p.test.data <- presettlement.set[-p.training.samples, ]

# Estimate preprocessing parameters
p.preproc.param <- p.train.data %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
p.train.transformed <- p.preproc.param %>% predict(p.train.data)
p.test.transformed <- p.preproc.param %>% predict(p.test.data)

# Fit the model
p.model <- fda(droplevels(NA_L2NAME)~., data = p.train.transformed)
# Make predictions
p.predicted.classes <- p.model %>% predict(p.test.transformed)
# Model accuracy
mean(p.predicted.classes == droplevels(p.test.transformed$NA_L2NAME))

######################## Try this again without lat long data
# Get only the columns we're worred about
surface.set2 <- s.ambart.props[,c(5:90,92,92)]
# Split surface samples into training set and test set
set.seed(123)
training.samples2 <- droplevels(surface.set2$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data2 <- surface.set2[training.samples2, ]
test.data2 <- surface.set2[-training.samples2, ]

# Estimate preprocessing parameters
preproc.param2 <- train.data2 %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed2 <- preproc.param2 %>% predict(train.data2)
test.transformed2 <- preproc.param2 %>% predict(test.data2)

# Fit the model
model2 <- fda(droplevels(NA_L2NAME)~., data = train.transformed2)
# Make predictions
predicted.classes2 <- model2 %>% predict(test.transformed2)
# Model accuracy
mean(predicted.classes2 == droplevels(test.transformed2$NA_L2NAME))
plot(model2)

######################## Try this again with only climate data
# Get only the columns we're worred about
surface.set3 <- s.pollen.regions2 %>%
  select("NA_L2NAME", "AnnP", "ambtoart")
surface.set3$ambtoart <- log(surface.set3$ambtoart)
# Split surface samples into training set and test set
set.seed(123)
training.samples3 <- droplevels(surface.set3$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data3 <- surface.set3[training.samples3, ]
test.data3 <- surface.set3[-training.samples3, ]

# Estimate preprocessing parameters
preproc.param3 <- train.data3 %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed3 <- preproc.param3 %>% predict(train.data3)
test.transformed3 <- preproc.param3 %>% predict(test.data3)

# Fit the model
model3 <- fda(droplevels(NA_L2NAME)~., data = train.transformed3, method = mars)
# Make predictions
predicted.classes3 <- model3 %>% predict(test.transformed3)
# Model accuracy
mean(predicted.classes3 == droplevels(test.transformed3$NA_L2NAME))
plot(model3)

# Run 1,000 times

surf.means <- NULL

surf.replicate <- function(df) {
  training.samples3 <- droplevels(df$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data3 <- surface.set3[training.samples3, ]
test.data3 <- surface.set3[-training.samples3, ]

# Estimate preprocessing parameters
preproc.param3 <- train.data3 %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed3 <- preproc.param3 %>% predict(train.data3)
test.transformed3 <- preproc.param3 %>% predict(test.data3)

# Fit the model
model3 <- fda(droplevels(NA_L2NAME)~., data = train.transformed3, method = mars)
# Make predictions
predicted.classes3 <- model3 %>% predict(test.transformed3)
# Model accuracy
mean(predicted.classes3 == droplevels(test.transformed3$NA_L2NAME))
}
surf.replicate(surface.set3)
tic()
surf.means <- replicate(1000, surf.replicate(surface.set3))
summary(surf.means)
toc()

surf.replicate2 <- function(df) {
  training.samples3 <- droplevels(df$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data3 <- df[training.samples3, ]
test.data3 <- df[-training.samples3, ]

# Estimate preprocessing parameters
preproc.param3 <- train.data3 %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed3 <- preproc.param3 %>% predict(train.data3)
test.transformed3 <- preproc.param3 %>% predict(test.data3)

# Fit the model
model3 <- fda(droplevels(NA_L2NAME)~., data = train.transformed3, method = mars)
# Make predictions
predicted.classes3 <- model3 %>% predict(test.transformed3)
# Model accuracy
mean(predicted.classes3 == droplevels(test.transformed3$NA_L2NAME))
}

surf.replicate2(s.pollen.minus.noart)
surf.means2 <- replicate(1000, surf.replicate(s.pollen.minus.noart))
summary(surf.means2)
histogram(surf.means2)

# Get only the columns we're worred about
surface.full <- s.ambart.props[,c(5:89)]
# Split surface samples into training set and test set
set.seed(123)
training.samples.full <- droplevels(surface.full$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data.full <- surface.full[training.samples.full, ]
test.data.full <- surface.full[-training.samples.full, ]

# Estimate preprocessing parameters
preproc.param.full <- train.data.full %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed.full <- preproc.param.full %>% predict(train.data.full)
test.transformed.full <- preproc.param.full %>% predict(test.data.full)

# Fit the model
model.full <- fda(droplevels(NA_L2NAME)~., data = train.transformed.full)
# Make predictions
predicted.classes.full <- model.full %>% predict(test.transformed.full)
# Model accuracy
mean(predicted.classes.full == droplevels(test.transformed.full$NA_L2NAME))
plot(model.full)


# Get only the columns we're worred about
surface.setclim <- s.pollen.regions2[,c(101,119)]
# Split surface samples into training set and test set
set.seed(123)
training.samplesclim <- droplevels(surface.setclim$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
train.dataclim <- surface.setclim[training.samplesclim, ]
test.dataclim <- surface.setclim[-training.samplesclim, ]

# Estimate preprocessing parameters
preproc.paramclim <- train.dataclim %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformedclim <- preproc.paramclim %>% predict(train.dataclim)
test.transformedclim <- preproc.paramclim %>% predict(test.dataclim)

# Fit the model
modelclim <- fda(droplevels(NA_L2NAME)~., data = train.transformedclim)
# Make predictions
predicted.classesclim <- modelclim %>% predict(test.transformedclim)
# Model accuracy
mean(predicted.classesclim == droplevels(test.transformedclim$NA_L2NAME))
plot.fda(modelclim)
perce
#######################################
# Regression



pres.pollen.regions2$category <- NULL
pres.pollen.regions2$category[pres.pollen.regions2$AMBROSIA == 0] <- "NoAmbrosia"
pres.pollen.regions2$category[pres.pollen.regions2$ARTEMISIA == 0] <- "NoArtemisia"
pres.pollen.regions2$category[pres.pollen.regions2$AMBROSIA != 0 & pres.pollen.regions2$ARTEMISIA !=0] <- "BothPresent"


ggplot(data=surface.pollen.regions.df, aes(x=(AnnP), y=log(ambtoart), group=NA_L2NAME, color=NA_L2NAME)) + geom_point() + geom_smooth() +
  facet_wrap(vars(NA_L2NAME), nrow=4) + theme(legend.position = "none")

aggregate(log(surface.pollen.corrected$ambtoart), list(surface.pollen.corrected$category), range)
aggregate(log(surface.pollen.corrected$ambtoart), list(surface.pollen.corrected$category), mean)
aggregate(log(surface.pollen.corrected$ambtoart), list(surface.pollen.corrected$category), summary)
aggregate((surface.pollen.corrected$AnnP), list(surface.pollen.corrected$category), range)
aggregate((surface.pollen.corrected$AnnP), list(surface.pollen.corrected$category), mean)
surface.pollen.corrected %>% 
  group_by(category) %>%
  summarise(no_rows = length(category))

surface.pollen.corrected %>% 
  group_by(NA_L2NAME) %>%
  summarise(no_rows = length(category))

pres.pollen.regions2 %>% 
  group_by(NA_L2NAME) %>%
  summarise(no_rows = length(category))

hist(surface.pollen.regions.df$AnnP[surface.pollen.regions.df$category == "NoAmbrosia"], probability=TRUE, ylim=c(0,0.015))
dens1 <- surface.pollen.regions.df$AnnP[surface.pollen.regions.df$category == "NoAmbrosia"]/sum(surface.pollen.regions.df$AnnP[surface.pollen.regions.df$category == "NoAmbrosia"])
lines(surface.pollen.regions.df$AnnP[surface.pollen.regions.df$category == "NoAmbrosia"],dens1,col="red")

plot()
#####################
# Surface pollen models
# Use surface.pollen.corrected

summary(gam(formula=log(ambtoart) ~ splines::bs(AnnP,3), method="ML",
            data=surface.pollen.corrected[grepl("BothPresent", surface.pollen.corrected$category),])) 
both.gam <- gam(formula=log(ambtoart) ~ splines::bs(AnnP,3), method="GCV.Cp",
                data=surface.pollen.corrected[grepl("BothPresent", surface.pollen.corrected$category),])

reverse.both.gam <- gam(formula=AnnP ~ splines::bs(log(ambtoart),3), method="ML",
                data=surface.pollen.corrected[grepl("BothPresent", surface.pollen.corrected$category),])
summary(reverse.both.gam)

pres.both <- pres.pollen.regions2[grepl("BothPresent", pres.pollen.regions2$category),]
pres.both$testpredict <- predict(reverse.both.gam, newdata = pres.both, interval = "confidence")

surf.both <- surface.pollen.corrected[grepl("BothPresent", surface.pollen.corrected$category),]

pres.noamb <- pres.pollen.regions2[grepl("NoAmbrosia", pres.pollen.regions2$category),]

surf.noamb <- surface.pollen.corrected[grepl("NoAmbrosia", surface.pollen.corrected$category),]

noamb.gam <- gam(formula=log(ambtoart) ~ AnnP, method="ML",
                data=surface.pollen.corrected[grepl("NoAmbrosia", surface.pollen.corrected$category),])
summary(noamb.gam)
reverse.noamb.gam <- gam(formula=AnnP ~ log(ambtoart), method="ML",
                data=surface.pollen.corrected[grepl("NoAmbrosia", surface.pollen.corrected$category),])
pres.noamb$testpredict <- predict(reverse.noamb.gam, newdata = pres.noamb, interval = "confidence")
plot(pres.noamb$AnnP, pres.noamb$testpredict)


pres.noart <- pres.pollen.regions2[grepl("NoArtemisia", pres.pollen.regions2$category),]

surf.noart <- surface.pollen.corrected[grepl("NoArtemisia", surface.pollen.corrected$category),]

noamb.lm <- lm(formula=log(ambtoart) ~ AnnP, data=surface.pollen.corrected[grepl("NoAmbrosia", surface.pollen.corrected$category),])
summary(noamb.lm)

noart.lm <- gam(formula=log(ambtoart) ~ AnnP, data=surface.pollen.corrected[grepl("NoArtemisia", surface.pollen.corrected$category),])
summary(noart.lm)
reverse.noart <- gam(formula=AnnP ~ log(ambtoart), data=surface.pollen.corrected[grepl("NoArtemisia", surface.pollen.corrected$category),])
pres.noart$testpredict <- predict(reverse.noart, newdata = pres.noart, interval = "confidence")

reverse.noamb <- gam(formula=AnnP ~ log(ambtoart), data=surface.pollen.corrected[grepl("NoAmbrosia", surface.pollen.corrected$category),])
summary(reverse.noamb)

testpredict <- as.data.frame(predict(reverse.noamb, newdata = pres.noamb, interval = "confidence"))
plot(testpredict)
pres.noamb <- merge(pres.noamb,testpredict)
plot(pres.noamb$fit, pres.noamb$AnnP)
histogram(pres.noamb$AnnP-pres.noamb$fit)

ggplot() + geom_polygon(data=region2, 
                           aes(x=long, y=lat, group=group), colour="grey78", fill="white") + 
  theme(legend.position = "bottom") +
  geom_point(data=pres.both, aes(x=long, y=lat, fill=(AnnP-testpredict), shape=category), size=4, color="grey60") + 
  geom_point(data=fortify(pres.noamb), aes(x=long, y=lat, fill=(AnnP-testpredict), shape=category), size=4, color="grey60") +
#  geom_point(data=pres.noart, aes(x=long, y=lat, fill=(AnnP-testpredict), shape=category),size=4, color="grey60") +
  scale_fill_gradient2(low="#0000FF", mid="#FFFFCC", high="#FF0000", midpoint=0,
                        name="Modern - Reconstructed \n(mm/year)") + 
  scale_shape_manual(name="Category", values=c(21,22), labels=c("Both Present", expression(paste("No", italic(" Artemisia"))))) +
  coord_map(xlim = c(-120, -80),ylim = c(30, 56)) + theme_bw(base_size = 14) +
  theme(legend.position = "bottom") + xlab("Longitude") + ylab("Latitude") +
            theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14)) + 
  guides(shape = guide_legend(nrow = 3, title.hjust = 0.5, title.vjust=1, vjust=1, hjust=0.5), 
         fill = guide_colourbar(title.position="left", title.hjust = 0.5, title.vjust=1, vjust=0.5, hjust=0.5, barwidth = 6.5)) +   geom_segment(aes(x=-81, y=31.6, xend=-81, yend=32.3), size = 2, arrow = arrow(length = unit(0.25, "cm"))) +
  annotate("text", x=-81, y=31, label="N") 

ggsave("modernvsreconstructed.png", dpi=600, scale=1.5)

ggplot() + geom_histogram(data=pres.both, aes(x=AnnP-testpredict), binwidth=50, color="grey50", fill=NA) +
           geom_histogram(data=pres.noamb, aes(x=AnnP-testpredict), binwidth=50, color="royalblue4", fill="grey60")

# Plot of modern observed and presettlement reconstructed
ggplot(data=surf.both, aes(x=AnnP, y=log(ambtoart), group=category, shape="Surface", color="Surface")) + geom_point() + 
  geom_smooth(data=surf.both, 
             method = lm, formula = y ~ splines::bs(x, 3), se = F) + 
  xlab("Mean Annual Precipitation (mm)") + ylab(expression(paste("Log ", italic("Ambrosia "), "to ", italic("Artemisia "), "Ratio"))) +
  theme_classic(base_size = 16) + 
  geom_point(data=pres.both, 
       aes(x=testpredict, y=log(ambtoart), group=category, shape="Presettlement", color="Presettlement"))  +
  geom_smooth(data=pres.both, 
             method = lm, formula = y ~ splines::bs(x, 3), se = F, color="red") +
  ylim(-12,12) + xlim(200,1500) + scale_shape_manual("Dataset", values=c(2,1)) +
  theme(legend.position = c(0.5,0.95), legend.direction = "horizontal") + scale_color_manual(values=c("gray50", "black"), guide=FALSE)

ggplot(data=surf.noamb, aes(x=AnnP, y=log(ambtoart), group=category, shape="Surface", color="Surface")) + geom_point() + 
  geom_smooth(data=surf.noamb, 
             method = lm, formula = y ~ x, se = F) + 
  xlab("Mean Annual Precipitation (mm)") + ylab(expression(paste("Log ", italic("Ambrosia "), "to ", italic("Artemisia "), "Ratio"))) +
  theme_classic(base_size = 16) + 
  geom_point(data=pres.noamb, 
       aes(x=testpredict, y=log(ambtoart), group=category, shape="Presettlement", color="Presettlement"))  +
  geom_smooth(data=pres.noamb, 
             method = lm, formula = y ~ x, color="red") +
  ylim(-12,12) + xlim(200,1500) + scale_shape_manual("Dataset", values=c(2,1)) +
  theme(legend.position = c(0.5,0.95), legend.direction = "horizontal") + scale_color_manual(values=c("gray50", "black"), guide=FALSE)


# Plot of modern vs presettlement reconstructed precip
ggplot() + geom_point(data=pres.both, aes(x=AnnP, y=testpredict, fill=(AnnP-testpredict)), size=3, shape=21) +
           geom_point(data=pres.noamb, aes(x=AnnP, y=testpredict, fill=(AnnP-testpredict)), size=3, shape=21) +
  xlim(200, 1200) + ylim(200, 1200) + geom_abline(intercept=0, slope=1) +
  xlab("Modern Annual Precipitation (mm)") + ylab("Modeled Presettlement Annual Precipitation (mm)") + coord_equal() + theme_bw(base_size=11) +
  scale_fill_gradient2(low="#0000FF", mid="#FFFFCC", high="#FF0000", midpoint=0,
                        name="Modern - Reconstructed \n(mm/year)") + theme(legend.position = "bottom") + 
#  annotate("text", x = 900, y = 1050, label = "1:1 line", size=6) +
  guides(fill = guide_colourbar(title.position="left", title.hjust = 0.5, title.vjust=1,hjust=0.5, barwidth = 10))

ggsave("precipreconstruction.png", dpi=600, scale=1.2)

ggplot(data=pres.noamb, aes(x=AnnP, y=testpredict, color=(AnnP-testpredict))) + geom_point(size=2) + 
  xlim(200, 1200) + ylim(200, 1200) + geom_abline(intercept=0, slope=1) +
  xlab("Modern Annual Precipitation (mm)") + ylab("Modeled Presettlement Annual Precipitation (mm)") + coord_equal() + theme_bw(base_size=14) +
  scale_color_viridis_c(name="Modern - Reconstructed Precip (mm/year)") + theme(legend.position = "right") + annotate("text", x = 950, y = 1000, label = "1:1 line")


# Both together
ggplot(data=pres.both, aes(x=AnnP, y=testpredict, color=(AnnP-testpredict))) + geom_point(size=2, shape=19) + 
  geom_point(data=pres.noamb, aes(x=AnnP, y=testpredict, color=(AnnP-testpredict)), size=2, shape=17) +
  xlim(200, 1200) + ylim(200, 1200) + geom_abline(intercept=0, slope=1) +
  xlab("Modern Annual Precipitation (mm)") + ylab("Modeled Presettlement Annual Precipitation (mm)") + coord_equal() + theme_bw(base_size=14) +
  scale_color_viridis_c(name="Modern - Reconstructed Precip (mm/year)") + theme(legend.position = "right") + annotate("text", x = 950, y = 1000, label = "1:1 line")

summary(pres.both$AnnP)
summary(pres.noamb$AnnP)
histogram((pres.noamb$AnnP-pres.noamb$testpredict))
histogram(pres.both$calendar.ybp)
ggsave("precipreconstruction.png", dpi=300, scale=2)
histogram((pres.both$AnnP-pres.both$testpredict))
summary(pres.both$testpredict)
summary((pres.both$AnnP-pres.both$testpredict))
summary(pres.noamb$testpredict)
summary((pres.noamb$AnnP-pres.noamb$testpredict))

summary((pres.noart$AnnP))
summary(pres.noart$testpredict)
summary((pres.noart$AnnP-pres.noart$testpredict))

#xval <- approx(x = formula$fitted, y = formula$x, xout = 0)$y 
#estprecip <- approx(x = noart.lm$fitted.values, y = noart.lm$xlevels, xout = 20)$y

par(mfrow = c(2,2))
gam.check(both.gam) # Surprisingly good
ggsave("bothpresentgamcheck.jpeg", dpi=600)


summary(gam(formula=log(ambtoart) ~ AnnP, data=surface.pollen.corrected[grepl("NoAmbrosia", surface.pollen.corrected$category),]))

summary(lm(formula=log(ambtoart) ~ AnnP, data=surface.pollen.corrected[grepl("NoAmbrosia", surface.pollen.corrected$category),]))

noart.lm <- lm(formula=log(ambtoart) ~ AnnP, data=surface.pollen.corrected[grepl("NoArtemisia", surface.pollen.corrected$category),])
plot(noart.lm)
summary(lm(formula=log(ambtoart) ~ AnnP, data=surface.pollen.corrected[grepl("NoArtemisia", surface.pollen.corrected$category),]))
summary(lm(formula=log(ambtoart) ~ TSeasonality, data=surface.pollen.corrected[grepl("NoArtemisia", surface.pollen.corrected$category),]))

 ggplot(data=surface.pollen.corrected, aes(x=AnnP, y=log(ambtoart), group=category, color=category, shape=category)) + geom_point(size=2) + 
  geom_smooth(data=surface.pollen.corrected[grepl("BothPresent",surface.pollen.corrected$category),], 
             method = lm, formula = y ~ splines::bs(x, 3), se = F) +
  geom_smooth(data=surface.pollen.corrected[grepl("NoArtemisia", surface.pollen.corrected$category),], 
              method = lm, formula = y ~ x, se=F) +
  geom_smooth(data=surface.pollen.corrected[grepl("NoAmbrosia", surface.pollen.corrected$category),], 
            method = lm, formula = y ~ x, se = F) +
  xlab("Mean Annual Precipitation (mm)") + ylab(expression(paste("Log ", italic("Ambrosia "), "to ", italic("Artemisia "), "Ratio"))) +
  annotate(y=-10.0779449, x=(682.2679+250), "text", size=5, label=expression(paste("Zero ", italic("Ambrosia")))) + 
  annotate(y=0.7665241, x=(1181.9250+250), "text", size=5, label=expression(paste("Both Present"))) + 
  annotate(y=(10.2), x=(1450), "text", size=5, label=expression(paste("Zero ", italic("Artemisia")))) +
  theme_classic(base_size = 14) + scale_color_manual(values=c("#450d54", "#33638d", "#20a387")) + theme(legend.position = "none") +
   scale_shape_manual(values=c(19,15,17))

 ggsave("modelplot.png", dpi=600, scale = 1.1)

 ggplot() + geom_point(data=surface.pollen.regions.df[!grepl("NoArtemisia", surface.pollen.regions.df$category),], aes(x=AnnP, y=log(ambtoart), shape=category), color="blue") +
           geom_point(data=pres.pollen.regions2[!grepl("NoArtemisia", pres.pollen.regions2$category),], aes(x=AnnP, y=log(ambtoart), shape=category), color="orange")

# Remove three low Zero Art values
ggplot(data=surface.pollen.regions.df, aes(x=AnnP, y=log(ambtoart), group=category)) + geom_point() + 
  geom_smooth(data=surface.pollen.regions.df[grepl("BothPresent", surface.pollen.regions.df$category),], 
             method = lm, formula = y ~ splines::bs(x, 3), se = F) +
  geom_smooth(data=surface.pollen.regions.df[grepl("NoArtemisia", surface.pollen.regions.df$category),],
              method = lm, formula = y ~ x, se = F, linetype=2, color="gray40") +
  geom_smooth(data=surface.pollen.regions.df[grepl("NoArtemisia", surface.pollen.regions.df$category) &
                                            surface.pollen.regions.df$AnnP > 600,],
              method = lm, formula = y ~ x, se = F) +
  geom_smooth(data=surface.pollen.regions.df[grepl("NoArtemisia", surface.pollen.regions.df$category) &
                                            surface.pollen.regions.df$AnnP > 1000 | surface.pollen.regions.df$AnnP < 800,],
              method = lm, formula = y ~ x, se = F, linetype=4, color="red4") +
  geom_smooth(data=surface.pollen.regions.df[grepl("NoAmbrosia", surface.pollen.regions.df$category),], 
            method = lm, formula = y ~ x, se = F) +
  xlab("Mean Annual Precipitation (mm)") + ylab(expression(paste("Log ", italic("Ambrosia "), "to ", italic("Artemisia "), "Ratio"))) +
  annotate(y=-10.0779449, x=(682.2679+250), "text", size=5, label=expression(paste("Zero ", italic("Ambrosia")))) + 
  annotate(y=0.7665241, x=(1181.9250+250), "text", size=5, label=expression(paste("Both Present"))) + 
  annotate(y=(10.950505-2), x=(1006.1716), "text", size=5, label=expression(paste("Zero ", italic("Artemisia")))) +
  theme_classic(base_size = 16)


ggplot(data=surface.pollen.regions.df, aes(x=(AnnP), y=log(ambtoart), group=category)) + geom_point() + geom_smooth(method="lm") +
  xlab("Mean Annual Precipitation (mm)") + ylab(expression(paste("Log ", italic("Ambrosia "), "to ", italic("Artemisia "), "Ratio"))) +
  theme_classic(base_size = 20) + facet_wrap(vars(droplevels(NA_L2NAME)), nrow=2)

ggsave("facetplot.jpeg",dpi=600,scale=3)
ggsave("surface_lmplot.truncated.jpeg",dpi=600)  



#Presettlement
# Create column to facet by

pres.pollen.regions2$category <- NULL
pres.pollen.regions2$category[pres.pollen.regions2$AMBROSIA == 0] <- "NoAmbrosia"
pres.pollen.regions2$category[pres.pollen.regions2$ARTEMISIA == 0] <- "NoArtemisia"
pres.pollen.regions2$category[pres.pollen.regions2$ARTEMISIA != 0 & pres.pollen.regions2$AMBROSIA != 0] <- "BothPresent"


ggplot(data=pres.pollen.regions2[!grepl("NoArtemisia", pres.pollen.regions2$category),], 
       aes(x=AnnP, y=log(ambtoart), group=category)) + geom_point() + 
  geom_smooth(data=pres.pollen.regions2[grepl("BothPresent", pres.pollen.regions2$category),], 
             method = lm, formula = y ~ splines::bs(x, 3), se = F) +
  geom_smooth(data=pres.pollen.regions2[grepl("NoArtemisia", pres.pollen.regions2$category),], 
              method = lm, formula = y ~ x, se = F) +
  geom_smooth(data=pres.pollen.regions2[grepl("NoAmbrosia", pres.pollen.regions2$category),], 
            method = lm, formula = y ~ x, se = F) +
  xlab("Mean Annual Precipitation (mm)") + ylab(expression(paste("Log ", italic("Ambrosia "), "to ", italic("Artemisia "), "Ratio"))) +
  annotate(y=-20.0779449, x=(682.2679+250), "text", size=5, label=expression(paste("Zero ", italic("Ambrosia")))) + 
  annotate(y=0.7665241, x=(1181.9250+250), "text", size=5, label=expression(paste("Both Present"))) + 
  annotate(y=(16.950505-2), x=(1006.1716), "text", size=5, label=expression(paste("Zero ", italic("Artemisia")))) +
  theme_classic(base_size = 16)

# presettlement and surface together

ggplot(data=surface.pollen.regions.df, aes(x=AnnP, y=log(ambtoart), group=category, shape="Surface", color="Surface")) + geom_point() + 
  geom_smooth(data=surface.pollen.regions.df[grepl("BothPresent", surface.pollen.regions.df$category),], 
             method = lm, formula = y ~ splines::bs(x, 3), se = T) +
  geom_smooth(data=surface.pollen.regions.df[grepl("NoArtemisia", surface.pollen.regions.df$category),], 
              method = lm, formula = y ~ x, se = T) +
  geom_smooth(data=surface.pollen.regions.df[grepl("NoAmbrosia", surface.pollen.regions.df$category),], 
            method = lm, formula = y ~ x, se = T) +
  xlab("Mean Annual Precipitation (mm)") + ylab(expression(paste("Log ", italic("Ambrosia "), "to ", italic("Artemisia "), "Ratio"))) +
  annotate(y=-10.0779449, x=(682.2679+250), "text", size=5, label=expression(paste("Zero ", italic("Ambrosia")))) + 
  annotate(y=0.7665241, x=(1181.9250+100), "text", size=5, label=expression(paste("Both Present"))) + 
  #annotate(y=(16.950505-2), x=(1006.1716), "text", size=5, label=expression(paste("Zero ", italic("Artemisia")))) +
  theme_classic(base_size = 16) + 
  geom_point(data=pres.pollen.regions2, 
       aes(x=AnnP, y=log(ambtoart), group=category, shape="Presettlement", color="Presettlement"))  +
  geom_smooth(data=pres.pollen.regions2[grepl("BothPresent", pres.pollen.regions2$category),], 
             method = lm, formula = y ~ splines::bs(x, 3), se = T, color="red") +
  geom_smooth(data=pres.pollen.regions2[grepl("NoArtemisia", pres.pollen.regions2$category),], 
              method = lm, formula = y ~ x, se = T, color="red") +
  geom_smooth(data=pres.pollen.regions2[grepl("NoAmbrosia", pres.pollen.regions2$category),], 
            method = lm, formula = y ~ x, se = T, color="red") + ylim(-12,12) + xlim(200,1500) + scale_shape_manual("Dataset", values=c(2,1)) +
  theme(legend.position = c(0.5,0.95), legend.direction = "horizontal") + scale_color_manual(values=c("gray50", "black"), guide=FALSE)

min(surface.pollen.regions.df$AnnP)
ggsave("ambtoartplot.png", dpi=300)
t.noart <- pres.pollen.regions2$ambtoart[pres.pollen.regions2$ARTEMISIA == 0]
ggplot(data=pres.pollen.regions2[pres.pollen.regions2$ARTEMISIA == 0,], aes(x=AnnP, y=ambtoart)) + geom_point()

aggregate(log(pres.pollen.regions2$ambtoart), list(pres.pollen.regions2$category), range)
aggregate(log(pres.pollen.regions2$ambtoart), list(pres.pollen.regions2$category), mean)
aggregate(log(pres.pollen.regions2$ambtoart), list(pres.pollen.regions2$category), summary)
aggregate((pres.pollen.regions2$AnnP), list(pres.pollen.regions2$category), range)
aggregate((pres.pollen.regions2$AnnP), list(pres.pollen.regions2$category), mean)
pres.pollen.regions2 %>% 
  group_by(category) %>%
  summarise(no_rows = length(category))

histogram(log(pres.pollen.regions2$ambtoart[pres.pollen.regions2$category == "BothPresent"]))

#####################
# pres pollen models

newdata <- pres.pollen.regions2[grepl("BothPresent", pres.pollen.regions2$category),]
newdata$ambtoart <- log(newdata$ambtoart)
newdata$predprecip <- predict.gam(reverse.fit, newdata)
histogram(newdata$predprecip)
histogram(newdata$AnnP)
ggplot(newdata, aes(x=log(AnnP), y=log(predprecip))) + geom_point()

summary(gam(formula=log(ambtoart) ~ splines::bs(AnnP,3), data=pres.pollen.regions2[grepl("BothPresent", pres.pollen.regions2$category),])) 
both.gam <- gam(formula=log(ambtoart) ~ splines::bs(AnnP,3), 
                data=pres.pollen.regions.df[grepl("BothPresent", pres.pollen.regions.df$category),])
both.gam

summary(gam(formula=log(ambtoart) ~ AnnP, data=pres.pollen.regions.df[grepl("NoAmbrosia", pres.pollen.regions.df$category),]))
summary(gam(formula=log(ambtoart) ~ splines::bs(AnnP,3), data=pres.pollen.regions.df[grepl("NoArtemisia", pres.pollen.regions.df$category),]))

summary(lm(formula=log(ambtoart) ~ AnnP, data=pres.pollen.regions.df[grepl("NoAmbrosia", pres.pollen.regions.df$category),]))

summary(lm(formula=log(ambtoart) ~ AnnP, data=pres.pollen.regions.df[grepl("NoArtemisia", pres.pollen.regions.df$category),]))
summary(lm(formula=log(ambtoart) ~ TSeasonality, data=pres.pollen.regions.df[grepl("NoArtemisia", pres.pollen.regions.df$category),]))


ggplot(data=pres.pollen.regions2, aes(x=(AnnP), y=log10(ambtoart), group=NA_L2NAME, color=NA_L2NAME)) + geom_point() + geom_smooth() +
  theme(legend.position = "none")

range(log(pres.pollen.regions2$ambtoart))
plot((surface.pollen.regions.df$AnnP), log(surface.pollen.regions.df$ambtoart))
histogram(log(pres.pollen.clean2$ambtoart))


buffermap + geom_point(data=pres.pollen.clean2, aes(x=long,y=lat)) +
                ggtitle("Presettlement Datasets Within 2 Degree Buffer of Great Plains")

ggsave("presettlementmap.png", dpi=600, scale=2)


###########################3
obs.surf <- count(surface.set3, droplevels(NA_L2NAME))
obs.pres <- count(presettlement.set, droplevels(NA_L2NAME))

#labels = c("5.2  Mixed Wood Shield",
 #                              "5.4  Boreal Plain",
  #                             "6.2  Western Cordillera",
   #                            "8.1  Mixed Wood Plain",
    #                           "8.2  Central USA Plains",
     #                          "8.3  Southeastern USA Plains",
      #                         "8.4  Ozark/Ouachita-Appalachian Forests",
       #                        "8.5  Mississippi Alluvial and Southeast USA Coastal Plains",
        #                       "9.2  Temperate Prairies",
         #                      "9.3  West-Central Semiarid Prairies",
          #                     "9.4  South Central Semiarid Prairies",
           #                    "9.6  Tamaulipas-Texa Semiarid Plain",
            #                   "10.1  Cold Deserts",
             #                  "10.2  Warm Deserts",
              #                 "13.1  Upper Gila Mountains")

scale_fill_manual(values = c("MIXED WOOD SHIELD" = "#aedee4",
                               "BOREAL PLAIN" = "#bed7c0",
                               "WESTERN CORDILLERA" = "#5ebc55",
                               "MIXED WOOD PLAINS" = "#a4cb4b",
                               "CENTRAL USA PLAINS" = "#e1f1e7",
                               "TEMPERATE PRAIRIES" = "#fee5ca",
                               "WEST-CENTRAL SEMIARID PRAIRIES" = "#fadca1",
                               "SOUTH CENTRAL SEMIARID PRAIRIES" = "#f8ce9e"),
                    breaks = c("MIXED WOOD SHIELD",
                               "BOREAL PLAIN",
                               "WESTERN CORDILLERA",
                               "MIXED WOOD PLAINS",
                               "CENTRAL USA PLAINS",
                               "TEMPERATE PRAIRIES",
                               "WEST-CENTRAL SEMIARID PRAIRIES",
                               "SOUTH CENTRAL SEMIARID PRAIRIES"),
                    labels = c("Mixed Wood Shield",
                               "Boreal Plain",
                               "Western Cordillera",
                               "Mixed Wood Plain",
                               "Central USA Plains",
                               "Temperate Prairies",
                               "West-Central Semiarid Prairies",
                               "South Central Semiarid Prairies"), 
                    name=NULL)

ggplot() + geom_boxplot(data=surface.set3, aes(x=NA_L2NAME, y= log(ambtoart), fill=droplevels(NA_L2NAME)), 
                        color="black") + 
  scale_fill_manual(values = c("MIXED WOOD SHIELD" = "#aedee4",
                               "BOREAL PLAIN" = "#bed7c0",
                               "WESTERN CORDILLERA" = "#5ebc55",
                               "MIXED WOOD PLAINS" = "#a4cb4b",
                               "CENTRAL USA PLAINS" = "#e1f1e7",
                               "TEMPERATE PRAIRIES" = "#fee5ca",
                               "WEST-CENTRAL SEMIARID PRAIRIES" = "#fadca1",
                               "SOUTH CENTRAL SEMIARID PRAIRIES" = "#f8ce9e"),
                    breaks = c("MIXED WOOD SHIELD",
                               "BOREAL PLAIN",
                               "WESTERN CORDILLERA",
                               "MIXED WOOD PLAINS",
                               "CENTRAL USA PLAINS",
                               "TEMPERATE PRAIRIES",
                               "WEST-CENTRAL SEMIARID PRAIRIES",
                               "SOUTH CENTRAL SEMIARID PRAIRIES"),
                    labels = c("Mixed Wood Shield",
                               "Boreal Plain",
                               "Western Cordillera",
                               "Mixed Wood Plain",
                               "Central USA Plains",
                               "Temperate Prairies",
                               "West-Central Semiarid Prairies",
                               "South Central Semiarid Prairies")) + 
  scale_x_discrete(breaks = c("MIXED WOOD SHIELD",
                               "BOREAL PLAIN",
                               "WESTERN CORDILLERA",
                               "MIXED WOOD PLAINS",
                               "CENTRAL USA PLAINS",
                               "TEMPERATE PRAIRIES",
                               "WEST-CENTRAL SEMIARID PRAIRIES",
                               "SOUTH CENTRAL SEMIARID PRAIRIES"),
                    labels = c("Mixed wood shield",
                               "Boreal plain",
                               "Western cordillera",
                               "Mixed wood plain",
                               "Central USA plains",
                               "Temperate prairies",
                               "West-central semiarid prairies",
                               "South central semiarid prairies"), name=NULL) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, family="sans", size=14)) +  
        labs(y= expression(paste(italic("Ambrosia")," to ", italic("Artemisia")," ratio")),
             x= NULL) + 
  guides(fill=FALSE) +
  theme(text=element_text(size=14))

ggplot() + geom_boxplot(data=s.pollen.regions2, aes(x=reorder(NA_L2NAME, log(ambtoart), FUN = median), y=log(ambtoart)), fill="gray90", 
                        color="black") + theme_classic() +
  scale_x_discrete(breaks = c("MIXED WOOD SHIELD",
                               "BOREAL PLAIN",
                               "MIXED WOOD PLAINS",
                               "CENTRAL USA PLAINS",
                               "TEMPERATE PRAIRIES",
                               "WEST-CENTRAL SEMIARID PRAIRIES",
                               "SOUTH CENTRAL SEMIARID PRAIRIES",
                               "MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS",
                               "OZARK/OUACHITA-APPALACHIAN FORESTS",
                               "SOUTHEASTERN USA PLAINS"),
                    labels = c("Mixed wood shield",
                               "Boreal plain",
                               "Mixed wood plains",
                               "Central USA plains",
                               "Temperate prairies",
                               "West-central semiarid prairies",
                               "South central semiarid prairies",
                               "MS alluvial and SE USA coastal plains",
                               "Ozark/Ouachita-Appalachian forests",
                               "Southeastern USA plains"), name=NULL) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, family="sans", size=12)) +  
        labs(y= expression(paste("Log", italic(" Ambrosia")," to ", italic("Artemisia")," ratio")),
             x= NULL) + 
  guides(fill=FALSE) +
  theme(text=element_text(size=12,  family="sans"))

ggsave("regionboxplot.png", dpi=600, scale=1.5)

surface.set3$NA_L2NAME <- factor(droplevels(surface.set3$NA_L2NAME), 
                     levels=c("MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS",
                              "OZARK/OUACHITA-APPALACHIAN FORESTS",
                              "SOUTHEASTERN USA PLAINS",
                              "MIXED WOOD PLAINS",
                              "CENTRAL USA PLAINS",
                              "SOUTH CENTRAL SEMIARID PRAIRIES",
                              "MIXED WOOD SHIELD",
                              "TEMPERATE PRAIRIES",
                              "WEST-CENTRAL SEMIARID PRAIRIES",
                              "WESTERN CORDILLERA",
                              "BOREAL PLAIN"
                              ), ordered = TRUE)
ggplot() + geom_boxplot(data=surface.set3, aes(x=NA_L2NAME, y=ambtoart), fill="gray90", 
                        color="black") + theme_classic() +
#  scale_fill_manual(values = p.colors) + 
  scale_x_discrete(breaks = c("MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS",
                              "OZARK/OUACHITA-APPALACHIAN FORESTS",
                              "SOUTHEASTERN USA PLAINS",
                              "MIXED WOOD PLAINS",
                              "CENTRAL USA PLAINS",
                              "SOUTH CENTRAL SEMIARID PRAIRIES",
                              "MIXED WOOD SHIELD",
                              "TEMPERATE PRAIRIES",
                              "WEST-CENTRAL SEMIARID PRAIRIES",
                              "WESTERN CORDILLERA",
                              "BOREAL PLAIN"
                              ),
                    labels = c("MS alluvial and SE USA coastal plains",
                               "Ozark/Ouachita-Appalachian forests",
                               "Southeastern USA plains",
                               "Mixed wood plains",
                               "Central USA plains",
                               "South central semiarid prairies",
                               "Mixed wood shield",
                               "Temperate prairies",
                               "West-central semiarid prairies",
                               "Western cordillera",
                               "Boreal plain" 
                               ), limits= rev(levels(surface.set3$NA_L2NAME)), name=NULL) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, family="sans", size=12)) +  
        labs(y= expression(paste("Log", italic(" Ambrosia")," to ", italic("Artemisia")," ratio")),
             x= NULL) + 
  guides(fill=FALSE) +
  theme(text=element_text(size=12,  family="sans"))

ggsave("surface_boxplot.png", dpi=400, scale=1.25)

pre.colors <- p.colors[c(1,3,4,6:8)]

ggplot() + geom_boxplot(data=presettlement.set, aes(x=NA_L2NAME, y=ambtoart, fill=droplevels(NA_L2NAME)), 
                        color="black")+ 
  scale_fill_manual(values = pre.colors) + 
  scale_x_discrete(breaks = c("MIXED WOOD SHIELD",
                               "BOREAL PLAIN",
                               "WESTERN CORDILLERA",
                               "MIXED WOOD PLAINS",
                               "TEMPERATE PRAIRIES",
                               "WEST-CENTRAL SEMIARID PRAIRIES"),
                    labels = c("Mixed wood shield",
                               "Boreal plain",
                               "Western cordillera",
                               "Mixed wood plain",
                               "Temperate prairies",
                               "West-central semiarid prairies"), name=NULL) + 
 theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, family="sans", size=14)) + 
        labs(y= expression(paste(italic("Ambrosia")," to ", italic("Artemisia")," ratio")),
             x= NULL) + 
  guides(fill=FALSE) + 
  theme(text=element_text(size=14,  family="sans"))
ggsave("presettlement_boxplot.jpeg", dpi=600, scale=1.25)
         
ggplot() + geom_boxplot(data=presettlement.set, aes(x=NA_L2NAME, y=ambtoart, fill=droplevels(NA_L2NAME)), 
                        color="black")+ 
  scale_fill_manual(values = c("MIXED WOOD SHIELD" = "#aedee4",
                               "BOREAL PLAIN" = "#bed7c0",
                               "WESTERN CORDILLERA" = "#5ebc55",
                               "MIXED WOOD PLAINS" = "#a4cb4b",
                               "TEMPERATE PRAIRIES" = "#fee5ca",
                               "WEST-CENTRAL SEMIARID PRAIRIES" = "#fadca1")) + 
  scale_x_discrete(breaks = c("MIXED WOOD SHIELD",
                               "BOREAL PLAIN",
                               "WESTERN CORDILLERA",
                               "MIXED WOOD PLAINS",
                               "TEMPERATE PRAIRIES",
                               "WEST-CENTRAL SEMIARID PRAIRIES",
                               "MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS",
                               "CENTRAL USA PLAINS"),
                    labels = c("Mixed wood shield",
                               "Boreal plain",
                               "Western cordillera",
                               "Mixed wood plain",
                               "Temperate prairies",
                               "West-central semiarid prairies"), name=NULL) + 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, family="sans", size=14)) + 
        labs(y= expression(paste(italic("Ambrosia")," to ", italic("Artemisia")," ratio")),
             x= NULL) + 
  guides(fill=FALSE) + 
  theme(text=element_text(size=14,  family="sans"))
ggsave("presettlement_boxplot_mapmatch.jpeg", dpi=600, scale=1.25)

# With climate data
aov.pres2 <- aov(ambtoart ~ droplevels(NA_L2NAME) + AnnP*MeanT, data=presettlement.set)
summary(aov.pres2)

aov.surf2 <- aov(ambtoart ~ droplevels(NA_L2NAME) + AnnP*PrecipWarmQ, data=surface.set3)
summary(aov.surf2)

######### Analogue

matches <- cbind(s.analogs.ambart, num.analogs = rowSums(s.analogs.ambart <= 0.205))

surf.no.ambart <- s.pollen.regions2[,c(2,6:95,101)]
rownames(surf.no.ambart) <- s.pollen.regions2[,8]
rownames(surf.no.ambart2) <- make.names(s.pollen.regions2[,101], unique=TRUE, allow_=TRUE)
s.noambart.props <- cbind(surf.no.ambart[c(1:4,92)], prop.table(as.matrix(surf.no.ambart[-c(1:4,92)]), margin = 1))
# Split surface samples into training set and test set
#set.seed(123)
training.surf.no.ambart <- droplevels(s.noambart.props$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
train.surf.no.ambart <- s.noambart.props[training.surf.no.ambart, ]
test.surf.no.ambart <- s.noambart.props[-training.surf.no.ambart, ]
noambart.analog <- analog(x=train.surf.no.ambart[,c(6:92)],  y=test.surf.no.ambart[,c(6:92)],
                      method="SQchord", keep.train = TRUE)
summary(noambart.analog)
cma(noambart.analog)


#########
surf.ambart <- s.pollen.regions2[,c(2, 6:8, 9:95,99:101)]
surf.ambart2 <- s.pollen.regions2[,c(97,98,101)]

rownames(surf.ambart2) <- make.names(s.pollen.regions2[,101], unique=TRUE, allow_=TRUE)
s.ambart.props <- cbind(surf.ambart[c(1:4,93)], prop.table(as.matrix(surf.ambart[-c(1:4,93)]), margin = 1))
# Split surface samples into training set and test set
#set.seed(123)
training.surf.ambart <- droplevels(s.ambart.props$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
train.surf.ambart <- s.ambart.props[training.surf.ambart, ]
test.surf.ambart <- s.ambart.props[-training.surf.ambart, ]
ambart.analog <- analog(x=train.surf.ambart[,c(6:93)],  y=test.surf.ambart[,c(6:93)],
                      method="SQchord", keep.train = TRUE)

training.surf.ambart2 <- droplevels(surf.ambart2$NA_L2NAME) %>%
  createDataPartition(p = 0.8, list = FALSE)
train.surf.ambart2 <- surf.ambart2[training.surf.ambart2, ]
test.surf.ambart2 <- surf.ambart2[-training.surf.ambart2, ]
ambart.analog2 <- analog(x=train.surf.ambart2[,c(1:2)],  y=test.surf.ambart2[,c(1:2)],
                      method="SQchord", keep.train = TRUE)


summary(ambart.analog2)
cma(ambart.analog)


ambart.analog.df <- as.data.frame(ambart.analog[["analogs"]])
ambart.analog.df2 <- which(as.numeric(rownames(s.pollen.regions2)) %in% ambart.analog.df)
testmds <- metaMDS(sqrt(s.ambart.props[,c(6:93)]), distance = "euclidean", plot=TRUE)
plot(testmds, type="text", display="sites")



######

# Define color scheme
geom_polygon(data = prairie, aes(x=long, y = lat, 
                group=group), fill="yellowgreen", colour = NA, alpha = 0.3) + geom_polygon(data = ozark, 
                aes(x=long, y = lat, group=group), fill="springgreen4", colour = NA, 
                alpha = 0.3) + geom_polygon(data = seplains, aes(x=long, y = lat, group=group),
                fill="royalblue4", colour = NA, alpha = 0.3) + geom_polygon(data = wcentral, 
                aes(x=long, y = lat, group=group), fill="bisque4", colour = NA, 
                alpha = 0.3) + geom_polygon(data = scentral, aes(x=long, y = lat, group=group),
                fill="peru", colour = NA, alpha = 0.3) + geom_polygon(data = mixed, aes(x=long, 
                y = lat, group=group), fill="azure4", colour = NA, alpha = 0.3) + geom_polygon(data = mixedplains, 
                aes(x=long, y = lat, group=group), fill="lightcoral", colour = NA, 
                alpha = 0.3) + geom_polygon(data = central, aes(x=long, y = lat, group=group), 
                fill="chocolate4", colour = NA, alpha = 0.3) + coord_map() + scale_fill_identity(name="Level II Ecoregions", 
                                  guide="legend", labels=c("Prairie"))



s.dissim2 <- dissim(ambart.analog2)
plot(s.dissim2)

clust <- hclust(as.dist(ambart.analog$train))
plot(clust, labels=NULL)
grps <- cutree(clust, k=5)
rect.hclust(clust, k = 5, border = 2:5)
g24 <- cutree(clust, k = c(4,6))
table(grp4 = g24[,"4"], grp6 = g24[,"6"])
s.ambart.roc <- roc(ambart.analog, groups=grps)

clust2 <- hclust(as.dist(ambart.analog2$train))
plot(clust2, labels=NULL)
grps2 <- cutree(clust2, k=5)
rect.hclust(clust2, k = 5, border = 2:5)
g24 <- cutree(clust, k = c(4,6))
table(grp4 = g24[,"4"], grp6 = g24[,"6"])
s.ambart.roc <- roc(ambart.analog, groups=grps)

color_by_region <- as.numeric(ambart.analog[,5])
colors_to_use
colors_to_use <- colors_to_use[order.dendrogram(clust)]
colors_to_use


#ROC curve of dissimilarities

#Discrimination for all groups:

#Optimal Dissimilarity = 0.205 

#AUC = 0.99, p-value: < 2.22e-16
#No. within: 348   No. outside: 3828 
s.ambart.bayes <- bayesF(s.ambart.roc)
plot(s.ambart.bayes) # 6 groups

# Noambart

naa.dissim <- dissim(noambart.analog)
plot(naa.dissim)
naa.clust <- hclust(as.dist(noambart.analog$train))
naa.grps <- cutree(naa.clust, k=5)
plot(naa.clust)
rect.hclust(naa.clust, k = 6, border = 1:6)
naa.s.ambart.roc <- roc(noambart.analog, groups=naa.grps)
naa.s.ambart.bayes <- bayesF(naa.s.ambart.roc)
plot(naa.s.ambart.bayes) # 6 groups

write.csv(ambart.analog[["analogs"]], file="sambartanalogs.csv", row.names = TRUE)

###############################################################################
# Datasets included in analysis

# Presettlement
presettlement.set.citations <- pres.pollen.regions2[!grepl("OZARK/OUACHITA-APPALACHIAN FORESTS", pres.pollen.regions2$NA_L2NAME),]
presettlement.set.citations <- presettlement.set.citations[!grepl("SOUTHEASTERN USA PLAINS", presettlement.set.citations$NA_L2NAME),]
surrounding.fossil <- get_dataset(loc = c(-120.352968, 32.603899, -87.130311, 57.457753), 
                      datasettype = 'pollen',
                      ageold = 700)

pres.list <- as.numeric(unique(pres.pollen.regions2$dataset))
pub <- NULL
for (i in 1:length(pres.list)){
pub[[i]] <- get_publication(datasetid = pres.list[i])
}

pub2 <- do.call("rbind", lapply(pub, '[[', 1))

pub3 <- rbindlist(pub2, fill=TRUE)

pub3 <- pub3 %>% unique()

write.csv(pub3, file="presettlementcitations.csv", row.names = FALSE)

surf.list <- as.numeric(unique(surface.pollen.corrected$dataset))
surfpub <- NULL
for (i in 1:length(surf.list)){
surfpub[[i]] <- get_publication(datasetid = surf.list[i])
}
spub2 <- do.call("rbind", lapply(surfpub, '[[', 1))
spub3 <- rbindlist(spub2, fill=TRUE)
spub3 <- spub3 %>% unique()

write.csv(spub3, file="surfacecitations.csv", row.names = FALSE)

#############################################################################
# Histograms of raw data

ggplot() + geom_density(data=pres.pollen.regions2, aes(ambrosia.proportion, fill="Presettlement"), alpha=0.6) +
           geom_density(data=surface.pollen.corrected, aes(ambrosia.proportion, fill="Surface"), alpha=0.6) +
           scale_fill_manual(values=c("lightblue2", "gray60"))

amb.boxplot <- ggplot() + geom_boxplot(data=pres.pollen.regions2, aes(x="Presettlemet", y=ambrosia.proportion, fill="Presettlement")) +
           geom_boxplot(data=surface.pollen.corrected, aes(x="Surface", y=ambrosia.proportion, fill="Surface")) +
           scale_fill_manual(values=c("white", "gray60")) + xlab("") +
           ylab("Proportion") + theme_classic() + theme(legend.position = "none") +
           ggtitle(expression(paste(italic("Ambrosia")," Proportion"))) +
           theme(axis.text=element_text(size=12, face="bold", color = "black"),
                 axis.title=element_text(size=12,face="bold", color = "black"), 
                 plot.title = element_text(hjust = 0.5, size=14, face="bold")) + ylim(0,0.6)

art.boxplot <- ggplot() + geom_boxplot(data=pres.pollen.regions2, aes(x="Presettlemet", y=artemisia.proportion, fill="Presettlement")) +
           geom_boxplot(data=surface.pollen.corrected, aes(x="Surface", y=artemisia.proportion, fill="Surface")) +
           scale_fill_manual(values=c("white", "gray60")) + xlab("") +
           ylab(NULL) + theme_classic() +theme(legend.position = "none") +  
            ggtitle(expression(paste(italic("Artemisia")," Proportion"))) +
  theme(axis.text=element_text(size=12, face="bold", color = "black"),
                 axis.title=element_text(size=12, face="bold", color = "black"), 
                 plot.title = element_text(hjust = 0.5, size=14, face="bold")) + ylim(0,0.6)
cowplot::plot_grid(amb.boxplot, art.boxplot)
ggsave("ambartboxplot.png", dpi=600, scale = 1.1)


ggplot() + geom_boxplot(data=surface.pollen.corrected, aes(x=category, y=log(ambtoart)), fill="gray90") +
           xlab("") + theme_classic() +
           theme(legend.position = "none") +
           theme(axis.text=element_text(size=12, colour = "black"),
                 axis.title=element_text(size=12, colour = "black")) +
       ylab(expression(paste("Log ", italic("Ambrosia "), "to ", italic("Artemisia "), "Ratio"))) + 
  scale_x_discrete(labels=c("BothPresent" = "Both Present", "NoAmbrosia" = expression(paste("No ", italic("Ambrosia"))),
                              "NoArtemisia" = expression(paste("No ", italic("Artemisia")))))
ggsave("surfacescenarios.png", dpi=600)
