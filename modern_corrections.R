pre.ratios <- pres.pollen.regions2[,c(11,12,115:117,119,137)]

pre.ratios$type <- "Presettlement"

surf.ratios <- surface.pollen.regions.df[ , which(names(surface.pollen.regions.df) %in% c("lat", "long", "ambrosia.proportion", "artemisia.proportion", "ambtoart", "NA_L2NAME", "AnnP"))]

surf.ratios$type <- "Surface"

pre.ratios.spdf <- SpatialPointsDataFrame(coords=pre.ratios[,c(2,1)], data=pre.ratios[,-c(2,1)], 
                                         proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
pre.ratios.spdf <- spTransform(pre.ratios.spdf, CRS("+proj=robin +datum=WGS84")) 

surf.ratios.spdf <- SpatialPointsDataFrame(coords=surf.ratios[,c(2,1)], data=surf.ratios[,-c(2,1)], 
                                         proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
surf.ratios.spdf <- spTransform(surf.ratios.spdf, CRS("+proj=robin +datum=WGS84")) 

nearest.neighbor <- gDistance(surf.ratios.spdf,pre.ratios.spdf,byid = TRUE)

min.d <- apply(nearest.neighbor, 1, function(x) order(x, decreasing=F)[1])

newdata <- cbind(pre.ratios, surf.ratios[min.d,], by="type")
colnames(newdata) <- c("p.lat","p.long","p.amb.prop","p.art.prop","p.ambart","p.NA_L2NAME","p.AnnP","p.type",
                       "s.lat","s.long","s.amb.prop","s.art.prop","s.ambart","s.NA_L2NAME","s.AnnP","s.type","by")

ggplot(newdata) + geom_point(aes(x=s.amb.prop/2, y=p.amb.prop), shape=4) + coord_equal(xlim = c(0,0.25), ylim=c(0,0.25))

ggplot(newdata) + geom_point(aes(x=s.art.prop, y=p.art.prop), shape=5) + coord_equal(xlim = c(0,0.6), ylim=c(0,0.6))
ggplot() + geom_histogram(data=pres.pollen.clean2, aes(x=log(ambrosia.proportion), fill="Presettlement"), binwidth = 0.01) +
           geom_histogram(data=surface.pollen.clean2, aes(x=log(ambrosia.proportion/2.959758), fill="Surface"), alpha=0.4, binwidth = 0.01) +
  xlim(-9,0)

histogram(pre.ratios$ambrosia.proportion)
histogram(surf.ratios$ambrosia.proportion)

# Ambrosia

# Max:
0.570261/0.23110

# Mean:
0.097672/0.03300

# Median:
0.048433/0.02403

# 1st Q:
0.003003/0.00906

summary(pres.pollen.regions2$ambrosia.proportion)
summary(surface.pollen.regions.df$ambrosia.proportion)


# Ambrosia has more than doubled since settlement - making modern climate appear artifically wet 
# without a correction (or, making past cliamte appear artificially dry without a correction)
median(s.pollen.regions2$ambrosia.proportion)/median(pres.pollen.regions2$ambrosia.proportion)
mean(s.pollen.regions2$ambrosia.proportion)/mean(pres.pollen.regions2$ambrosia.proportion)

# Artemisia
summary(pres.pollen.regions2$artemisia.proportion)
summary(surface.pollen.regions.df$artemisia.proportion)


ggplot() + #geom_histogram(data=pres.pollen.clean2, aes(x=log(artemisia.proportion), fill="Presettlement"), binwidth = 0.01) +
           geom_histogram(data=surface.pollen.regions.df, aes(x=log(artemisia.proportion), fill="Surface"), alpha=0.4, binwidth = 0.01)


