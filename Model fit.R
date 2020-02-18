# Randomly remove 20% of the surface data
# Split surface samples into training set and test set
library(boot)
library(rioja)

twocol <- surface.pollen.corrected[grepl("BothPresent", surface.pollen.corrected$category),] %>%
  select("AnnP", "ambtoart")
twocol$ambtoart <- log(twocol$ambtoart)

both.gam <- gam(formula=log(ambtoart) ~ splines::bs(AnnP,3), method="ML",
                data=surface.pollen.corrected[grepl("BothPresent", surface.pollen.corrected$category),])

rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
}

results <- boot(data=twocol, statistic=rsq, 
   R=1000, formula=ambtoart ~ splines::bs(AnnP,3))
results
plot(results)
boot.ci(results, type="bca")

png(filename="bootstrap.png", 
    type="cairo",
    units="in", 
    width=6, 
    height=5, 
    pointsize=12, 
    res=600)
plot(results)
dev.off()

png(filename="gamcheck.png", 
    type="cairo",
    units="in", 
    width=6.5, 
    height=7, 
    pointsize=12, 
    res=600)
par(mfrow = c(2,2))
gam.check(both.gam)
dev.off()


Ambart <- surface.pollen.regions.df[grepl("BothPresent", surface.pollen.regions.df$category), c(97:99)] 
Ambart$ambtoart <- log(Ambart$ambtoart)
Ambart$ambrosia.proportion <- Ambart$ambrosia.proportion*100
Ambart$artemisia.proportion <- Ambart$artemisia.proportion*100
Ambart <- Ambart[-166,]
AnnP <- surface.pollen.regions.df[grepl("BothPresent", surface.pollen.regions.df$category), c(101,119)] 
AnnP <- AnnP[-166,]
mrfit <- MR(y=Ambart, x=AnnP$AnnP)
summary(mrfit)

crossvaltest <- rioja::crossval(mrfit)
plot(AnnP$AnnP, crossvaltest$predicted)
performance(mrfit)

##############
library(caret)


# define training control
train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(ambtoart~., data=twocol, trControl=train_control, method="nb")
# summarize results
print(model)

# By ecoregion
print(both.gam)
both.gam$gcv.ubre

library(gamclass)
crossvalboth <- gamclass::CVgam(surface.pollen.regions.df[grepl("BothPresent", surface.pollen.regions.df$category),], both.gam, 
                          nfold = nrow(surface.pollen.regions.df[grepl("BothPresent", surface.pollen.regions.df$category),]))
gam.check(both.gam)
crossvalboth$delta
plot(crossvalboth)


threecol <- surface.pollen.regions.df[grepl("BothPresent", surface.pollen.regions.df$category), c(99,119,101)]
threecol$ambtoart <- log(threecol$ambtoart)

ggplot(data=threecol, aes())
ecoregion.both.gam <- gam(formula=log(ambtoart) ~ splines::bs(AnnP,3), method="ML",
                data=surface.pollen.regions.df[grepl("BothPresent", surface.pollen.regions.df$category),])

