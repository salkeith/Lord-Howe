# SK 21/08/2014

##########################################################################
##########################################################################

## LORD HOWE VS GBR BREAK
## DOES LORD HOWE HAVE A SIGNIFICANTLY DIFFERENT TRAIT COMPOSITION FROM 
## GBR ASSEMBLAGES?

## COLLABORATORS: ANDREW BAIRD, ERIKA WOOLSEY, MARIA BYRNE

## INCLUDES BOOTSTRAP TO ACCOUNT FOR DIFFERENT SAMPLING EFFORT #########

#########################################################################
#########################################################################

rm(list=ls())

library(vegan)
library(lme4)
library(visreg)
library(reshape)
library(MuMIn)
library(effects)
library(MASS)

###################################################
###################################################

load("LongFormatCoverData.RData")
head(cover.data)

# Lizard, crest, North Point site has 13 transects
# One Tree, crest, Long Bank has 8 transects
# One Tree, crest, Third Lagoon has 8 transects
# One Tree, crest, Second Lagoon has 6 transects
# All other sites have 12 transects

# BOOTSTRAP SAMPLES TO ACCOUNT FOR DIFFERENT SAMPLING EFFORT
# For sites >6 transects, randomly select 6 transects and calculate mean relative cover across those transects
# Do this 1000 times and get confidence intervals around means

# Create unique identifier for each transect/habitat/site/island combination
id <- paste(cover.data$island,cover.data$site,cover.data$habitat,cover.data$transect,sep="_")
d <- cbind(id,cover.data)
head(d)

# Get the unique identifiers and separate by site so can select from these lists for bootstrapping
lizc1id <- unique(id)[37:48]
otc1id <- unique(id)[110:117]
lhc1id <- unique(id)[132:143]
lizc2id <- unique(id)[49:61]
otc2id <- unique(id)[118:123]
lhc2id <- unique(id)[144:155]
lizc3id <- unique(id)[62:73]
otc3id <- unique(id)[124:131]
lhc3id <- unique(id)[156:167]
id.vec <- c("lizc1id","otc1id","lhc1id","lizc2id","otc2id","lhc2id","lizc3id","otc3id","lhc3id")


bootstrap.n <- 1000  # number of bootstraps
sample.n <- 6  # number of samples
boot.id <- array(NA,dim=c(length(id.vec),sample.n,bootstrap.n))  # to hold sampled transect names
boot.d <- list()  # list because it doesn't seem to work if I try to save results into an array 

# generate sample of sample.n transects from each island/site on crest
for(z in 1:bootstrap.n){
   for(i in 1:length(id.vec)){
      # list of random selected transects in island/site
      boot.id[i,,z] <- get(id.vec[i])[sample(1:length(get(id.vec[i])),sample.n)]   
   }
   # subset cover data to only the 6 samples for each island/site
   # selects the unique id if it matches an id found in boot.id
   boot.d[[z]] <- subset(d,d$id %in% boot.id[,,z])
}

# order transects alphabetically so ensure they match up in other bits of code
boot.d <- lapply(boot.d,function(x) x[order(x$id),])
# simplify boot.id into 2D array
boot.id2 <- matrix(NA,54,bootstrap.n)
for(z in 1:bootstrap.n){
   boot.id2[,z] <- sort(as.vector(t(boot.id[,,z])))
}
boot.id <- t(boot.id2)

save(boot.d,boot.id,file="BootstrappedData.RData")



#######################################################
#######################################################

## CALCULATE RELATIVE ABUNDANCE

load("BootstrappedData.RData")

boot.mat <- list()
for(z in 1:bootstrap.n){
   # create empty vector to hold result of loop
   cover.relative.bs <- c()
   # calculate relative cover (i.e., proportion of total coral cover)
      for(tr in 1:(sample.n*9)){
         x <- subset(boot.d[[z]],boot.d[[z]]$id==boot.id[z,tr])
         cover.relative.bs <- c(cover.relative.bs,x$cover.cm/sum(x$cover.cm))
      }
   boot.d[[z]] <- cbind(boot.d[[z]],cover.relative.bs)

}
# check data looks right
head(boot.d[[1]])



#######################################################
#######################################################

## JOIN BOOTSTRAPPED CREST DATA BACK WITH ORIGINAL 
## LAGOON DATA & FORMAT INTO SITE BY SPECIES MATRIX


# create list to hold output
boot.join <- list()
# loop through bootstrapped data frames
for(z in 1:bootstrap.n){
   # add new column so it will bind with bootstrapped data
   d$cover.relative.bs <- d$cover.relative
   boot.join[[z]] <- rbind(subset(d,d$habitat=="lagoon"),boot.d[[z]])
   # format as site by species matrix
   boot.mat[[z]] <- cast(boot.join[[z]],species~id,value="cover.relative.bs")
   rownames(boot.mat[[z]]) <- boot.mat[[z]][,1]
   boot.mat[[z]] <- boot.mat[[z]][,-1]
}
# check data looks right
boot.mat[[1]][1:30,1:10]
dim(boot.mat[[1]])


save(boot.mat,boot.join,boot.d,boot.id,file="BootstrappedData.RData")


#######################################################
#######################################################

## GET SIMILARITY MEASURES FOR BOOTSTRAPPED DATA
## AND THEN USE THE MEAN SIMILARITY BETWEEN EACH SITE

load("BootstrappedData.RData")

boot.dist <- list()
bootstrap.n <- 1000

for(z in 1:bootstrap.n){
   boot.dist[[z]] <- vegdist(t(boot.mat[[z]]),method="jaccard")
}

boot.dist.mat <- sapply(boot.dist,as.numeric)
boot.dist.rowmean <- rowMeans(boot.dist.mat)
boot.dist.mean <- as.dist(boot.dist.rowmean)

MDS <- monoMDS(boot.dist.mean,maxit=1,k=2)
par(mfcol=c(2,1))
stressplot(MDS)
plot(MDS)

par(mfcol=c(1,1))
# plot MDS with each island a different colour
plot(MDS,display="sites",col="white",cex.lab=2,cex.axis=1.5,font.lab=2,ylim=c(-1,1),xlim=c(-1,1))
points(MDS, display="sites",pch=1,col=1,select=is== "lizard")
points(MDS, display="sites",pch=19,col="grey",select=is== "one.tree")
points(MDS, display="sites",pch=19,col=1,select=is== "lord.howe")
# add to plot 95% confidence ellipses for each island
ordiellipse(MDS,is,kind="sd",conf=0.95,lty=2)
legend("topright",c("Lizard","One Tree","Lord Howe"),pch=c(1,19,19),col=c(1,"grey",1),cex=2.5)


#######################################################
#######################################################

## RUN ANOSIM & NMDS ON BOOTSTRAPPED DATA (i.e., x 1000)
## GET MEAN AND CONFIDENCE INTERVALS OF AXIS SCORES

habitat.short <- rep(c("lagoon","crest"),c(108,54))
island.short <- rep(c("lizard","one.tree","lord.howe","lizard","one.tree","lord.howe"),c(36,36,36,18,18,18))
region <- rep(c("GBR","lord.howe","GBR","lord.howe"),c(72,36,36,18))

empty.sp <- c()  # empty vector to record absent species 
boot.anosim <- matrix(NA,bootstrap.n,6)
for(z in 1:bootstrap.n){
   if(z %in% seq(50,1000,50)) print(z)
   # remove species with no records (to make ANOSIM work)
   x <- subset(boot.mat[[z]],rowSums(boot.mat[[z]])>0)
   # record number of species absent in each data set
   empty.sp[z] <- nrow(boot.mat[[z]]-nrow(x))
   # ANOSIM
   x <- t(x)
   a.island <- anosim(x,island.short,permutations=1000,distance="bray")
   boot.anosim[z,1] <- a.island$statistic
   boot.anosim[z,2] <- a.island$signif
   a.region <- anosim(x,region,permutations=1000,distance="bray")
   boot.anosim[z,3] <- a.region$statistic
   boot.anosim[z,4] <- a.region$signif
   a.habitat <- anosim(x,habitat.short,permutations=1000,distance="bray")
   boot.anosim[z,5] <- a.habitat$statistic
   boot.anosim[z,6] <- a.habitat$signif
}
# check the data look right
head(boot.anosim)
summary(empty.sp)
colnames(boot.anosim) <- c("island.R","island.p","region.R","region.p","habitat.R","habitat.p")
colMeans(boot.anosim)  # mean values across all bootstrapped data sets

# might need to use mean rank of score rather actual score
# because axes might look totally different for each MDS
boot.mds.sites <- list()
boot.mds.sp <- list()
mds.stress <- c()  # record stress value for each MDS

for(z in 1:bootstrap.n){
   if(z %in% seq(50,1000,50)) print(z)
   x <- subset(boot.mat[[z]],rowSums(boot.mat[[z]])>0)
   x <- t(x)
   mds <- metaMDS(x,trymax=50)
   mds.stress[z] <- round(mds$stress,3)   
   boot.mds.sites[[z]] <- mds$points
   boot.mds.sp[[z]] <- mds$species
}


save(boot.anosim,empty.sp,mds.stress,boot.mds.sites,boot.mds.sp,
        file="BootstrappedANOSIMandMDS.RData")


par(mfcol=c(3,3))
for(i in 110:112){
   bms <- boot.mds.sites[[i]]
   plot(bms,col="white",main="island")   
   points(bms[which(island.short=="lizard"),],col=1)
   points(bms[which(island.short=="one.tree"),],col=2)
   points(bms[which(island.short=="lord.howe"),],col=3)
}
for(i in 110:112){
   bms <- boot.mds.sites[[i]]
   plot(bms,col="white",main="habitat")   
   points(bms[which(habitat.short=="lagoon"),],col=1)
   points(bms[which(habitat.short=="crest"),],col=2)
}
for(i in 110:112){
   bms <- boot.mds.sites[[i]]
   plot(bms,col="white",main="region")   
   points(bms[which(region=="GBR"),],col=1)
   points(bms[which(region=="lord.howe"),],col=2)
}



#########################################################
## SK 25/08/2014
## USE LAGOON SITES ONLY

load("WideFormatRelCoverData.RData")

spbysite <- rel.cover.mat[,-1]
sitebysp <- t(spbysite)
colnames(sitebysp) <- rel.cover.mat[,1]
island.short <- rep(c("lizard","one.tree","lord.howe"),c(73,58,72))
region <- rep(c("GBR","lord.howe"),c(73+58,72))
habitat.short <- rep(c("lagoon","crest","lagoon","crest","crest","lagoon"),c(36,37,36,22,36,36))

sbsp <- sitebysp[which(habitat.short=="lagoon"),]
is <- island.short[which(habitat.short=="lagoon")]
r <- region[which(habitat.short=="lagoon")]
anosim(sbsp,is,permutations=1000,distance="bray")
anosim(sbsp,r,permutations=1000,distance="bray")
# ANOSIM between each island pair
anosim(sbsp[which(is=="lizard"|is=="one.tree"),],is[which(is=="lizard"|is=="one.tree")],permutations=1000,distance="bray")
anosim(sbsp[which(is=="lizard"|is=="lord.howe"),],is[which(is=="lizard"|is=="lord.howe")],permutations=1000,distance="bray")
anosim(sbsp[which(is=="lord.howe"|is=="one.tree"),],is[which(is=="lord.howe"|is=="one.tree")],permutations=1000,distance="bray")


MDS <- metaMDS(sbsp,trymax=50)
plot(MDS)

par(mfcol=c(1,1))
# plot MDS with each island a different colour
plot(MDS,display="sites",col="white",cex.lab=2,cex.axis=1.5,font.lab=2,ylim=c(-1,1),xlim=c(-1,1))
points(MDS, display="sites",pch=1,col=1,select=is== "lizard")
points(MDS, display="sites",pch=19,col="grey",select=is== "one.tree")
points(MDS, display="sites",pch=19,col=1,select=is== "lord.howe")
# add to plot 95% confidence ellipses for each island
ordiellipse(MDS,is,kind="sd",conf=0.95,lty=2)
legend("topright",c("Lizard","One Tree","Lord Howe"),pch=c(1,19,19),col=c(1,"grey",1),cex=2.5)

# isolate MDS axis scores for use in GLMM
MDSsp <- MDS$species
MDSsp <- as.data.frame(MDSsp[,1:2])
MDSsp <- cbind(colnames(sbsp),MDSsp)
colnames(MDSsp) <- c("species","ax1","ax2")
rownames(MDSsp) <- NULL
# remove all the species not present in lagoon
MDSsp <- MDSsp[-c(6,7,10,21,27,31,44,52,54,55,66,89,92,96,102,104,117),]

traits <- read.csv("Traits_LHI_AB26Aug14.csv")
# Dummy variables in data that I changed in csv file
# water clarity: turbid = -1; both = 0; clear = 1
# wave exposure: exposed = -1; broad = 0; protected = 1
# Which traits do we want to include?
colnames(traits)
# convert all columns to numeric format
for(i in 2:ncol(traits)){
   traits[,i] <- as.numeric(traits[,i])
}
str(traits)
# merge species traits with NMDS axis scores
sptr <- merge(MDSsp,traits,by="species",all.x=T)
save(sptr,file="SpeciesAxisScoresTraitsLagoon.RData")

#### ACROPORA CUNEATA AND A. PALIFERA ARE ISOPORA


####################################################
####################################################

## REGRESSION MODEL
## HOW WELL CAN TRAITS PREDICT THE NMDS AXIS SCORE?

####################################################
####################################################

load("SpeciesAxisScoresTraitsLagoon.RData")

###################################################
## REGRESSION DIAGNOSTICS

# normalize trait data (mean = 0, sd = 1)
t.norm <- scale(sptr[,4:ncol(sptr)])

# visualise distribution of values for each normalized variable
par(mfcol=c(3,3))
for(z in 1:9){
   hist(t.norm[,z], main = colnames(t.norm)[z])
}
par(mfcol=c(3,2))
for(z in 10:ncol(t.norm)){
   hist(t.norm[,z], main = colnames(t.norm)[z])
}

# check collinearity of traits
pairs(t.norm)
multcol <- cor(t.norm,use="complete.obs")
multcol< (-0.6) | multcol>0.6
# CORRELATIONS:
# depth range & lower depth
# larval development mode & development rate
# ACTION:
# remove lower depth because less information than depth range
# remove development rate
t.norm <- t.norm[,-c(4,10)]

# put data together into one data frame ready for regression
d <- cbind(sptr[,1:3],t.norm)
colnames(d) <- c(colnames(d)[1:3],"colony.size","valley.size","genus.age","upp.depth","depth.range",
                 "water.clarity","wave.exposure","larval.mode","sex","zoox")
summary(d)

# check quadratics for traits (if delta AICc >3, use quadratic)
for(i in 4:length(d)){
   par(mfcol=c(2,2))
   print(colnames(d)[i])
   lin <- lm(ax1~d[,i],data=d)
   plot(lin)
   quad <- lm(ax1~d[,i]+I(d[,i]^2),data=d)  
   print(AICc(lin,quad))
}

# remove rows with NAs to ensure models all use same underlying data (no na.omit)
d <- d[-c(50,93,96),]

# recheck regression diagnostics
for(i in 4:length(d)){
   par(mfcol=c(2,2))
   print(colnames(d)[i])
   lin <- lm(ax1~d[,i],data=d)
   plot(lin)
}

# save file so don't have to repeat steps above
save(d,file="DataForRegressionLHIlagoon.RData")



###################################################
###################################################
## LINEAR REGRESSION AXIS 1

load(d,file="DataForRegressionLHI.RData")

# Full model including interactions that make biological sense and genus as a random effect
mod <- lm(ax1~colony.size+valley.size+upp.depth+depth.range+wave.exposure+water.clarity+genus.age+sex+
               larval.mode+zoox+upp.depth*depth.range+wave.exposure*water.clarity+larval.mode*sex,
               data=d,na.action=na.fail)
summary(mod)
# remove non-significant interactions
mod <- lm(ax1~colony.size+valley.size+upp.depth+depth.range+wave.exposure+water.clarity+genus.age+sex+
               larval.mode+zoox+upp.depth*depth.range,data=d,na.action=na.fail)
summary(mod)

# take a look at the partial coefficients
par(mfcol=c(3,3))
visreg(mod,partial=F)



#############################################
# MODEL SELECTION AND AVERAGING
# model selection by delta < 3 (Bolker 2009) with polynomials
dmod <- dredge(mod)
mod.sel <- summary(model.avg(get.models(subset(dmod,delta < 3))))
mod.sel

mod.best <- lm(ax1~larval.mode+genus.age,data=d,na.action=na.fail)
summary(mod.best)

save(mod.sel,mod.best,file="ModelAveragedandBestModResultLagoon.RData")

write.csv(coefTable(mod.sel),"ModelAveragedCoefsLagoon.csv")

#############################################
# MODEL-AVERAGED CONFIDENCE INTERVALS
cimod<- confint(mod.sel)
# full model-averaged coefficients not needed because all models in the set contain
# both interaction terms and main effects, or cubic term
comod <- coef(mod.sel)
CImod <- cbind(comod,cimod) 
CImod
# return variables that do not have coefficients that overlap zero
# (i.e., those that are significant)
# N.B. doesn't apply to interactions or cubic 
CImod[!(CImod[,2]<0 & CImod[,3]>0),]

write.csv(CImod,"ModelAveragedCoefs95CILagoon.csv")


#######################################################
# PARTIAL COEFFICIENT PLOTS for significant variables 
# manual plotting required, visreg doesn't work on model-averged object

pdf("PartialCoefsLordHoweTraitsLagoonBest.pdf")

par(mfrow=c(2,2))

# LARVAL MODE
x <- seq(min(d$larval.mode),max(d$larval.mode),length=100)
larval.pred <- with(d,(predict(mod.best,interval="confidence",newdata=data.frame(larval.mode=x,
                 upp.depth=mean(upp.depth),water.clarity=mean(water.clarity)))))
plot(ax1~larval.mode,data=d,col="lightgray",pch=16,ylab="NMDS Axis 1",
     xlab="Larval mode",cex.lab=1.5)
points(larval.pred[,1]~x,type="l",col=1,lwd=2)
lines(larval.pred[,2]~x,lty=2,type="l")
lines(larval.pred[,3]~x,lty=2,type="l")

# UPPER DEPTH
x <- seq(min(d$upp.depth),max(d$upp.depth),length=100)
upp.pred <- with(d,(predict(mod.best,interval="confidence",newdata=data.frame(upp.depth=x,
                  larval.mode=mean(larval.mode),water.clarity=mean(water.clarity)))))
plot(ax1~upp.depth,data=d,col="lightgray",pch=20,ylab="NMDS Axis 1",
     xlab="Upper depth",cex.lab=1.5)
points(upp.pred[,1]~x,type="l",col=1,lwd=2)
lines(upp.pred[,2]~x,lty=2,type="l")
lines(upp.pred[,3]~x,lty=2,type="l")

# WATER CLARITY
x <- seq(min(d$water.clarity),max(d$water.clarity),length=100)
wc.pred <- with(d,(predict(mod.best,interval="confidence",newdata=data.frame(water.clarity=x,
                  larval.mode=mean(larval.mode),upp.depth=mean(upp.depth)))))
plot(ax1~water.clarity,data=d,col="lightgray",pch=20,ylab="NMDS Axis 1",
     xlab="Water clarity",cex.lab=1.5)
points(wc.pred[,1]~x,type="l",col=1,lwd=2)
lines(wc.pred[,2]~x,lty=2,type="l")
lines(wc.pred[,3]~x,lty=2,type="l")

dev.off()



###################################################
###################################################
## LINEAR REGRESSION AXIS 2

load(d,file="DataForRegressionLHILagoon.RData")

# check quadratics for traits
for(i in 5:length(d)){
   par(mfcol=c(2,2))
   print(colnames(d)[i])
   lin <- lm(ax2~d[,i],data=d)
   plot(lin)
   quad <- lm(ax2~d[,i]+I(d[,i]^2),data=d)  
   print(AICc(lin,quad))
}

# Full model including interactions that make biological sense and genus as a random effect
mod <- lmer(ax2~valley.size+upp.depth+depth.range+wave.exposure+water.clarity+genus.age+sex+
               larval.mode+zoox+upp.depth*depth.range+wave.exposure*water.clarity+larval.mode*sex+
               (1|genus),data=d,na.action=na.fail)
summary(mod)
dotplot(ranef(mod, postVar=TRUE))
# mixed model not required - very low variance accounted for and overlapping confidence intervals in dotplot
mod <- lm(ax2~valley.size+upp.depth+depth.range+wave.exposure+water.clarity+genus.age+sex+
             larval.mode+zoox+upp.depth*depth.range+wave.exposure*water.clarity+larval.mode*sex,data=d,na.action=na.fail)
summary(mod)
# remove non-significant interactions
mod <- lm(ax2~valley.size+upp.depth+depth.range+wave.exposure+water.clarity+genus.age+sex+
             larval.mode+zoox,data=d,na.action=na.fail)
summary(mod)

# take a look at the partial coefficients
par(mfcol=c(3,3))
visreg(mod,partial=F)



#############################################
# MODEL SELECTION AND AVERAGING
# model selection by delta < 3 (Bolker 2009) with polynomials
dmod <- dredge(mod)
mod.sel <- summary(model.avg(get.models(subset(dmod,delta < 3))))
mod.sel

mod.best <- lm(ax2~upp.depth+wave.exposure+genus.age,data=d,na.action=na.fail)
summary(mod.best)

save(mod.sel,mod.best,file="ModelAveragedandBestModResultLagoon2.RData")

write.csv(coefTable(mod.sel),"ModelAveragedCoefsLagoon2.csv")

#############################################
# MODEL-AVERAGED CONFIDENCE INTERVALS
cimod<- confint(mod.sel)
# full model-averaged coefficients not needed because all models in the set contain
# both interaction terms and main effects, or cubic term
comod <- coef(mod.sel)
CImod <- cbind(comod,cimod) 
CImod
# return variables that do not have coefficients that overlap zero
# (i.e., those that are significant)
# N.B. doesn't apply to interactions or cubic 
CImod[!(CImod[,2]<0 & CImod[,3]>0),]

write.csv(CImod,"ModelAveragedCoefs95CILagoon2.csv")


#######################################################
# PARTIAL COEFFICIENT PLOTS for significant variables 
# manual plotting required, visreg doesn't work on model-averged object

pdf("PartialCoefsLordHoweTraitsLagoonBest2.pdf")

par(mfrow=c(2,2))

# GENUS AGE
x <- seq(min(d$genus.age),max(d$genus.age),length=100)
ga.pred <- with(d,(predict(mod.best,interval="confidence",newdata=data.frame(genus.age=x,
                    upp.depth=mean(upp.depth),wave.exposure=mean(wave.exposure)))))
plot(ax1~genus.age,data=d,col="lightgray",pch=16,ylab="NMDS Axis 2",
     xlab="Genus age",cex.lab=1.5)
points(ga.pred[,1]~x,type="l",col=1,lwd=2)
lines(ga.pred[,2]~x,lty=2,type="l")
lines(ga.pred[,3]~x,lty=2,type="l")

# UPPER DEPTH
x <- seq(min(d$upp.depth),max(d$upp.depth),length=100)
ud.pred <- with(d,(predict(mod.best,interval="confidence",newdata=data.frame(upp.depth=x,
                      genus.age=mean(genus.age),wave.exposure=mean(wave.exposure)))))
plot(ax1~upp.depth,data=d,col="lightgray",pch=16,ylab="NMDS Axis 2",
     xlab="Upper depth",cex.lab=1.5)
points(ud.pred[,1]~x,type="l",col=1,lwd=2)
lines(ud.pred[,2]~x,lty=2,type="l")
lines(ud.pred[,3]~x,lty=2,type="l")

# WAVE EXPOSURE
x <- seq(min(d$wave.exposure),max(d$wave.exposure),length=100)
we.pred <- with(d,(predict(mod.best,interval="confidence",newdata=data.frame(wave.exposure=x,
                       upp.depth=mean(upp.depth),genus.age=mean(genus.age)))))
plot(ax1~wave.exposure,data=d,col="lightgray",pch=16,ylab="NMDS Axis 2",
     xlab="Wave exposure",cex.lab=1.5)
points(we.pred[,1]~x,type="l",col=1,lwd=2)
lines(we.pred[,2]~x,lty=2,type="l")
lines(we.pred[,3]~x,lty=2,type="l")

dev.off()





