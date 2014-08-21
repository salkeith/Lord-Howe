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


save(boot.anisim,empty.sp,mds.stress,boot.mds.sites,boot.mds.sp,
        file="BootstrappedANOSIMandMDS.RData")


par(mfcol=c(3,3))
for(i in 10:12){
   bms <- boot.mds.sites[[i]]
   plot(bms,col="white",main="island")   
   points(bms[which(island.short=="lizard"),],col=1)
   points(bms[which(island.short=="one.tree"),],col=2)
   points(bms[which(island.short=="lord.howe"),],col=3)
}
for(i in 10:12){
   bms <- boot.mds.sites[[i]]
   plot(bms,col="white",main="habitat")   
   points(bms[which(habitat.short=="lagoon"),],col=1)
   points(bms[which(habitat.short=="crest"),],col=2)
}
for(i in 10:12){
   bms <- boot.mds.sites[[i]]
   plot(bms,col="white",main="region")   
   points(bms[which(region=="GBR"),],col=1)
   points(bms[which(region=="lord.howe"),],col=2)
}






