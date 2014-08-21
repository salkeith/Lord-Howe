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









