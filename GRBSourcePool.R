# SK 07/07/2014

##########################################################################

## GBR SOURCE POOL ANALYSIS

## COLLABORATORS: ANDREW BAIRD, 

#########################################################################

library(vegan)

load("WideFormatRelCoverData.RData")
load("LongFormatCoverData.RData")
load("species by province matrix, published data.RData")
traits <- read.csv("CoralTraitsAug2014SK.csv")

## PREP DATA FOR NULL MODEL OF SOURCE POOL

spbysite <- rel.cover.mat[,-1]
rownames(spbysite) <- rel.cover.mat[,1]

# LHI has some temperate corals so where do we draw the line around the source pool?
# Perhaps we can test with entire E. coast Aus? Then only GBR?
# Or I could use the Australian province? Probably most easily defensible

aus <- mat.pd[,3]
# subset to only those species that are present
aus.pres <- subset(aus,aus==1)
sum(aus.pres)  # how many species?
aus.sp <- names(aus.pres)  # list of species present in Aus province

# how many species did we observe at each island?
liz <- sum(rowSums(spbysite[,1:72])>0)
ot <- sum(rowSums(spbysite[,73:131])>0)
lhi <- sum(rowSums(spbysite[,132:203])>0)
# subset matrix so there's a separate one for each island
liz.all <- spbysite[,1:72]
ot.all <- spbysite[,73:131]
lhi.all <- spbysite[,132:203]
# get species list for each island
liz.sp <- as.data.frame(rownames(subset(liz.all,rowSums(liz.all)>0)))
ot.sp <- as.data.frame(rownames(subset(ot.all,rowSums(ot.all)>0)))
lhi.sp <- as.data.frame(rownames(subset(lhi.all,rowSums(lhi.all)>0)))
colnames(liz.sp) <- "species"
colnames(ot.sp) <- "species"
colnames(lhi.sp) <- "species"

# NULL MODELS ONLY MAKE SENSE IN THE CONTEXT OF TRAITS
# proportion of traits in null assemblage vs in LIZ/OT/LHI assemblages
# NOTE this is on presence data only

# Which traits do we want to include?
colnames(traits)
traits <- traits[,c(1,13,14:17,20,21,23,30:34)]
# convert all columns to numeric format
for(i in 2:ncol(traits)){
   traits[,i] <- as.numeric(traits[,i])
}
str(traits)

# find observed trait space within sites
lizt <- merge(liz.sp,traits,by="species")
ott <- merge(ot.sp,traits,by="species")
lhit <- merge(lhi.sp,traits,by="species")
lizt.mean <- colMeans(lizt[,2:ncol(lizt)],na.rm=T)
ott.mean <- colMeans(ott[,2:ncol(lizt)],na.rm=T)
lhit.mean <- colMeans(lhit[,2:ncol(lizt)],na.rm=T)


############################################
############################################
# SK 05/08/2014
# GENERATE NULL MODEL AND OUTPUT RESULTS
# AS PDF FILE OF HISTOGRAMS

source("RandomDrawFunction.txt")
source("PlotTraitDistributionFunction.txt")
# trait.space <- function(island.name,sp.pool,ndraw=1000,weights=NULL)
# plot.trait <- function(island.name,island.trait.mean,trait.space.res)

ndraw <- 10000  # set number of random draws

# Lizard Island random draws
liz.ts <- trait.space(liz,aus.sp,ndraw)
plot.trait("Lizard",lizt.mean,liz.ts)
# One Tree Island random draws
ot.ts <- trait.space(ot,aus.sp,ndraw)
plot.trait("One Tree",ott.mean,ot.ts)
# Lord Howe Island random draws
lhi.ts <- trait.space(lhi,aus.sp,ndraw)
plot.trait("Lord Howe",lhit.mean,lhi.ts)


###########################################
###########################################
## CALCULATE 95% CONFIDENCE INTERVALS

source("TraitConfIntFunction.txt")
source("SignificanceCI.txt")
# trait.CI <- function(island.ts,island.trait.mean,quantile.lower,quantile.upper)
# sig.table <- function(island.CI1,island.CI2,island.CI3)

CIlow <- 0.025
CIupp <- 0.975
lizCI <- trait.CI(liz.ts,lizt.mean,CIlow,CIupp)
otCI <- trait.CI(ot.ts,ott.mean,CIlow,CIupp)
lhiCI <- trait.CI(lhi.ts,lhit.mean,CIlow,CIupp)

sig.res <- sig.table(lizCI,otCI,lhiCI)
sig.res 


##############################################
##############################################

## INCIDENCE WEIGHTED RANDOM DRAWS

# Load data from AB with abundance of GBR species
GBRab <- read.csv("LIT_GBRNorth&South_ahb_abundance.csv")
head(GBRab)

# should incidence be calculated on number of transects (nested within sites),
# sites (n = 39, nested in reefs, most reefs have 2 sites but some only 1) 
# or number of reefs (n = 22, lower resolution)? Try all
num.transects <- as.data.frame(table(GBRab$species))
num.sites1 <- as.data.frame(table(GBRab$species,GBRab$site))
num.sites2 <- subset(num.sites1,num.sites1$Freq>0)
num.sites <- as.data.frame(table(num.sites2$Var1))
num.reefs1 <- as.data.frame(table(GBRab$species,GBRab$reef))
num.reefs2 <- subset(num.reefs1,num.reefs1$Freq>0)
num.reefs <- as.data.frame(table(num.reefs2$Var1))
inc <- cbind(num.transects,num.sites[,2],num.reefs[,2])
colnames(inc) <- c("species","transects","sites","reefs") 
# save data frame with species incidence at transect, site and reef level
save(inc,file="SpeciesIncidenceGBR.RData")

# SK 06/08/2014
# how many of the species from the islands are represented in AB data?
length(merge(lizt,inc,by="species")[,1])  # 70 out of 84
length(merge(ott,inc,by="species")[,1])   # 49 out of 65
length(merge(lhit,inc,by="species")[,1])  # 25 out of 30
# we don't know the incidence of the missing species so use a probability
# of sampling = 0.5 for those. Creates some issues though!
# Merge and keep all species present on the islands (i.e., all.x=TRUE)
lizt.inc <- merge(lizt,inc,by="species",all.x=T) 
ott.inc <- merge(ott,inc,by="species",all.x=T)
lhit.inc <- merge(lhit,inc,by="species",all.x=T)

# convert incidence to 'probability', divide each incidence by column sum
# replace NA values for incidence with ??? Median? Lowest? <lowest? half the lowest?
# divide lowest by some number... CHECK SENSITIVITY to number later!!!
source("WeightIncidenceFunction.txt")
# weight.inc <- function(island.inc,cols,na.divide)
lizt.inc.prob <- weight.inc(lizt.inc,c(15,16,17),2)
ott.inc.prob <- weight.inc(ott.inc,c(15,16,17),2)
lhit.inc.prob <- weight.inc(lhit.inc,c(15,16,17),2)



##############################################
##############################################

## ABUNDANCE WEIGHTED RANDOM DRAWS

