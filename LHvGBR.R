# SK 07/07/2014

##########################################################################

## LORD HOWE VS GBR BREAK
## DOES LORD HOWE HAVE A SIGNIFICANTLY DIFFERENT TRAIT COMPOSITION FROM 
## NULL MODEL AND FROM GBR ASSEMBLAGES?

## COLLABORATORS: ANDREW BAIRD, ERIKA WOOLSEY, MARIA BYRNE

#########################################################################


###################################################
###################################################

## REFORMAT ORIGINAL DATA FROM E. WOOLSEY SO IT 
## CAN BE EASILY USED IN R

d <- read.csv("Cover_Matrix2.csv")
head(d)
d[is.na(d)] <- 0

rownames(d) <- d[,1]
d <- d[,-1]
d[1:10,1:10]

# how many species?
nsp <- nrow(spd)

library(reshape)
d2 <- melt(d)
d3 <- cbind(rep(rownames(d),ncol(d)),d2)
head(d3)
colnames(d3) <- c("species","transect","cover.cm")
save(d3,file="LongFormatCoverData.RData")
table(d3[,2])

# add other columns for site data
island <- rep(c("lizard","one.tree","lord.howe"),c(120*73,120*58,120*72))
habitat <- rep(c("lagoon","crest","lagoon","crest","crest","lagoon"),c(120*36,120*37,120*36,120*22,120*36,120*36))
site <- rep(c("horseshoe","reef","north.of.trimodal","lizard.head","north.point","trimodal","second.lagoon",
              "shark.alley","the.gutter","long.bank","second.lagoon","third.lagoon","blackburn.island",
              "jetty","north.bay","comets.hole","north.bay","sylphs.hole"),
              c(120*12,120*12,120*12,120*12,120*13,120*12,120*12,120*12,120*12,120*8,120*6,120*8,
                120*12,120*12,120*12,120*12,120*12,120*12))
transect <- rep(c(1:203),rep(120,203))
d4 <- cbind(island,habitat,site,transect,d3)
head(d4)
cover.data <- d4[,c(1,2,3,4,5,7)]

# calculate proportional cover (i.e., proportion of total transect length)
cover.data <- cbind(cover.data,cover.data[,6]/1000)
colnames(cover.data)[7] <- "cover.prop"
# add record number
cover.data <- cbind(1:nrow(cover.data),cover.data)
colnames(cover.data)[1] <- "record"

# create empty vector to hold result of loop
cover.relative <- c()
# calculate relative cover (i.e., proportion of total coral cover)
for(i in 1:length(unique(transect))){
      x <- subset(cover.data,cover.data$transect==i)
      cover.relative <- c(cover.relative,x$cover.cm/sum(x$cover.cm))
}

cover.data <- cbind(cover.data,cover.relative)
rel.cover.mat <- cast(cover.data,species~transect,value="cover.relative")

save(rel.cover.mat,file="WideFormatRelCoverData.RData")
save(cover.data,file="LongFormatCoverData.RData")

## PHEW!


############################################################
############################################################

## LOAD PREPARED DATA AND RUN NMDS & ANOSIM

load("LongFormatCoverData.RData")
head(cover.data)

library(vegan)

spbysite <- rel.cover.mat[,-1]
rownames(spbysite) <- rel.cover.mat[,1]

island.short <- rep(c("lizard","one.tree","lord.howe"),c(73,58,72))
region <- rep(c("GBR","lord.howe"),c(73+58,72))
habitat.short <- rep(c("lagoon","crest","lagoon","crest","crest","lagoon"),c(36,37,36,22,36,36))
anosim(t(spbysite),island.short,permutations=1000,distance="bray")
anosim(t(spbysite),region,permutations=1000,distance="bray")
anosim(t(spbysite),habitat.short,permutations=1000,distance="bray")

# ANOSIM between each island pair
island.LizOT <- rep(c("lizard","one.tree"),c(73,58))
anosim(t(spbysite[,1:131]),island.LizOT,permutations=1000,distance="bray")
island.LizLH <- rep(c("lizard","lord.howe"),c(73,72))
anosim(t(spbysite[,c(1:73,132:203)]),island.LizLH,permutations=1000,distance="bray")
island.LHOT <- rep(c("one.tree","lord.howe"),c(58,72))
anosim(t(spbysite[,74:203]),island.LHOT,permutations=1000,distance="bray")

MDS <- metaMDS(t(spbysite),trymax=50)
plot(MDS)
MDS$points

par(mfcol=c(2,1))
# plot MDS with each island a different colour
plot(MDS,display="sites",col="white",cex.lab=1,font.lab=2,main="MDS transects")
points(MDS, display="sites",pch=20,col=1,select=island.short== "lizard")
points(MDS, display="sites",pch=20,col=2,select=island.short== "one.tree")
points(MDS, display="sites",pch=20,col=3,select=island.short== "lord.howe")
# add to plot 95% confidence ellipses for each island
ordiellipse(MDS,island.short,kind="sd",conf=0.90,lty=2)
legend("topright",c("lizard","one tree","lord howe"),pch=20,col=c(1,2,3))
# plot MDS with each habitat a different colour
plot(MDS,display="sites",col="white",cex.lab=1,font.lab=2,main="MDS transects")
points(MDS, display="sites",pch=20,col=1,select=habitat.short== "lagoon")
points(MDS, display="sites",pch=20,col=2,select=habitat.short== "crest")
# Lord Howe separates out on axis 1, lagoon/crest on axis 2
# add to plot 95% confidence ellipses for each habitat
ordiellipse(MDS,habitat.short,kind="sd",conf=0.90,lty=2)
legend("topright",c("lagoon","crest"),pch=20,col=c(1,2))

# isolate MDS axis scores for use in GLMM
MDSax <- MDS$points
MDSax1 <- MDSax[,1]
MDSax2 <- MDSax[,2]



#############################################
#############################################

## PREP DATA FOR NULL MODEL OF SOURCE POOL

traits <- read.csv("CoralTraitsAug2014SK.csv")
load("species by province matrix, published data.RData")

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

transects.w <- round(lizt.inc[,15]/sum(lizt.inc[,15],na.rm=T),3)
median(transects.w,na.rm=T)

colnames(lizt.inc) <- c(colnames(lizt.inc)[1:17],"transects.w","sites.w","reefs.w")



liz.ts <- trait.space(liz,ndraw,weights=  )


##############################################
##############################################

## ABUNDANCE WEIGHTED RANDOM DRAWS

