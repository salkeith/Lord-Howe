# SK 07/07/2014
# SK SUBSTANTIALLY MODIFIED 06/08/2014

##########################################################################

## LORD HOWE VS GBR BREAK
## DOES LORD HOWE HAVE A SIGNIFICANTLY DIFFERENT TRAIT COMPOSITION FROM 
## GBR ASSEMBLAGES?

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

## PREP DATA FOR MIXED EFFECTS MODEL
## BASED ON NMDS AXIS SCORES

# Which traits do we want to include?
colnames(traits)
traits <- traits[,c(1,13,14:17,20,21,23,30:34)]
# convert all columns to numeric format
for(i in 2:ncol(traits)){
   traits[,i] <- as.numeric(traits[,i])
}
str(traits)

lizt <- merge(liz.sp,traits,by="species")
ott <- merge(ot.sp,traits,by="species")
lhit <- merge(lhi.sp,traits,by="species")


