# SK 07/07/2014

########################################

## LORD HOWE VS GBR BREAK

########################################

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

## LOAD PREPARED DATA

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

## NULL MODEL OF SOURCE POOL

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


ndraw <- 10000  # set number of random draws

# Function for generating trait space (distributions) in random assemblages
# i.e., mean of each trait for each random draw
trait.space <- function(island.name,ndraw=1000){
   # prepare data frame to hold result
   trait.space.res <- as.data.frame(matrix(nrow=ndraw,ncol=ncol(traits)-1))
   colnames(trait.space.res) <- colnames(traits)[2:ncol(traits)]
   # loop through random draws
   for(i in 1:ndraw){
      rd <- as.data.frame(sample(aus.sp,island.name,replace=F))  
      colnames(rd) <- "species"
      # merge randomly drawn species with traits
      rdt <- merge(rd,traits,by="species")
      trait.space.res[i,] <- colMeans(rdt[,2:ncol(rdt)],na.rm=T)
   }
   return(trait.space.res)
}

# Function to plot histograms of trait distributions from random draws
# with vertical line for observed mean trait value
plot.trait <- function(island.name,trait.mean,trait.space.res){
   pdf(paste(island.name,"trait histograms random draw.pdf"))
   par(mfcol=c(3,2))
   # plot variables with -1, 0, 1 values
   for(z in c(6,7,9:11)){
      hist(trait.space.res[,z],main=colnames(trait.space.res)[z],xlim=c(-1,1)) # random values
      abline(v=trait.mean[z],col=2)  # observed value
   }
   # plot variables with continuous values
   par(mfcol=c(3,3))
   for(z in c(1:5,8,12,13)){
      hist(trait.space.res[,z],main=colnames(trait.space.res)[z]) # random values
      abline(v=trait.mean[z],col=2) # observed value
   }
   dev.off()   
}

# Lizard Island random draws
liz.ts <- trait.space(liz,ndraw)
plot.trait("Lizard",lizt.mean,liz.ts)
# One Tree Island random draws
ot.ts <- trait.space(ot,ndraw)
plot.trait("One Tree",ott.mean,ot.ts)
# Lord Howe Island random draws
lhi.ts <- trait.space(lhi,ndraw)
plot.trait("Lord Howe",lhit.mean,lhi.ts)





##### CALCULATE 95% CONFIDENCE INTERVALS
cil.liz <- apply(trait.space.liz,2,function(x) mean(x)-(1.96*(sd(x)/sqrt(nrow(trait.space.liz)))))
ciu.liz <- apply(trait.space.liz,2,function(x) mean(x)+(1.96*(sd(x)/sqrt(nrow(trait.space.liz)))))
liz.ci <- cbind(lizt.mean,cil.liz,ciu.liz)
sig.liz <- as.data.frame(ifelse(liz.ci[,1]>cil.liz & liz.ci[,1]<ciu.liz,0,1))
colnames(sig.liz) <- "lizard"

cil.ot <- apply(trait.space.ot,2,function(x) mean(x)-(1.96*(sd(x)/sqrt(nrow(trait.space.ot)))))
ciu.ot <- apply(trait.space.ot,2,function(x) mean(x)+(1.96*(sd(x)/sqrt(nrow(trait.space.ot)))))
ot.ci <- cbind(ott.mean,cil.ot,ciu.ot)
sig.ot <- as.data.frame(ifelse(ot.ci[,1]>cil.ot & ot.ci[,1]<ciu.ot,0,1))
colnames(sig.ot) <- "one.tree"

cil.lhi <- apply(trait.space.lh,2,function(x) mean(x)-(1.96*(sd(x)/sqrt(nrow(trait.space.lh)))))
ciu.lhi <- apply(trait.space.lh,2,function(x) mean(x)+(1.96*(sd(x)/sqrt(nrow(trait.space.lh)))))
lhi.ci <- cbind(lhit.mean,cil.lhi,ciu.lhi)
sig.lhi <- as.data.frame(ifelse(lhi.ci[,1]>cil.lhi & lhi.ci[,1]<ciu.lhi,0,1))
colnames(sig.lhi) <- "lord.howe"

sig <- cbind(sig.liz,sig.ot,sig.lhi)

######### AH! BUT THIS IS CONFIDENCE INTERVAL OF THE MEAN
## not what I want to find

summary(trait.space.lh)

liz95 <- t(as.data.frame(apply(trait.space.liz,2,function(x) quantile(x,probs=c(.025,.975)))))
ot95 <- t(as.data.frame(apply(trait.space.ot,2,function(x) quantile(x,probs=c(.025,.975)))))
lh95 <- t(as.data.frame(apply(trait.space.lh,2,function(x) quantile(x,probs=c(.025,.975)))))
liz95 <- cbind(liz95,lizt.mean)
ot95 <- cbind(ot95,ott.mean)
lh95 <- cbind(lh95,lhit.mean)

sig.liz <- as.data.frame(ifelse(liz95[,3]>liz95[,1] & liz95[,3]<liz95[,2],0,1))
colnames(sig.liz) <- "lizard"
sig.ot <- as.data.frame(ifelse(ot95[,3]>ot95[,1] & ot95[,3]<ot95[,2],0,1))
colnames(sig.ot) <- "one.tree"
sig.lhi <- as.data.frame(ifelse(lh95[,3]>lh95[,1] & lh95[,3]<lh95[,2],0,1))
colnames(sig.lhi) <- "lord.howe"
sig <- cbind(sig.liz,sig.ot,sig.lhi)
sig


##############################################
##############################################

# SK 05/08/2014
## INCIDENCE WEIGHTED RANDOM DRAWS




##############################################
##############################################

## ABUNDANCE INCIDENCE WEIGHTED RANDOM DRAWS

