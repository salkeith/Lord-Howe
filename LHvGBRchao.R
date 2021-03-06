# SK 04/09/2014

##########################################################################

## LORD HOWE VS GBR BREAK
## DOES LORD HOWE HAVE A SIGNIFICANTLY DIFFERENT TRAIT COMPOSITION FROM 
## GBR ASSEMBLAGES?

## COLLABORATORS: ANDREW BAIRD, ERIKA WOOLSEY, MARIA BYRNE

#########################################################################

library(vegan)
library(lme4)
library(visreg)
library(reshape)
library(MuMIn)
library(effects)

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


d2 <- melt(d)
d3 <- cbind(rep(rownames(d),ncol(d)),d2)
head(d3)
colnames(d3) <- c("species","transect","cover.cm")
save(d3,file="LongFormatCoverData.RData")
table(d3[,2])

# add other columns for site data
island <- rep(c("lizard","one.tree","lord.howe"),c(120*73,120*58,120*72))
habitat <- rep(c("lagoon","crest","lagoon","crest","crest","lagoon"),c(120*36,120*37,120*36,120*22,120*36,120*36))
site <- rep(c("horseshoe","lagoon.reef","north.of.trimodal","lizard.head","north.point","trimodal","second.lagoon",
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

levels(rel.cover.mat[,1])[9] <- "Isopora cuneata"
levels(rel.cover.mat[,1])[28] <- "Isopora palifera"
levels(cover.data$species)[9] <- "Isopora cuneata"
levels(cover.data$species)[28] <- "Isopora palifera"

save(rel.cover.mat,file="WideFormatRelCoverData.RData")
save(cover.data,file="LongFormatCoverData.RData")

## PHEW!


############################################################
############################################################

## LOAD PREPARED DATA AND RUN NMDS & ANOSIM

load("LongFormatCoverData.RData")
load("WideFormatRelCoverData.RData")
head(cover.data)

# merge Porites annae & P. rus because Id issues
Porites.massive <- rel.cover.mat[105,2:nrow(rel.cover.mat)] + rel.cover.mat[113,2:nrow(rel.cover.mat)]
# shuffle rows so still alphabetical
rel.cover.mat[113,] <- rel.cover.mat[112,]
rel.cover.mat[112,] <- Porites.massive
# add name in as a new factor
rel.cover.mat[,1] <- factor(rel.cover.mat[,1], levels=c(levels(rel.cover.mat[,1]), "Porites massive"))
rel.cover.mat[112,1] <- "Porites massive"
rel.cover.mat <- rel.cover.mat[-105,]


sitebysp <- t(rel.cover.mat[,-1])
colnames(sitebysp) <- rel.cover.mat[,1]

# try as integer data by multiplying by 100 (from proportion to percentage)
# so can use Chao distance metric, which is good for unequal sampling
sitebysp2 <- as.integer(round(sitebysp*100,0))
sitebysp3 <- matrix(sitebysp2,nrow=203,ncol=119,byrow=F)
colnames(sitebysp3) <- colnames(sitebysp)

# create vectors with group labels
island.short <- rep(c("lizard","one.tree","lord.howe"),c(73,58,72))
region <- rep(c("GBR","lord.howe"),c(73+58,72))
habitat.short <- rep(c("lagoon","crest","lagoon","crest","crest","lagoon"),c(36,37,36,22,36,36))

as.is <- anosim(sitebysp3,island.short,permutations=1000,distance="chao")
as.r <- anosim(sitebysp3,region,permutations=1000,distance="chao")
as.h <- anosim(sitebysp3,habitat.short,permutations=1000,distance="chao")
as.is; as.r; as.h
par(mfcol=c(3,2))
hist(as.is$perm,xlim=c(-0.1,0.5),ylim=c(0,300),main="(a) Islands",xlab="Global R",cex.lab=1.5,cex.main=1.5,cex.axis=1.5); abline(v=as.is$statistic,lty=2)
hist(as.r$perm,xlim=c(-0.1,0.5),ylim=c(0,300),main="(b) Region",xlab="Global R",cex.lab=1.5,cex.main=1.5,cex.axis=1.5); abline(v=as.r$statistic,lty=2)
hist(as.h$perm,xlim=c(-0.1,0.5),ylim=c(0,300),main="(c) Habitat",xlab="Global R",cex.lab=1.5,cex.main=1.5,cex.axis=1.5); abline(v=as.h$statistic,lty=2)

# ANOSIM between each island pair
island.LizOT <- rep(c("lizard","one.tree"),c(73,58))
as.lizot <- anosim(sitebysp[1:131,],island.LizOT,permutations=1000,distance="chao")
island.LizLH <- rep(c("lizard","lord.howe"),c(73,72))
as.lizlh <- anosim(sitebysp[c(1:73,132:203),],island.LizLH,permutations=1000,distance="chao")
island.LHOT <- rep(c("one.tree","lord.howe"),c(58,72))
as.lhot <- anosim(sitebysp[74:203,],island.LHOT,permutations=1000,distance="chao")
as.lizot; as.lizlh; as.lhot
hist(as.lizot$perm,xlim=c(-0.1,0.5),ylim=c(0,300),main="(d) Lizard - One Tree",xlab="Global R",cex.lab=1.5,cex.main=1.5,cex.axis=1.5); abline(v=as.lizot$statistic,lty=2)
hist(as.lizlh$perm,xlim=c(-0.1,0.5),ylim=c(0,300),main="(e) Lizard - Lord Howe",xlab="Global R",cex.lab=1.5,cex.main=1.5,cex.axis=1.5); abline(v=as.lizlh$statistic,lty=2)
hist(as.lhot$perm,xlim=c(-0.1,0.5),ylim=c(0,300),main="(f) One Tree - Lord Howe",xlab="Global R",cex.lab=1.5,cex.main=1.5,cex.axis=1.5); abline(v=as.lhot$statistic,lty=2)


## MDS
MDS <- metaMDS(sitebysp3,trymax=50,distance="chao")
plot(MDS)
MDS$points

par(mfrow=c(2,2),mar=c(5,5,4,2))
# plot MDS with each island a different colour
plot(MDS,display="sites",col="white",cex.lab=1.5,cex.axis=1.5,font.lab=1.5,ylab="NMDS Axis 2",xlab="NMDS Axis 1",
     main="(a) NMDS Site Scores by Island",cex.main=1.5)
points(MDS, display="sites",pch=1,col=1,select=island.short== "lizard")
points(MDS, display="sites",pch=19,col="grey",select=island.short== "one.tree")
points(MDS, display="sites",pch=19,col=1,select=island.short== "lord.howe")
# add to plot 90% confidence ellipses for each island
ordiellipse(MDS,island.short,kind="sd",conf=0.90,lty=2)
legend("topright",c("Lizard","One Tree","Lord Howe"),pch=c(1,19,19),col=c(1,"grey",1),cex=1.5,bty="n")
# plot MDS with each habitat a different colour
plot(MDS,display="sites",col="white",cex.lab=1,font.lab=2,ylab="NMDS Axis 2",xlab="NMDS Axis 1",cex.lab=1.5,
     cex.axis=1.5,font.lab=1.5,main="(b) NMDS Site Scores by Habitat",cex.main=1.5)
points(MDS, display="sites",pch=1,col=1,select=habitat.short== "lagoon")
points(MDS, display="sites",pch=19,col=1,select=habitat.short== "crest")
# Lord Howe separates out on axis 1, lagoon/crest on axis 2
# add to plot 95% confidence ellipses for each habitat
ordiellipse(MDS,habitat.short,kind="sd",conf=0.90,lty=2)
legend("topright",c("Lagoon","Crest"),pch=c(1,19),col=c(1,1),cex=1.5,bty="n")

# isolate MDS axis scores for use in GLMM
MDSsp <- MDS$species
MDSsp <- as.data.frame(MDSsp[,1:2])
MDSsp <- cbind(colnames(sitebysp3),MDSsp)
colnames(MDSsp) <- c("species","ax1","ax2")
rownames(MDSsp) <- NULL

save(MDSsp,MDS,file="LordHoweMDSChao.RData")



#############################################
#############################################

## PREP DATA FOR MIXED EFFECTS MODEL
## BASED ON NMDS AXIS SCORES

load("LordHoweMDSChao.RData")
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

# add in genus
genus <- read.csv("CoralTraitsAug2014SK.csv")
genus <- genus[,c(1,7)]
sptr <- merge(sptr,genus,by="species")
sptr <- sptr[,c(1,16,2,3,4:15)]
colnames(sptr)[2] <- "genus"
save(sptr,file="SpeciesAxisScoresTraits.RData")




####################################################
####################################################

## REGRESSION MODEL FOR NMDS AXIS 1
## HOW WELL CAN TRAITS PREDICT THE NMDS AXIS SCORE?

####################################################
####################################################

load("SpeciesAxisScoresTraits.RData")

###################################################
## REGRESSION DIAGNOSTICS

# normalize trait data (mean = 0, sd = 1)
t.norm <- scale(sptr[,5:10])
t.norm <- cbind(t.norm,sptr[,11:16])

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
# larval development mode & egg size class
# ACTION:
# remove lower depth because less information than depth range
# remove development rate - doesn't fit as clearly into our hypotheses
t.norm <- t.norm[,-c(4,10)]

# put data together into one data frame ready for regression
d <- cbind(sptr[,1:4],t.norm)
colnames(d) <- c(colnames(d)[1:4],"colony.size","valley.size","genus.age","upp.depth","depth.range",
                 "water.clarity","wave.exposure","larval.mode","sex","zoox")
summary(d)

# check higher order terms for traits
for(i in 5:length(d)){
   par(mfcol=c(2,2))
   print(colnames(d)[i])
   lin <- lm(ax1~d[,i],data=d)
   plot(lin)
   quad <- lm(ax1~d[,i]+I(d[,i]^2),data=d)
   print(AICc(lin,quad))
}

# valley size as quadratic

# remove rows with NAs to ensure models all use same underlying data (no na.omit)
d <- d[-c(58,108,111),]

# looks good, ready to go
# save file so don't have to repeat steps above
save(d,file="DataForRegressionLHI.RData")


###################################################
## LINEAR REGRESSION 

load("DataForRegressionLHI.RData")

# Full model, include genus as a fixed effect
mod <- lm(ax1~genus+valley.size+I(valley.size^2)+upp.depth+depth.range+wave.exposure+water.clarity+genus.age+sex+
               larval.mode+zoox+upp.depth*depth.range+wave.exposure*water.clarity+larval.mode*sex,data=d,na.action=na.fail)
summary(mod)
# Full model including interactions that make biological sense and genus as a random effect
mod <- lmer(ax1~valley.size+I(valley.size^2)+upp.depth+depth.range+wave.exposure+water.clarity+genus.age+sex+
              larval.mode+zoox+upp.depth*depth.range+wave.exposure*water.clarity+larval.mode*sex+
              (1|genus),data=d,na.action=na.fail)
summary(mod)
dotplot(ranef(mod, condVar=TRUE))
# Variance partition coefficient 
# VPC = variance for random effect/(variance for random effect + 3.29)
# The 3.29 is because it uses a binomial distribution
tax <- 0.008/(0.008+3.29)
tax
# mixed model not required - very low variance accounted for and overlapping confidence intervals in dotplot

# Full model, no genus
mod <- lm(ax1~valley.size+I(valley.size^2)+upp.depth+depth.range+wave.exposure+water.clarity+genus.age+sex+
             larval.mode+zoox+upp.depth*depth.range+wave.exposure*water.clarity+larval.mode*sex,data=d,na.action=na.fail)
summary(mod)

# remove non-significant interactions and higher order terms
mod <- lm(ax1~valley.size+upp.depth+depth.range+wave.exposure+water.clarity+
             genus.age+sex+larval.mode+zoox+upp.depth*depth.range,data=d,na.action=na.fail)
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

save(mod.sel,file="ModelAveragedResult.RData")

write.csv(coefTable(mod.sel),"ModelAveragedCoefs.csv")

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

write.csv(CImod,"ModelAveragedCoefs95CI.csv")


#######################################################
# PARTIAL COEFFICIENT PLOTS for significant variables 

# PLOT INTERACTION BETWEEN DEPTH VARIABLES
# N.B. the function above doesnt work with model averaged result
best.mod <- lm(ax1~upp.depth+depth.range+wave.exposure+larval.mode+upp.depth*depth.range,data=d,na.action=na.fail)
# plot all partial coefficients on the same screen - cool :)
z <- allEffects(mod=best.mod,default.levels=4)
pdf("AllEffectsBestModel.pdf")
plot(z,rug=F)
dev.off()

#######################################################





####################################################
####################################################

## REGRESSION MODEL FOR NMDS AXIS 2

####################################################
####################################################

load("SpeciesAxisScoresTraits.RData")

###################################################
## REGRESSION DIAGNOSTICS

# check higher order terms for traits
for(i in 5:length(d)){
   par(mfcol=c(2,2))
   print(colnames(d)[i])
   lin <- lm(ax2~d[,i],data=d)
   plot(lin)
   quad <- lm(ax2~d[,i]+I(d[,i]^2),data=d)
   print(AICc(lin,quad))
}

# remove wave exposure and water clarity because they are based
# on the habitat a coral is found in and therefore would be circular

d2 <- d[,-c(10,11)]

# looks good, ready to go
# save file so don't have to repeat steps above
save(d2,file="DataForRegressionLHIAx2.RData")


###################################################
## LINEAR REGRESSION 

load("DataForRegressionLHIAx2.RData")

# Full model, include genus as a fixed effect
mod <- lm(ax2~genus+valley.size+upp.depth+depth.range+genus.age+sex+
             larval.mode+zoox+upp.depth*depth.range+larval.mode*sex,data=d2,na.action=na.fail)
summary(mod)
# Full model including interactions that make biological sense and genus as a random effect
mod <- lmer(ax2~valley.size+upp.depth+depth.range+genus.age+sex+
               larval.mode+zoox+upp.depth*depth.range+larval.mode*sex+
               (1|genus),data=d2,na.action=na.fail)
summary(mod)
dotplot(ranef(mod, condVar=TRUE))
# Variance partition coefficient 
# VPC = variance for random effect/(variance for random effect + 3.29)
# The 3.29 is because it uses a binomial distribution
tax <- 0.007/(0.007+3.29)
tax
# mixed model not required - very low variance accounted for and overlapping confidence intervals in dotplot

# Full model, no genus
mod <- lm(ax2~valley.size+upp.depth+depth.range+genus.age+sex+
             larval.mode+zoox+upp.depth*depth.range+larval.mode*sex,data=d2,na.action=na.fail)
summary(mod)

# remove non-significant interactions and higher order terms
mod <- lm(ax2~valley.size+upp.depth+depth.range+
             genus.age+sex+larval.mode+zoox,data=d2,na.action=na.fail)
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

save(mod.sel,file="ModelAveragedResultAx2.RData")

write.csv(coefTable(mod.sel),"ModelAveragedCoefsAx2.csv")

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

write.csv(CImod,"ModelAveragedCoefs95CIAx2.csv")


#######################################################
# PARTIAL COEFFICIENT PLOTS for significant variables 

# PLOT INTERACTION BETWEEN DEPTH VARIABLES
# N.B. the function above doesnt work with model averaged result
best.mod <- lm(ax2~genus.age+sex,data=d2,na.action=na.fail)
# plot all partial coefficients on the same screen - cool :)
z <- allEffects(mod=best.mod,default.levels=4)
par(mfrow=c(2,2))
pdf("AllEffectsBestModelAx2.pdf")
plot(z,rug=F)
dev.off()

#######################################################







####################################################
####################################################

## SIZE DISTRIBUTIONS

####################################################
####################################################

library(dplyr)
library(tidyr)

y <- read.csv("LIT_Keithetal_2014.csv")
head(y)

# use tidyr to put all intercepts into one column
y2 <- gather(y,key=colony,value=intercept,intercept_1:intercept_18)
# merge with traits so we know if species are brooders or spawners
load("SpeciesAxisScoresTraits.RData")
bs <- sptr[,c(1,2,13)]
colnames(bs) <- c("species","genus","spawn")

y3 <- merge(y2,bs,by="species")
y3 <- y3[,c(1,4,9,10,11)]
# remove NAs using dplyr package
y4 <- filter(y3,value>0)
colnames(y4) <- c(colnames(y4)[1:2],"intercept",colnames(y4)[4:5])

# remove extra data frames
rm(y,y2,y3)

# use dplyr to subset by island
y.liz <- filter(y4,Reef=="Lizard Island")
y.ot <- filter(y4,Reef=="One Tree Island")
y.lh <- filter(y4,Reef=="Lord Howe Island")

# Plot log size distributions for each island and brooder/spawner
bandwidth <- 0.175
plot(density(log((filter(y.liz,spawn==1))$intercept),bw=bandwidth),ylim=c(0,1),main="Brooder & Spawner Size Distribution",xlab="log(Size distribution)")
lines(density(log((filter(y.liz,spawn==0))$intercept),bw=bandwidth),lty=2)
lines(density(log((filter(y.ot,spawn==1))$intercept),bw=bandwidth),lty=1,col=3)
lines(density(log((filter(y.ot,spawn==0))$intercept),bw=bandwidth),lty=2,col=3)
lines(density(log((filter(y.lh,spawn==1))$intercept),bw=bandwidth),lty=1,col=4)
lines(density(log((filter(y.lh,spawn==0))$intercept),bw=bandwidth),lty=2,col=4)
legend("topright",c("Lizard S","Lizard B","One Tree S","One Tree B","Lord Howe S","Lord Howe B"),lty=c(1,2,1,2,1,2),col=c(1,1,3,3,4,4))

# Plot size distributions (without logging) for each island and brooder/spawner
plot(density((filter(y.liz,spawn==1))$intercept,bw=2),ylim=c(0,0.08),xlim=c(0,200))
lines(density((filter(y.liz,spawn==0))$intercept,bw=2),lty=2)
lines(density((filter(y.ot,spawn==1))$intercept,bw=2),lty=1,col=3)
lines(density((filter(y.ot,spawn==0))$intercept,bw=2),lty=2,col=3)
lines(density((filter(y.lh,spawn==1))$intercept,bw=2),lty=1,col=4)
lines(density((filter(y.lh,spawn==0))$intercept,bw=2),lty=2,col=4)
legend("topright",c("Lizard S","Lizard B","One Tree S","One Tree B","Lord Howe S","Lord Howe B"),lty=c(1,2,1,2,1,2),col=c(1,1,3,3,4,4))


plot(density(log((filter(y.liz,spawn==1,genus=="Pocillopora"))$intercept)),ylim=c(0,2),main="Brooder & Spawner Size Distribution",xlab="log(Size distribution)")
lines(density(log((filter(y.liz,spawn==0,genus=="Pocillopora"))$intercept)),lty=2)
lines(density(log((filter(y.ot,spawn==1,genus=="Pocillopora"))$intercept),bw=0.15),lty=1,col=3)
lines(density(log((filter(y.ot,spawn==0,genus=="Pocillopora"))$intercept)),lty=2,col=3)
lines(density(log((filter(y.lh,spawn==1,genus=="Pocillopora"))$intercept)),lty=1,col=4)
lines(density(log((filter(y.lh,spawn==0,genus=="Pocillopora"))$intercept)),lty=2,col=4)
legend("topright",c("Lizard S","Lizard B","One Tree S","One Tree B","Lord Howe S","Lord Howe B"),lty=c(1,2,1,2,1,2),col=c(1,1,3,3,4,4))





####################################################

## SK 24/09/2014
## CREATE STAR PLOTS OF MOST ABUNDANT SPECIES

####################################################

# uses a base function in R

load("WideFormatRelCoverData.RData")
rel.cover.mat[1:10,1:10]
island <- rep(c("lizard","one.tree","lord.howe"),c(73,58,72))
habitat <- rep(c("lagoon","crest","lagoon","crest","crest","lagoon"),c(36,37,36,22,36,36))

liz.mean.rel.ab <- rowMeans(rel.cover.mat[,2:74])
ot.mean.rel.ab <- rowMeans(rel.cover.mat[,75:132])
lh.mean.rel.ab <- rowMeans(rel.cover.mat[,133:204])
dstar <- rbind(liz.mean.rel.ab,ot.mean.rel.ab,lh.mean.rel.ab)
rownames(dstar) <- c("Lizard Island","One Tree Island","Lord Howe Island")
colnames(dstar) <- rel.cover.mat[,1]
stars(dstar)

# Select 10 most abundant species to focus on and lump others under "other species"
plot(rev(sort(colSums(dstar))))  # distribution of summed mean abundance
top10 <- rev(sort(colSums(dstar)))[1:12]  # top 10 species
top10sp <- attributes(top10)$names

ds <- dstar[,which(colnames(dstar)%in%top10sp==T)] # top 10 species
dstar2 <- dstar[,which(colnames(dstar)%in%top10sp==F)] # all other species
other.species <- as.data.frame(rowSums(dstar2))
colnames(other.species) <- "Other species"
dstar.top <- cbind(ds,other.species)
stars(dstar.top)

# looks unclear. Normalise to vary between zero and 1?
dstar.top/max(dstar.top)
# still looks odd
# plot without the other species
stars(ds)


stars(ds,scale=T,key.labels=top10sp,draw.segments=T,col.segments=terrain.colors(12))
# order so that brooders are on bottom half and one colour scheme, 
# spawners on top half with a different colour scheme

# have a look at which colour is which....
plot(1:12,1:12,col=terrain.colors(12),pch=19)
# create a vector in the same order as the top 10 species with brooder=0,spawner=1
colnames(ds)
brsp <- c(0,1,1,1,0,1,1,0,0,1,0,0)
ds2 <- ds[,c(1,5,8,9,11,12,2,3,4,6,7,10)]

star.col <- terrain.colors(14)

stars(ds2,scale=F,draw.segments=T,col.segments=star.col[c(1:6,9:14)],flip.labels=T)

# scale of the bars is not very good
stars(sqrt(ds2),scale=F,draw.segments=T,col.segments=star.col[c(1:6,9:14)],flip.labels=T)
stars(sqrt(cbind(ds2,other.species)),scale=F,draw.segments=T,col.segments=star.col[c(1:6,8:14)],flip.labels=T)
# specify where the star plots should go on the page
# offset so labels dont overlap each other
par(oma=c(3,0,0,0))
stars(sqrt(cbind(ds2,other.species)),scale=F,draw.segments=T,col.segments=star.col[c(1:6,8:14)],
       location=cbind(c(0.5,2,3.5),c(3,2,1)),cex=1.5)

# try stretching the values between 0 and 1
ds3 <- cbind(ds2,other.species)
stretch <- function(x){(x-min(x))/(max(x)-min(x))}
stretch(ds3)

par(oma=c(3,0,0,0))
stars(ds3,scale=F,draw.segments=T,col.segments=star.col[c(1:6,8:14)],
      location=cbind(c(0.5,2,3.5),c(3,2,1)),cex=1.5)

ds4 <- stretch(ds2)
par(oma=c(3,0,0,0))
stars(ds4,scale=F,draw.segments=T,col.segments=star.col[c(1:6,8:14)],
      location=cbind(c(0.5,2,3.5),c(3,2,1)),cex=1.5)


# AB suggests other species should be split into brooders/spawners too
# good idea andrew!

# here goes...

load("SpeciesAxisScoresTraits.RData")
b.sp <- sptr[which(sptr[,13]==0),1]  # list of brooding species
s.sp <- sptr[which(sptr[,13]==1),1]  # list of spawning species
dstar2


