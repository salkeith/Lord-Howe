# SK 07/07/2014
# SK SUBSTANTIALLY MODIFIED 06/08/2014

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

par(mfcol=c(1,1))
# plot MDS with each island a different colour
plot(MDS,display="sites",col="white",cex.lab=2,cex.axis=1.5,font.lab=2)
points(MDS, display="sites",pch=1,col=1,select=island.short== "lizard")
points(MDS, display="sites",pch=19,col="grey",select=island.short== "one.tree")
points(MDS, display="sites",pch=19,col=1,select=island.short== "lord.howe")
# add to plot 95% confidence ellipses for each island
ordiellipse(MDS,island.short,kind="sd",conf=0.95,lty=2)
legend("topright",c("Lizard","One Tree","Lord Howe"),pch=c(1,19,19),col=c(1,"grey",1),cex=2.5)
# plot MDS with each habitat a different colour
plot(MDS,display="sites",col="white",cex.lab=1,font.lab=2,main="MDS transects")
points(MDS, display="sites",pch=20,col=1,select=habitat.short== "lagoon")
points(MDS, display="sites",pch=20,col=2,select=habitat.short== "crest")
# Lord Howe separates out on axis 1, lagoon/crest on axis 2
# add to plot 95% confidence ellipses for each habitat
ordiellipse(MDS,habitat.short,kind="sd",conf=0.90,lty=2)
legend("topright",c("lagoon","crest"),pch=20,col=c(1,2))

# isolate MDS axis scores for use in GLMM
MDSsp <- MDS$species
MDSsp <- as.data.frame(MDSsp[,1:2])
MDSsp <- cbind(rownames(spbysite2),MDSsp)
colnames(MDSsp) <- c("species","ax1","ax2")
rownames(MDSsp) <- NULL


##############################################################################################
## TRY REMOVING OUTLIER site 140
spbysite2 <- spbysite[,-140]
island.short2 <- island.short[-140]
region2 <- region[-140]
habitat.short2 <- habitat.short[-140]
MDS <- metaMDS(t(spbysite2),trymax=50)
plot(MDS)
par(mfcol=c(1,1))
# plot MDS with each island a different colour
plot(MDS,display="sites",col="white",cex.lab=2,cex.axis=1.5,font.lab=2)
points(MDS, display="sites",pch=1,col=1,select=island.short2== "lizard")
points(MDS, display="sites",pch=19,col="grey",select=island.short2== "one.tree")
points(MDS, display="sites",pch=19,col=1,select=island.short2== "lord.howe")
# add to plot 95% confidence ellipses for each island
ordiellipse(MDS,island.short2,kind="sd",conf=0.95,lty=2)
legend("topright",c("Lizard","One Tree","Lord Howe"),pch=c(1,19,19),col=c(1,"grey",1),cex=2.5)
anosim(t(spbysite2),island.short2,permutations=1000,distance="bray")
anosim(t(spbysite2),region2,permutations=1000,distance="bray")
anosim(t(spbysite2),habitat.short2,permutations=1000,distance="bray")
# ANOSIM between each island pair
island.LizOT <- rep(c("lizard","one.tree"),c(73,58))
anosim(t(spbysite2[,1:131]),island.LizOT,permutations=1000,distance="bray")
island.LizLH <- rep(c("lizard","lord.howe"),c(73,71))
anosim(t(spbysite2[,c(1:73,132:202)]),island.LizLH,permutations=1000,distance="bray")
island.LHOT <- rep(c("one.tree","lord.howe"),c(58,71))
anosim(t(spbysite2[,74:202]),island.LHOT,permutations=1000,distance="bray")

rownames(spbysite2) <- rel.cover.mat[,1]
##############################################################################################
# doesn't affect ANOSIM or model result so leave in. Ellipses capture central region.
##############################################################################################


#############################################
#############################################

# SK 06/08/2014
## PREP DATA FOR MIXED EFFECTS MODEL
## BASED ON NMDS AXIS SCORES

traits <- read.csv("CoralTraitsAug2014SK.csv")
# Dummy variables in data that I changed in csv file
# water clarity: turbid = -1; both = 0; clear = 1
# wave exposure: exposed = -1; broad = 0; protected = 1

# Which traits do we want to include?
colnames(traits)
traits <- traits[,c(1,7,13,14:17,20,21,23,30:34)]
# convert all columns to numeric format
for(i in 3:ncol(traits)){
   traits[,i] <- as.numeric(traits[,i])
}
str(traits)

# merge species traits with NMDS axis scores
sptr <- merge(MDSsp,traits,by="species")
save(sptr,file="SpeciesAxisScoresTraits.RData")




####################################################
####################################################

## REGRESSION MODEL
## HOW WELL CAN TRAITS PREDICT THE NMDS AXIS SCORE?

####################################################
####################################################

load("SpeciesAxisScoresTraits.RData")

###################################################
## REGRESSION DIAGNOSTICS

# normalize trait data (mean = 0, sd = 1)
t.norm <- scale(sptr[,5:ncol(sptr)])

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
# remove egg size class - doesn't fit as clearly into our hypotheses
t.norm <- t.norm[,-c(3,13)]

# put data together into one data frame ready for regression
d <- cbind(sptr[,1:4],t.norm)
colnames(d) <- c(colnames(d)[1:3],"genus","valley.size","colony.size","upp.depth","depth.range",
                 "wave.exposure","water.clarity","genus.age","sex","larval.mode","zoox","egg.diam")
summary(d)
# colony size is NA for 12 species. Remove trait.
d <- d[,-6]

# check higher order terms for traits
for(i in 5:length(d)){
   par(mfcol=c(2,2))
   print(colnames(d)[i])
   lin <- lm(ax1~d[,i],data=d)
   plot(lin)
   quad <- lm(ax1~d[,i]+I(d[,i]^2),data=d)
   cub <- lm(ax1~d[,i]+I(d[,i]^2)+I(d[,i]^3),data=d)
   print(AICc(lin,quad,cub))
}

# valley size is best as cubic
# upper depth is best as quadratic

# Regression diagnostics show approximately normal distribution of residuals
# row 46 is a strong outlier - sqrt(standardised residuals) >3
d[46,]  # large axis one score
# remove it for now but note that this may skew result!
d <- d[-46,]

# remove rows with NAs to ensure models all use same underlying data (no na.omit)
d <- d[-c(59,72,74,76,108,111,113,119),]

# with outlier removed, upper depth can now be linear
# valley size still cubic

# recheck regression diagnostics
for(i in 5:length(d)){
   par(mfcol=c(2,2))
   print(colnames(d)[i])
   lin <- lm(ax1~d[,i],data=d)
   plot(lin)
}
par(mfcol=c(2,2))
cub <- lm(ax1~d[,4]+I(d[,4]^2)+I(d[,4]^3),data=d)
plot(quad)

# looks good, ready to go
# save file so don't have to repeat steps above
save(d,file="DataForRegressionLHI.RData")


###################################################
## LINEAR REGRESSION 

load(d,file="DataForRegressionLHI.RData")

# Full model including interactions that make biological sense and genus as a random effect
mod <- lmer(ax1~valley.size+I(valley.size^2)+I(valley.size^3)+upp.depth+depth.range+wave.exposure+water.clarity+genus.age+sex+
              larval.mode+zoox+egg.diam+upp.depth*depth.range+wave.exposure*water.clarity+larval.mode*sex+larval.mode*egg.diam+
              sex*egg.diam+(1|genus),data=d,na.action=na.fail)
summary(mod)
dotplot(ranef(mod, postVar=TRUE))
# mixed model not required - very low variance accounted for and overlapping confidence intervals in dotplot

# remove non-significant interactions
mod <- lm(ax1~valley.size+I(valley.size^2)+I(valley.size^3)+upp.depth+depth.range+wave.exposure+water.clarity+
             genus.age+sex+larval.mode+zoox+egg.diam+upp.depth*depth.range,
             data=d,na.action=na.fail)
summary(mod)

# take a look at the partial coefficients
par(mfcol=c(3,4))
visreg(mod,partial=F)


#############################################
# MODEL SELECTION AND AVERAGING
# expression so that cubic term cannot be included unless linear & quadratic term present
# and interaction terms cannot be included without the main effects
msubset <- expression(c(valley.size|!`I(valley.size^2)`),(`I(valley.size^2)`|!`I(valley.size^3)`),
                      (upp.depth|!depth.range*upp.depth),(depth.range|!depth.range*upp.depth))
# model selection by delta < 3 (Bolker 2009) with polynomials
dmod <- dredge(mod,subset=msubset)
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
# manual plotting required, visreg doesn't work on model-averged object

pdf("PartialCoefsLordHoweTraits.pdf")

par(mfcol=c(2,2))

# LARVAL MODE
x <- seq(min(d$larval.mode),max(d$larval.mode),length=100)
larval.pred <- with(d,(predict(mod.sel,se.fit=T,newdata=data.frame(larval.mode=x,
               depth.range=mean(depth.range),upp.depth=mean(upp.depth),wave.exposure=mean(wave.exposure),
               water.clarity=mean(water.clarity),valley.size=mean(valley.size),genus.age=mean(genus.age),
               sex=mean(sex),egg.diam=mean(egg.diam),zoox=mean(zoox)))))
plot(ax1~larval.mode,data=d,col="lightgray",pch=16,ylab="NMDS Axis 1",
     xlab="Larval mode",cex.lab=1.5)
points(larval.pred$fit~x,type="l",col=1,lwd=2)
lines(larval.pred$fit+larval.pred$se.fit~x,lty=2,type="l")
lines(larval.pred$fit-larval.pred$se.fit~x,lty=2,type="l")

# WAVE EXPOSURE
x <- seq(min(d$wave.exposure),max(d$wave.exposure),length=100)
wave.pred <- with(d,(predict(mod.sel,se.fit=T,newdata=data.frame(wave.exposure=x,
              depth.range=mean(depth.range),upp.depth=mean(upp.depth),larval.mode=mean(larval.mode),
              water.clarity=mean(water.clarity),valley.size=mean(valley.size),genus.age=mean(genus.age),
              sex=mean(sex),egg.diam=mean(egg.diam),zoox=mean(zoox)))))
plot(ax1~wave.exposure,data=d,col="lightgray",pch=16,ylab="NMDS Axis 1",
     xlab="Wave exposure",cex.lab=1.5)
points(wave.pred$fit~x,type="l",col=1,lwd=2)
lines(wave.pred$fit+wave.pred$se.fit~x,lty=2,type="l")
lines(wave.pred$fit-wave.pred$se.fit~x,lty=2,type="l")

# VALLEY SIZE
x <- seq(min(d$valley.size),max(d$valley.size),length=100)
valley.pred <- with(d,(predict(mod.sel,se.fit=T,newdata=data.frame(valley.size=x,
             depth.range=mean(depth.range),upp.depth=mean(upp.depth),larval.mode=mean(larval.mode),
             water.clarity=mean(water.clarity),wave.exposure=mean(wave.exposure),genus.age=mean(genus.age),
             sex=mean(sex),egg.diam=mean(egg.diam),zoox=mean(zoox)))))
plot(ax1~valley.size,data=d,col="lightgray",pch=16,ylab="NMDS Axis 1",
     xlab="Valley size",cex.lab=1.5)
points(valley.pred$fit~x,type="l",col=1,lwd=2)
lines(valley.pred$fit+valley.pred$se.fit~x,lty=2,type="l")
lines(valley.pred$fit-valley.pred$se.fit~x,lty=2,type="l")

# plot interaction 
# plot depth range as binned variable, upper depth as continuous
# UPPER DEPTH * DEPTH RANGE
par(mfrow=c(2,2))
drbins <- c(-1.5,0,1.5,3)
for(i in 1:4){
   drbin <- drbins[i]
   x <- seq(min(d$upp.depth),max(d$upp.depth),length=100)
   udepth.pred <- with(d,(predict(mod.sel,se.fit=T,newdata=data.frame(upp.depth=x,
                        depth.range=drbin,valley.size=mean(valley.size),larval.mode=mean(larval.mode),
                        water.clarity=mean(water.clarity),wave.exposure=mean(wave.exposure),genus.age=mean(genus.age),
                        sex=mean(sex),egg.diam=mean(egg.diam),zoox=mean(zoox)))))
   plot(ax1~upp.depth,data=d,col="lightgray",pch=16,ylab="NMDS Axis 1",
        xlab="Upper depth",cex.lab=1.5,main=paste("Depth range bin ", drbin,sep=""))
   points(udepth.pred$fit~x,type="l",col=1,lwd=2)
   lines(udepth.pred$fit+udepth.pred$se.fit~x,lty=2,type="l")
   lines(udepth.pred$fit-udepth.pred$se.fit~x,lty=2,type="l")
}

dev.off()



############################################################
## PLOT PARTIAL COEFFICIENTS WITH 95% CONFIDENCE INTERVALS

# LARVAL MODE
x<-seq(min(d$larval.mode),max(d$larval.mode),length=100)
plot(d$larval.mode,d$ax1,pch=20,cex=0.5,col="grey") 
y3 <- CImod[1,1] + CImod[3,1]*x 
lines(x,y3,lwd=2) 
# Confidence intervals 
y3 <- CImod[1,2] + CImod[3,2]*x 
y4 <- CImod[1,3] + CImod[3,3]*x 
# create polygon so the CI region is filled in
xx <- cbind(x,rev(x))
yy <- cbind(y3,rev(y4))
# for rgb the last number is transparency
polygon(xx, yy, col=rgb(0.05,0.05,0.05,0.2),border=NA)
  
# WAVE EXPOSURE
x<-seq(min(d$wave.exposure),max(d$wave.exposure),length=100)
plot(d$wave.exposure,d$ax1,pch=20,cex=0.5,col="grey") 
y3 <- CImod[1,1] + CImod[8,1]*x 
lines(x,y3,lwd=2) 
# Confidence intervals 
y3 <- CImod[1,2] + CImod[8,2]*x 
y4 <- CImod[1,3] + CImod[8,3]*x 
# create polygon so the CI region is filled in
xx <- cbind(x,rev(x))
yy <- cbind(y3,rev(y4))
# for rgb the last number is transparency
polygon(xx, yy, col=rgb(0.05,0.05,0.05,0.2),border=NA)

# VALLEY SIZE
x<-seq(min(d$valley.size),max(d$valley.size),length=100)
plot(d$valley.size,d$ax1,pch=20,cex=0.5,col="grey") 
y3 <- CImod[1,1] + CImod[5,1]*x + CImod[6,1]*x^2 + CImod[7,1]*x^3 
lines(x,y3,lwd=2) 
# Confidence intervals 
y3 <- CImod[1,2] + CImod[5,2]*x + CImod[6,2]*x^2 + CImod[7,2]*x^3  
y4 <- CImod[1,3] + CImod[5,3]*x + CImod[6,3]*x^2 + CImod[7,3]*x^3 
# create polygon so the CI region is filled in
xx <- cbind(x,rev(x))
yy <- cbind(y3,rev(y4))
# for rgb the last number is transparency
polygon(xx, yy, col=rgb(0.05,0.05,0.05,0.2),border=NA)

## WTF?!



#################################################################
#################################################################

## I'VE READ SOME STUFF THAT SAYS DON'T USE CUBIC SO TRY WITHOUT

#################################################################
#################################################################

# SK 07/08/2014

# linear has lower AICc than quadratic for valley size
mod <- lm(ax1~valley.size+upp.depth+depth.range+wave.exposure+water.clarity+
             genus.age+sex+larval.mode+zoox+egg.diam+upp.depth*depth.range,
          data=d,na.action=na.fail)
summary(mod)

# take a look at the partial coefficients
par(mfcol=c(3,4))
visreg(mod,partial=F)

# MODEL SELECTION AND AVERAGING
# expression so that cubic term cannot be included unless linear & quadratic term present
# and interaction terms cannot be included without the main effects
msubset <- expression(c(upp.depth|!depth.range*upp.depth),(depth.range|!depth.range*upp.depth))
# model selection by delta < 3 (Bolker 2009) with polynomials
dmod <- dredge(mod,subset=msubset)
mod.sel <- summary(model.avg(get.models(subset(dmod,delta < 3))))
mod.sel

best.mod <- lm(ax1~upp.depth+depth.range+wave.exposure+water.clarity+
                  larval.mode+upp.depth*depth.range,data=d,na.action=na.fail)
summary(best.mod)

# PLOT INTERACTION BETWEEN DEPTH VARIABLES
# N.B. the function above doesnt work with model averaged result
plot(effect(term="upp.depth:depth.range",mod=best.mod,default.levels=5),multiline=TRUE)
z <- effect(term="upp.depth:depth.range",mod=best.mod,default.levels=5)
# gives coefficients & CIs for each depth bin
summary(z)  
# plot all partial coefficients on the same screen - cool :)
z <- allEffects(mod=best.mod,default.levels=4)
pdf("AllEffectsBestModel.pdf")
plot(z)
dev.off()


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

write.csv(CImod,"ModelAveragedCoefs95CIValleySizeLinear.csv")


pdf("PartialCoefsLordHoweTraitsCI.pdf")

# LARVAL MODE
x<-seq(min(d$larval.mode),max(d$larval.mode),length=100)
plot(d$larval.mode,d$ax1,pch=20,cex=0.5,col="grey") 
y3 <- CImod[1,1] + CImod[3,1]*x 
lines(x,y3,lwd=2) 
# Confidence intervals 
y3 <- CImod[1,2] + CImod[3,2]*x 
y4 <- CImod[1,3] + CImod[3,3]*x 
# create polygon so the CI region is filled in
xx <- cbind(x,rev(x))
yy <- cbind(y3,rev(y4))
# for rgb the last number is transparency
polygon(xx, yy, col=rgb(0.05,0.05,0.05,0.2),border=NA)

# WAVE EXPOSURE
x<-seq(min(d$wave.exposure),max(d$wave.exposure),length=100)
plot(d$wave.exposure,d$ax1,pch=20,cex=0.5,col="grey") 
y3 <- CImod[1,1] + CImod[6,1]*x 
lines(x,y3,lwd=2) 
# Confidence intervals 
y3 <- CImod[1,2] + CImod[6,2]*x 
y4 <- CImod[1,3] + CImod[6,3]*x 
# create polygon so the CI region is filled in
xx <- cbind(x,rev(x))
yy <- cbind(y3,rev(y4))
# for rgb the last number is transparency
polygon(xx, yy, col=rgb(0.05,0.05,0.05,0.2),border=NA)

# WATER CLARITY
x<-seq(min(d$water.clarity),max(d$water.clarity),length=100)
plot(d$water.clarity,d$ax1,pch=20,cex=0.5,col="grey") 
y3 <- CImod[1,1] + CImod[5,1]*x  
lines(x,y3,lwd=2) 
# Confidence intervals 
y3 <- CImod[1,2] + CImod[5,2]*x  
y4 <- CImod[1,3] + CImod[5,3]*x 
# create polygon so the CI region is filled in
xx <- cbind(x,rev(x))
yy <- cbind(y3,rev(y4))
# for rgb the last number is transparency
polygon(xx, yy, col=rgb(0.05,0.05,0.05,0.2),border=NA)


dev.off()



#######################################################
# PARTIAL COEFFICIENT PLOTS for significant variables 
# manual plotting required, visreg doesn't work on model-averged object

pdf("PartialCoefsLordHoweTraitsSE.pdf")

par(mfcol=c(2,2))

# LARVAL MODE
x <- seq(min(d$larval.mode),max(d$larval.mode),length=100)
larval.pred <- with(d,(predict(mod.sel,se.fit=T,newdata=data.frame(larval.mode=x,
                        depth.range=mean(depth.range),upp.depth=mean(upp.depth),wave.exposure=mean(wave.exposure),
                        water.clarity=mean(water.clarity),valley.size=mean(valley.size),genus.age=mean(genus.age),
                        sex=mean(sex),egg.diam=mean(egg.diam),zoox=mean(zoox)))))
plot(ax1~larval.mode,data=d,col="lightgray",pch=16,ylab="NMDS Axis 1",
     xlab="Larval mode",cex.lab=1.5)
points(larval.pred$fit~x,type="l",col=1,lwd=2)
lines(larval.pred$fit+larval.pred$se.fit~x,lty=2,type="l")
lines(larval.pred$fit-larval.pred$se.fit~x,lty=2,type="l")

# WAVE EXPOSURE
x <- seq(min(d$wave.exposure),max(d$wave.exposure),length=100)
wave.pred <- with(d,(predict(mod.sel,se.fit=T,newdata=data.frame(wave.exposure=x,
                        depth.range=mean(depth.range),upp.depth=mean(upp.depth),larval.mode=mean(larval.mode),
                        water.clarity=mean(water.clarity),valley.size=mean(valley.size),genus.age=mean(genus.age),
                        sex=mean(sex),egg.diam=mean(egg.diam),zoox=mean(zoox)))))
plot(ax1~wave.exposure,data=d,col="lightgray",pch=16,ylab="NMDS Axis 1",
     xlab="Wave exposure",cex.lab=1.5)
points(wave.pred$fit~x,type="l",col=1,lwd=2)
lines(wave.pred$fit+wave.pred$se.fit~x,lty=2,type="l")
lines(wave.pred$fit-wave.pred$se.fit~x,lty=2,type="l")

# VALLEY SIZE
x <- seq(min(d$valley.size),max(d$valley.size),length=100)
valley.pred <- with(d,(predict(mod.sel,se.fit=T,newdata=data.frame(valley.size=x,
                        depth.range=mean(depth.range),upp.depth=mean(upp.depth),larval.mode=mean(larval.mode),
                        water.clarity=mean(water.clarity),wave.exposure=mean(wave.exposure),genus.age=mean(genus.age),
                        sex=mean(sex),egg.diam=mean(egg.diam),zoox=mean(zoox)))))
plot(ax1~valley.size,data=d,col="lightgray",pch=16,ylab="NMDS Axis 1",
     xlab="Water clarity",cex.lab=1.5)
points(valley.pred$fit~x,type="l",col=1,lwd=2)
lines(valley.pred$fit+valley.pred$se.fit~x,lty=2,type="l")
lines(valley.pred$fit-valley.pred$se.fit~x,lty=2,type="l")

# plot interaction 
# plot depth range as binned variable, upper depth as continuous
# UPPER DEPTH * DEPTH RANGE
par(mfrow=c(2,2))
drbins <- c(-1.5,0,1.5,3)
for(i in 1:4){
   drbin <- drbins[i]
   x <- seq(min(d$upp.depth),max(d$upp.depth),length=100)
   udepth.pred <- with(d,(predict(mod.sel,se.fit=T,newdata=data.frame(upp.depth=x,
                        depth.range=drbin,valley.size=mean(valley.size),larval.mode=mean(larval.mode),
                        water.clarity=mean(water.clarity),wave.exposure=mean(wave.exposure),genus.age=mean(genus.age),
                        sex=mean(sex),egg.diam=mean(egg.diam),zoox=mean(zoox)))))
   plot(ax1~upp.depth,data=d,col="lightgray",pch=16,ylab="NMDS Axis 1",
        xlab="Upper depth",cex.lab=1.5,main=paste("Depth range bin ", drbin,sep=""))
   points(udepth.pred$fit~x,type="l",col=1,lwd=2)
   lines(udepth.pred$fit+udepth.pred$se.fit~x,lty=2,type="l")
   lines(udepth.pred$fit-udepth.pred$se.fit~x,lty=2,type="l")
}

dev.off()


#############################
## BACKTRANSFORM

tmean <- apply(sptr[,5:ncol(sptr)],2,function(x) mean(x,na.rm=T))
tsd <- apply(sptr[,5:ncol(sptr)],2,function(x) sd(x,na.rm=T))

# wave exposure
we <- rep()
x <- -1*tsd[6]
x+tmean[6]


##################################################################################
##################################################################################
## EXTRA CRAP THAT TURNED OUT NOT TO BE SO HELPFUL

## TRY OUT A GAM
mod.gam <- gam(ax1~valley.size+upp.depth+depth.range+wave.exposure+water.clarity+
             genus.age+sex+larval.mode+zoox+egg.diam+upp.depth*depth.range,
          data=d,family="gaussian",na.action=na.fail)
summary(mod.gam)
1-pchisq(mod.gam$deviance,mod.gam$df.resid)


## EXPLORE PREDICT 95% CI ERROR MESSAGE WITH MORE SIMPLE MODEL
mod.test <- lm(ax1~valley.size+upp.depth,data=d,na.action=na.fail)
# model selection by delta < 3 (Bolker 2009) with polynomials
dmod <- dredge(mod.test)
mod.sel.test <- summary(model.avg(get.models(subset(dmod,delta < 3))))
mod.sel.test
x <- seq(min(d$valley.size),max(d$valley.size),length=100)
valley.pred <- with(d,(predict(mod.sel.test,interval="c",level=0.95,se.fit=T,newdata=data.frame(valley.size=x,
                        upp.depth=mean(upp.depth)))))




