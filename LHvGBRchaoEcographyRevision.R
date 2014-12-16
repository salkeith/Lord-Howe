# SK 07/12/2014

##########################################################################

## LORD HOWE VS GBR BREAK
## DOES LORD HOWE HAVE A SIGNIFICANTLY DIFFERENT TRAIT COMPOSITION FROM 
## GBR ASSEMBLAGES?

## COLLABORATORS: ANDREW BAIRD, JOSH MADIN, ERIKA WOOLSEY, MARIA BYRNE

## REVISIONS FOR ECOGRAPHY - ADD GROWTH FORM (PANDOLFI REVIEW)

#########################################################################

rm(list=ls())

library(vegan)
library(lme4)
library(visreg)
library(reshape)
library(MuMIn)
library(effects)


load("LordHoweMDSChao.RData")
traits <- read.csv("Traits_LHI_AB26Aug14.csv")
gf <- read.csv("ctdbGrowthForm/species_gf.csv")
traits <- merge(traits,gf,by="species",all.x=T)
# Dummy variables in data that I changed in csv file
# water clarity: turbid = -1; both = 0; clear = 1
# wave exposure: exposed = -1; broad = 0; protected = 1
# Which traits do we want to include?
colnames(traits)
# convert all columns to numeric format
for(i in c(2:7,11)){
   traits[,i] <- as.numeric(traits[,i])
}
str(traits)

for(i in c(8:10,12:14)){
   traits[,i] <- as.factor(traits[,i])
}
str(traits)

# set reference category for water clarity and wave exposure as 0
traits[,8] <- relevel(traits[,8],"0")
traits[,9] <- relevel(traits[,9],"0")

# merge species traits with NMDS axis scores
sptr <- merge(MDSsp,traits,by="species",all.x=T)

# add in genus
genus <- read.csv("CoralTraitsAug2014SK.csv")
genus <- genus[,c(1,7)]
sptr <- merge(sptr,genus,by="species")
sptr <- sptr[,c(1,17,2,3,4:16)]
colnames(sptr)[2] <- "genus"
save(sptr,file="SpeciesAxisScoresTraitsDec2014.RData")


###################################################
## REGRESSION DIAGNOSTICS

# normalize trait data (mean = 0, sd = 1)
t.norm <- scale(sptr[,5:10])
t.norm <- cbind(t.norm,sptr[,11:17])

# visualise distribution of values for each normalized variable
par(mfcol=c(3,2))
for(z in 1:6){
   hist(t.norm[,z], main = colnames(t.norm)[z])
}

# check collinearity of numeric traits
pairs(t.norm)
multcol <- cor(t.norm[,1:6],use="complete.obs")
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
                 "water.clarity","wave.exposure","larval.mode","sex","zoox","growth.form")
summary(d)


# RE-DO CORRELATIONS TO ALSO HAVE BINARY (POINT BISERIAL) AND 1,0,-1 AS SPEARMAN
# normalize trait data (mean = 0, sd = 1)
t.norm2 <- cbind(scale(sptr[,8]),sptr[,14]) # add lower depth back in
colnames(t.norm2) <- c("lower.depth","dev.rate")
# make factors numeric so they can be arguments in the cor function
num.wc <- as.numeric(levels(d[,10])[d[,10]])
num.we <- as.numeric(levels(d[,11])[d[,11]])
num.larval <- as.numeric(levels(d[,12])[d[,12]])
num.sex <- as.numeric(levels(d[,13])[d[,13]])
num.zoox <- as.numeric(levels(d[,14])[d[,14]])
dcor <- cbind(d[,5:9],t.norm2,num.wc,num.we,num.larval,num.sex,num.zoox)
# write the results straight to csv file
write.csv(rbind(cor(dcor,use="complete.obs",method="pearson"),cor(dcor,use="complete.obs",method="spearman")),"PearsonSpearmanDec2014.csv")


# check higher order terms for numeric traits
for(i in 5:10){
   par(mfcol=c(2,2))
   print(colnames(d)[i])
   lin <- lm(ax1~d[,i],data=d)
   plot(lin)
   quad <- lm(ax1~d[,i]+I(d[,i]^2),data=d)
   print(AICc(lin,quad))
}

# valley size as quadratic

# remove rows with NAs to ensure models all use same underlying data (no na.omit)
d <- d[-c(28,58,108,111),]

# looks good, ready to go
# save file so don't have to repeat steps above
save(d,file="DataForRegressionLHIDec2014.RData")




## LINEAR REGRESSION 

load("DataForRegressionLHIDec2014.RData")



# Full model, include genus as a fixed effect
mod <- lm(ax1~genus+valley.size+I(valley.size^2)+upp.depth+depth.range+wave.exposure+water.clarity+genus.age+sex+
             larval.mode+zoox+growth.form+wave.exposure*water.clarity+larval.mode*sex+upp.depth*depth.range,data=d,na.action=na.fail)
summary(mod)
# Full model including interactions that make biological sense and genus as a random effect
mod <- lmer(ax1~valley.size+I(valley.size^2)+upp.depth+depth.range+wave.exposure+water.clarity+genus.age+sex+
               larval.mode+zoox+wave.exposure*water.clarity+larval.mode*sex+upp.depth*depth.range+
               (1|genus),data=d,na.action=na.fail)
summary(mod)
dotplot(ranef(mod, condVar=TRUE))
# Variance partition coefficient 
# VPC = variance for random effect/(variance for random effect + 3.29)
# The 3.29 is because it uses a binomial distribution
tax <- 0.006/(0.006+3.29)
tax
# mixed model not required - very low variance accounted for and overlapping confidence intervals in dotplot

# Full model, no genus
mod <- lm(ax1~valley.size+I(valley.size^2)+upp.depth+depth.range+wave.exposure+water.clarity+genus.age+sex+
             larval.mode+zoox+growth.form+larval.mode*sex+upp.depth*depth.range,data=d,na.action=na.fail)
summary(mod)

# remove non-significant interactions and higher order terms
mod <- lm(ax1~valley.size+depth.range+upp.depth+wave.exposure+water.clarity+genus.age+
             sex+larval.mode+zoox+growth.form+upp.depth*depth.range,data=d,na.action=na.fail)
summary(mod)

# take a look at the partial coefficients
par(mfcol=c(3,3))
visreg(mod,partial=F)


#############################################
# MODEL SELECTION AND AVERAGING
# model selection by delta < 3 (Bolker 2009) with polynomials
dmod <- dredge(mod)
mod.sel <- summary(model.avg(dmod,subset=delta<3))
mod.sel

save(mod.sel,file="factorModelAveragedResultDec2014.RData")

write.csv(coefTable(mod.sel),"factorModelAveragedCoefsDec2014.csv")

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

write.csv(CImod,"factorModelAveragedCoefs95CI.csv")


#######################################################
# PARTIAL COEFFICIENT PLOTS for significant variables 

# PLOT INTERACTION BETWEEN DEPTH VARIABLES
# N.B. the function above doesnt work with model averaged result
best.mod <- lm(ax1~upp.depth+depth.range+upp.depth*depth.range+larval.mode,data=d,na.action=na.fail)
# plot all partial coefficients on the same screen - cool :)
z <- allEffects(mod=best.mod,default.levels=4)
pdf("factorAllEffectsBestModel.pdf")
plot(z,rug=F)
dev.off()

# BACKTRANSFORM scaled depth range
# Depth range is held constant at 4 different values through
# the allEffects function to visualise the interaction
# Get original values so can find mean and sd
orig.dr <- sptr[-c(28,58,108,111),10]
# pull levels out of allEffects result & backtransform
z[[2]][[4]][[2]][[3]] * sd(orig.dr) + mean(orig.dr)



######################################################


# TURNOVER BETWEEN SITES
# REQUESTED BY REVIEWER 1

load("dstar.RData")
dpa <- t(dstar)>0
dpa[,1] <- as.numeric(dpa[,1])
dpa <- as.data.frame(dpa)
dpa <- dpa[-120,]
rownames(dpa) <- MDSsp[,1]

write.csv(dpa,"PresAbsIslands.csv")

