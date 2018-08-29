#List of pacakages to load for subsequent analyses (not all may be used; general set of packages I include at start of all analyses)
library(ade4)
library(vegan)
library(MASS)
library(car)
library(lme4)
library(lmerTest)
library(multcomp)
library(gclus)
library(ape)
library(rgl)
library(vegan3d)
library(labdsv)
library(agricolae)
library(fitdistrplus)
library(lavaan)
library(lavaanPlot)

#Script to run Permanova analyses credit Dr. Pedro Martinez Arbizu, taken from https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.10] <-'.'
  sig[p.adjusted <= 0.05] <-'*'
  sig[p.adjusted <= 0.01] <-'**'
  sig[p.adjusted <= 0.001] <-'***'
  sig[p.adjusted <= 0.0001] <-'****'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0.0001**** 0.001*** 0.01** 0.05* 0.1.")
  return(pairw.res)
  
} 

##############################################################################
#Loading vegetation data######################################################
##############################################################################
vegdata<-read.csv("Bakken 2016 veg cover data.csv")
names(vegdata)
vegdata<-vegdata[1:90,]#Loads extra (empty) rows in .csv file
vegdata[which(vegdata$Site=="8221"),]=NA#removing 8221, which was on broken lands

site.full<-vegdata$Site
trans.full<-vegdata$Transect
part<-vegdata$Part#This was part of a sampling method that we ultimately decided to not use
veg.A<-vegdata[1:45,4:112]
veg.B<-vegdata[46:90,4:112]
veg.both<-(veg.A+veg.B)#Creates total cover for all 15 quadrats per transect
veg.cover<-veg.both/15#generates average cover values per transect
veg.cover<-veg.cover[-(7:9),]#Deletes 8221 from dataset
trans<-trans.full[1:45]
trans<-as.vector(trans)
trans<-as.factor(trans)
trans<-trans[-(7:9)]
site<-site.full[1:45]
site<-as.vector(site)
site<-as.factor(site)
site<-site[-(7:9)]

#Separates plants into exotic or native data frames
exotics.all<-data.frame(veg.cover[,1:10],veg.cover$MEOF)
exotics.common<-data.frame(veg.cover[,1:8],veg.cover$MEOF)
natives.common<-data.frame(veg.cover[,11:43],veg.cover$ANTSP,veg.cover$DAPU5,veg.cover$FLAX,veg.cover$FLAX_SP,veg.cover$GETR,veg.cover$LYJU,veg.cover$OPPO,veg.cover$PEAR,veg.cover$PLPA,veg.cover$SORI)
natives.rare<-veg.cover[,50:96]#Includes MEOF, which is exotic
natives.rare<-natives.rare[,-22]#removes MEOF
natives.all<-cbind(natives.common,natives.rare)

natives.sum<-rowSums(natives.all)#Total native plant cover
exotics.sum<-rowSums(exotics.all)#Total exotic plant cover
veg.sum<-exotics.sum+natives.sum#Total overall plant cover

natives.H<-diversity(natives.all)#Shannon diversity
natives.S<-specnumber(natives.all)#Species richness
natives.J<-natives.H/log(natives.S)#Pielou's evenness index

exotics.H<-diversity(exotics.all)
exotics.S<-specnumber(exotics.all)
exotics.J<-exotics.H/log(exotics.S)

veg.all<-cbind(exotics.all,natives.all)
veg.H<-diversity(veg.all)
veg.S<-specnumber(veg.all)
veg.J<-veg.H/log(veg.S)

exot.prop<-exotics.sum/(exotics.sum+natives.sum)#Proportion of plant cover represented by exotic species
natives.prop<-natives.sum/(exotics.sum+natives.sum)

grass<-veg.cover$AGCR+veg.cover$ARPU+veg.cover$BOCU+veg.cover$BOGR+veg.cover$BRIN+veg.cover$BRJA+veg.cover$CALO+veg.cover$DISP+veg.cover$ELSU+veg.cover$ELTR+veg.cover$HECO+veg.cover$HOJU+veg.cover$KOMA+veg.cover$MUCU+veg.cover$NAVI+veg.cover$PASM+veg.cover$POCO+veg.cover$POPR+veg.cover$SCPA+veg.cover$SCPA_SP+veg.cover$SCSC+veg.cover$SPCR+veg.cover$UNK_G+veg.cover$KOMA_SP
forb<-veg.cover$ACMI+veg.cover$ANCA8+veg.cover$ANT_SP+veg.cover$ANTSP+veg.cover$ASPU+veg.cover$ASOB+veg.cover$CIAR+veg.cover$CIUN+veg.cover$COUM+veg.cover$ECAN+veg.cover$GETR+veg.cover$GRSQ+veg.cover$HEPA+veg.cover$HYFI+veg.cover$KOSC+veg.cover$KRLA+veg.cover$LIPU+veg.cover$FLAX_SP+veg.cover$FLAX+veg.cover$LYJU+veg.cover$OECA+veg.cover$PHHO+veg.cover$PLPA+veg.cover$POAL+veg.cover$RACO+veg.cover$SORI+veg.cover$SPCO+veg.cover$TAOF+veg.cover$TRDU+veg.cover$UNK_ASTER_1
legume<-veg.cover$DAPU5+veg.cover$MELU+veg.cover$MEOF+veg.cover$PEAR+veg.cover$UNKASTR+veg.cover$UNK_ASTRAG+veg.cover$UNK_ASTRG+veg.cover$UNK_FORB+veg.cover$UNK_OXYT+veg.cover$UNKOXY
shrub<-veg.cover$ARCA13+veg.cover$ATGA+veg.cover$SYSP+veg.cover$JUHO
subshrub<-veg.cover$ARCAC+veg.cover$ARFR+veg.cover$GUSA+veg.cover$ROAR3
sedge<-veg.cover$CADU+veg.cover$CAFI
cactus<-veg.cover$OPFR+veg.cover$OPPO

veg.guilds<-data.frame(grass,forb,legume,shrub,subshrub,sedge,cactus)#will have to use [,X] to pull these out, based on order above^

ruderals<-natives.common$GRSQ+exotics.common$MELU+exotics.common$TAOF+exotics.common$TRDU
invasives<-exotics.common$AGCR+exotics.common$BRIN+exotics.common$BRJA+exotics.common$POCO+exotics.common$POPR+exotics.common$veg.cover.MEOF
common.other.cover<-veg.cover[,44:49]
rare.other.cover<-veg.cover[,97:109]
other.cover<-data.frame(common.other.cover,rare.other.cover)
other.sum<-rowSums(other.cover)

plants.ord<-cbind(exotics.common,natives.common)#Data frame with most common species to use for ordination, permanova, indicator species etc.

reclaim.year<-c(rep(1997,3),rep(1991,3),rep(1997,3),rep(1992,3),rep(1990,3),rep(1998,3),rep(1994,3),rep(1990,3),rep(1998,3),rep(1993,3),rep(1983,3),rep(1991,3),rep(2007,3),rep(2014,3))
sample.year<-2016
since.reclaim<-sample.year-reclaim.year#Time between reclamation and sampling

###########################################################################################
#Hypothesis 1: Increase in undesirable plant cover and bare ground, decreased S on reclaims
###########################################################################################
veg.S.test<-lmer(veg.S~trans+(1|site))
summary(veg.S.test)#S lower on pad vs. away |t| = 2.68 P = 0.012
anova(veg.S.test,ddf="Kenward-Roger")#F2,26=7.28 P = 0.0031
summary(glht(veg.S.test,linfct=mcp(trans="Tukey")),test=adjusted("holm"))#50 = 150 > Pad
boxplot(veg.S~trans,main="Vegetation species richness")

natives.S.test<-lmer(natives.S~trans+(1|site))
summary(natives.S.test)#pad significantly lower |t|=3.81 P = 0.0007
anova(natives.S.test,ddf="Kenward-Roger")#F2,26=13.65 P < 0.0001
summary(glht(natives.S.test,linfct=mcp(trans="Tukey")),test=adjusted("holm"))#50 = 150 > Pad
boxplot(natives.S~trans,ylab="Native species richness")

exotics.S.test<-lmer(exotics.S~trans+(1|site))
summary(exotics.S.test)#Pad significantly higher |t|=3.59 P = 0.001
anova(exotics.S.test,ddf="Kenward-Roger")#F2,26=10.98 P = 0.0003
summary(glht(exotics.S.test,linfct=mcp(trans="Tukey")),test=adjusted("holm"))#Pad > 50 = 150
boxplot(exotics.S~trans,ylab="Exotic species richness")

nativesvtrans<-lmer(natives.sum~trans+(1|site))
summary(nativesvtrans)#150=50>Pad
anova(nativesvtrans,ddf="Kenward-Roger")#F2,26=21.41, P < 0.0001
boxplot(natives.sum~trans)
summary(glht(nativesvtrans,linfct=mcp(trans="Tukey")),test=adjusted("holm"))#Pad < 150 = 50

exot.test<-lmer(exot.prop~trans+(1|site))#Smaller portion of vegetation cover, this controls for decreases in exotics simply because less vegetation cover overall
summary(exot.test)
anova(exot.test,ddf="Kenward-Roger")#F2,26=6.04 P = 0.007
summary(glht(exot.test,linfct=mcp(trans="Tukey")),test=adjusted("holm"))#Pad > 50 = 150
boxplot(exot.prop~trans,ylab="Proportion of exotic cover")

ruderals.trans<-lmer(log1p(ruderals)~trans+(1|site))
ruderals.trans<-glmer(ruderals~trans+(1|site),family=poisson)
summary(ruderals.trans)
anova(ruderals.trans,ddf="Kenward-Roger")#F2,26 = 11.47, P = 0.0003
summary(glht(ruderals.trans,linfct=mcp(trans="Tukey")),test=adjusted("holm"))#Pad > 150 = 50

invasives.trans<-lmer(log1p(invasives)~trans+(1|site))
summary(invasives.trans)
anova(invasives.trans,ddf="Kenward-Roger")#F2,26 = 2.14, P = 0.14
summary(glht(invasives.trans,linfct=mcp(trans="Tukey")),test=adjusted("holm"))#Pad > 150 = 50
plot(fitdist((invasives+0.01),"lnorm","mme"))

grass.test<-lmer(grass~trans+(1|site))
summary(grass.test)#Significant difference on Pad (lower cover)
anova(grass.test,ddf="Kenward-Roger")#F2,26=4.37 P = 0.023
summary(glht(grass.test,linfct=mcp(trans="Tukey")),test=adjusted("holm"))#50=150>pad

forb.test<-lmer(log1p(forb)~trans+(1|site))
summary(forb.test)
anova(forb.test,ddf="Kenward-Roger")#F2,26 = 1.45 P = 0.25
summary(glht(forb.test,linfct=mcp(trans="Tukey")),test=adjusted("holm"))#50 equal to 150 but greater than pad, 150 equal to pad

########Indicator species analysis to identify individual species to test
plantgroups<-c(rep(1:3,14))#This sets up predetermined groups where 1 = 150m, 2 = 50m and 3 = reclaim
transtype<-c(rep(c(1,1,2),14))#This sets up predetermined groups where 1 = native prairie and 2 = reclaim
IS.plants<-indval(plants.ord,plantgroups)#This is for the common species only
summary(IS.plants)#Species over 0.5 (plus only 150 m species): 150: CAFI; 50: BOGR; reclaim: ELTR, AGCR, MELU, DISP

cafi.test<-lmer(plants.ord$CAFI~trans+(1|site))
summary(cafi.test)
anova(cafi.test,ddf="Kenward-Roger")#P = 0.001, F = 8.75
summary(glht(cafi.test,linfct=mcp(trans="Tukey"),test=adjusted("holm")))#50 = 150 > Pad

bogr.test<-lmer(plants.ord$BOGR~trans+(1|site))
summary(bogr.test)
anova(bogr.test,ddf="Kenward-Roger")#P = 0.0001 F = 12.93
summary(glht(bogr.test,linfct=mcp(trans="Tukey"),test=adjusted("holm")))#50 = 150 > Pad

eltr.test<-lmer(log1p(ELTR)~trans+(1|site))
summary(eltr.test)
anova(eltr.test,ddf="Kenward-Roger")#F = 8.36 P = 0.002
summary(glht(eltr.test,linfct=mcp(trans="Tukey"),test=adjusted("holm")))

agcr.test<-lmer(log1p(AGCR)~trans+(1|site))
summary(agcr.test)
anova(agcr.test,ddf="Kenward-Roger")#P = 0.0004 F = 10.56
summary(glht(agcr.test,linfct=mcp(trans="Tukey"),test=adjusted("holm")))#50 = 150 > Pad

melu.test<-lmer(log1p(plants.ord$MELU)~trans+(1|site))
summary(melu.test)
anova(melu.test,ddf="Kenward-Roger")#P = 0.001 F = 9.04
summary(glht(melu.test,linfct=mcp(trans="Tukey"),test=adjusted("holm")))#50 = 150 < Pad

disp.test<-lmer(log1p(DISP)~trans+(1|site))
summary(disp.test)
anova(disp.test,ddf="Kenward-Roger")#P = 0.001 F = 8.48
summary(glht(disp.test,linfct=mcp(trans="Tukey"),test=adjusted("holm")))#50 = 150 < Pad

################Testing non-plant cover by transect type
bare.test<-lmer(log1p(common.other.cover$BARE)~trans+(1|site))
summary(bare.test)#Bare cover higher on pads, |t|=3.4 P = 0.002; F2,26=5.5 P = 0.01
anova(bare.test,ddf="Kenward-Roger")#P = 0.007, F = 5.97

lit.test<-lmer(common.other.cover$LITTER~trans+(1|site))
summary(lit.test)
anova(lit.test,ddf="Kenward-Roger")#NS
summary(glht(lit.test,linfct=mcp(trans="Tukey"),test=adjusted("holm")))

#################NMS and PERMANOVA analyses of plant community
two.trans<-c(rep(c("native rangeland","native rangeland","reclaimed pad"),14))

#use previous.best=plant.nmds a few times to make sure lowest stress value isn't a local minimum
plant.nmds<-metaMDS(plants.ord,maxit=9999,trymax=9999,distance="bray",k=3)
plant.nmds
plant.nmds$stress
#Shepard plot
stressplot(plant.nmds,main="Shepard Plot")
#goodness of fit plot
plot(plant.nmds,main="Goodness of fit",display="sites",type="t")
points(plant.nmds,display="sites",cex=goodness(plant.nmds)*200)

####Run the pairwise.adonis function from R functions script file
pairwise.adonis(plants.ord,trans)#50 and 150 equivalent, pad different; P = 0.003

########NMS and PERMANOVA on plant guilds by transect type
plant.guilds.nms<-metaMDS(veg.guilds,maxit=9999,trymax=9999,distance="bray",k=3)
stressplot(plant.guilds.nms,main="Shepard Plot")
plot(plant.guilds.nms,main="Goodness of fit",display="sites",type="t")
points(plant.guilds.nms,display="sites",cex=goodness(plant.guilds.nms)*200)

pairwise.adonis(veg.guilds,trans)#Pads different from other transects

############################################################################
#Hypothesis 2: Change with time across transects
############################################################################
dist.regress<-c(rep(c(150,50,0),14))#Distance from reclaim in meters
CAFI<-veg.all$CAFI
BOGR<-veg.all$BOGR
ELTR<-veg.all$ELTR
MELU<-veg.all$MELU
AGCR<-veg.all$AGCR
DISP<-veg.all$DISP

veg.S.age<-lm(veg.S~dist.regress)
summary(veg.S.age)#P = 0.15, R2 = 0.05

natives.S.age<-lm(natives.S~dist.regress)
summary(natives.S.age)#P = 0.04, R2 = 0.1

exotics.S.age<-lm(exotics.S~dist.regress*since.reclaim)
summary(exotics.S.age)#P = 0.0019 R2 = 0.27 interaction significant P = 0.013
exotics.S.pad<-lm(exotics.S~since.reclaim,data=vegage.data[vegage.data$trans=="PAD",])
summary(exotics.S.pad)#NS P = 0.09 R2 = 0.22
exotics.S.50<-lm(exotics.S~since.reclaim,data=vegage.data[vegage.data$trans=="50",])
summary(exotics.S.50)#NS P = 0.997
exotics.S.150<-lm(exotics.S~since.reclaim,data=vegage.data[vegage.data$trans=="150",])
summary(exotics.S.150)#P = 0.018 R2 = 0.39

natives.age<-lm(natives.sum~+dist.regress,data=vegage.data)
summary(natives.age)#P = 0.0035 R2 = 0.19 increases away from pad
plot(natives.sum~dist.regress)
abline(natives.age)
natives.bin<-lm(natives.sum~reclaim.bin,data=pad.only)
summary(natives.bin)#NS Intercept = 25.15 Slope = -1.127

exotics.age<-lm((exot.prop)~dist.regress,data=vegage.data)
summary(exotics.age)#P = 0.03 R2 = 0.11 decrease with distance
plot(exot.prop~dist.regress)
abline(exotics.age)

ruderals.age<-lm(log1p(ruderals)~dist.regress)
summary(ruderals.age)#P = 0.004, R2 = 0.19

invasives.age<-lm(log1p(invasives)~dist.regress)
summary(invasives.age)#P = 0.1 R2 = 0.06

grass.age<-lm(grass~dist.regress)
summary(grass.age)#NS P = 0.076

forb.age<-lm(log1p(forb)~since.reclaim)
summary(forb.age)#P = 0.94, R2 = 0.0002 

cafi.age<-lm(plants.ord$CAFI~dist.regress)
summary(cafi.age)#P = 0.009 R2 = 0.16

bogr.age<-lm(BOGR~dist.regress)
summary(bogr.age)#P = 0.034 R2 = 0.11 Increase with distance
plot(BOGR~dist.regress)
abline(bogr.age)

eltr.age<-lm(log1p(ELTR)~dist.regress)
summary(eltr.age)#P = 0.01 R2 = 0.15

agcr.age<-lm(log1p(AGCR)~dist.regress)
summary(agcr.age)#P = 0.004 R2 = 0.19

melu.age<-lm(log1p(MELU)~dist.regress)
summary(melu.age)#P = 0.003, R2 = 0.20 decrease with distance
plot(log1p(MELU)~dist.regress)
abline(melu.age)

disp.age<-lm(log1p(DISP)~dist.regress)
summary(disp.age)#P = 0.009 R2 = 0.16

bare.age<-lm(log1p(common.other.cover$BARE)~dist.regress)
summary(bare.age)#P = 0.04 R2 = 0.1

############################################################################
#Hypothesis 3: Change in soil factors across transect types
############################################################################

#Loading soil factor data###################################################

soil.props<-read.csv("Soils by site.csv")
soil.vars<-soil.props[,4:21]
soil.A<-soil.vars[1:45,]
soil.B<-soil.vars[46:90,]
soil.env<-(soil.A+soil.B)/2
soil.env<-soil.env[-(7:9),]

salt.trans<-lmer(soil.env$Soluble_salts~trans+(1|site))
summary(salt.trans)#
anova(salt.trans,ddf="Kenward-Roger")#F = 5.16 P = 0.013
summary(glht(salt.trans,linfct=mcp(trans="Tukey"),test=adjusted("holm")))#50 = 150 < pad

Ca.trans<-lmer(soil.env$Calcium~trans+(1|site))
summary(Ca.trans)#higher on pad
anova(Ca.trans,ddf="Kenward-Roger")#F2,26=8.39, P=0.0015
summary(glht(Ca.trans,linfct=mcp(transect="Tukey")),test=adjusted("holm"))#50 = 150 < Pad

sodium.trans<-lmer(log1p(soil.env$Sodium)~trans+(1|site))
summary(sodium.trans)#Higher on pad
anova(sodium.trans,ddf="Kenward-Roger")#F2,26=5.04, P=0.014
summary(glht(sodium.trans,linfct=mcp(transect="Tukey")),test=adjusted("holm"))#Pad > 50 = 150

pH.trans<-lmer(soil.env$pH~trans+(1|site))
summary(pH.trans)#
anova(pH.trans,ddf="Kenward-Roger")#F2,26 = 7.88 Higher on Pad P = 0.002
summary(glht(pH.trans,linfct=mcp(transect="Tukey")),test=adjusted("holm"))#Pad > 50 = 150

CEC.trans<-lmer(soil.env$CEC~trans+(1|nmds.labels$site))
summary(CEC.trans)#higher on pad
anova(CEC.trans,ddf="Kenward-Roger")#F2,26=13.4, P<0.0001
summary(glht(CEC.trans,linfct=mcp(transect="Tukey")),test=adjusted("holm"))#Pad > 50 = 150

orgmat.trans<-lmer((soil.env$PercentOM_LOI)~trans+(1|site))
summary(orgmat.trans)#suggestively lower on pad
anova(orgmat.trans,ddf="Kenward-Roger")#F2,26=2.01, P=0.15 (NS)
summary(glht(orgmat.trans,linfct=mcp(trans="Tukey"),test=adjusted("holm")))

n.trans<-lmer(log1p(soil.env$Nitrate)~trans+(1|site))
summary(n.trans)# F = 0.91 P = 0.42
anova(n.trans,ddf="Kenward-Roger")#NS
boxplot(soil.env$Nitrate~nmds.labels$trans)

sand.trans<-lmer(soil.env$Sand~trans+(1|site))
summary(sand.trans)#No difference F = 0.35 P = 0.71
anova(sand.trans,ddf="Kenward-Roger")#
boxplot(soil.env$Sand~nmds.labels$trans)

silt.trans<-lmer(soil.env$Silt~trans+(1|nmds.labels$site))
summary(silt.trans)#Lower on pad
anova(silt.trans,ddf="Kenward-Roger")#F2,26=10.57, P=0.0004
summary(glht(silt.trans,linfct=mcp(transect="Tukey")),test=adjusted("holm"))#Pad < 50 = 150

####Reading in soil compaction data
compact<-read.csv("compaction.csv")

compact.three<-aggregate(compact$three,list(compact$site_transect),mean,na.rm=TRUE)
compact.six<-aggregate(compact$six,list(compact$site_transect),mean,na.rm=TRUE)

compaction<-data.frame(nmds.labels,compact.three$x,compact.six$x,since.reclaim)

three.trans<-lmer(compaction$compact.three.x~trans+(1|site))
summary(three.trans)#
anova(three.trans,ddf="Kenward-Roger")#NS
boxplot(compaction$compact.three.x~compaction$trans)
summary(glht(three.trans,linfct=mcp(transect="Tukey")),test=adjusted("holm"))

six.trans<-lmer(compaction$compact.six.x~trans+(1|site))
summary(six.trans)
anova(six.trans,ddf="Kenward-Roger")#F2,26 = 6.354 P = 0.0057
boxplot(compaction$compact.six.x~compaction$trans)
summary(glht(six.trans,linfct=mcp(transect="Tukey")),test=adjusted("holm"))#higher on pads

three.age<-lm(compact.three.x~since.reclaim,data=compaction[compaction$trans=="PAD",])
summary(three.age)
plot(compact.three.x~since.reclaim,data=compaction[compaction$trans=="PAD",])
abline(three.age)

six.age<-lm(compact.six.x~since.reclaim,data=compaction[compaction$trans=="PAD",])
summary(six.age)
plot(compact.three.x~since.reclaim,data=compaction[compaction$trans=="PAD",])
abline(three.age)

############################################################################
#Hypothesis 4: Change in nematode communities across transect types
############################################################################

#Loading nematode data###################################################
nemdata<-read.csv("BakkenNems2016.csv")
nemdata[which(nemdata$Site=="8221"),]=NA

extract<-nemdata$ExtractMass
prop<-nemdata$ProportionCounted
dry<-nemdata$DrySoil
wet<-nemdata$WetSoil
SM<-(nemdata$WetSoil-nemdata$DrySoil)/nemdata$DrySoil
SMtrans<-(SM[1:45]+SM[46:90])/2

#Code below standardizes nematode abundances to # individuals per kg dry soil
bf.all<-1000*((nemdata$Bacterivorous/prop)/(extract*(dry/wet)))
#We originally had each sampling transect subdivided by two to provide a bit finer resolution to the landscape-scale
#comparisons, but scrapped this for single data points per transect. This code averages across both half-transect samples
#to create that single value.
bf<-(bf.all[1:45]+bf.all[46:90])/2

ff.all<-1000*((nemdata$Fungivorous/prop)/(extract*(dry/wet)))
ff<-(ff.all[1:45]+ff.all[46:90])/2

pp.all<-1000*((nemdata$Herbivorous/prop)/(extract*(dry/wet)))
pp<-(pp.all[1:45]+pp.all[46:90])/2

om.all<-1000*((nemdata$Omnivorous/prop)/(extract*(dry/wet)))
om<-(om.all[1:45]+om.all[46:90])/2

pr.all<-1000*((nemdata$Predatory/prop)/(extract*(dry/wet)))
pr<-(pr.all[1:45]+pr.all[46:90])/2

tot.all<-1000*((nemdata$Total/prop)/(extract*(dry/wet)))
tot<-(tot.all[1:45]+tot.all[46:90])/2

bf<-bf[-(7:9)]#deleting data from the 8221 site we discarded upon realizing it didn't fit our site selection requirements
ff<-ff[-(7:9)]
pp<-pp[-(7:9)]
om<-om[-(7:9)]
pr<-pr[-(7:9)]
tot<-tot[-(7:9)]

omnicarn<-om+pr#Creates a new, unified predatory and omnivorous feeding guild

#Testing nematode feeding guilds against transect types#################
bf.test<-lmer(bf~trans+(1|site))
summary(bf.test)#no effect
anova(bf.test,ddf="Kenward-Roger")#F = 0.78 P = 0.47

ff.test<-lmer(log1p(ff)~trans+(1|site))
summary(ff.test)
anova(ff.test,ddf="Kenward-Roger")# F2,26=1.96 P = 0.16 NS

pp.test<-lmer(pp~trans+(1|site))
summary(pp.test)
anova(pp.test,ddf="Kenward-Roger")#F2,26=2.42 P = 0.11
summary(glht(pp.test,linfct=mcp(trans2="Tukey")),test=adjusted("holm"))

tot.test<-lmer(tot~trans+(1|site))
summary(tot.test)
anova(tot.test,ddf="Kenward-Roger")#Pad different F2,26=3.4 P = 0.049
summary(glht(tot.test,linfct=mcp(trans="Tukey")),test=adjusted("holm"))#pad lower but 50 intermediate

omnicarn.test<-lmer(omnicarn~trans+(1|site))
summary(omnicarn.test)
anova(omnicarn.test,ddf="Kenward-Roger")#F2,26 = 4.98 P = 0.015
summary(glht(omnicarn.test,linfct=mcp(trans="Tukey")),test=adjusted("holm"))#50 intermediate

#Testing feeding guilds for interactions between time and distance from reclaims
bf.regress<-lm(bf~dist.regress)
summary(bf.regress)#P = 0.85 R2 = 0.0009

ff.regress<-lm(log1p(ff)~dist.regress)
summary(ff.regress)#P = 0.25, R2 = 0.03

pp.regress<-lm(pp~dist.regress)
summary(pp.regress)#P = 0.08 R2 = 0.07

tot.regress<-lm(tot~dist.regress)
summary(tot.regress)#P = 0.07 R2 = 0.08

omnicarn.regress<-lm(omnicarn~dist.regress)
summary(omnicarn.regress)#P = 0.02 R2 = 0.12

#NMS and Permanova analysis
nema.ord<-data.frame(bf,ff,pp,om,pr)

nema.nmds<-metaMDS(nema.ord,maxit=9999,trymax=9999,distance="bray",k=3)
nema.nmds
nema.nmds$stress
#Shepard plot
stressplot(nema.nmds,main="Shepard Plot")
#goodness of fit plot
plot(nema.nmds,main="Goodness of fit",display="sites",type="t")
points(nema.nmds,display="sites",cex=goodness(nema.nmds)*200)

pairwise.adonis(nema.ord,nmds.labels$trans)#NS

#Testing nematode groups against soil factors
bf.mult<-lm(bf~soil.env$PercentOM_LOI)
summary(bf.mult)#P = 0.012, R2 = 0.15 OM+

ff.mult<-lm(ff~soil.env$PercentOM_LOI)
summary(ff.mult)#P = 0.0058 R2 = 0.18 OM+

pp.mult<-lm(pp~veg.cover$BARE+soil.env$pH)
#pp.mult2<-lm(pp~veg.cover$BARE+soil.env$Silt)
AIC(pp.mult,pp.mult2)
summary(pp.mult)#P = 0.006 R2 = 0.19 pH+, Bare-

tot.mult<-lm(tot~veg.cover$BARE+soil.env$pH)
summary(tot.mult)#P = 0.008 R2 = 0.22 Bare-, pH+

omnicar.mult<-lm(omnicar~veg.cover$CAFI)
summary(omnicar.mult)#P = 0.004 R2 = 0.19 CAFI+

############################################################################
#Examining plant-soil-nematode interactions
############################################################################

#Plant regressions against soil factors to determine what may be driving reclamation impacts
veg.S.mult<-lm(veg.S~soil.env$Soluble_salts)
summary(veg.S.mult)#NS P = 0.25 R2 = 0.03

natives.S.mult<-lm(natives.S~exot.prop)
summary(natives.S.mult)#R2 = 0.30, P = 0.0001

exotics.S.mult<-lm(exotics.S~soil.env$Soluble_salts+exot.prop)
summary(exotics.S.mult)#P = 0.0009 R2 = 0.27 

natives.sum.mult<-lm(natives.sum~soil.env$Soluble_salts)
summary(natives.sum.mult)#NS

exot.prop.mult<-lm(exot.prop~soil.env$PercentOM_LOI)
summary(exot.prop.mult)#Close, suggestive positive

ruderals.mult<-lm(log1p(ruderals)~common.other.cover$BARE)
summary(ruderals.mult)#P = 0.001 R2 = 0.24

invasives.mult<-lm(log1p(invasives)~soil.env$PercentOM_LOI)
summary(invasives.mult)#P = 0.01, R2 = 0.15

grass.mult<-lm(grass~soil.env$Soluble_salts+soil.env$PercentOM_LOI)
summary(grass.mult)#P < 0.0001 R2 = 0.45

forb.mult<-lm(log1p(forb)~soil.env$Silt)
summary(forb.mult)#P = 0.019 R2 = 0.13

cafi.mult<-lm(CAFI~soil.env$Sodium)
summary(cafi.mult)#

bogr.mult<-lm(BOGR~soil.env$CEC)
summary(bogr.mult)#P = 0.04 R2 = 0.1

eltr.mult<-lm(log1p(ELTR)~soil.env$PercentOM_LOI)
summary(eltr.mult)#P = 0.045 R2 = 0.1

agcr.mult<-lm(log1p(veg.cover$AGCR)~soil.env$PercentOM_LOI)
summary(agcr.mult)#NS

melu.mult<-lm(MELU~soil.env$PercentOM_LOI)
summary(melu.mult)#P = 0.0008 R2 = 0.25

disp.mult<-lm(DISP~soil.env$Soluble_salts)
summary(disp.mult)#P = 0.001 R2 = 0.23

bare.mult<-lm(log1p(veg.cover$BARE)~soil.env$PercentOM_LOI+soil.env$pH+soil.env$Sodium)
summary(bare.mult)#P < 0.0001 R2 = 0.54 OM-, pH+, Na+

################################################
#SEM for veg/nem/soil
################################################

####Running PCA of nematode communities to take first axis as summary variable
nem.summary<-rda(nema.ord,scale=TRUE)

plot(nem.summary)
summary(nem.summary)#axis 1 = 37% data explained

nem.community<-scores(nem.summary,choices=1,display="sites")
nem.community<-nem.community[1:42]

########Creating a single dataframe with all SEM variables included to facilitate running model
soil.om<-soil.env$PercentOM_LOI
soil.na<-soil.env$Sodium/1000
soil.salts<-soil.env$Soluble_salts
silt<-soil.env$Silt
litter.sem<-common.other.cover$LITTER/100
bare.sem<-common.other.cover$BARE
perennials.native<-(rowSums(natives.common[,-15]))/100
perennials.exotic<-(exotics.common$AGCR+exotics.common$BRIN+exotics.common$BRJA+exotics.common$POCO+exotics.common$POPR+exotics.common$veg.cover.MEOF)/100
ruderals<-natives.common$GRSQ+exotics.common$MELU+exotics.common$TAOF+exotics.common$TRDU
nema.veg<-data.frame(bare.sem,nem.community,perennials.native,perennials.exotic,ruderals,soil.om,soil.na,soil.salts,silt,litter.sem)

nema.veg<-within(nema.veg,{
  bare.sem<-lm(bare.sem~site,data=nema.veg)$residuals
  soil.om<-lm(soil.om~site,data=nema.veg)$residuals
  soil.na<-lm(soil.na~site,data=nema.veg)$residuals
  soil.salts<-lm(soil.salts~site,data=nema.veg)$residuals
  perennials.native<-lm(perennials.native~site,data=nema.veg)$residuals
  perennials.exotic<-lm(perennials.exotic~site,data=nema.veg)$residuals
  ruderals<-lm(ruderals~site,data=nema.veg)$residuals
  nem.community<-lm(nem.community~site,data=nema.veg)$residuals
  litter.sem<-lm(litter.sem~site,data=nema.veg)$residuals
})

#Initial Hypothesis model
soil.veg.model<-'
nem.community~soil.om+litter.sem+perennials.native+perennials.exotic+bare.sem
ruderals~perennials.native+perennials.exotic+soil.salts+bare.sem+soil.om
perennials.native~perennials.exotic+soil.om+soil.salts
perennials.exotic~soil.om+bare.sem
litter.sem~perennials.native+perennials.exotic
'

#Final model
soil.veg.model<-'
nem.community~bare.sem
ruderals~soil.salts+bare.sem
perennials.native~perennials.exotic+soil.om+bare.sem
perennials.exotic~soil.om+bare.sem
'

soil.veg.fit<-sem(soil.veg.model,data=nema.veg,estimator="MLM")
#This provides the fit statistics and multiple regression table (Supplementary Table S5) for the SEM
summary(soil.veg.fit)

#This provides the effect size and significance for the path coefficients of the SEM
standardizedSolution(soil.veg.fit)

#This provides the R2 values for each component of the SEM
inspect(soil.veg.fit,"r2")

#Quick plot to look at structure of SEM to facilitate creating box-and-arrow diagram
lavaanPlot(soil.veg.fit,node_options=list(shape="box",fontname=""),edge_options=list(color="black"))

###################################################################
#Code for figures
###################################################################

#####################
#Code for Figure 1
#####################
par(mfrow=c(2,2))
boxplot(plants.ord$CAFI~nmds.labels$trans,ylab="% Cover",main="Carex filifolia",frame="F")#Pad lower P < 0.0001
boxplot(plants.ord$BOGR~nmds.labels$trans,ylab="% Cover",main="Bouteloua gracilis",frame="F")#50 > 150 > pad P < 0.0001
boxplot(plants.ord$AGCR~nmds.labels$trans,ylab="% Cover",main="Agropyron cristatum",frame="F")#Pad greater P = 0.035
boxplot(plants.ord$DISP~nmds.labels$trans,ylab="% Cover",main="Distichlis spicata",frame="F")#Pad greater P < 0.0001
par(mfrow=c(1,1))

#####################
#Code for Figure 2
#####################
par(mfrow=c(2,3))
boxplot(natives.S~trans,ylab="Species richness",main="Native plant richness",frame="F",axes="T")
boxplot(exotics.S~trans,ylab="Species richness",main="Exotic plant richness",frame="F",axes="T")
boxplot(natives.prop~trans,main="Native plant species cover",ylab="Relative cover",frame="F",axes="T")
boxplot(bare~trans,ylab="% Cover",main="Bare ground by transect",frame="F",axes="T")
boxplot(ruderals~trans,ylab="% Cover",main="Ruderal cover by transect",frame="F",axes="T")
boxplot(grass~trans,ylab="% Cover",main="Grass cover by transect",frame="F",axes="T")
par(mfrow=c(1,1))

#################################
#Code for Figure 3
#################################
plot(plant.nmds,type="n",xlim=c(-1,1.5),ylim=c(-1,1),main=paste("Plant Species NMDS/Bray - stress=",round(plant.nmds$stress,3)),frame="F",axes="F")
points(plant.nmds,display="sites",pch=transvec,cex=1.25)
#legend("bottomleft",legend=c("150 m","50 m","Reclaim"),bty="n",pch=transleg,cex=1.25)
axis(1)
axis(2)

#################################
#Code for Figure 4
#################################
par(mfrow=c(2,3))
boxplot(soil.env$Soluble_salts~trans,main="Soluble salts",ylab="Parts per million",frame="F")
boxplot(soil.env$Calcium~nmds.labels$trans,main="Calcium",ylab="Parts per million",frame="F")
boxplot(soil.env$Sodium~nmds.labels$trans,main="Sodium",ylab="Parts per million",frame="F")
boxplot(soil.env$pH~nmds.labels$trans,main="pH",ylab="pH",frame="F")
boxplot(soil.env$Silt~nmds.labels$trans,main="Silt",ylab="% Composition",frame="F")
boxplot(soil.env$CEC~nmds.labels$trans,main="Cation exchange capacity",ylab="Cation meq/100g",frame="F")
par(mfrow=c(1,1))

####################################
#Code for Figure 5
####################################
par(mfrow=c(1,2))
boxplot(tot~trans,main="Total nematodes",ylab="Individuals per kg dry soil",frame="F")
boxplot(omnicarn~trans,main="Omni-carnivorous nematodes",ylab="Individuals per kg dry soil",frame="F")
par(mfrow=c(1,1))

###################################################################
#Code for Supplementary Material
###################################################################
for.tables<-data.frame(nmds.labels,bf,ff,pp,omnicar,soil.env)

#######Supplementary Table S1
soils.u<-aggregate(for.tables[,7:24],list(for.tables$trans),mean)
soils.sd<-aggregate(for.tables[,7:24],list(for.tables$trans),sd)
soils.se<-soils.sd
soils.se[,2:19]<-soils.se[,2:19]/sqrt(14)#calculates standard errors

#######Supplementary Table S4
nem.u<-aggregate(for.tables[,3:6],list(for.tables$trans),mean)
nem.sd<-aggregate(for.tables[,3:6],list(for.tables$trans),sd)
nem.se<-nem.sd
nem.se[,2:5]<-nem.se[,2:5]/sqrt(14)
nemtot<-bf+ff+pp+omnicar
tot.u<-aggregate(nemtot,list(for.tables$trans),mean)
tot.sd<-aggregate(nemtot,list(for.tables$trans),sd)
tot.sd<-tot.se
tot.se[,2]<-tot.se[,2]/sqrt(14)

#######Supplementary Figure S1
hist(since.reclaim,xlab="Years since reclamation",main="")
