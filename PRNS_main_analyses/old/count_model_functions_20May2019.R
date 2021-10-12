#Models for PRNS elk project (Count survey data). Created by Kevin Shoemaker and edited by Lacey Hughey on May 20, 2019.


#############
# Clear workspace
#############

rm(list=ls())

setwd ("/Users/Lacey/Box Sync/PhD_Project/R.Analyses/Elk/Raw data")

#############
# Install Development versions for glmmTMB and DHARMa (KTS: I did this but I don't think it's needed actually)

# devtools::install_github("glmmTMB/glmmTMB",build_vignettes=FALSE,dependencies=TRUE,subdir="glmmTMB")
# devtools::install_github("glmmTMB/glmmTMB/glmmTMB")

# devtools::install_github(repo = "florianhartig/DHARMa", subdir = "DHARMa", dependencies = T, build_vignettes = T)

# devtools::install_github("bbolker/broom.mixed")

#############
# Load packages
#############

library(bbmle)
library(glmmTMB)
library(DHARMa)
library(car)
library(lme4)
library(buildmer)


#############
# Load custom functions
#############

VisualizeRelation <- function(data=df,model=full,predvar="std.slope",allvars=names(df)[-c(1,6:13)],season){
  len <- 100
  
  predvar2 <- strsplit(predvar,"\\.")[[1]][2]
  
  dataclasses <- sapply(data,class)
  
  dim <- data[,predvar]
  range <- seq(min(dim),max(dim),length=len)
  
  realmean <- mean(data[[predvar2]])
  realsd <- sd(data[[predvar2]])
  
  newdata <- data.frame(temp=range)
  names(newdata) <- c(predvar)
  
  othervars <- allvars[!allvars%in%c(predvar,"TotalElk")]
  
  
  var = othervars[5]
  for(var in othervars){
    thisvar <- data[[var]]
    if(is.factor(thisvar)){
      tab <- table(thisvar)
      vals <- names(tab)
      levs <- levels(thisvar)
      mostcom <- vals[which.max(tab)]
      newvec <- factor(rep(mostcom,times=nrow(newdata)),levels=levs)
      newdata[,var] <- newvec
    }else{
      newdata[,var] <- 0 #mean(thisvar)
    }
  }
  
  
  pred_pz <- predict(model,newdata,type="zprob",se.fit=T)
  pred_ct <- predict(model,newdata,type="conditional",se.fit=T)
  #pred <- predict(model,newdata,type="response",se.fit=T)
  
  layout(matrix(c(1:2),nrow=2))
  par(mai=c(0.9,0.8,0.1,0.1))
  plot(range,pred_pz$fit,xlab=predvar2,ylab="Prob of zero",type="l",lwd=2,xaxt="n",ylim=c(0,1))
  points(range,pred_pz$fit+pred_pz$se.fit,type="l",lty=2)
  points(range,pred_pz$fit-pred_pz$se.fit,type="l",lty=2)
  ats <- seq(min(range),max(range),length=6)
  axis(1,ats,labels = round(realmean+ats*realsd,2))
  rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  
  plot(range,pred_ct$fit,xlab=predvar2,ylab="Exp Count",type="l",lwd=2,xaxt="n",ylim=c(0,20))
  points(range,pred_ct$fit+pred_ct$se.fit,type="l",lty=2)
  points(range,pred_ct$fit-pred_ct$se.fit,type="l",lty=2)
  ats <- seq(min(range),max(range),length=6)
  axis(1,ats,labels = round(realmean+ats*realsd,2))
  rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  
  svg(sprintf("pdplots_%s_%s.svg",predvar2,season),width=4.5,height = 6)
  
  layout(matrix(c(1:2),nrow=2))
  par(mai=c(0.9,0.8,0.1,0.1))
  plot(range,pred_pz$fit,xlab=predvar2,ylab="Prob of zero",type="l",lwd=2,xaxt="n",ylim=c(0,1))
  points(range,pred_pz$fit+pred_pz$se.fit,type="l",lty=2)
  points(range,pred_pz$fit-pred_pz$se.fit,type="l",lty=2)
  ats <- seq(min(range),max(range),length=6)
  axis(1,ats,labels = round(realmean+ats*realsd,2))
  rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  
  plot(range,pred_ct$fit,xlab=predvar2,ylab="Exp Count",type="l",lwd=2,xaxt="n",ylim=c(0,20))
  points(range,pred_ct$fit+pred_ct$se.fit,type="l",lty=2)
  points(range,pred_ct$fit-pred_ct$se.fit,type="l",lty=2)
  ats <- seq(min(range),max(range),length=6)
  axis(1,ats,labels = round(realmean+ats*realsd,2))
  rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  
  dev.off()
}

VisualizeInteraction <- function(data=df,model=bestmod2,var1="std.StockingMed",var2="std.ndvi",allvars=names(df)[-c(1,6:13)],season){
  len <- 30     # increase this for higher-res figures
  
  dataclasses <- sapply(data,class)
  
  var1_2 <- strsplit(var1,"\\.")[[1]][2]
  var2_2 <- strsplit(var2,"\\.")[[1]][2]
  
  realmean1 <- mean(data[[var1_2]])
  realsd1 <- sd(data[[var1_2]])
  realmean2 <- mean(data[[var2_2]])
  realsd2 <- sd(data[[var2_2]])
  standvar1 <- var1
  standvar2 <- var2
  
  
  firstdim <- data[,standvar1]
  seconddim <- data[,standvar2]
  range1 <- seq(min(firstdim),max(firstdim),length=len)
  range2 <- seq(min(seconddim),max(seconddim),length=len)
  newdata <- expand.grid(range1,range2)
  # head(newdata,50)
  names(newdata) <- c(standvar1,standvar2)
  
  othervars <- allvars[!allvars%in%c(standvar1,standvar2,"TotalElk")]
  
  var = othervars[2]
  for(var in othervars){
    thisvar <- data[[var]]
    if(is.factor(thisvar)){
      tab <- table(thisvar)
      vals <- names(tab)
      levs <- levels(thisvar)
      mostcom <- vals[which.max(tab)]
      newvec <- factor(rep(mostcom,times=nrow(newdata)),levels=levs)
      newdata[[var]] <- newvec
    }else{
      newdata[[var]] <- mean(thisvar)
    }
  }
  
  pred_pz <- predict(model,newdata,type="zprob")
  pred_ct <- predict(model,newdata,type="conditional")
  
  predmat_pz <-  matrix(pred_pz,nrow=len,ncol=len)
  predmat_ct <-  matrix(pred_ct,nrow=len,ncol=len)
  
  
  
  svg(sprintf("intplots_%s_%s_%s.svg",var1_2,var2_2,season),width=4.5,height = 6)
  
  layout(matrix(c(1:2),nrow=2))
  par(mai=c(0,0,0,0))
  
  persp(realmean1+realsd1*range1,realmean2+realsd2*range2,predmat_pz,xlab=var1_2,ylab=var2_2,theta = 55, phi = 40, r = sqrt(10), d = 3, 
        ticktype = "detailed", mgp = c(4, 1, 0))
  
  persp(realmean1+realsd1*range1,realmean2+realsd2*range2,predmat_ct,xlab=var1_2,ylab=var2_2,theta = 55, phi = 40, r = sqrt(10), d = 3, 
        ticktype = "detailed", mgp = c(4, 1, 0))
  
  dev.off() 
}


#############
# Load data
#############

####rename files to comply with requirements of script. only do this once, then call these rds files in the future.
# r.temp <- readRDS("rut.obs.19May2019.rds")
# p.temp <- readRDS("part.obs.19May2019.rds")
# s.temp <- readRDS("summer.obs.19May2019.rds")
# w.temp <- readRDS("winter.obs.19May2019.rds")
# 
# #rename variables to comply with requirements of script
# colnames(r.temp)[17] <- "stockingMeanNonz"
# colnames(r.temp)[18] <- "stockingMean"
# colnames(r.temp)[19] <- "stockingMed"
# 
# colnames(p.temp)[17] <- "stockingMeanNonz"
# colnames(p.temp)[18] <- "stockingMean"
# colnames(p.temp)[19] <- "stockingMed"
# 
# colnames(s.temp)[17] <- "stockingMeanNonz"
# colnames(s.temp)[18] <- "stockingMean"
# colnames(s.temp)[19] <- "stockingMed"
# 
# colnames(w.temp)[17] <- "stockingMeanNonz"
# colnames(w.temp)[18] <- "stockingMean"
# colnames(w.temp)[19] <- "stockingMed"
# 
# saveRDS(r.temp, "rut.obs.27May2019.rds")
# saveRDS(p.temp, "part.obs.27May2019.rds")
# saveRDS(s.temp, "summer.obs.27May2019.rds")
# saveRDS(w.temp, "winter.obs.27May2019.rds")

#load files

alldata<-list()

alldata$rut <- readRDS("rut.obs.27May2019.rds")
alldata$parturition <- readRDS("part.obs.27May2019.rds")
alldata$summer <- readRDS("summer.obs.27May2019.rds")
alldata$winter <- readRDS("winter.obs.27May2019.rds")

##############
# Process data

#replace with subset of columns of interest for appropriate season
alldata2 <- lapply(alldata,function(t) t[c(1, 2, 9, 12, 14, 16:26)])


#standardize variables...

standardze <- function(df){
  df$std.pasture <- ((df$pasture)-(mean(df$pasture)))/(2*sd(df$pasture))
  df$std.stockingMeanNonz <- ((df$stockingMeanNonz )-(mean(df$stockingMeanNonz , na.rm=TRUE)))/(2*sd(df$stockingMeanNonz , na.rm=TRUE))
  df$std.stockingMean <- ((df$stockingMean)-(mean(df$stockingMean, na.rm=TRUE)))/(2*sd(df$stockingMean, na.rm=TRUE))
  df$std.stockingMed <- ((df$stockingMed)-(mean(df$stockingMed, na.rm=TRUE)))/(2*sd(df$stockingMed, na.rm=TRUE))
  df$std.elevation <- ((df$elevation)-(mean(df$elevation)))/(2*sd(df$elevation))
  df$std.ndvi <- ((df$ndvi)-(mean(df$ndvi, na.rm=TRUE)))/(2*sd(df$ndvi, na.rm=TRUE))
  df$std.ponds <- ((df$ponds)-(mean(df$ponds)))/(2*sd(df$ponds))
  df$std.scrub <- ((df$scrub)-(mean(df$scrub)))/(2*sd(df$scrub))
  df$std.grassland <- ((df$grassland)-(mean(df$grassland)))/(2*sd(df$grassland))
  df$std.slope <- ((df$slope)-(mean(df$slope)))/(2*sd(df$slope))
  df$std.aspect <- ((df$aspect)-(mean(df$aspect)))/(2*sd(df$aspect))
  return(df)
}

data.std <- lapply(alldata2,function(t) standardze(t))

sapply(data.std,function(t) nrow(t)-length(which(complete.cases(t))))

final.data <- lapply(data.std,function(t) t[complete.cases(t), ])


tofac <- function(df){
  df$year <- as.factor(df$year)
  df$month <- as.factor(df$month)
  df$ID <- as.factor(df$ID)
  return(df)
}

final.data <- lapply(final.data,function(t) tofac(t))


# ############
# # Select a season to focus on for initial selection of model structure- then repeat for all seasons
# ############
# 
# season <- "parturition"
# 
# df <- final.data[[season]]
# 
# ############
# # Perform initial diagnoses of appropriate count distributions
# ############
# 
# # test Poisson
# fit <- vcd::goodfit(df$TotalElk,type="poisson")               # test fails for Poisson count model every time
# summary(fit)
# vcd::rootogram(fit)                                           # terrible fit
# vcd::Ord_plot(df$TotalElk)                                    # this diagnoses that neg binom (or no) distribution is appropriate!
# 
# # test NegBin
# fit <- vcd::goodfit(df$TotalElk, type="nbinom",method="ML")   # goodness of fit test for neg binom
# summary(fit)                                                  # terrible fit
# vcd::rootogram(fit)                                           # terrible fit
# vcd::distplot(df$TotalElk, type="nbinom")                     # Neg binom not good either- but observations (open circles) fall within the CI
# 
# 
# hist(df$TotalElk)   
# hist(log(df$TotalElk))
# 
# 
# ###########
# # Explore relationships visually
# ###########
# 
# plot(df$stocking.mean.nonz,df$TotalElk)
# plot(df$pasture,df$TotalElk)
# plot(df$elevation,df$TotalElk)    # possibility for quadratic effect? 
# plot(df$slope,df$TotalElk)   
# plot(df$ndvi,df$TotalElk)         # quadratic?
# plot(df$scrub,df$TotalElk)        # quadratic?
# plot(df$ponds,df$TotalElk)
# plot(df$slope,df$TotalElk)
# 
# table(df$year)
# base::tapply(df$TotalElk,df$year,sum)


# ##############
# # STEP 1: determine which distribution to use
# ##############
# 
# ### try poisson regression with no zero component
# 
# test1.mod <- glmmTMB(TotalElk ~ std.pasture + std.stockingMeanNonz + std.ponds 
#                      + std.elevation + std.ndvi + std.scrub + std.grassland + std.slope 
#                      + std.aspect + std.ndvi:std.scrub + std.ndvi:std.grassland 
#                      + std.stockingMeanNonz:std.ndvi
#                      + std.slope:std.aspect
#                      + (1|ID) + (1|year),
#                      df,
#                      ziformula = ~0 ,
#                      family= poisson(link = "log"))
# 
# summary(test1.mod)
# 
# 
# test1.res <- DHARMa::simulateResiduals(test1.mod,n=300)
# plot(test1.res)                        # terrible fit
# DHARMa::testUniformity(test1.res)      # fail
# DHARMa::testResiduals(test1.res)       # fails uniformity
# DHARMa::testZeroInflation(test1.res)   # fail
# 
# 
# ### try negative binomial regression with no zero component
# 
# test2.mod <- glmmTMB(TotalElk ~ std.pasture + std.stockingMeanNonz + std.ponds 
#                      + std.elevation + std.ndvi + std.scrub + std.grassland + std.slope 
#                      + std.aspect + std.ndvi:std.scrub + std.ndvi:std.grassland 
#                      + std.stockingMeanNonz:std.ndvi
#                      + std.slope:std.aspect
#                      + (1|ID) + (1|year),
#                      df,
#                      ziformula = ~0 ,
#                      dispformula = ~1,                     # should we look at alternative dispersion formulas? Why just month?
#                      family= nbinom1(link = "log"))
# 
# summary(test2.mod)
# 
# 
# test2.res <- DHARMa::simulateResiduals(test2.mod,n=300)
# plot(test2.res)      
# DHARMa::testUniformity(test2.res)         # pass
# DHARMa::testResiduals(test2.res)          # passes all tests (but glmer.nb did not)
# DHARMa::testZeroInflation(test2.res)      # passes zero inflation test
# 
# 
# ### try hurdle model with Poisson distribution
# 
# test3.mod <- glmmTMB(TotalElk ~ std.pasture + std.stockingMeanNonz + std.ponds 
#                      + std.elevation + std.ndvi + std.scrub + std.grassland + std.slope 
#                      + std.aspect + std.ndvi:std.scrub + std.ndvi:std.grassland 
#                      + std.stockingMeanNonz:std.ndvi
#                      + std.slope:std.aspect
#                      + (1|ID) + (1|year), 
#                      df,
#                      ziformula = ~. ,
#                      dispformula = ~1,                     # should we look at alternative dispersion formulas? Why just month?
#                      family= truncated_poisson(link = "log"))
# 
# summary(test3.mod)     #looks good
# 
# test3.res <- DHARMa::simulateResiduals(test3.mod,n=300)
# plot(test3.res)     
# DHARMa::testUniformity(test3.res)        # pass
# DHARMa::testResiduals(test3.res)         # fails outlier test
# DHARMa::testZeroInflation(test3.res)     # pass
# 
# 
# ### try hurdle model with NegBinom distribution
# 
# test4.mod <- glmmTMB(TotalElk ~ std.pasture + std.stockingMeanNonz + std.ponds 
#                      + std.elevation + std.ndvi + std.scrub + std.grassland + std.slope 
#                      + std.aspect + std.ndvi:std.scrub + std.ndvi:std.grassland 
#                      + std.stockingMeanNonz:std.ndvi
#                      + std.slope:std.aspect
#                      + (1|ID) + (1|year),     
#                      df,
#                      ziformula = ~. ,
#                      dispformula = ~1,                     # should we look at alternative dispersion formulas? Why just month?
#                      family= truncated_nbinom1(link = "log"))
# 
# summary(test4.mod)
# test4.res <- DHARMa::simulateResiduals(test4.mod,n=300)
# plot(test4.res)            
# DHARMa::testUniformity(test4.res)          # pass
# DHARMa::testResiduals(test4.res)           # passes all tests
# DHARMa::testZeroInflation(test4.res)       # pass      
# 
# ### try hurdle model with NegBinom distribution and without dispersion formula
# 
# test5.mod <- glmmTMB(TotalElk ~ std.pasture + std.stockingMeanNonz + std.ponds 
#                      + std.elevation + std.ndvi + std.scrub + std.grassland + std.slope 
#                      + std.aspect + std.ndvi:std.scrub + std.ndvi:std.grassland 
#                      + std.stockingMeanNonz:std.ndvi
#                      + std.slope:std.aspect
#                      + (1|ID) + (1|year),     
#                      df,
#                      ziformula = ~. ,
#                      family= truncated_nbinom1(link = "log"))
# 
# summary(test5.mod)
# test5.res <- DHARMa::simulateResiduals(test5.mod,n=300)
# plot(test5.res)            
# DHARMa::testUniformity(test5.res)          # pass
# DHARMa::testResiduals(test5.res)           # passes all tests
# DHARMa::testZeroInflation(test5.res)       # pass      
# 
# # AIC model selection table
# bbmle::AICtab(test1.mod,test2.mod,test3.mod,test4.mod, test5.mod,weights=TRUE,mnames=c("Poisson","NegBin","Hurdle Pois","Hurdle NegBin", "Hudle NegBin NoDisp") )
# 
# 
# ###################
# # Conclusions from structural model selection procedure (applies to Rut, Summer, Winter, Part:   
# # Hurdle NegBin model gets all the AIC weight!
# # Hurdle NegBin model passes the GOF tests
# # Perform remaining analyses using hurdle NegBin structure
# # NegBin model (without zero model) fits well, but explains much less of the total variance
# 
# 
# 
# ##################
# # STEP 2: CHOOSE TOP CATTLE STOCKING VARIABLE 
# ##################
# 
# ###Run through each season one at a time
# 
# season <- "summer"
# 
# df <- final.data[[season]]
# 
# cow.mean.nonz <- glmmTMB(TotalElk ~ std.pasture + std.stockingMeanNonz + std.ponds
#                      + std.elevation + std.ndvi + std.scrub + std.grassland + std.slope
#                      + std.aspect + std.ndvi:std.scrub + std.ndvi:std.grassland
#                      + std.stockingMeanNonz:std.ndvi
#                      + std.slope:std.aspect
#                      + (1|ID) + (1|year),
#                      df,
#                      ziformula = ~. ,
#                      family= truncated_nbinom1(link = "log"))
# 
# cow.mean <- glmmTMB(TotalElk ~ std.pasture + std.stockingMean + std.ponds
#                     + std.elevation + std.ndvi + std.scrub + std.grassland + std.slope
#                     + std.aspect + std.ndvi:std.scrub + std.ndvi:std.grassland
#                     + std.stockingMean:std.ndvi
#                     + std.slope:std.aspect
#                     + (1|ID) + (1|year),
#                     df,
#                     ziformula = ~. ,
#                     family= truncated_nbinom1(link = "log"))
# 
# cow.med <- glmmTMB(TotalElk ~ std.pasture + std.stockingMed + std.ponds
#                      + std.elevation + std.ndvi + std.scrub + std.grassland + std.slope
#                      + std.aspect + std.ndvi:std.scrub + std.ndvi:std.grassland
#                      + std.stockingMed:std.ndvi
#                      + std.slope:std.aspect
#                      + (1|ID) + (1|year),
#                      df,
#                      ziformula = ~. ,
#                      family= truncated_nbinom1(link = "log"))
# 
# bbmle::AICtab(cow.mean.nonz, cow.mean, cow.med,weights=TRUE,mnames=c("Mean Nonz", "Mean", "Median") )

# ######Results from AIC test:
# 
# ##Parturition:
# # dAIC df weight
# # Median     0.0 33 0.43  
# # Mean       0.2 33 0.38  
# # Mean Nonz  1.6 33 0.19  
# 
# ##Summer:
# # dAIC df weight
# # Median     0.0 33 0.659 
# # Mean Nonz  1.7 33 0.282 
# # Mean       4.8 33 0.059 
# 
# #Winter:
# # Median     0.0 33 0.9939
# # Mean      10.9 33 0.0044
# # Mean Nonz 12.7 33 0.0017
# 
# #Rut:
# # dAIC df weight
# # Mean Nonz  0.0 33 0.83  
# # Mean       3.5 33 0.14  
# # Median     7.4 33 0.02  

###########
# STEP 3: AUTOMATE THE REMAINDER OF MODEL SELECTION AND VISUALIZATIONS
###########

############
# Loop through seasons (this will take a bit of time)
############

##finished:


allseasons <- names(final.data)
bestmods <- list()


season <- allseasons[3]
for(season in allseasons){
  df <- final.data[[season]]
  
  bestmods[[season]] <- buildmer::buildglmmTMB(TotalElk ~ std.pasture + std.stockingMed 
                                               + std.ponds + std.elevation + std.ndvi 
                                               + std.scrub + std.grassland + std.slope 
                                               + std.aspect + std.ndvi:std.scrub 
                                               + std.ndvi:std.grassland 
                                               + std.stockingMed:std.ndvi
                                               + std.slope:std.aspect
                                               + (1|ID) + (1|year),ziformula = ~.,data=df, dispformula = ~1, direction="backward",crit="AIC",
                                               reduce.random=F,family=truncated_nbinom1 (link="log"))
  bestmod2 <- bestmods[[season]]@model
  #summary(bestmod2)
  tmp <- rownames(coef(summary(bestmod2))$cond)[-1]
  topvars <- tmp[!grepl(":",tmp)]
  ints <- strsplit(tmp[grep(":",tmp)],":")
  if(length(topvars)>0) tmp <- sapply(1:length(topvars), function(t) VisualizeRelation(df,bestmod2,topvars[t],season=season) )     # run and save partial dependence plots for all top variables
  
  if(length(ints)>0 ) tmp <- sapply(1:length(ints), function(t) VisualizeInteraction(df,bestmod2,ints[[t]][1],ints[[t]][2],season=season) )
  
  
  
  
}


############
# Visualize coefficients


temp_count <- lapply(bestmods, function(t) coef(summary(t))$cond )    # count model coefs
temp_zero <- lapply(bestmods, function(t) coef(summary(t))$zi )     # zero model coefs

allcoefnames <- unique(unlist(sapply(temp_count,function(t) rownames(t) )))
allmains <- sort(allcoefnames[!grepl(":",allcoefnames)])
allints <- sort(allcoefnames[grepl(":",allcoefnames)])
allcoefnames <- c(allmains,allints)

coef_ndx <- lapply(temp_count,function(t) match(allcoefnames,rownames(t)) )

temp <- rep(NA,times=length(allcoefnames))
allcoefs <- lapply(1:length(temp_count), function(t) temp_count[[t]][,'Estimate'][coef_ndx[[t]]])
names(allcoefs) <- allseasons

allses <- lapply(1:length(temp_count), function(t) temp_count[[t]][,'Std. Error'][coef_ndx[[t]]])
names(allses) <- allseasons

allcoefsz <- lapply(1:length(temp_zero), function(t) temp_zero[[t]][,'Estimate'][coef_ndx[[t]]])
names(allcoefs) <- allseasons

allsesz <- lapply(1:length(temp_zero), function(t) temp_zero[[t]][,'Std. Error'][coef_ndx[[t]]])
names(allsesz) <- allseasons


svg("coefplot_zanb.svg",8,8)

layout(matrix(1:9,nrow=3,byrow = T))
var_ndx <- 2
for(var_ndx in c(1:9)){
  par(mai=c(0.1,0.1,0.5,0))
  plot(15,15,xlim=c(0,15),ylim=c(0,15),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main=allcoefnames[var_ndx],pch="")    # blank plot
  xcents <- c(7,11)
  #abline(v=xcents,lwd=2)
  ycents <- seq(14,2,length=4)
  segments(xcents,rep(range(ycents)[1],2)-1,xcents,rep(range(ycents)[2],2)+1,lwd=2)
  text(rep(1,4),ycents, tools::toTitleCase(allseasons),adj=0)
  points(sapply(allcoefs,function(t) t[var_ndx])+xcents[1],ycents,pch=20,cex=1.5)
  arrows(sapply(allcoefs,function(t) t[var_ndx])+xcents[1]-sapply(allses,function(t) t[var_ndx]),
         ycents,sapply(allcoefs,function(t) t[var_ndx])+xcents[1]+sapply(allses,function(t) t[var_ndx]),
         ycents,angle = 90,code=3,length=0.05)
  points(sapply(allcoefsz,function(t) t[var_ndx])+xcents[2],ycents,pch=20,cex=1.5)
  arrows(sapply(allcoefsz,function(t) t[var_ndx])+xcents[2]-sapply(allsesz,function(t) t[var_ndx]),
         ycents,sapply(allcoefsz,function(t) t[var_ndx])+xcents[2]+sapply(allsesz,function(t) t[var_ndx]),
         ycents,angle = 90,code=3,length=0.05)
  text(xcents,c(0,0),c("Count","Zero"),cex=1.2)
  
}

dev.off()

graphics.off()

############
# End of script
############




########
# TODO
########

# Add visualization of coefficients for both the zero and count models- for each season









#############
# OLD CODE:

# #########
# # STEP 3: Select best model from the "global model" using backward stepwise selection
# 
# # full <- test4.mod
# 
# full <- glmmTMB(TotalElk ~ std.pasture + std.stocking + std.ponds + std.elevation 
#                 + std.insolation + std.ndvi + std.scrub + std.slope 
#                 + std.scrub:std.ndvi + std.stocking:std.ndvi
#                 + (1|ID) + (1|year),  
#                 df,
#                 ziformula = ~.,
#                 dispformula = ~month,
#                 family= truncated_nbinom1 (link="log"))
# 
# summary(full)
# 
# 
# bestmod <- buildmer::buildglmmTMB(TotalElk ~ std.pasture + std.stocking + std.ponds + std.elevation 
#                                   + std.insolation + std.ndvi + std.scrub + std.slope 
#                                   + std.scrub:std.ndvi + std.stocking:std.ndvi
#                                   + (1|ID) + (1|year),ziformula = ~.,data=df, dispformula = ~month,direction="backward",crit="AIC",
#                                   reduce.random=F,family=truncated_nbinom1 (link="log"))
# 
# bestmod2 <- bestmod@model
# summary(bestmod2)
# 
# #######
# # partial dependence plots for full model (visualize key relationships)
# 
# #VisualizeRelation(df,full,'std.pasture',season=season)
# VisualizeRelation(df,bestmod2,'std.stocking',season=season)
# VisualizeRelation(df,bestmod2,'std.ponds',season=season)
# VisualizeRelation(df,bestmod2,'std.elevation',season=season)
# VisualizeRelation(df,bestmod2,'std.insolation',season=season)
# VisualizeRelation(df,bestmod2,'std.ndvi',season=season)
# VisualizeRelation(df,bestmod2,'std.scrub',season=season)
# #VisualizeRelation(df,full,'std.slope',season=season)
# 
# VisualizeInteraction(df,bestmod2,"std.stocking","std.ndvi",season=season)












