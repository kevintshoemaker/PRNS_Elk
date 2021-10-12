#Models for PRNS elk project (Count survey data and collar data). 
# Created by Kevin Shoemaker and edited by Lacey Hughey on May 20, 2019.


#############
# Clear workspace
#############

rm(list=ls())

#setwd ("/Users/Lacey/Box Sync/PhD_Project/R.Analyses/Elk/Raw data")

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
  plot(range,1-pred_pz$fit,xlab=predvar2,ylab="P(1+)",type="l",lwd=2,xaxt="n",ylim=c(0,1))
  points(range,1-(pred_pz$fit+pred_pz$se.fit),type="l",lty=2)
  points(range,1-(pred_pz$fit-pred_pz$se.fit),type="l",lty=2)
  ats <- seq(min(range),max(range),length=6)
  axis(1,ats,labels = round(realmean+ats*2*realsd,2))
  rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  
  plot(range,pred_ct$fit,xlab=predvar2,ylab="Exp Count",type="l",lwd=2,xaxt="n",ylim=c(0,20))
  points(range,pred_ct$fit+pred_ct$se.fit,type="l",lty=2)
  points(range,pred_ct$fit-pred_ct$se.fit,type="l",lty=2)
  ats <- seq(min(range),max(range),length=6)
  axis(1,ats,labels = round(realmean+ats*2*realsd,2))
  rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  
  svg(sprintf("pdplots_%s_%s.svg",predvar2,season),width=4.5,height = 6)
  
  layout(matrix(c(1:2),nrow=2))
  par(mai=c(0.9,0.8,0.1,0.1))
  plot(range,1-pred_pz$fit,xlab=predvar2,ylab="P(1+)",type="l",lwd=2,xaxt="n",ylim=c(0,1))
  points(range,1-(pred_pz$fit+pred_pz$se.fit),type="l",lty=2)
  points(range,1-(pred_pz$fit-pred_pz$se.fit),type="l",lty=2)
  ats <- seq(min(range),max(range),length=6)
  axis(1,ats,labels = round(realmean+ats*2*realsd,2))
  rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  
  plot(range,pred_ct$fit,xlab=predvar2,ylab="Exp Count",type="l",lwd=2,xaxt="n",ylim=c(0,20))
  points(range,pred_ct$fit+pred_ct$se.fit,type="l",lty=2)
  points(range,pred_ct$fit-pred_ct$se.fit,type="l",lty=2)
  ats <- seq(min(range),max(range),length=6)
  axis(1,ats,labels = round(realmean+ats*2*realsd,2))
  rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  
  dev.off()
}

VisualizeInteraction <- function(data=df,model=bestmod2,var1=ints[[t]][1],var2=ints[[t]][2],allvars=names(df)[-c(1,6:13)],season=season){
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
  
  persp(realmean1+realsd1*2*range1,realmean2+2*realsd2*range2,1-predmat_pz,
        xlab=sprintf("\n%s",var1_2),ylab=sprintf("\n%s",var2_2),zlab = "\n\nP(1+)",
        theta = 55, phi = 40, r = sqrt(10), d = 3, 
        ticktype = "detailed", mgp = c(4, 1, 0)) -> res
  
  points(trans3d(data[[var1_2]], data[[var2_2]], max(1-predmat_pz)-0.03, pmat = res), col = gray(0.8), pch = 19, cex=0.5)
  
  persp(realmean1+realsd1*2*range1,realmean2+2*realsd2*range2,predmat_ct,
        xlab=sprintf("\n%s",var1_2),ylab=sprintf("\n%s",var2_2),zlab = "\n\nExp_cnt",
        theta = 55, phi = 40, r = sqrt(10), d = 3, 
        ticktype = "detailed", mgp = c(4, 1, 0)) -> res

  points(trans3d(data[[var1_2]], data[[var2_2]], max(predmat_ct)-0.03, pmat = res), col = gray(0.8), pch = 19, cex=0.5)
  
  dev.off() 
}

## this function takes an individual name and runs a logistic regression. If the model throws a warning due to complete separation, a
# firth-corrected (penalized likelihood) version is run.

build_indmodels <- function(thisind=indnames[2],season="rut"){
  
  formula <- as.formula(sprintf("Used ~ std.pasture + std.%s + std.ponds + std.elevation + std.ndvi + std.scrub + std.grassland + std.slope + std.aspect + std.ndvi:std.scrub + std.ndvi:std.grassland + std.%s:std.ndvi + std.slope:std.aspect",best.stocking[season],best.stocking[season]))
  
  dat=allinds2[[thisind]]
  
  full <- tryCatch( {     
    glm(formula,
        data= dat, family=binomial(link= "logit"))
  }
  , warning = function(w) { 
    
    f <- logistf::logistf(formula,
                          data= dat, family=binomial(link= "logit"))
    return(f)
  })
  
  return(full)
}


#allvar_c <- c('Used','std.pasture','std.stockingMeanNonz','std.ponds','std.elevation',
      #        'std.ndvi','std.scrub','std.grassland', 'std.slope', 'std.aspect')


VisualizeRelation_col <- function(data=allcollar[[season]],meancoefs=thisres,predvar="std.ndvi",allvars=allvar_c,season){
  len <- 100
  
  predvar2 <- strsplit(predvar,"\\.")[[1]][2]
  
  dataclasses <- sapply(data,class)
  
  dim <- data[,predvar]
  range <- seq(min(dim),max(dim),length=len)
  
  realmean <- mean(data[[predvar2]])
  realsd <- sd(data[[predvar2]])
  
  newdata <- data.frame(temp=range)
  names(newdata) <- c(predvar)
  
  othervars <- allvars[!allvars%in%c(predvar,"Used")]
  
  
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
  
  ndx <- match(names(newdata),meancoefs$names)
  
  int <- meancoefs$meancoefs[meancoefs$names=="(Intercept)"]
  coef <- meancoefs$meancoefs[meancoefs$names==predvar]
  
  pred <- sapply(1:nrow(newdata), function(t)              # make prediction plots based on population-averaged model
    plogis(int + as.matrix(newdata[t,]) %*% as.matrix(meancoefs$meancoefs[ndx]) ) )
  
  m <- c(int=int,coef=coef)
  v <- diag(c(meancoefs$secoefs[meancoefs$names=="(Intercept)"],meancoefs$secoefs[meancoefs$names==predvar]))
  
  se_pred <- sapply(1:nrow(newdata), function(t)
    as.numeric(car::deltaMethod(m,sprintf("1/(1+exp(-1*(int+%s*coef)))",newdata[t,predvar]),vcov.=v)["SE"]) )
  
  svg(sprintf("pdplots_col_%s_%s.svg",predvar2,season),width=4.5,height = 4.5)
  
  tcrit <- qnorm(0.975) #qt(0.95,nrow(meancoefs))
  
  par(mai=c(0.9,0.8,0.1,0.1))
  plot(range,pred,xlab=predvar2,ylab="Selection propensity",type="l",lwd=2,xaxt="n",ylim=c(0,min(1,max(pred)*1.4)))
  points(range,pred+tcrit*se_pred,type="l",lty=2)
  points(range,pred-tcrit*se_pred,type="l",lty=2)
  ats <- seq(min(range),max(range),length=6)
  axis(1,ats,labels = round(realmean+ats*2*realsd,2))
  rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  
  
  dev.off()
}

VisualizeInteraction_col <- function(data=df,meancoefs=thisres_f,var1="std.stockingMeanNonz",var2="std.ndvi",allvars=allvar_c,season){
  len <- 50     # increase this for higher-res figures
  
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
  
  othervars <- allvars[!allvars%in%c(standvar1,standvar2,"Used")]
  
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
  
  int <- meancoefs$meancoefs[meancoefs$names=="(Intercept)"]
  inter <- meancoefs$meancoefs[grepl(sprintf("%s:%s",var1,var2),meancoefs$names)]
  
  ndx <- match(names(newdata),meancoefs$names)
  
  pred <- sapply(1:nrow(newdata), function(s)              # make prediction plots based on population-averaged model
    plogis(int + as.matrix(newdata[s,]) %*% as.matrix(meancoefs$meancoefs[ndx]) + inter*(newdata[s,var1]*newdata[s,var2])) )
  
  predmat <-  matrix(pred,nrow=len,ncol=len)
  
  
  
  svg(sprintf("intplots_col_%s_%s_%s.svg",var1_2,var2_2,season),width=4.5,height = 6)
  
  par(mai=c(0.5,0.5,0,0))
  
  persp(realmean1+realsd1*2*range1,realmean2+realsd2*2*range2,predmat,xlab=var1_2,ylab=var2_2,theta = 55, phi = 40, r = sqrt(10), d = 3, 
        ticktype = "detailed", mgp = c(4, 1, 0)) -> res
  
  points(trans3d(data[[var1_2]], data[[var2_2]], max(predmat)-0.02, pmat = res), col = gray(0.8), pch = 19, cex=0.5)
  
  
  dev.off() 
}

#############
# Load data
#############

getwd()

#load files (count data)

alldata<-list()

alldata$rut <- readRDS("rut.obs.25Jun2019.rds")              
alldata$parturition <- readRDS("part.obs.25Jun2019.rds")
alldata$summer <- readRDS("summer.obs.25Jun2019.rds")
alldata$winter <- readRDS("winter.obs.25Jun2019.rds")

allcollar <- list()

#load data
allcollar$rut <- readRDS("all.extracts.collars.rut.25Jun2019.rds")
allcollar$parturition <- readRDS("all.extracts.collars.part.25Jun2019.rds")
allcollar$summer <- readRDS("all.extracts.collars.summer.25Jun2019.rds")
allcollar$winter <- readRDS("all.extracts.collars.winter.25Jun2019.rds")


##############
# Process data

    #replace with subset of columns of interest for appropriate season
alldata2 <- lapply(alldata,function(t) t[c(1, 2, 9, 12, 14, 16:27)])

#replace with subset of columns of interest for appropriate season 
allcollar2 <- lapply(allcollar,function(t) t[c(3:5, 9:14, 17:29)])


#standardize variables...

standardze <- function(df){
  df$std.pasture <- ((df$pasture)-(mean(df$pasture)))/(2*sd(df$pasture))
  df$std.stockingMeanNonz <- ((df$stockingMeanNonz )-(mean(df$stockingMeanNonz , na.rm=TRUE)))/(2*sd(df$stockingMeanNonz , na.rm=TRUE))
  df$std.stockingMean <- ((df$stockingMean)-(mean(df$stockingMean, na.rm=TRUE)))/(2*sd(df$stockingMean, na.rm=TRUE))
  df$std.stockingMed <- ((df$stockingMed)-(mean(df$stockingMed, na.rm=TRUE)))/(2*sd(df$stockingMed, na.rm=TRUE))
  df$std.stockingMax <- ((df$stockingMax)-(mean(df$stockingMax, na.rm=TRUE)))/(2*sd(df$stockingMax, na.rm=TRUE))
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


standardze2 <- function(df){
  df$std.pasture <- ((df$pasture)-(mean(df$pasture)))/(2*sd(df$pasture))
  df$std.stockingMeanNonz <- ((df$stockingMeanNonz )-(mean(df$stockingMeanNonz , na.rm=TRUE)))/(2*sd(df$stockingMeanNonz , na.rm=TRUE))
  df$std.stockingMean <- ((df$stockingMean)-(mean(df$stockingMean, na.rm=TRUE)))/(2*sd(df$stockingMean, na.rm=TRUE))
  df$std.stockingMed <- ((df$stockingMed)-(mean(df$stockingMed, na.rm=TRUE)))/(2*sd(df$stockingMed, na.rm=TRUE))
  df$std.stockingMax <- ((df$stockingMax)-(mean(df$stockingMax, na.rm=TRUE)))/(2*sd(df$stockingMax, na.rm=TRUE))
  df$std.elevation <- ((df$elevation)-(mean(df$elevation)))/(2*sd(df$elevation))
  df$std.ndvi <- ((df$ndvi)-(mean(df$ndvi, na.rm=TRUE)))/(2*sd(df$ndvi, na.rm=TRUE))
  df$std.ponds <- ((df$ponds)-(mean(df$ponds)))/(2*sd(df$ponds))
  df$std.scrub <- ((df$scrub)-(mean(df$scrub)))/(2*sd(df$scrub))
  df$std.grassland <- ((df$grassland)-(mean(df$grassland)))/(2*sd(df$grassland))
  df$std.slope <- ((df$slope)-(mean(df$slope)))/(2*sd(df$slope))
  df$std.aspect <- ((df$aspect)-(mean(df$aspect)))/(2*sd(df$aspect))
  return(df)
}

collar.std <- lapply(allcollar2,function(t) standardze2(t))

final.collar <- lapply(collar.std,function(t) t[complete.cases(t), ])

#subset available points to be re-combined with "used" subset for each ID above
avails <- lapply(final.collar, function(t)  t[t$sex %in% 'NA',])


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
##################
# STEP 2: CHOOSE TOP CATTLE STOCKING VARIABLE
##################

###Run through each season one at a time

# season <- "parturition"
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
# cow.max <- glmmTMB(TotalElk ~ std.pasture + std.stockingMax + std.ponds
#                    + std.elevation + std.ndvi + std.scrub + std.grassland + std.slope
#                    + std.aspect + std.ndvi:std.scrub + std.ndvi:std.grassland
#                    + std.stockingMax:std.ndvi
#                    + std.slope:std.aspect
#                    + (1|ID) + (1|year),
#                    df,
#                    ziformula = ~. ,
#                    family= truncated_nbinom1(link = "log"))
# 
# bbmle::AICtab(cow.mean.nonz, cow.mean, cow.med,weights=TRUE,mnames=c("Mean Nonz", "Mean", "Median", "Max") )

# ######Results from AIC test:
# 
# ####25 m radius below
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

######50 m radius below
# #Summer:
# dAIC df weight
# Mean Nonz  0.0 33 0.708 
# Median     2.6 33 0.194 
# Mean       4.0 33 0.098 

# #Rut:
# Mean       0.0 33 0.802 
# Mean Nonz  4.2 33 0.099 
# Median     4.2 33 0.099 

# #Winter:
# dAIC df weight
# Median     0.0 33 0.916 
# Mean       5.1 33 0.072 
# Mean Nonz  8.7 33 0.012 

# # Parturition:
# dAIC df weight
# Median     0.0 33 0.48  
# Mean       0.5 33 0.37  
# Mean Nonz  2.3 33 0.15  

stockingVars <- c("stockingMeanNonz", "stockingMean", "stockingMed", "stockingMax")
best.stocking <- c(
  rut=stockingVars[1],
  parturition=stockingVars[1],
  summer=stockingVars[1],
  winter=stockingVars[3]
)
best.stocking

###########
# STEP 3: AUTOMATE THE REMAINDER OF MODEL SELECTION AND VISUALIZATIONS
###########

############
# Loop through seasons (this will take a bit of time)
############

##finished:

formulas <- list()

allseasons <- names(final.data)
bestmods <- list()

season <- allseasons[3]
for(season in allseasons){
  df <- final.data[[season]]
  
  formulas[[season]] <- as.formula(sprintf("TotalElk ~ std.pasture + std.%s  + std.ponds + std.elevation + std.ndvi + std.scrub + std.grassland + std.slope + std.aspect + std.ndvi:std.scrub  + std.ndvi:std.grassland + std.%s:std.ndvi + std.slope:std.aspect + (1|ID) + (1|year)",best.stocking[season],best.stocking[season] ))
  
  bestmods[[season]] <- buildmer::buildglmmTMB(formulas[[season]],
             ziformula = ~.,data=df, dispformula = ~1, direction="backward",crit="AIC",
             reduce.random=F,family=truncated_nbinom1 (link="log"))
  bestmod2 <- bestmods[[season]]@model
  #summary(bestmod2)
  tmp <- rownames(coef(summary(bestmod2))$cond)[-1]
  topvars <- tmp[!grepl(":",tmp)]
  ints <- strsplit(tmp[grep(":",tmp)],":")
  if(length(topvars)>0) tmp <- sapply(1:length(topvars), function(t) VisualizeRelation(df,bestmod2,topvars[t],season=season) )     # run and save partial dependence plots for all top variables
  
  if(length(ints)>0 ) tmp <- sapply(1:length(ints), function(t) VisualizeInteraction(df,bestmod2,ints[[t]][1],ints[[t]][2],season=season) )
  
}

#######
# .. and for collar data

allseasons <- names(final.collar)

allmeans <- list()

thisres <- list()

season <- allseasons[2]

for(season in allseasons){
  
  df <- final.collar[[season]]
  
  
  table(df$ID)
  
  ############
  #subset by ID for indvidual level analysis
  ############
  
  indnames <- c('24034B','31710A','31710B','31711A','31711B','31713A','31713B','31739A')
  indsex <- rep(c("M","F"),times=c(3,5))
  names(indsex) <- indnames
  
  allinds <- list()
  test <- lapply(indnames,function(t) allinds[[t]] <<- subset(df,ID==t) )
  
  
  # remove individuals with too few observations
  indnames <- indnames[sapply(allinds,nrow)>100]
  allinds <- allinds[indnames]
  
  nboot <- 100
  coefs <- list()
  
  ### loop through boot samples
  b=1
  for(b in 1:nboot){
    nobs <- sapply(allinds,nrow)
    
    nbackground <- nrow(avails[[season]])
    t <- indnames[1]
    avails2 <- lapply(indnames, function(t) 
      avails[[season]][sample(1:nrow(avails[[season]]),
      size=min(nrow(avails[[season]]),nobs[t]*5),replace=TRUE),])  # scale the number of background points to the number of observations
    names(avails2) <- indnames
    
    allinds2 <- lapply(indnames,function(t) rbind(allinds[[t]][sample(1:nrow(allinds[[t]]),nrow(allinds[[t]]),
                  replace = TRUE),], avails2[[t]]) )
    names(allinds2) <- indnames
    
    allmods <- lapply(indnames,build_indmodels,season=season)    # run models
    names(allmods) <- indnames
    
    #summary(allmods[[1]])
    
    coefs[[b]] <- sapply(allmods,coefficients)    
  }
  
  allmeans[[season]] <- sapply(coefs,function(t)  apply(t,1,mean))  
  thisres[[season]] = data.frame(
    meancoefs = apply(allmeans[[season]],1,mean),
    secoefs = apply(allmeans[[season]],1,sd),                       
    uppr = apply(allmeans[[season]],1,function(t) quantile(t,0.95)),
    lowr = apply(allmeans[[season]],1,function(t) quantile(t,0.05)),     
    names=rownames(allmeans[[season]])
  )
  
  # svg(sprintf("collar_coefs_all_%s.svg",season),width = 4,height=5)
  # 
  # p<- ggplot(thisres,aes(x=meancoefs, y=names)) + 
  #   geom_point() +
  #   ggtitle(season) +
  #   ylab("") +
  #   xlab("Mean coefficient estimate") +
  #   geom_vline(xintercept = 0,linetype=1, size=1) +
  #   geom_errorbarh(aes(xmin=lowr, xmax=uppr),height=0.2)
  # print(p)
  # # ggsave(p, filename = paste("coef.plot", season, "pdf", sep = "."), height = 11, width = 8)
  # 
  # dev.off()
  
  #save coefficients table
  #  write.csv(thisres, file = paste("coef.table", season, "csv", sep = "."))
  
  #save AIC
  #  sink(file = paste("aic", season, "txt", sep = "."))
  
  #  AIC(mdl)
  #  sink()
  
  
#   ###########
#   # females only            # KTS: do we need to keep this- didn't seem like it was showing anything too interesting, and it's not as analogous with the count models this way
#   ###########
#   
#   ndx <-  intersect(indnames, names(which(indsex=="F")))
#   
#   allmeans2 <- sapply(coefs,function(t)  apply(t[,ndx],1,mean))
#   
#   # thisres_f <- data.frame(
#   #   meancoefs=sapply(1:nrow(coefs),function(t) wtd.mean(coefs[t,][ndx],1/(ses[t,][ndx])^2)),
#   #   secoefs=sapply(1:nrow(coefs),function(t) sqrt(wtd.var(coefs[t,][ndx],1/(ses[t,][ndx])^2))/sqrt(length(indnames[ndx]))),
#   #   tcrit=qt(0.95,df=length(indnames)),
#   #   names=rownames(coefs)
#   # # )
#   # thisres$up <- thisres$meancoefs+thisres$secoefs*thisres$tcrit
#   # thisres$low <- thisres$meancoefs-thisres$secoefs*thisres$tcrit
#   
#   thisres_f = data.frame(
#     meancoefs = apply(allmeans2,1,mean),
#     secoefs = apply(allmeans,1,sd),
#     uppr = apply(allmeans2,1,function(t) quantile(t,0.95)),
#     lowr = apply(allmeans2,1,function(t) quantile(t,0.05)),
#     names=rownames(allmeans2)
#   )
#   
#   svg(sprintf("collar_coefs_fem_%s.svg",season),width = 4,height=5)
#   
#   p<- ggplot(thisres_f,aes(x=meancoefs, y=names)) + 
#     geom_point() +
#     ggtitle(season) +
#     ylab("") +
#     xlab("Mean coefficient estimate") +
#     geom_vline(xintercept = 0,linetype=1, size=1) +
#     geom_errorbarh(aes(xmin=lowr, xmax=uppr),height=0.2)
#   print(p)
#   
#   dev.off()
#   
  tmp <- as.character(thisres$names[grepl("std",thisres$names)])
  allmains <- tmp[!grepl(":",tmp)]
  allints <- tmp[grepl(":",tmp)]
  tmp <- sapply(1:length(allmains), function(t) VisualizeRelation_col(data=df,meancoefs=thisres,predvar = allmains[t],allvars=allmains,season=season) )     # run and save partial dependence plots for all top variables

  ints <- strsplit(allints[grep(":",allints)],":")
  tmp <- sapply(1:length(ints), function(t) VisualizeInteraction_col(data=df,meancoefs=thisres,var1=ints[[t]][1],var2=ints[[t]][2],allvars=allmains,season=season) )
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
names(allcoefsz) <- allseasons

allsesz <- lapply(1:length(temp_zero), function(t) temp_zero[[t]][,'Std. Error'][coef_ndx[[t]]])
names(allsesz) <- allseasons

best.stocking

varlist <- data.frame(
  var = as.character(thisres[[1]]$names)[-1],
  readable = c("Pasture","Stocking","Ponds","Elevation","NDVI","Scrub","Grassland","Slope","Cos Aspect",
               "NDVI*Scrub","NDVI*Grassland","NDVI*Stocking","Slope*Aspect"),
  stringsAsFactors = F
) 

varlist2 <- list()

temp <- lapply(allseasons,function(t) varlist2[[t]] <<- varlist)

temp <- lapply(allseasons,function(t) varlist2[[t]]$var[2] <<- sprintf("std.%s",best.stocking[t]) )
temp <- lapply(allseasons,function(t) varlist2[[t]]$var[12] <<- sprintf("std.%s:std.ndvi",best.stocking[t]) )

varlist2

graphics.off()

svg("coefplot_final.svg",7.5,4.5)

layout(matrix(1:12,nrow=3,byrow = T),widths = c(1,2,2,2))          
var_ndx <- 1
onleft <- c(1,5,9)
for(panel in c(1:12)){
  par(mai=c(0.1,0.1,0.5,0.01))
  if(panel %in% onleft){
    plot(20,20,xlim=c(1,4),ylim=c(-0.5,16),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    text(rep(1,4),ycents, tools::toTitleCase(allseasons),adj=0,cex=1.3)
  }else{
    plot(20,20,xlim=c(8,25),ylim=c(-0.5,16),xlab="",ylab="",xaxt="n",yaxt="n",main=varlist2$rut$readable[var_ndx],pch="")    # blank plot
    xcents <- c(11,17,23)
    #abline(v=xcents,lwd=2)
    ycents <- seq(15,3,length=4)
    yshift <- 0
    
    polygon(c(xcents[1]-2,xcents[1]+2,xcents[1]+2,xcents[1]-2),c(rep(range(ycents)[1],2)-1-yshift,
                                                                 rep(range(ycents)[2],2)+1-yshift),border=NA,col=gray(0.8))
    polygon(c(xcents[2]-2,xcents[2]+2,xcents[2]+2,xcents[2]-2),c(rep(range(ycents)[1],2)-1+yshift,
                                                                 rep(range(ycents)[2],2)+1+yshift),border=NA,col=gray(0.8))
    polygon(c(xcents[3]-2,xcents[3]+2,xcents[3]+2,xcents[3]-2),c(rep(range(ycents)[1],2)-1-yshift,
                                                                 rep(range(ycents)[2],2)+1-yshift),border=NA,col=gray(0.8))
    
    segments(xcents,rep(range(ycents)[1],3)-1+c(-yshift,yshift,-yshift),xcents,rep(range(ycents)[2],3)+1+c(-yshift,yshift,-yshift),lwd=1)
    segments(xcents[2],range(ycents)[1]-1+yshift,xcents[2],range(ycents)[2]+1+yshift,lwd=2,col=gray(0.4))
    
    
    
    meanvalz <- sapply(allseasons,function(t) thisres[[t]]$meancoefs[as.character(thisres[[t]]$names)==varlist2[[t]]$var[var_ndx]] )
    points(meanvalz+xcents[1],ycents-yshift,pch=20,cex=1.5)
    lower <- sapply(allseasons,function(t) thisres[[t]]$lowr[as.character(thisres[[t]]$names)==varlist2[[t]]$var[var_ndx]] )
    upper <- sapply(allseasons,function(t) thisres[[t]]$uppr[as.character(thisres[[t]]$names)==varlist2[[t]]$var[var_ndx]] )
    arrows(lower+xcents[1],ycents-yshift,upper+xcents[1],ycents-yshift,angle = 90,code=3,length=0.05)
    
    points(sapply(allseasons,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])+xcents[2],ycents+yshift,pch=20,cex=1.5,col=gray(0.4))
    arrows(sapply(allseasons,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])+xcents[2]-1.66*sapply(allseasons,function(t) allses[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents+yshift,sapply(allseasons,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])+xcents[2]+1.66*sapply(allseasons,function(t) allses[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents+yshift,angle = 90,code=3,length=0.05,col=gray(0.4))
    
    points((sapply(allseasons,function(t) allcoefsz[[t]][varlist2[[t]]$var[var_ndx]]))*-1+xcents[3],ycents-yshift,pch=20,cex=1.5)
    arrows((sapply(allseasons,function(t) allcoefsz[[t]][varlist2[[t]]$var[var_ndx]]))*-1+xcents[3]-1.66*sapply(allseasons,function(t) allsesz[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents-yshift,(sapply(allseasons,function(t) allcoefsz[[t]][varlist2[[t]]$var[var_ndx]]))*-1+xcents[3]+1.66*sapply(allseasons,function(t) allsesz[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents-yshift,angle = 90,code=3,length=0.05)
    
    text(xcents,c(0,0,0),c("Collar","Count","P(1+)"),cex=1.2,col=c("black",gray(0.4),"black"))
    abline(h=c(ycents[1]+2,ycents[-length(ycents)]-2),lty=2,lwd=0.5)
    var_ndx <- var_ndx+1
  }
  
  
}

dev.off()

#### TODO: update this code to look like the univariate coefficients plot... 

svg("intcoefs_final.svg",6,3.5)


layout(matrix(1:6,nrow=2,byrow = T),widths = c(0.7,2,2))          
var_ndx <- 10
onleft <- c(1,4)
for(panel in c(1:6)){
  par(mai=c(0.1,0.1,0.5,0.01))
  if(panel %in% onleft){
    plot(20,20,xlim=c(1,4),ylim=c(-0.5,16),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    text(rep(1,4),ycents, tools::toTitleCase(allseasons),adj=0,cex=1.3)
  }else{
    plot(20,20,xlim=c(9,32),ylim=c(-0.5,16),xlab="",ylab="",xaxt="n",yaxt="n",main=varlist2$rut$readable[var_ndx],pch="")    # blank plot
    xcents <- c(14,22,30)
    #abline(v=xcents,lwd=2)
    ycents <- seq(15,3,length=4)
    yshift <- 0
    
    polygon(c(xcents[1]-2,xcents[1]+2,xcents[1]+2,xcents[1]-2),c(rep(range(ycents)[1],2)-1-yshift,
                                                                 rep(range(ycents)[2],2)+1-yshift),border=NA,col=gray(0.8))
    polygon(c(xcents[2]-2,xcents[2]+2,xcents[2]+2,xcents[2]-2),c(rep(range(ycents)[1],2)-1+yshift,
                                                                 rep(range(ycents)[2],2)+1+yshift),border=NA,col=gray(0.8))
    polygon(c(xcents[3]-2,xcents[3]+2,xcents[3]+2,xcents[3]-2),c(rep(range(ycents)[1],2)-1-yshift,
                                                                 rep(range(ycents)[2],2)+1-yshift),border=NA,col=gray(0.8))
    
    segments(xcents,rep(range(ycents)[1],3)-1+c(-yshift,yshift,-yshift),xcents,rep(range(ycents)[2],3)+1+c(-yshift,yshift,-yshift),lwd=1)
    segments(xcents[2],range(ycents)[1]-1+yshift,xcents[2],range(ycents)[2]+1+yshift,lwd=2,col=gray(0.4))
    
    
    
    meanvalz <- sapply(allseasons,function(t) thisres[[t]]$meancoefs[as.character(thisres[[t]]$names)==varlist2[[t]]$var[var_ndx]] )
    points(meanvalz+xcents[1],ycents-yshift,pch=20,cex=1.5)
    lower <- sapply(allseasons,function(t) thisres[[t]]$lowr[as.character(thisres[[t]]$names)==varlist2[[t]]$var[var_ndx]] )
    upper <- sapply(allseasons,function(t) thisres[[t]]$uppr[as.character(thisres[[t]]$names)==varlist2[[t]]$var[var_ndx]] )
    arrows(lower+xcents[1],ycents-yshift,upper+xcents[1],ycents-yshift,angle = 90,code=3,length=0.05)
    
    if(var_ndx==12){
      points(sapply(allseasons,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])/c(5,1,1,1)+xcents[2],ycents+yshift,pch=20,cex=1.5,col=gray(0.4))
      arrows(sapply(allseasons,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])/c(5,1,1,1)+xcents[2]-1.66*sapply(allseasons,function(t) allses[[t]][varlist2[[t]]$var[var_ndx]])/c(5,1,1,1),
             ycents+yshift,sapply(allseasons,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])/c(5,1,1,1)+xcents[2]+1.66*sapply(allseasons,function(t) allses[[t]][varlist2[[t]]$var[var_ndx]])/c(5,1,1,1),
             ycents+yshift,angle = 90,code=3,length=0.05,col=gray(0.4))
      text(xcents[2]-2,ycents[1]-1,"*",cex=1.8)
    }else{
      points(sapply(allseasons,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])+xcents[2],ycents+yshift,pch=20,cex=1.5,col=gray(0.4))
      arrows(sapply(allseasons,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])+xcents[2]-1.66*sapply(allseasons,function(t) allses[[t]][varlist2[[t]]$var[var_ndx]]),
             ycents+yshift,sapply(allseasons,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])+xcents[2]+1.66*sapply(allseasons,function(t) allses[[t]][varlist2[[t]]$var[var_ndx]]),
             ycents+yshift,angle = 90,code=3,length=0.05,col=gray(0.4))
    }
    
    
    points((sapply(allseasons,function(t) allcoefsz[[t]][varlist2[[t]]$var[var_ndx]]))*-1+xcents[3],ycents-yshift,pch=20,cex=1.5)
    arrows((sapply(allseasons,function(t) allcoefsz[[t]][varlist2[[t]]$var[var_ndx]]))*-1+xcents[3]-1.66*sapply(allseasons,function(t) allsesz[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents-yshift,(sapply(allseasons,function(t) allcoefsz[[t]][varlist2[[t]]$var[var_ndx]]))*-1+xcents[3]+1.66*sapply(allseasons,function(t) allsesz[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents-yshift,angle = 90,code=3,length=0.05)
    
    text(xcents,c(0,0,0),c("Collar","Count","P(1+)"),cex=1.2,col=c("black",gray(0.4),"black"))
    abline(h=c(ycents[1]+2,ycents[-length(ycents)]-2),lty=2,lwd=0.5)
    var_ndx <- var_ndx+1
  }
  
  
}

dev.off()

graphics.off()

############
# End of script
############




########
# TODO
########

# Add visualization of coefficients for both the collar, zero and count models- for each season
# make sure that I convert Pzero to P(1+)









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



###########################################################################
##IF JUST RE-RUNNING WITH EXISTING RDATA FILE TO CHANGE FIGURES, START HERE
###########################################################################
rm(list=ls())

#setwd ("/Users/Lacey/Box Sync/PhD_Project/R.Analyses/Elk/Raw data")

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
# Load data
#############
load("/Users/Lacey/Box Sync/PhD_Project/PORE Project/Manuscripts/Results_ms/Figures/Summary plots/July2/27Jun2019_collars_surveys.RData")








