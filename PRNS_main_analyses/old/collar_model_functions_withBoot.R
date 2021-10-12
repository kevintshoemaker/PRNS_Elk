#GLM's for PRNS elk project (Collar data). Created by Kevin Shoemaker 18 May, 2019. Updated by Lacey Hughey May 20, 2019


###########
# Clear workspace
###########

rm(list=ls())        #clear workspace if needed


###########
# set working directory
###########

setwd ("/Users/Lacey/Box Sync/PhD_Project/R.Analyses/Elk/Raw data")


############
# Load packages
############

library(Hmisc)
library(ggplot2)

############
# Custom functions
############


## this function takes an individual name and runs a logistic regression. If the model throws a warning due to complete separation, a
# firth-corrected (penalized likelihood) version is run.

build_indmodels <- function(thisind=indnames[2]){
  dat=allinds2[[thisind]]
  
  full <- tryCatch( {     
    glm(Used ~ std.pasture + std.stockingMeanNonz + std.ponds + std.elevation +
          std.ndvi + std.scrub + std.grassland + std.slope + std.aspect +
          std.ndvi:std.scrub + std.ndvi:std.grassland + std.stockingMeanNonz:std.ndvi + 
          std.slope:std.aspect,
        data= dat, family=binomial(link= "logit"))
  }
  , warning = function(w) { 
    
    f <- logistf::logistf(Used ~ std.pasture + std.stockingMeanNonz + std.ponds + std.elevation +
                            std.ndvi + std.scrub + std.grassland + std.slope + std.aspect +
                            std.ndvi:std.scrub + std.ndvi:std.grassland + std.stockingMeanNonz:std.ndvi + 
                            std.slope:std.aspect,
                          data= dat, family=binomial(link= "logit"))
    return(f)
  })
  
  return(full)
}


allvar_c <- c('Used','std.pasture','std.stockingMeanNonz','std.ponds','std.elevation',
              'std.ndvi','std.scrub','std.grassland', 'std.slope', 'std.aspect')


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
  axis(1,ats,labels = round(realmean+ats*realsd,2))
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
  
  persp(realmean1+realsd1*range1,realmean2+realsd2*range2,predmat,xlab=var1_2,ylab=var2_2,theta = 55, phi = 40, r = sqrt(10), d = 3, 
        ticktype = "detailed", mgp = c(4, 1, 0))
  
  
  dev.off() 
}


############
# Load data
############

allcollar <- list()

#load data
allcollar$rut <- readRDS("all.extracts.rut.27May2019.rds")
allcollar$part <- readRDS("all.extracts.part.27May2019.rds")
allcollar$summer <- readRDS("all.extracts.summer.27May2019.rds")
allcollar$winter <- readRDS("all.extracts.winter.27May2019.rds")

##############
# Process data

#replace with subset of columns of interest for appropriate season 
allcollar2 <- lapply(allcollar,function(t) t[c(3:5, 9:14, 17:28)])

#standardize variables...

standardze2 <- function(df){
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

collar.std <- lapply(allcollar2,function(t) standardze2(t))

final.collar <- lapply(collar.std,function(t) t[complete.cases(t), ])

#subset available points to be re-combined with "used" subset for each ID above
avails <- lapply(final.collar, function(t)  t[t$sex %in% 'NA',])

############
# Loop through seasons
############

allseasons <- names(final.collar)

season <- allseasons[3]

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
  
  for(b in 1:nboot){
    nobs <- sapply(allinds,nrow)
    
    nbackground <- nrow(avails[[season]])
    t <- indnames[1]
    avails2 <- lapply(indnames, function(t) 
      avails[[season]][sample(1:nrow(avails[[season]]),size=min(nrow(avails[[season]]),nobs[t]*5),replace=TRUE),])  # scale the number of background points to the number of observations
    names(avails2) <- indnames
    
    allinds2 <- lapply(indnames,function(t) rbind(allinds[[t]][sample(1:nrow(allinds[[t]]),nrow(allinds[[t]]),replace = TRUE),], avails2[[t]]) )
    names(allinds2) <- indnames
    
    allmods <- lapply(indnames,build_indmodels)    # run models
    names(allmods) <- indnames
    
    summary(allmods[[1]])
    
    coefs[[b]] <- sapply(allmods,coefficients)
  }
  
  
  # mdl=allmods[[1]]
  # sefunc <- function(mdl){
  #   sqrt(diag(vcov(mdl)))    # extract standard error
  # }
  # ses <- sapply(allmods,sefunc) #coef(summary(t))[, 2])
  # rownames(ses) <- rownames(coefs)
  # thisres <- data.frame(
  #   meancoefs=sapply(1:nrow(coefs),function(t) wtd.mean(coefs[t,],1/(ses[t,])^2)),
  #   secoefs=sapply(1:nrow(coefs),function(t) sqrt(wtd.var(coefs[t,],1/(ses[t,])^2))/sqrt(length(indnames))),
  #   tcrit=qt(0.95,df=length(indnames)),
  #   names=rownames(coefs)
  # )
  # thisres$up <- thisres$meancoefs+thisres$secoefs*thisres$tcrit
  # thisres$low <- thisres$meancoefs-thisres$secoefs*thisres$tcrit
  allmeans <- sapply(coefs,function(t)  apply(t,1,mean))
  thisres = data.frame(
    meancoefs = apply(allmeans,1,mean),
    secoefs = apply(allmeans,1,sd),
    uppr = apply(allmeans,1,function(t) quantile(t,0.95)),
    lowr = apply(allmeans,1,function(t) quantile(t,0.05)),
    names=rownames(allmeans)
  )
  
  svg(sprintf("collar_coefs_all_%s.svg",season),width = 4,height=5)
  
  p<- ggplot(thisres,aes(x=meancoefs, y=names)) + 
    geom_point() +
    ggtitle(season) +
    ylab("") +
    xlab("Mean coefficient estimate") +
    geom_vline(xintercept = 0,linetype=1, size=1) +
    geom_errorbarh(aes(xmin=lowr, xmax=uppr),height=0.2)
  print(p)
 # ggsave(p, filename = paste("coef.plot", season, "pdf", sep = "."), height = 11, width = 8)
  
  dev.off()
  
  #save coefficients table
#  write.csv(thisres, file = paste("coef.table", season, "csv", sep = "."))
  
  #save AIC
#  sink(file = paste("aic", season, "txt", sep = "."))
  
#  AIC(mdl)
#  sink()
  
  
  ###########
  # females only
  ###########
  
  ndx <-  intersect(indnames, names(which(indsex=="F")))
  
  allmeans2 <- sapply(coefs,function(t)  apply(t[,ndx],1,mean))
  
  # thisres_f <- data.frame(
  #   meancoefs=sapply(1:nrow(coefs),function(t) wtd.mean(coefs[t,][ndx],1/(ses[t,][ndx])^2)),
  #   secoefs=sapply(1:nrow(coefs),function(t) sqrt(wtd.var(coefs[t,][ndx],1/(ses[t,][ndx])^2))/sqrt(length(indnames[ndx]))),
  #   tcrit=qt(0.95,df=length(indnames)),
  #   names=rownames(coefs)
  # # )
  # thisres$up <- thisres$meancoefs+thisres$secoefs*thisres$tcrit
  # thisres$low <- thisres$meancoefs-thisres$secoefs*thisres$tcrit
  
  thisres_f = data.frame(
    meancoefs = apply(allmeans2,1,mean),
    secoefs = apply(allmeans,1,sd),
    uppr = apply(allmeans2,1,function(t) quantile(t,0.95)),
    lowr = apply(allmeans2,1,function(t) quantile(t,0.05)),
    names=rownames(allmeans2)
  )
  
  svg(sprintf("collar_coefs_fem_%s.svg",season),width = 4,height=5)
  
  p<- ggplot(thisres_f,aes(x=meancoefs, y=names)) + 
    geom_point() +
    ggtitle(season) +
    ylab("") +
    xlab("Mean coefficient estimate") +
    geom_vline(xintercept = 0,linetype=1, size=1) +
    geom_errorbarh(aes(xmin=lowr, xmax=uppr),height=0.2)
  print(p)
  
  dev.off()
  
  tmp <- as.character(thisres$names[grepl("std",thisres$names)])
  allmains <- tmp[!grepl(":",tmp)]
  allints <- tmp[grepl(":",tmp)]
  tmp <- sapply(1:length(allmains), function(t) VisualizeRelation_col(data=df,meancoefs=thisres,predvar = allmains[t],allvars=allmains,season=season) )     # run and save partial dependence plots for all top variables
  
  ints <- strsplit(allints[grep(":",allints)],":")
  tmp <- sapply(1:length(ints), function(t) VisualizeInteraction_col(data=df,meancoefs=thisres,var1=ints[[t]][1],var2=ints[[t]][2],allvars=allmains,season=season) )
}


#######
# run gof test
#######


gof <- DHARMa::simulateResiduals(allmods[[1]])
plot(gof)
DHARMa::testResiduals(gof)



############
# TODO
###########

# visualize relationships

plot(p)
