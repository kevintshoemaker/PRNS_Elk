#Models for PRNS elk project (Count survey data and collar data). 
# Created by Kevin Shoemaker and edited by Lacey Hughey on May 20, 2019.


###########
# TODO: 
############

# rerun collar data models
# run again to enable prediction??? [why can't I predict from these models?]
# make univariate plot and interaction figures




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
#Load previous workspace (if just plotting)
#load("15Jul2019_workspace_all.RData")

#############
# Load custom functions
#############

VisualizeRelation <- function(data=df,model=full,predvar="std.slope",varname="Slope",allvars=names(df)[-c(1,6:13)],season,svg=F){
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
  
  newdata$ID <- NA # data$ID[1]
  newdata$year <- NA
  
  pred_pz <- predict(model,newdata,type="zprob",se.fit=T)
  pred_ct <- predict(model,newdata,type="conditional",se.fit=T)
  #pred <- predict(model,newdata,type="response",se.fit=T)
  
  if(svg){
    layout(matrix(c(1:2),nrow=2))
    par(mai=c(0.9,0.8,0.1,0.1))
  }else{
      
    plot(range,pred_ct$fit,xlab=varname,ylab="Exp Count",type="l",lwd=2,xaxt="n",ylim=c(0,20))
    points(range,pred_ct$fit+pred_ct$se.fit,type="l",lty=2)
    points(range,pred_ct$fit-pred_ct$se.fit,type="l",lty=2)
    ats <- seq(min(range),max(range),length=6)
    axis(1,ats,labels = round(realmean+ats*2*realsd,2))
    rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  
    plot(range,1-pred_pz$fit,xlab=varname,ylab="P(1+)",type="l",lwd=2,xaxt="n",ylim=c(0,0.1))
    points(range,1-(pred_pz$fit+pred_pz$se.fit),type="l",lty=2)
    points(range,1-(pred_pz$fit-pred_pz$se.fit),type="l",lty=2)
    ats <- seq(min(range),max(range),length=6)
    axis(1,ats,labels = round(realmean+ats*2*realsd,2))
    rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
    
  } 
  
  if(svg){
    svg(sprintf("pdplots_%s_%s.svg",predvar2,season),width=4.5,height = 6)
    
    layout(matrix(c(1:2),nrow=2))
    par(mai=c(0.9,0.8,0.1,0.1))
    plot(range,1-pred_pz$fit,xlab=predvar2,ylab="P(1+)",type="l",lwd=2,xaxt="n",ylim=c(0,.1))
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
  
}

VisualizeInteraction <- function(data=final.data$summer,model=bestmods$summer,var2="std.stockingMeanNonz",var1="std.ndvi",allvars=names(df)[-c(1,6:13)],season="summer",svg=F){
  len <- 20     # increase this for higher-res figures
  
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
  
  newdata$ID <- NA # data$ID[1]
  newdata$year <- NA
  
  pred_pz <- predict(model,newdata,type="zprob")
  pred_ct <- predict(model,newdata,type="conditional")
  
  predmat_pz <-  1-matrix(pred_pz,nrow=len,ncol=len)
  predmat_ct <-  matrix(pred_ct,nrow=len,ncol=len)
  
  
  
  if(svg){
    svg(sprintf("intplots_%s_%s_%s.svg",var1_2,var2_2,season),width=4.5,height = 6)
    
    layout(matrix(c(1:2),nrow=2))
    par(mai=c(0,0,0,0))
  }else{
    var1_2 <- varlist2[[season]]$readable[which(varlist2[[season]]$var==var1)] #strsplit(var1,"\\.")[[1]][2]
    var2_2 <- varlist2[[season]]$readable[which(varlist2[[season]]$var==var2)] #strsplit(var2,"\\.")[[1]][2]
  }
  
  x.axis <- realmean1+realsd1*2*range1
  min.x <- min(x.axis)
  max.x <- max(x.axis)
  xrange <- max.x-min.x
  x.axis2 <- seq(min.x,max.x,length=5)[c(2:4)]
  y.axis <- realmean2+realsd2*2*range2
  min.y <- min(y.axis)
  max.y <- max(y.axis)
  yrange <- max.y-min.y
  y.axis2 <- seq(min.y,max.y,length=5)[c(2:4)]
  z.axis <- seq(min(predmat_ct), min(20,max(predmat_ct)), length=5)[c(2:4)]
  min.z <- min(predmat_ct)
  max.z <- min(20,max(predmat_ct))
  zrange <- max.z-min.z
  predmat_ct[which(predmat_ct>max.z,arr.ind = T)] <- max.z
  xlab=paste("",var1_2,sep="")
  ylab=paste("",var2_2,sep="")
  zlab="Count"
  persp(x.axis,y.axis,predmat_ct,zlim=c(min.z,max.z),
        axes=F,
        theta = 55, phi = 40, r = sqrt(10), d = 3, 
        ticktype = "detailed", mgp = c(4, 1, 0)) -> res
  # lines(trans3d(x.axis, min.y, min.z, res) , col="black",lwd=2)
  # lines(trans3d(max.x, y.axis, min.z, res) , col="black",lwd=2)
  # lines(trans3d(min.x, min.y, z.axis, res) , col="black",lwd=2)
  
  labels <- sprintf("%.1f",x.axis2)
  label.pos <- trans3d(x.axis2, (min.y - yrange*0.3), min.z, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=30, cex=1)
  
  labels <- xlab
  label.pos <- trans3d(x.axis2[1]-xrange*0.3, (min.y - yrange*0.4), min.z, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=30+100+180, cex=1)
  
  labels <- sprintf("%.1f",y.axis2)
  label.pos <- trans3d((max.x + xrange*0.15), y.axis2, min.z, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), cex=1) 
  
  labels <- ylab
  label.pos <- trans3d((max.x + xrange*0.4), y.axis2[2], min.z, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), cex=1,srt=25)
  
  labels <- sprintf("%.0f",z.axis)
  label.pos <- trans3d(min.x, (min.y - yrange*0.1), z.axis, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(1, NA), cex=1)
  
  labels <- zlab
  label.pos <- trans3d((min.x), (min.y -5.3), z.axis[3], res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(1, NA), cex=1,srt=100)
  
  
  var1_2 <- strsplit(var1,"\\.")[[1]][2]
  var2_2 <- strsplit(var2,"\\.")[[1]][2]
  subst <- data[sample(1:nrow(data),1000,replace = F),]
  points(trans3d(subst[[var1_2]], subst[[var2_2]], max(predmat_ct)+0, pmat = res), col = gray(0.8), pch = 19, cex=0.5)
  
  if(!svg){
    var1_2 <- varlist2[[season]]$readable[which(varlist2[[season]]$var==var1)] #strsplit(var1,"\\.")[[1]][2]
    var2_2 <- varlist2[[season]]$readable[which(varlist2[[season]]$var==var2)] #strsplit(var2,"\\.")[[1]][2]
  }
  
  x.axis <- realmean1+realsd1*2*range1
  min.x <- min(x.axis)
  max.x <- max(x.axis)
  xrange <- max.x-min.x
  x.axis2 <- seq(min.x,max.x,length=5)[c(2:4)]
  y.axis <- realmean2+realsd2*2*range2
  min.y <- min(y.axis)
  max.y <- max(y.axis)
  yrange <- max.y-min.y
  y.axis2 <- seq(min.y,max.y,length=5)[c(2:4)]
  z.axis <- seq(min(predmat_pz), max(predmat_pz), length=5)[c(2:4)]
  min.z <- min(predmat_pz)
  max.z <-max(predmat_pz)
  zrange <- max.z-min.z
  xlab=paste("\n\n",var1_2,sep="")
  ylab=paste("\n\n",var2_2,sep="")
  zlab="\n\nP(1+)"
  persp(x.axis,y.axis,predmat_pz,zlim=c(min.z,max.z),
        axes=F,
        theta = 55, phi = 40, r = sqrt(10), d = 3, 
        ticktype = "detailed", mgp = c(4, 1, 0)) -> res
  # lines(trans3d(x.axis, min.y, min.z, res) , col="black",lwd=2)
  # lines(trans3d(max.x, y.axis, min.z, res) , col="black",lwd=2)
  # lines(trans3d(min.x, min.y, z.axis, res) , col="black",lwd=2)
  labels <- sprintf("%.1f",x.axis2)
  label.pos <- trans3d(x.axis2, (min.y - yrange*0.3), min.z, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=30, cex=1)
  
  labels <- xlab
  label.pos <- trans3d(x.axis2[1]-xrange*0.3, (min.y - yrange*0.25), min.z, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=30+100+180, cex=1)
  
  labels <- sprintf("%.1f",y.axis2)
  label.pos <- trans3d((max.x + xrange*0.12), y.axis2, min.z, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), cex=1) 
  
  labels <- ylab
  label.pos <- trans3d((max.x + xrange*0.2), y.axis2[3]-yrange*0.1, min.z, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), cex=1,srt=25)
  
  labels <- sprintf("%.2f",z.axis)
  label.pos <- trans3d(min.x, (min.y - yrange*0.1), z.axis, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(1, NA), cex=1)
  
  labels <- zlab
  label.pos <- trans3d((min.x), (min.y -5.3), z.axis[3], res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(1, NA), cex=1,srt=100)
  
  var1_2 <- strsplit(var1,"\\.")[[1]][2]
  var2_2 <- strsplit(var2,"\\.")[[1]][2]
 # subst <- data[sample(1:nrow(data),1000,replace = F),]
  points(trans3d(subst[[var1_2]], subst[[var2_2]], max(predmat_pz)-0, pmat = res), col = gray(0.8), pch = 19, cex=0.5)
  
  if(svg) dev.off() 
}


#allvar_c <- c('Used','std.pasture','std.stockingMeanNonz','std.ponds','std.elevation',
      #        'std.ndvi','std.scrub','std.grassland', 'std.slope', 'std.aspect')


VisualizeRelation_col <- function(data=df,model=full,predvar="std.slope",varname="Slope",allvars=names(df)[-c(1,6:13)],season,svg=F){
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
  
  newdata$ID <- NA #data$ID[1]
  newdata$weights = 1
  
  pred <- predict(model,newdata,type="conditional",se.fit=T)
  
  if(svg){
    layout(matrix(c(1:2),nrow=2))
    par(mai=c(0.9,0.8,0.1,0.1))
  }else{
    plot(range,pred$fit*1000,xlab=varname,ylab="Sel. Intensity",type="l",lwd=2,xaxt="n",ylim=c(0,max(pred$fit)*1000+max(pred$se.fit*1000)))
    points(range,(pred$fit*1000+pred$se.fit*1000),type="l",lty=2)
    points(range,(pred$fit*1000-pred$se.fit*1000),type="l",lty=2)
    ats <- seq(min(range),max(range),length=6)
    axis(1,ats,labels = round(realmean+ats*2*realsd,2))
    rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  } 
  
  if(svg){
    svg(sprintf("pdplots_%s_%s.svg",predvar2,season),width=4.5,height = 6)
    
    layout(matrix(c(1:2),nrow=2))
    par(mai=c(0.9,0.8,0.1,0.1))
    plot(range,pred$fit*1000,xlab=varname,ylab="Sel. Intensity",type="l",lwd=2,xaxt="n",ylim=c(0,max(pred$fit)*1000+max(pred$se.fit*1000)))
    points(range,(pred$fit*1000+pred$se.fit*1000),type="l",lty=2)
    points(range,(pred$fit*1000-pred$se.fit*1000),type="l",lty=2)
    ats <- seq(min(range),max(range),length=6)
    axis(1,ats,labels = round(realmean+ats*2*realsd,2))
    rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
    
    dev.off()
  }
  
}


VisualizeInteraction_col <- function(data=final.data$mating,model=bestmods$mating,var1=ints[[t]][1],var2=ints[[t]][2],allvars=names(df)[-c(1,6:13)],season=season,svg=T){
  len <- 20     # increase this for higher-res figures
  
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
  
  newdata$ID <- NA # data$ID[1]
  newdata$weights <- 1
  
  pred  <- predict(model,newdata,type="response",allow.new.levels=TRUE)*1000
  
  predmat <-  matrix(pred,nrow=len,ncol=len)
  #predmat_ct <-  matrix(pred_ct,nrow=len,ncol=len)
  
  
  
  if(svg){
    svg(sprintf("intplots_%s_%s_%s.svg",var1_2,var2_2,season),width=4.5,height = 6)
    
    layout(matrix(c(1:2),nrow=2))
    par(mai=c(0,0,0,0))
  }else{
    var1_2 <- varlist2[[season]]$readable[which(varlist2[[season]]$var==var1)] #strsplit(var1,"\\.")[[1]][2]
    var2_2 <- varlist2[[season]]$readable[which(varlist2[[season]]$var==var2)] #strsplit(var2,"\\.")[[1]][2]
  }
  
  x.axis <- realmean1+realsd1*2*range1
  min.x <- min(x.axis)
  max.x <- max(x.axis)
  xrange <- max.x-min.x
  x.axis2 <- seq(min.x,max.x,length=5)[c(2:4)]
  y.axis <- realmean2+realsd2*2*range2
  min.y <- min(y.axis)
  max.y <- max(y.axis)
  yrange <- max.y-min.y
  y.axis2 <- seq(min.y,max.y,length=5)[c(2:4)]
  z.axis <- seq(min(predmat), min(10,max(predmat)), length=5)[c(2:4)]
  min.z <- min(predmat)
  max.z <- min(10,max(predmat))
  zrange <- max.z-min.z
  predmat[which(predmat>max.z,arr.ind = T)] <- max.z
  xlab=paste("",var1_2,sep="")
  ylab=paste("",var2_2,sep="")
  zlab="Sel. intens."
  persp(x.axis,y.axis,predmat,zlim=c(min.z,max.z),
        axes=F,
        theta = 55, phi = 40, r = sqrt(10), d = 3, 
        ticktype = "detailed", mgp = c(4, 1, 0)) -> res
  # lines(trans3d(x.axis, min.y, min.z, res) , col="black",lwd=2)
  # lines(trans3d(max.x, y.axis, min.z, res) , col="black",lwd=2)
  # lines(trans3d(min.x, min.y, z.axis, res) , col="black",lwd=2)
  
  labels <- sprintf("%.1f",x.axis2)
  label.pos <- trans3d(x.axis2, (min.y - yrange*0.3), min.z, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=30, cex=1)
  
  labels <- xlab
  label.pos <- trans3d(x.axis2[1]-xrange*0.3, (min.y - yrange*0.4), min.z, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=30+100+180, cex=1)
  
  labels <- sprintf("%.1f",y.axis2)
  label.pos <- trans3d((max.x + xrange*0.15), y.axis2, min.z, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), cex=1) 
  
  labels <- ylab
  label.pos <- trans3d((max.x + xrange*0.4), y.axis2[2], min.z, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), cex=1,srt=25)
  
  labels <- sprintf("%.1f",z.axis)
  label.pos <- trans3d(min.x, (min.y - yrange*0.1), z.axis, res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(1, NA), cex=1)
  
  labels <- zlab
  label.pos <- trans3d((min.x), (min.y -5.5), z.axis[3], res)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(1, NA), cex=1,srt=100)
  
  
  var1_2 <- strsplit(var1,"\\.")[[1]][2]
  var2_2 <- strsplit(var2,"\\.")[[1]][2]
  subst <- data[sample(1:nrow(data),1000,replace = F),]
  points(trans3d(subst[[var1_2]], subst[[var2_2]], max(predmat)+0, pmat = res), col = gray(0.8), pch = 19, cex=0.5)
  
  if(!svg){
    var1_2 <- varlist2[[season]]$readable[which(varlist2[[season]]$var==var1)] #strsplit(var1,"\\.")[[1]][2]
    var2_2 <- varlist2[[season]]$readable[which(varlist2[[season]]$var==var2)] #strsplit(var2,"\\.")[[1]][2]
  }
  
  if(svg) dev.off() 
}




#############
# Load data
#############

getwd()

#load files (count data)

alldata<-list()

alldata$mating <- readRDS("rut.obs.01Jan2020.rds")              
alldata$parturition <- readRDS("part.obs.01Jan2020.rds")
alldata$summer <- readRDS("summer.obs.01Jan2020.rds")
alldata$winter <- readRDS("winter.obs.01Jan2020.rds")

allcollar <- list()

#load data
allcollar$mating <- readRDS("all.extracts.collars.rut.02Jan2020.rds")
allcollar$parturition <- readRDS("all.extracts.collars.part.02Jan2020.rds")
allcollar$summer <- readRDS("all.extracts.collars.summer.02Jan2020.rds")
allcollar$winter <- readRDS("all.extracts.collars.winter.02Jan2020.rds")

allseasons <- c("mating","parturition","summer","winter" )

##############
# Process data

names(alldata[[1]])

    #replace with relevant subset of columns 
alldata2 <- lapply(alldata,function(t) t[c(1, 2, 9, 12, 14, 16:30)])

names(allcollar[[1]])
    #replace with relevant subset of columns 
allcollar2 <- lapply(allcollar,function(t) t[c(3:5, 9:14, 17:32)])

head(allcollar2[[1]])

############
# Check for zero-sum constraint

summary(alldata$mating)

graphics.off()
sapply(allseasons,function(t) hist(alldata2[[t]]$grassland,main=t))
sapply(allseasons,function(t) hist(alldata2[[t]]$scrub,main=t))
sapply(allseasons,function(t) hist(alldata2[[t]]$pasture,main=t))

t="mating"
layout(matrix(1:4,nrow=2))
sapply(allseasons,function(t) hist(alldata2[[t]]$grassland+alldata2[[t]]$scrub,main=t,xlab="sum of grassland and scrub"))
sapply(allseasons,function(t) hist(allcollar2[[t]]$grassland+allcollar2[[t]]$scrub,main=t,xlab="sum of grassland and scrub"))
        # grassland and scrub don't sum to 1- sometimes add up to 200%

sapply(allseasons,function(t) hist(alldata2[[t]]$wetgrass+alldata2[[t]]$scrub,main=t,xlab="sum of wetgrass and scrub"))
        # but wetgrass looks good

names(alldata2[[1]])

sapply(allseasons,function(t) hist(alldata2[[t]]$wetgrass+alldata2[[t]]$drygrass+alldata2[[t]]$ag,main=t,xlab="sum of wetgrass and dry"))
sapply(allseasons,function(t) hist(allcollar2[[t]]$wetgrass+allcollar2[[t]]$drygrass+allcollar2[[t]]$ag,main=t,xlab="sum of wetgrass and dry"))
# dry and wet grass DO sum to 1



##########
# confirm that collar data comes rarefied?

allcollar$mating


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
  df$std.wetgrass <- ((df$wetgrass)-(mean(df$wetgrass)))/(2*sd(df$wetgrass))
  df$std.drygrass <- ((df$drygrass)-(mean(df$drygrass)))/(2*sd(df$drygrass))
  df$std.ag <- ((df$ag)-(mean(df$ag)))/(2*sd(df$ag))
  return(df)
}

data.std <- lapply(alldata2,function(t) standardze(t))

sapply(data.std,function(t) nrow(t)-length(which(complete.cases(t))))

final.data <- lapply(data.std,function(t) t[complete.cases(t), ])

sapply(final.data,nrow)    # lots of data

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
  df$std.wetgrass <- ((df$wetgrass)-(mean(df$wetgrass)))/(2*sd(df$wetgrass))
  df$std.drygrass <- ((df$drygrass)-(mean(df$drygrass)))/(2*sd(df$drygrass))
  df$std.ag <- ((df$ag)-(mean(df$ag)))/(2*sd(df$ag))
  return(df)
}

collar.std <- lapply(allcollar2,function(t) standardze2(t))

final.collar <- lapply(collar.std,function(t) t[complete.cases(t), ])

sapply(final.collar,nrow)


#subset available points to be re-combined with "used" subset for each ID above
avails <- lapply(final.collar, function(t)  t[t$sex %in% 'NA',])
sapply(avails,nrow)

useds <- lapply(final.collar, function(t)  t[!t$sex %in% 'NA',])
sapply(useds,nrow)


rarefy <- min(sapply(final.data,nrow))   
#t=final.data[[1]]
final.data2 <- lapply(final.data, function(t) t[ceiling(seq(0.1,nrow(t),length=rarefy)),]   )



tofac <- function(df){
  df$year <- as.factor(df$year)
  df$month <- as.factor(df$month)
  df$ID <- as.factor(df$ID)
  return(df)
}

final.data <- lapply(final.data2,function(t) tofac(t))

sapply(final.data,nrow)

rm(allcollar,allcollar2,alldata,alldata2,final.data2,rarefy,t,data.std,collar.std)


#############
# get data summaries for ms
#############

final.data2 <- do.call('rbind',final.data)

mean(apply(table(final.data2$year,final.data2$ID),2,mean))   # average of 43 surveys per year

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
# # Conclusions from structural model selection procedure (applies to mating, Summer, Winter, Part:   
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

# ###Run through each season one at a time
# 
# season <- "winter"
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
# bbmle::AICtab(cow.mean.nonz, cow.mean, cow.med, cow.max,weights=TRUE,mnames=c("Mean Nonz", "Mean", "Median", "Max") )

# ######Results from current AIC test: mean for all except summer, which is mean nonz
# Parturition
# dAIC df weight
# Max        0.0 33 0.44  
# Mean Nonz  0.7 33 0.31  
# Mean       1.6 33 0.20  
# Median     4.4 33 0.05  

# Summer
# dAIC df weight
# Mean Nonz  0.0 33 0.865 
# Mean       4.1 33 0.113 
# Median     7.3 33 0.022 
# Max         NA 33 NA   

# mating
# dAIC df weight
# Mean       0.0 33 0.9433
# Median     6.9 33 0.0305
# Mean Nonz  7.3 33 0.0248
# Max       13.0 33 0.0014

# Winter
# dAIC df weight
# Median     0.0 33 0.731 
# Mean       2.3 33 0.236 
# Mean Nonz  6.2 33 0.033 
# Max         NA 33 NA    




# ######Results from previous AIC test:
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
# #mating:
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

# #mating:
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
  mating=stockingVars[2],
  parturition=stockingVars[2],
  summer=stockingVars[1],
  winter=stockingVars[2]
)
best.stocking


###########
# Correlation analysis
###########

names(final.data[[1]])
cor(final.data[[1]][,c(6,8,11,12:20)])

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

season <- allseasons[4]
for(season in allseasons){
  df <- final.data[[season]]
  
  formulas[[season]] <- as.formula(sprintf("TotalElk ~ std.pasture + std.%s  + std.ponds + std.elevation + std.ndvi + std.scrub + std.wetgrass + std.ag + std.slope + std.aspect + std.ndvi:std.scrub  + std.ndvi:std.wetgrass + std.%s:std.ndvi + std.slope:std.aspect + (1|ID) + (1|year)",best.stocking[season],best.stocking[season] ))    #  
  
  #temp <- glmmTMB(formulas[[season]],
  #        ziformula = ~.,data=df, dispformula = ~1, family=truncated_nbinom1(link="log"))
  
  bestmods[[season]] <- buildmer::buildglmmTMB(formulas[[season]],
             ziformula = ~.,data=df, dispformula = ~1, direction="backward",crit="AIC",  # use LRT instead
             family=truncated_nbinom1(link="log"))    # reduce.random=T,
  bestmod2 <- bestmods[[season]]@model

  cat(sprintf("done with transect model for season %s\n",season))
  save(bestmod2,file = sprintf("bestmod_%s.RData",season))
  
  
    # #summary(bestmod2)       # move visualizations to a different part of code
  # tmp <- rownames(coef(summary(bestmod2))$cond)[-1]
  # topvars <- tmp[!grepl(":",tmp)]
  # ints <- strsplit(tmp[grep(":",tmp)],":")
  # if(length(topvars)>0) tmp <- sapply(1:length(topvars), function(t) VisualizeRelation(df,bestmod2,topvars[t],season=season) )     # run and save partial dependence plots for all top variables
  # 
  # if(length(ints)>0 ) tmp <- sapply(1:length(ints), function(t) VisualizeInteraction(df,bestmod2,ints[[t]][1],ints[[t]][2],season=season) )
  
}

save(bestmods,file="bestmods_allseasons.RData")

load("bestmods_allseasons.RData")

lapply(bestmods,summary)



##############
# Collar data model (single model for all data- one per season) 
##############

season=allseasons[1]
sapply(useds,nrow)
sapply(avails,nrow)
nused <- 2500
perind <- 4
navail <- nused*perind

bestmods_collar <- list()
bestmods_collar2 <- list()    # to enable prediction
for(season in allseasons){   
  # Set up the model, but do not yet fit it
  dfu <- useds[[season]]
  dfa <- avails[[season]]
  ndx_u <- ceiling(seq(0.1,nrow(dfu),length=nused))
  dfu <- dfu[ndx_u,]
  ids_a <- rep(as.character(dfu$ID),each=perind)
  ndx_a <- ceiling(seq(0.1,nrow(dfa),length=navail))
  dfa <- dfa[ndx_a,]
  dfa$ID <- ids_a    # randomly assign random points to individuals
  
  thisdat <- rbind(dfu,dfa)
  thisdat$weights = as.numeric(as.character(thisdat$Used))
  thisdat$weights[thisdat$weights==0] <- 1000
  
  rm(dfu, dfa, ndx_u,ndx_a,ids_a)
  
  full_formula <- as.formula(gsub("[\r\n]","",sprintf("Used ~ std.pasture + std.%s + std.ponds 
                             + std.elevation + std.ndvi + std.scrub + std.wetgrass + std.ag + std.slope + std.aspect 
                             + std.ndvi:std.scrub + std.ndvi:std.wetgrass + std.%s:std.ndvi + std.slope:std.aspect
                             + (1|ID) + (0 + std.pasture | ID) + (0 + std.ponds | ID) + (0 + std.elevation | ID) + (0 + std.ag | ID) 
                             + (0 + std.ndvi | ID) + (0 + std.%s | ID) + (0 + std.scrub | ID) 
                             + (0 + std.wetgrass | ID)+(0 + std.slope | ID)+(0 + std.aspect | ID) +
                             + (0 + std.ndvi:std.scrub | ID) + (0 + std.ndvi:std.wetgrass | ID) + (0 + std.%s:std.ndvi | ID) + (0 + std.slope:std.aspect | ID)",best.stocking[season],best.stocking[season],best.stocking[season],best.stocking[season])))
  

  temp_bestmod <- buildmer::buildglmmTMB(full_formula, weights=thisdat$weights,
                         data=thisdat, direction="backward",crit="AIC",
                         family=binomial(link="logit"),REML=FALSE)   #reduce.random=T,

  red_formula <- temp_bestmod@model$call$formula
  
  numre <- length(grep("std.",colnames(ranef( temp_bestmod@model)$cond$ID)))
  
  
  ## now fit the best model to get proper parameter estimates!  
  TMBStruc = glmmTMB(red_formula,
                     weights = weights,
                     family=binomial,
                     data=thisdat,
                     doFit=FALSE)
  
  
  # Fix the standard deviation of the first random term, which is the (1|id) component
  # in the above model equation
  TMBStruc$parameters$theta[1] = log(1e3)
  
  # Tell glmmTMB not to change the first entry of the vector of variances,
  # and give all other variances another indicator to make sure they can be freely estimated
  TMBStruc$mapArg = list(theta=factor(c(NA,1:numre)))
  
  bestmods_collar[[season]] = glmmTMB:::fitTMB(TMBStruc)   
  # summary(bestmods_collar[[season]])
  
  
       # enable predictions!
  bestmods_collar2[[season]] = glmmTMB::glmmTMB(red_formula,weights=weights,
                                                family=binomial,data=thisdat)

  
}

save(bestmods_collar,bestmods_collar2, file="collarmods2_allseasons.RData")
load("collarmods2_allseasons.RData")

summary(bestmods_collar$mating)
summary(bestmods_collar$summer)
summary(bestmods_collar$parturition)
summary(bestmods_collar$winter)

# reorder bestmods_collar in same season order as everything else 

##########################
############
# Visualize coefficients

temp_count <- lapply(bestmods, function(t) coef(summary(t))$cond )    # count model coefs
temp_zero <- lapply(bestmods, function(t) coef(summary(t))$zi )     # zero model coefs
temp_rscoefs <- lapply(bestmods_collar, function(t) coef(summary(t))$cond )     # resource selection coefs

allcoefnames <- unique(unlist(sapply(temp_count,function(t) rownames(t) )))
allmains <- sort(allcoefnames[!grepl(":",allcoefnames)])
allints <- sort(allcoefnames[grepl(":",allcoefnames)])
allcoefnames <- c(allmains,allints)

coef_ndx <- lapply(temp_count,function(t) match(allcoefnames,rownames(t)) )
coef_ndx_rsf <- lapply(temp_rscoefs,function(t) match(allcoefnames,rownames(t)) )

temp <- rep(NA,times=length(allcoefnames))

allcoefs <- lapply(1:length(temp_count), function(t) temp_count[[t]][,'Estimate'][coef_ndx[[t]]])
names(allcoefs) <- allseasons

allses <- lapply(1:length(temp_count), function(t) temp_count[[t]][,'Std. Error'][coef_ndx[[t]]])
names(allses) <- allseasons

allcoefsz <- lapply(1:length(temp_zero), function(t) temp_zero[[t]][,'Estimate'][coef_ndx[[t]]])
names(allcoefsz) <- allseasons

allsesz <- lapply(1:length(temp_zero), function(t) temp_zero[[t]][,'Std. Error'][coef_ndx[[t]]])
names(allsesz) <- allseasons

allcoefrsf <- lapply(1:length(temp_rscoefs), function(t) temp_rscoefs[[t]][,'Estimate'][coef_ndx_rsf[[t]]])
names(allcoefrsf) <- allseasons

allsesrsf <- lapply(1:length(temp_rscoefs), function(t) temp_rscoefs[[t]][,'Std. Error'][coef_ndx_rsf[[t]]])
names(allsesrsf ) <- allseasons

best.stocking

varlist <- data.frame(
  var = allcoefnames[-c(1,11,17)],    # as.character(thisres[[1]]$names)
  readable = c("High Intensity Ag","Cos Aspect","Elevation","NDVI","Pasture","Dist to Pond",
               "% Scrub","Slope","Cattle Dens","% Wet Grassland",
               "NDVI*Scrub","NDVI*WetGrass","Slope*CosAspect","NDVI*CattleDens"),
  stringsAsFactors = F
) 

varlist <- varlist[c(5,9,6,1,10,4,7,3,8,2,14,12,11,13),]

varlist2 <- list()

temp <- lapply(allseasons,function(t) varlist2[[t]] <<- varlist)

best.stocking2 <- best.stocking
#best.stocking2[3] <- "stockingMean"
temp <- lapply(allseasons,function(t) varlist2[[t]]$var[2] <<- sprintf("std.%s",best.stocking[t]) )
temp <- lapply(allseasons,function(t) varlist2[[t]]$var[11] <<- sprintf("std.%s:std.ndvi",best.stocking[t]) )

varlist2
#varlist2 <- lapply(varlist2,function(t) t[!t$var=="std.slope",])

allseasons2 <- allseasons[c(4,2,3,1)]    # for visualization
#allseasons2[1] <- "breeding"

graphics.off()

##### make sure the plot looks okay (shrink really wide conf ints)
v=2
flags <- list()
for(v in 1:nrow(varlist)){
  thisvar <- varlist[v,]$var
  flags[[thisvar]] <- list()
  season = allseasons2[1]
  for(season in allseasons2){
    flags[[thisvar]][[season]] <- numeric(3)
    thisvar2 <- varlist2[[season]][v,]$var
    
    thiscoefrsf <- allcoefrsf[[season]][thisvar2]
    thissersf <- allsesrsf[[season]][thisvar2]
    if(!any(c(is.na(thiscoefrsf),is.na(thissersf)))){
      if((abs(thiscoefrsf)>1.9)|(thissersf>1.3)){
        allcoefrsf[[season]][thisvar2] <- allcoefrsf[[season]][thisvar2]/5
        allsesrsf[[season]][thisvar2] <- allsesrsf[[season]][thisvar2]/5
        flags[[thisvar]][[season]][1] <- 1
      } 
    }
    
    thiscoefz <- allcoefsz[[season]][thisvar2]
    thissez <- allsesz[[season]][thisvar2]
    if(!any(c(is.na(thiscoefz),is.na(thissez)))){
      if((abs(thiscoefz)>1.9)|(thissez>1.3)){
        allcoefsz[[season]][thisvar2] <- allcoefsz[[season]][thisvar2]/5
        allsesz[[season]][thisvar2] <- allsesz[[season]][thisvar2]/5
        flags[[thisvar]][[season]][3] <- 1
      }
    }
    
    
    thiscoef <- allcoefs[[season]][thisvar2]
    thisse <- allses[[season]][thisvar2]
    if(!any(c(is.na(thiscoef),is.na(thisse)))){
      if((abs(thiscoef)>1.9)|(thisse>1.3)){
        allcoefs[[season]][thisvar2] <- allcoefs[[season]][thisvar2]/5
        allses[[season]][thisvar2] <- allses[[season]][thisvar2]/5
        flags[[thisvar]][[season]][2] <- 1
      }
    }
   
  }
}

flags



svg("coefplot_final3.svg",7.5,6.8)

layout(matrix(1:16,nrow=4,byrow = T),widths = c(0.7,2,2,2))          
var_ndx <- 1
onleft <- c(1,5,9,13)
panel=2
for(panel in c(1:14)){
  par(mai=c(0.1,0.1,0.3,0.01))
  xcents <- c(11,17,23)
  ycents <- seq(15,3,length=4)
  if(panel %in% onleft){
    plot(20,20,xlim=c(1,4),ylim=c(-2.5,16),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    text(rep(1,4),ycents, tools::toTitleCase(allseasons2),adj=0,cex=1.3)
  }else{
    thisvar = varlist[var_ndx,]$var
    plot(20,20,xlim=c(8,25),ylim=c(-2.5,16),xlab="",ylab="",xaxt="n",yaxt="n",main=varlist2$mating$readable[var_ndx],pch="")    # blank plot
    yshift <- 0
    
    polygon(c(xcents[1]-2,xcents[1]+2,xcents[1]+2,xcents[1]-2),c(rep(range(ycents)[1],2)-1.8-yshift,
                                                                 rep(range(ycents)[2],2)+1.5-yshift),border=NA,col=gray(0.8))
    polygon(c(xcents[2]-2,xcents[2]+2,xcents[2]+2,xcents[2]-2),c(rep(range(ycents)[1],2)-1.8+yshift,
                                                                 rep(range(ycents)[2],2)+1.5+yshift),border=NA,col=gray(0.8))
    polygon(c(xcents[3]-2,xcents[3]+2,xcents[3]+2,xcents[3]-2),c(rep(range(ycents)[1],2)-1.8-yshift,
                                                                 rep(range(ycents)[2],2)+1.5-yshift),border=NA,col=gray(0.8))
    
    segments(xcents,rep(range(ycents)[1],3)-1+c(-yshift,yshift,-yshift),xcents,rep(range(ycents)[2],3)+1+c(-yshift,yshift,-yshift),lwd=1)
    segments(xcents[2],range(ycents)[1]-1+yshift,xcents[2],range(ycents)[2]+1+yshift,lwd=2,col=gray(0.4))
    
    
    
    # meanvalz <- sapply(allseasons,function(t) thisres[[t]]$meancoefs[as.character(thisres[[t]]$names)==varlist2[[t]]$var[var_ndx]] )
    # points(meanvalz+xcents[1],ycents-yshift,pch=20,cex=1.5)
    # lower <- sapply(allseasons,function(t) thisres[[t]]$lowr[as.character(thisres[[t]]$names)==varlist2[[t]]$var[var_ndx]] )
    # upper <- sapply(allseasons,function(t) thisres[[t]]$uppr[as.character(thisres[[t]]$names)==varlist2[[t]]$var[var_ndx]] )
    # arrows(lower+xcents[1],ycents-yshift,upper+xcents[1],ycents-yshift,angle = 90,code=3,length=0.05)
    #
    points(sapply(allseasons2,function(t) allcoefrsf[[t]][varlist2[[t]]$var[var_ndx]])+xcents[1],ycents+yshift,pch=20,cex=1.5)  # col=gray(0.4)
    arrows(sapply(allseasons2,function(t) allcoefrsf[[t]][varlist2[[t]]$var[var_ndx]])+xcents[1]-1.66*sapply(allseasons2,function(t) allsesrsf[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents+yshift,sapply(allseasons2,function(t) allcoefrsf[[t]][varlist2[[t]]$var[var_ndx]])+xcents[1]+1.66*sapply(allseasons2,function(t) allsesrsf[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents+yshift,angle = 90,code=3,length=0.05,col=gray(0.4))
    shrunk <- as.logical(sapply(allseasons2, function(t) flags[[thisvar]][[t]][1]))
    text(xcents[1]-2.5,ycents[shrunk],"*",cex=1.8)
    
    points(sapply(allseasons2,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])+xcents[2],ycents+yshift,pch=20,cex=1.5,col=gray(0.4))
    arrows(sapply(allseasons2,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])+xcents[2]-1.66*sapply(allseasons2,function(t) allses[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents+yshift,sapply(allseasons2,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])+xcents[2]+1.66*sapply(allseasons2,function(t) allses[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents+yshift,angle = 90,code=3,length=0.05,col=gray(0.4))
    shrunk <- as.logical(sapply(allseasons2, function(t) flags[[thisvar]][[t]][2]))
    text(xcents[2]-2.5,ycents[shrunk],"*",cex=1.8)
    
    points((sapply(allseasons2,function(t) allcoefsz[[t]][varlist2[[t]]$var[var_ndx]]))*-1+xcents[3],ycents-yshift,pch=20,cex=1.5)
    arrows((sapply(allseasons2,function(t) allcoefsz[[t]][varlist2[[t]]$var[var_ndx]]))*-1+xcents[3]-1.66*sapply(allseasons,function(t) allsesz[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents-yshift,(sapply(allseasons2,function(t) allcoefsz[[t]][varlist2[[t]]$var[var_ndx]]))*-1+xcents[3]+1.66*sapply(allseasons2,function(t) allsesz[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents-yshift,angle = 90,code=3,length=0.05)
    shrunk <- as.logical(sapply(allseasons2, function(t) flags[[thisvar]][[t]][3]))
    text(xcents[3]-2.5,ycents[shrunk],"*",cex=1.8)
    
    abline(h=range(ycents)[1]-2,lwd=1.5)
    text(xcents[2],-2,"Effect size (standardized)",cex=1.2)
    text(xcents,c(-0,-0,-0),c("Collar","Count","P(1+)"),cex=1,col=c("black",gray(0.4),"black"))
    abline(h=c(ycents[1]+2,ycents[-length(ycents)]-2),lty=2,lwd=0.5)
    var_ndx <- var_ndx+1
  }

}

dev.off()


### make 'coefficients table' for supplement

# allcoefrsf, allcoefs, allcoefsz
# allsesrsf, allses, allsesz
# flags
# varlist

temp <- varlist[,-1,drop=F]
names(temp) <- "Variable"

coefdf <- NULL

allseasons2
s=1
for(s in 1:length(allseasons2)){
  temp2 <- temp
  temp2$Season <- allseasons2[s]
  temp2$Estimate_GPS <- allcoefrsf[[allseasons2[s]]][match(varlist$var,names(allcoefrsf[[allseasons2[s]]]))]
  temp2$SE_GPS <- allsesrsf[[allseasons2[s]]][match(varlist$var,names(allcoefrsf[[allseasons2[s]]]))]
  temp2$Estimate_Dens <- allcoefs[[allseasons2[s]]][match(varlist$var,names(allcoefrsf[[allseasons2[s]]]))]
  temp2$SE_Dens <- allses[[allseasons2[s]]][match(varlist$var,names(allcoefrsf[[allseasons2[s]]]))]
  temp2$Estimate_Zero <- allcoefsz[[allseasons2[s]]][match(varlist$var,names(allcoefrsf[[allseasons2[s]]]))]
  temp2$SE_Zero <- allsesz[[allseasons2[s]]][match(varlist$var,names(allcoefrsf[[allseasons2[s]]]))]
  coefdf <- rbind(coefdf,temp2)
}
coefdf

flags
v=1
for(v in 1:nrow(varlist)){
  thisvar <- varlist$readable[v]
  thisvar2 <- varlist$var[v]
  ndx1 <- which(coefdf$Variable==thisvar)
  names(ndx1) <- allseasons2
  thisflags <- flags[[thisvar2]]
  s=1
  for(s in 1:length(allseasons2)){
    thisseas <- allseasons2[s]
    thisflags2 <- thisflags[[thisseas]]
    if(thisflags2[1]==1){
      coefdf[ndx1[thisseas],]$Estimate_GPS <- coefdf[ndx1[thisseas],]$Estimate_GPS*5
      coefdf[ndx1[thisseas],]$SE_GPS <- coefdf[ndx1[thisseas],]$SE_GPS*5 
    }
    if(thisflags2[2]==1){
      coefdf[ndx1[thisseas],]$Estimate_Dens <- coefdf[ndx1[thisseas],]$Estimate_Dens*5
      coefdf[ndx1[thisseas],]$SE_Dens <- coefdf[ndx1[thisseas],]$SE_Dens*5 
    }
    if(thisflags2[3]==1){
      coefdf[ndx1[thisseas],]$Estimate_Zero <- coefdf[ndx1[thisseas],]$Estimate_Zero*5
      coefdf[ndx1[thisseas],]$SE_Zero <- coefdf[ndx1[thisseas],]$SE_Zero*5 
    }
  }
}
coefdf

getwd()
write.csv(coefdf,"coef_table.csv")


svg("intcoefs_final3.svg",6,3.5)


layout(matrix(1:6,nrow=2,byrow = T),widths = c(0.7,2,2))          
var_ndx <- 11
onleft <- c(1,4)

for(panel in c(1:6)){
  par(mai=c(0.1,0.1,0.3,0.01))
  xcents <- c(14,22,30)
  ycents <- seq(15,3,length=4)
  if(panel %in% onleft){
    plot(20,20,xlim=c(1,4),ylim=c(-2.5,16),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    text(rep(1,4),ycents, tools::toTitleCase(allseasons2),adj=0,cex=1.3)
  }else{
    plot(20,20,xlim=c(9,32),ylim=c(-2.5,16),xlab="",ylab="",xaxt="n",yaxt="n",main=varlist2$mating$readable[var_ndx],pch="")    # blank plot
    yshift <- 0
    
    polygon(c(xcents[1]-2,xcents[1]+2,xcents[1]+2,xcents[1]-2),c(rep(range(ycents)[1],2)-1.8-yshift,
                                                                 rep(range(ycents)[2],2)+1.5-yshift),border=NA,col=gray(0.8))
    polygon(c(xcents[2]-2,xcents[2]+2,xcents[2]+2,xcents[2]-2),c(rep(range(ycents)[1],2)-1.8+yshift,
                                                                 rep(range(ycents)[2],2)+1.5+yshift),border=NA,col=gray(0.8))
    polygon(c(xcents[3]-2,xcents[3]+2,xcents[3]+2,xcents[3]-2),c(rep(range(ycents)[1],2)-1.8-yshift,
                                                                 rep(range(ycents)[2],2)+1.5-yshift),border=NA,col=gray(0.8))
    
    segments(xcents,rep(range(ycents)[1],3)-1+c(-yshift,yshift,-yshift),xcents,rep(range(ycents)[2],3)+1+c(-yshift,yshift,-yshift),lwd=1)
    segments(xcents[2],range(ycents)[1]-1+yshift,xcents[2],range(ycents)[2]+1+yshift,lwd=2,col=gray(0.4))
    
    
    
    # meanvalz <- sapply(allseasons,function(t) thisres[[t]]$meancoefs[as.character(thisres[[t]]$names)==varlist2[[t]]$var[var_ndx]] )
    # points(meanvalz+xcents[1],ycents-yshift,pch=20,cex=1.5)
    # lower <- sapply(allseasons,function(t) thisres[[t]]$lowr[as.character(thisres[[t]]$names)==varlist2[[t]]$var[var_ndx]] )
    # upper <- sapply(allseasons,function(t) thisres[[t]]$uppr[as.character(thisres[[t]]$names)==varlist2[[t]]$var[var_ndx]] )
    # arrows(lower+xcents[1],ycents-yshift,upper+xcents[1],ycents-yshift,angle = 90,code=3,length=0.05)
    # 
    points(sapply(allseasons2,function(t) allcoefrsf[[t]][varlist2[[t]]$var[var_ndx]])+xcents[1],ycents+yshift,pch=20,cex=1.5)  # col=gray(0.4)
    arrows(sapply(allseasons2,function(t) allcoefrsf[[t]][varlist2[[t]]$var[var_ndx]])+xcents[1]-1.66*sapply(allseasons2,function(t) allsesrsf[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents+yshift,sapply(allseasons2,function(t) allcoefrsf[[t]][varlist2[[t]]$var[var_ndx]])+xcents[1]+1.66*sapply(allseasons2,function(t) allsesrsf[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents+yshift,angle = 90,code=3,length=0.05)   # col=gray(0.4)
    shrunk <- as.logical(sapply(allseasons2, function(t) flags[[thisvar]][[t]][1]))
    text(xcents[1]-2.5,ycents[shrunk],"*",cex=1.8)
    
      # points(sapply(allseasons,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])/c(5,1,1,1)+xcents[2],ycents+yshift,pch=20,cex=1.5,col=gray(0.4))
      # arrows(sapply(allseasons,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])/c(5,1,1,1)+xcents[2]-1.66*sapply(allseasons,function(t) allses[[t]][varlist2[[t]]$var[var_ndx]])/c(5,1,1,1),
      #        ycents+yshift,sapply(allseasons,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])/c(5,1,1,1)+xcents[2]+1.66*sapply(allseasons,function(t) allses[[t]][varlist2[[t]]$var[var_ndx]])/c(5,1,1,1),
      #        ycents+yshift,angle = 90,code=3,length=0.05,col=gray(0.4))
      # text(xcents[2]-2,ycents[1]-1,"*",cex=1.8)
      # 
    points(sapply(allseasons2,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])+xcents[2],ycents+yshift,pch=20,cex=1.5,col=gray(0.4))   #
    arrows(sapply(allseasons2,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])+xcents[2]-1.66*sapply(allseasons2,function(t) allses[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents+yshift,sapply(allseasons2,function(t) allcoefs[[t]][varlist2[[t]]$var[var_ndx]])+xcents[2]+1.66*sapply(allseasons2,function(t) allses[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents+yshift,angle = 90,code=3,length=0.05,col=gray(0.4))
    shrunk <- as.logical(sapply(allseasons2, function(t) flags[[thisvar]][[t]][2]))
    text(xcents[2]-2.5,ycents[shrunk],"*",cex=1.8)
  
    
    points((sapply(allseasons2,function(t) allcoefsz[[t]][varlist2[[t]]$var[var_ndx]]))*-1+xcents[3],ycents-yshift,pch=20,cex=1.5)
    arrows((sapply(allseasons2,function(t) allcoefsz[[t]][varlist2[[t]]$var[var_ndx]]))*-1+xcents[3]-1.66*sapply(allseasons2,function(t) allsesz[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents-yshift,(sapply(allseasons2,function(t) allcoefsz[[t]][varlist2[[t]]$var[var_ndx]]))*-1+xcents[3]+1.66*sapply(allseasons2,function(t) allsesz[[t]][varlist2[[t]]$var[var_ndx]]),
           ycents-yshift,angle = 90,code=3,length=0.05)
    shrunk <- as.logical(sapply(allseasons2, function(t) flags[[thisvar]][[t]][3]))
    text(xcents[3]-2.5,ycents[shrunk],"*",cex=1.8)
    
    abline(h=range(ycents)[1]-2,lwd=1.5)
    text(xcents[2],-2,"Effect size (standardized)",cex=1.2)
    text(xcents,c(-0,-0,-0),c("Collar","Count","P(1+)"),cex=1,col=c("black",gray(0.4),"black"))
    abline(h=c(ycents[1]+2,ycents[-length(ycents)]-2),lty=2,lwd=0.5)
    var_ndx <- var_ndx+1
  }
  
}


dev.off()

graphics.off()


##########
# Visualize key interactions
##########

  ## specificially, highlight the ndvi by stocking interactions

ordr <- data.frame(
  season=rep(allseasons2,each=2),
  dat=rep(c(2,1),times=4),
  stringsAsFactors = F
)

dfs <- list()
dfs[[1]] <- lapply(allseasons2, function(t) final.data[[t]])
dfs[[2]] <- lapply(allseasons2, function(t) final.collar[[t]])  #final.collar[[season]]
names(dfs[[1]]) <- allseasons2
names(dfs[[2]]) <- allseasons2
dfs[[1]]$mating

graphics.off()
base <- par()
base$mai

svg("int_dens_ndvi4.svg",7,7)

ordr2 <- ordr
#ordr2 <- ordr2[!ordr2[,"season"]=="summer",]

layout(matrix(1:20,nrow=5,byrow = T),widths = c(2,4,4,4),heights=c(1,4,4,4,4))          
onleft <- c(5,8,11,14)
ontop <- c(1,2,3,4)
counter=1
counter2=1
counter3=1
panel=0
panel=panel+1
for(panel in c(1:16)){
  if(panel %in% onleft){
    par(mai=c(0,0,0,0))
    plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    text(0,0.5, tools::toTitleCase(ordr2$season[counter3]),adj=0,cex=1.5)
    counter=counter+1
  }else if (panel %in% ontop){
    if(panel==1){
      par(mai=c(0,0,0,0))
      plot(0,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    }else{
      par(mai=c(0,0,0,0))
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5, c("Collar","Count","P(1+)")[counter2],adj=0.5,cex=1.5)
      counter2=counter2+1
    }
  } else{
    #par(mai=base$mai)
    par(mai=c(0.1,0.1,0.0,0.0))
    season <- ordr2$season[counter3]
    thisdf <- dfs[[ordr2$dat[counter3]]][[season]]
    tmp <- names(fixef(bestmods_collar2[[season]])$cond)[-1]
    allmains <- tmp[!grepl(":",tmp)]
    #if(season!="summer"){
      if(ordr2$dat[counter3]==2){
        if(!is.na(allcoefrsf[[season]][sprintf("std.%s:std.ndvi",best.stocking[season])])){
          VisualizeInteraction_col(data=thisdf,model=bestmods_collar2[[season]],  #meancoefs=thisres[[season]],
                                   var2=sprintf("std.%s",best.stocking[season]),
                                   var1="std.ndvi",
                                   allvars=allmains,
                                   season=season,svg=F)
        }else{
          plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
        }
        counter3 <- counter3+1
         
      }else{
        #if(season!="parturition"){
        if(!all(is.na(c(allcoefs[[season]][sprintf("std.%s:std.ndvi",best.stocking[season])],
                    allcoefsz[[season]][sprintf("std.%s:std.ndvi",best.stocking[season])])))){
          VisualizeInteraction(data=thisdf,model=bestmods[[season]]@model,
                             var2=sprintf("std.%s",best.stocking[season]),
                             var1="std.ndvi",
                             allvars=allmains,  #names(thisdf)[-c(1,2,3,4,5,6,7,8,9)]
                             season=season,svg=F)
          }else{
            plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
            plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
          }
        counter3 <- counter3+1
        
        #}else{
        #  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
        #  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
        #}
      } 
   
    # }
    

    #counter3 <- counter3+1
  }
}

dev.off()
  
#par()







svg("int_asp_slope4.svg",7,7)

ordr2 <- ordr
#ordr2 <- ordr2[!ordr2[,"season"]=="summer",]

layout(matrix(1:20,nrow=5,byrow = T),widths = c(2,4,4,4),heights=c(1,4,4,4,4))          
onleft <- c(5,8,11,14)
ontop <- c(1,2,3,4)
counter=1
counter2=1
counter3=1
panel=0
panel=panel+1
for(panel in c(1:16)){
  if(panel %in% onleft){
    par(mai=c(0,0,0,0))
    plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    text(0,0.5, tools::toTitleCase(ordr2$season[counter3]),adj=0,cex=1.5)
    counter=counter+1
  }else if (panel %in% ontop){
    if(panel==1){
      par(mai=c(0,0,0,0))
      plot(0,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    }else{
      par(mai=c(0,0,0,0))
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5, c("Collar","Count","P(1+)")[counter2],adj=0.5,cex=1.5)
      counter2=counter2+1
    }
  } else{
    #par(mai=base$mai)
    par(mai=c(0.1,0.1,0.0,0.0))
    season <- ordr2$season[counter3]
    thisdf <- dfs[[ordr2$dat[counter3]]][[season]]
    tmp <- names(fixef(bestmods_collar2[[season]])$cond)[-1]
    allmains <- tmp[!grepl(":",tmp)]
    #if(season!="summer"){
    if(ordr2$dat[counter3]==2){
      if(!is.na(allcoefrsf[[season]]["std.slope:std.aspect"])){
        VisualizeInteraction_col(data=thisdf,model=bestmods_collar2[[season]],  #meancoefs=thisres[[season]],
                                 var1="std.slope",
                                 var2="std.aspect",
                                 allvars=allmains,
                                 season=season,svg=F)
      }else{
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      }
      counter3 <- counter3+1
      
    }else{
      #if(season!="parturition"){
      if(!all(is.na(c(allcoefs[[season]]["std.slope:std.aspect"],
                      allcoefsz[[season]]["std.slope:std.aspect"])))){
        VisualizeInteraction(data=thisdf,model=bestmods[[season]]@model,
                             var1="std.slope",
                             var2="std.aspect",
                             allvars=allmains,  #names(thisdf)[-c(1,2,3,4,5,6,7,8,9)]
                             season=season,svg=F)
      }else{
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      }
      counter3 <- counter3+1
      
      #}else{
      #  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      #  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      #}
    } 
    
    # }
    
    
    #counter3 <- counter3+1
  }
}

dev.off()



svg("int_scrub_ndvi4.svg",7,7)

ordr2 <- ordr
#ordr2 <- ordr2[!ordr2[,"season"]=="summer",]

layout(matrix(1:20,nrow=5,byrow = T),widths = c(2,4,4,4),heights=c(1,4,4,4,4))          
onleft <- c(5,8,11,14)
ontop <- c(1,2,3,4)
counter=1
counter2=1
counter3=1
panel=0
panel=panel+1
for(panel in c(1:16)){
  if(panel %in% onleft){
    par(mai=c(0,0,0,0))
    plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    text(0,0.5, tools::toTitleCase(ordr2$season[counter3]),adj=0,cex=1.5)
    counter=counter+1
  }else if (panel %in% ontop){
    if(panel==1){
      par(mai=c(0,0,0,0))
      plot(0,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    }else{
      par(mai=c(0,0,0,0))
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5, c("Collar","Count","P(1+)")[counter2],adj=0.5,cex=1.5)
      counter2=counter2+1
    }
  } else{
    #par(mai=base$mai)
    par(mai=c(0.1,0.1,0.0,0.0))
    season <- ordr2$season[counter3]
    thisdf <- dfs[[ordr2$dat[counter3]]][[season]]
    tmp <- names(fixef(bestmods_collar2[[season]])$cond)[-1]
    allmains <- tmp[!grepl(":",tmp)]
    #if(season!="summer"){
    if(ordr2$dat[counter3]==2){
      if(!is.na(allcoefrsf[[season]]["std.ndvi:std.scrub"])){
        VisualizeInteraction_col(data=thisdf,model=bestmods_collar2[[season]],  #meancoefs=thisres[[season]],
                                 var1="std.scrub",
                                 var2="std.ndvi",
                                 allvars=allmains,
                                 season=season,svg=F)
      }else{
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      }
      counter3 <- counter3+1
      
    }else{
      #if(season!="parturition"){
      if(!all(is.na(c(allcoefs[[season]]["std.ndvi:std.scrub"],
                      allcoefsz[[season]]["std.ndvi:std.scrub"])))){
        VisualizeInteraction(data=thisdf,model=bestmods[[season]]@model,
                             var1="std.scrub",
                             var2="std.ndvi",
                             allvars=allmains,  #names(thisdf)[-c(1,2,3,4,5,6,7,8,9)]
                             season=season,svg=F)
      }else{
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      }
      counter3 <- counter3+1
      
      #}else{
      #  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      #  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      #}
    } 
    
    # }
    
    
    #counter3 <- counter3+1
  }
}

dev.off()



svg("int_grass_ndvi4.svg",7,7)

ordr2 <- ordr
#ordr2 <- ordr2[!ordr2[,"season"]=="summer",]

layout(matrix(1:20,nrow=5,byrow = T),widths = c(2,4,4,4),heights=c(1,4,4,4,4))          
onleft <- c(5,8,11,14)
ontop <- c(1,2,3,4)
counter=1
counter2=1
counter3=1
panel=0
panel=panel+1
for(panel in c(1:16)){
  if(panel %in% onleft){
    par(mai=c(0,0,0,0))
    plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    text(0,0.5, tools::toTitleCase(ordr2$season[counter3]),adj=0,cex=1.5)
    counter=counter+1
  }else if (panel %in% ontop){
    if(panel==1){
      par(mai=c(0,0,0,0))
      plot(0,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    }else{
      par(mai=c(0,0,0,0))
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5, c("Collar","Count","P(1+)")[counter2],adj=0.5,cex=1.5)
      counter2=counter2+1
    }
  } else{
    #par(mai=base$mai)
    par(mai=c(0.1,0.1,0.0,0.0))
    season <- ordr2$season[counter3]
    thisdf <- dfs[[ordr2$dat[counter3]]][[season]]
    tmp <- names(fixef(bestmods_collar2[[season]])$cond)[-1]
    allmains <- tmp[!grepl(":",tmp)]
    #if(season!="summer"){
    if(ordr2$dat[counter3]==2){
      if(!is.na(allcoefrsf[[season]]["std.ndvi:std.wetgrass"])){
        VisualizeInteraction_col(data=thisdf,model=bestmods_collar2[[season]],  #meancoefs=thisres[[season]],
                                 var1="std.wetgrass",
                                 var2="std.ndvi",
                                 allvars=allmains,
                                 season=season,svg=F)
      }else{
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      }
      counter3 <- counter3+1
      
    }else{
      #if(season!="parturition"){
      if(!all(is.na(c(allcoefs[[season]]["std.ndvi:std.wetgrass"],
                      allcoefsz[[season]]["std.ndvi:std.wetgrass"])))){
        VisualizeInteraction(data=thisdf,model=bestmods[[season]]@model,
                             var1="std.wetgrass",
                             var2="std.ndvi",
                             allvars=allmains,  #names(thisdf)[-c(1,2,3,4,5,6,7,8,9)]
                             season=season,svg=F)
      }else{
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      }
      counter3 <- counter3+1
      
      #}else{
      #  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      #  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      #}
    } 
    
    # }
    
    
    #counter3 <- counter3+1
  }
}

dev.off()


##################
# Basic univariate plots
##################


svg("pdp_pasture4.svg",8.5,8)

varndx = 1
var = varlist2[[season]]$var[varndx]
var2 = varlist2[[season]]$readable[varndx]
#mod = bestmods[[season]]@model

ordr2 <- ordr
#ordr2 <- ordr2[!ordr2[,"season"]=="summer",]

laymat <- matrix(1:30,nrow=6,byrow = T)
laymat[6,] = 40
laymat[2:5,2] = 39
laymat[] <- as.numeric(as.factor(laymat))
laymat[6,3:5] <- 24

layout(laymat,widths = c(1.5,1.2,4,4,4),heights=c(1,4,4,4,4,1))          
#layout.show(24)
onleft <- c(6,10,14,18)
ontop <- c(1,2,3,4,5,23)
axlabs <- c(22,24)
counter=1
counter2=1
counter3=1
panel=0
panel=panel+1
# panel=22
#sapply(1:21,function(t) plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch=""))
while(panel < 25){
  if(panel %in% onleft){
    par(mai=c(0,0,0,0))
    plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    text(0,0.5, tools::toTitleCase(ordr2$season[counter3]),adj=0,cex=1.5)
    panel=panel+1
  }else if (panel %in% ontop){
    if(panel%in%c(1,2,23)){
      par(mai=c(0,0,0,0))
      plot(0,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      panel=panel+1
    }else{
      par(mai=c(0,0,0,0))
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5, c("Collar","Count","P(1+)")[counter2],adj=0.5,cex=1.5)
      counter2=counter2+1
      panel=panel+1
    }
  }else if(!panel%in%axlabs){
    #par(mai=base$mai)
    par(mai=c(0.3,0.4,0.0,0.0))
    season <- ordr2$season[counter3]
    thisdf <- dfs[[ordr2$dat[counter3]]][[season]]
    tmp <- varlist2[[season]]$var   #   as.character(thisres[[season]]$names[grepl("std",thisres[[season]]$names)])
    allmains <- tmp[!grepl(":",tmp)]
    #if(season!="summer"){
    if(ordr2$dat[counter3]==2){
      if(!is.na(allcoefrsf[[season]][var])){
        VisualizeRelation_col(data=thisdf,model=bestmods_collar2[[season]],predvar=var,varname=var2,allvars = allmains, season=season,svg=F)      # run and save partial dependence plots for all top variable
      }else{
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      }
      counter3 <- counter3+1
      panel=panel+1
      
    }else{
      #if(season!="parturition"){
      if(!all(is.na(c(allcoefs[[season]][var],
                      allcoefsz[[season]][var])))){
        VisualizeRelation(data=thisdf,model=bestmods[[season]]@model,predvar=var,varname=var2,allvars = allmains, season=season,svg=F)      # run and save partial dependence plots for all top variable
        
      }else{
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      }
      counter3 <- counter3+1
      panel=panel+2
      
      #}else{
      #  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      #  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      #}
    }
      
  }else{
    if(panel==22){
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5,"Predicted outcome (see caption for details)",srt=90,adj=0.5,cex=1.5)
      panel=panel+1
    }else{
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5,"Pasture (% Grazed Land)",adj=0.5,cex=1.5)
      panel=panel+1
    }
  }
}

dev.off()




svg("pdp_cowdens4.svg",8.5,8)

varndx = 2
var = varlist2[[1]]$var[varndx]
var2 = varlist2[[1]]$readable[varndx]
#mod = bestmods[[season]]@model

ordr2 <- ordr
#ordr2 <- ordr2[!ordr2[,"season"]=="summer",]

laymat <- matrix(1:30,nrow=6,byrow = T)
laymat[6,] = 40
laymat[2:5,2] = 39
laymat[] <- as.numeric(as.factor(laymat))
laymat[6,3:5] <- 24

layout(laymat,widths = c(1.5,1.2,4,4,4),heights=c(1,4,4,4,4,1))          
#layout.show(24)
onleft <- c(6,10,14,18)
ontop <- c(1,2,3,4,5,23)
axlabs <- c(22,24)
counter=1
counter2=1
counter3=1
panel=0
panel=panel+1
# panel=22
#sapply(1:21,function(t) plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch=""))
while(panel < 25){
  if(panel %in% onleft){
    par(mai=c(0,0,0,0))
    plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    text(0,0.5, tools::toTitleCase(ordr2$season[counter3]),adj=0,cex=1.5)
    panel=panel+1
  }else if (panel %in% ontop){
    if(panel%in%c(1,2,23)){
      par(mai=c(0,0,0,0))
      plot(0,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      panel=panel+1
    }else{
      par(mai=c(0,0,0,0))
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5, c("Collar","Count","P(1+)")[counter2],adj=0.5,cex=1.5)
      counter2=counter2+1
      panel=panel+1
    }
  }else if(!panel%in%axlabs){
    #par(mai=base$mai)
    par(mai=c(0.3,0.4,0.0,0.0))
    season <- ordr2$season[counter3]
    var = varlist2[[season]]$var[varndx]
    var2 = varlist2[[season]]$readable[varndx]
    thisdf <- dfs[[ordr2$dat[counter3]]][[season]]
    tmp <- varlist2[[season]]$var   #   as.character(thisres[[season]]$names[grepl("std",thisres[[season]]$names)])
    allmains <- tmp[!grepl(":",tmp)]
    #if(season!="summer"){
    if(ordr2$dat[counter3]==2){
      if(!is.na(allcoefrsf[[season]][var])){
        VisualizeRelation_col(data=thisdf,model=bestmods_collar2[[season]],predvar=var,varname=var2,allvars = allmains, season=season,svg=F)      # run and save partial dependence plots for all top variable
      }else{
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      }
      counter3 <- counter3+1
      panel=panel+1
      
    }else{
      #if(season!="parturition"){
      if(!all(is.na(c(allcoefs[[season]][var],
                      allcoefsz[[season]][var])))){
        VisualizeRelation(data=thisdf,model=bestmods[[season]]@model,predvar=var,varname=var2,allvars = allmains, season=season,svg=F)      # run and save partial dependence plots for all top variable
        
      }else{
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
        plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      }
      counter3 <- counter3+1
      panel=panel+2
      
      #}else{
      #  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      #  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      #}
    }
    
  }else{
    if(panel==22){
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5,"Predicted outcome (see caption for details)",srt=90,adj=0.5,cex=1.5)
      panel=panel+1
    }else{
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5,"Cattle Density",adj=0.5,cex=1.5)
      panel=panel+1
    }
  }
}

dev.off()




#########
# practice code


#summary(bestmod2)       # move visualizations to a different part of code

season = allseasons2[3]
varndx = 2
var = varlist2[[season]]$var[varndx]
var2 = varlist2[[season]]$readable[varndx]
#mod = bestmods[[season]]@model
mod = bestmods_collar2[[season]]
data = dfs[[2]][[season]]

tmp <- rownames(coef(summary(mod))$cond)[-1]
topvars <- tmp[!grepl(":",tmp)]
ints <- strsplit(tmp[grep(":",tmp)],":")

if(!is.na(allcoefrsf[[season]][var])){
  VisualizeRelation_col(data=data,model=mod,predvar=var,varname=var2,allvars = topvars, season=season,svg=F)      # run and save partial dependence plots for all top variable
}else{
  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
}

if(!is.na(allcoefs[[season]][var])){
  VisualizeRelation(data=data,model=mod,predvar=var,varname=var2,allvars = topvars, season=season,svg=F)      # run and save partial dependence plots for all top variable
}else{
  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
}




############
# End of script
############





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
#test



# #######
# # .. and for collar data
# 
# allseasons <- names(useds)
# 
# allmeans <- list()
# 
# thisres <- list()
# 
# season <- allseasons[2]
# 
# rarefy <- min(sapply(useds,nrow))
# 
# 
# for(season in allseasons){
#   
#   df <- useds[[season]]
#   
#   #  sapply(final.collar,function(t) sapply(t,class))
#   table(df$ID)
#   
#   ############
#   #subset by ID for indvidual level analysis
#   ############
#   
#   indnames <- c('24034B','31710A','31710B','31711A','31711B','31713A','31713B','31739A')
#   indsex <- rep(c("M","F"),times=c(3,5))
#   names(indsex) <- indnames
#   
#   allinds <- list()
#   test <- lapply(indnames,function(t) allinds[[t]] <<- subset(df,ID==t) )
#   
#   
#   # remove individuals with too few observations
#   indnames <- indnames[sapply(allinds,nrow)>100]
#   allinds <- allinds[indnames]
#   table(allinds$`24034B`$ID)
#   
#   
#   nboot <- 50
#   coefs <- list()
#   
#   ### loop through boot samples
#   b=1
#   for(b in 1:nboot){
#     nobs <- sapply(allinds,nrow)
#     
#     nbackground <- nrow(avails[[season]])
#     #t <- indnames[1]
#     avails2 <- lapply(indnames, function(t) 
#       avails[[season]][sample(1:nrow(avails[[season]]),
#       size=min(nrow(avails[[season]]),nobs[t]*5),replace=FALSE),])  # scale the number of background points to the number of observations  
#     names(avails2) <- indnames
#     
#          # bootstrap step
#     t=indnames[1]
#     allinds2 <- lapply(indnames,function(t) rbind(allinds[[t]][sample(1:nrow(allinds[[t]]),nrow(allinds[[t]]),    # 
#                   replace = TRUE),], avails2[[t]]) )
#     names(allinds2) <- indnames
#     
#     allmods <- lapply(indnames,build_indmodels,season=season,df=allinds2)    # run models
#     names(allmods) <- indnames
#     
#     #summary(allmods[[1]])
#     
#     coefs[[b]] <- sapply(allmods,coefficients) 
#     cat(sprintf("bootstrap %s\n",b))
#   }
#   
#   allmeans[[season]] <- sapply(coefs,function(t)  apply(t,1,mean))  
#   thisres[[season]] = data.frame(
#     meancoefs = apply(allmeans[[season]],1,mean),
#     secoefs = apply(allmeans[[season]],1,sd),                       
#     uppr = apply(allmeans[[season]],1,function(t) quantile(t,0.95)),
#     lowr = apply(allmeans[[season]],1,function(t) quantile(t,0.05)),     
#     names=rownames(allmeans[[season]])
#   )
#   
#   # svg(sprintf("collar_coefs_all_%s.svg",season),width = 4,height=5)
#   # 
#   # p<- ggplot(thisres,aes(x=meancoefs, y=names)) + 
#   #   geom_point() +
#   #   ggtitle(season) +
#   #   ylab("") +
#   #   xlab("Mean coefficient estimate") +
#   #   geom_vline(xintercept = 0,linetype=1, size=1) +
#   #   geom_errorbarh(aes(xmin=lowr, xmax=uppr),height=0.2)
#   # print(p)
#   # # ggsave(p, filename = paste("coef.plot", season, "pdf", sep = "."), height = 11, width = 8)
#   # 
#   # dev.off()
#   
#   #save coefficients table
#   #  write.csv(thisres, file = paste("coef.table", season, "csv", sep = "."))
#   
#   #save AIC
#   #  sink(file = paste("aic", season, "txt", sep = "."))
#   
#   #  AIC(mdl)
#   #  sink()
#   
#   
# #   ###########
# #   # females only            # KTS: do we need to keep this- didn't seem like it was showing anything too interesting, and it's not as analogous with the count models this way
# #   ###########
# #   
# #   ndx <-  intersect(indnames, names(which(indsex=="F")))
# #   
# #   allmeans2 <- sapply(coefs,function(t)  apply(t[,ndx],1,mean))
# #   
# #   # thisres_f <- data.frame(
# #   #   meancoefs=sapply(1:nrow(coefs),function(t) wtd.mean(coefs[t,][ndx],1/(ses[t,][ndx])^2)),
# #   #   secoefs=sapply(1:nrow(coefs),function(t) sqrt(wtd.var(coefs[t,][ndx],1/(ses[t,][ndx])^2))/sqrt(length(indnames[ndx]))),
# #   #   tcrit=qt(0.95,df=length(indnames)),
# #   #   names=rownames(coefs)
# #   # # )
# #   # thisres$up <- thisres$meancoefs+thisres$secoefs*thisres$tcrit
# #   # thisres$low <- thisres$meancoefs-thisres$secoefs*thisres$tcrit
# #   
# #   thisres_f = data.frame(
# #     meancoefs = apply(allmeans2,1,mean),
# #     secoefs = apply(allmeans,1,sd),
# #     uppr = apply(allmeans2,1,function(t) quantile(t,0.95)),
# #     lowr = apply(allmeans2,1,function(t) quantile(t,0.05)),
# #     names=rownames(allmeans2)
# #   )
# #   
# #   svg(sprintf("collar_coefs_fem_%s.svg",season),width = 4,height=5)
# #   
# #   p<- ggplot(thisres_f,aes(x=meancoefs, y=names)) + 
# #     geom_point() +
# #     ggtitle(season) +
# #     ylab("") +
# #     xlab("Mean coefficient estimate") +
# #     geom_vline(xintercept = 0,linetype=1, size=1) +
# #     geom_errorbarh(aes(xmin=lowr, xmax=uppr),height=0.2)
# #   print(p)
# #   
# #   dev.off()
# #                  # move visualizations
#   # tmp <- as.character(thisres[[season]]$names[grepl("std",thisres[[season]]$names)])
#   # allmains <- tmp[!grepl(":",tmp)]
#   # allints <- tmp[grepl(":",tmp)]
#   # tmp <- sapply(1:length(allmains), function(t) VisualizeRelation_col(data=df,meancoefs=thisres[[season]],predvar = allmains[t],allvars=allmains,season=season) )     # run and save partial dependence plots for all top variables
#   # 
#   # ints <- strsplit(allints[grep(":",allints)],":")
#   # tmp <- sapply(1:length(ints), function(t) VisualizeInteraction_col(data=df,meancoefs=thisres[[season]],var1=ints[[t]][1],var2=ints[[t]][2],allvars=allmains,season=season) )
#   cat(sprintf("done with season %s\n",season))
# }
# 
# save(allmeans,thisres,file="collarmods_allseasons.RData")



# ## this function takes an individual name and runs a logistic regression. If the model throws a warning due to complete separation, a
# # firth-corrected (penalized likelihood) version is run.
# 
# build_indmodels <- function(thisind=indnames[2],season="mating",df=allinds2){
#   
#   formula <- as.formula(sprintf("Used ~ std.pasture + std.%s + std.ponds + std.elevation + std.ndvi + std.scrub + std.wetgrass + std.ag + std.slope + std.aspect + std.ndvi:std.scrub + std.ndvi:std.wetgrass + std.%s:std.ndvi + std.slope:std.aspect",best.stocking[season],best.stocking[season])) # 
#   
#   dat=df[[thisind]]
#   
#   dat$wts <- 1000
#   dat$wts[dat$Used==1] <- 1
#   
#   full <- tryCatch( {     
#     glm(formula,
#         data= dat, family=binomial(link= "logit"),weights = wts)
#   }
#   , warning = function(w) { 
#     
#     f <- logistf::logistf(formula,
#                           data= dat, family=binomial(link= "logit"),weights = wts)
#     return(f)
#   })
#   
#   return(full)
# }
# 


# 
# VisualizeInteraction_col <- function(data=dfs[[2]],meancoefs=thisres$mating,var1=best.stocking["mating"],var2="std.ndvi",allvars=allvars,season,svg=T){
#   len <- 20     # increase this for higher-res figures
#   
#   dataclasses <- sapply(data,class)
#   
#   
#   var1_2 <- strsplit(var1,"\\.")[[1]][2]
#   var2_2 <- strsplit(var2,"\\.")[[1]][2]
#   
#   
#   realmean1 <- mean(data[[var1_2]])
#   realsd1 <- sd(data[[var1_2]])
#   realmean2 <- mean(data[[var2_2]])
#   realsd2 <- sd(data[[var2_2]])
#   standvar1 <- var1
#   standvar2 <- var2
#   
#   
#   firstdim <- data[,standvar1]
#   seconddim <- data[,standvar2]
#   range1 <- seq(min(firstdim),max(firstdim),length=len)
#   range2 <- seq(min(seconddim),max(seconddim),length=len)
#   newdata <- expand.grid(range1,range2)
#   # head(newdata,50)
#   names(newdata) <- c(standvar1,standvar2)
#   
#   othervars <- allvars[!allvars%in%c(standvar1,standvar2,"Used")]
#   
#   var = othervars[1]
#   for(var in othervars){
#     thisvar <- data[[var]]
#     if(is.factor(thisvar)){
#       tab <- table(thisvar)
#       vals <- names(tab)
#       levs <- levels(thisvar)
#       mostcom <- vals[which.max(tab)]
#       newvec <- factor(rep(mostcom,times=nrow(newdata)),levels=levs)
#       newdata[[var]] <- newvec
#     }else{
#       newdata[[var]] <- mean(thisvar)
#     }
#   }
#   
#   int <- meancoefs$meancoefs[meancoefs$names=="(Intercept)"]
#   inter <- meancoefs$meancoefs[grepl(sprintf("%s:%s",var1,var2),meancoefs$names)]
#   
#   ndx <- match(names(newdata),meancoefs$names)
#   
#   pred <- sapply(1:nrow(newdata), function(s)              # make prediction plots based on population-averaged model
#     plogis(int + as.matrix(newdata[s,]) %*% as.matrix(meancoefs$meancoefs[ndx]) + inter*(newdata[s,var1]*newdata[s,var2])) )
#   
#   predmat <-  matrix(pred,nrow=len,ncol=len)
#   
#   
#   if(svg){
#     svg(sprintf("intplots_col_%s_%s_%s.svg",var1_2,var2_2,season),width=4.5,height = 6)
#     
#     par(mai=c(0.5,0.5,0,0))
#   }else{
#     var1_2 <- varlist2[[season]]$readable[which(varlist2[[season]]$var==var1)] #strsplit(var1,"\\.")[[1]][2]
#     var2_2 <- varlist2[[season]]$readable[which(varlist2[[season]]$var==var2)] #strsplit(var2,"\\.")[[1]][2]
#     
#   }
#   
#   x.axis <- realmean1+realsd1*2*range1
#   min.x <- min(x.axis)
#   max.x <- max(x.axis)
#   xrange <- max.x-min.x
#   x.axis2 <- seq(min.x,max.x,length=5)[c(2:4)]
#   y.axis <- realmean2+realsd2*2*range2
#   min.y <- min(y.axis)
#   max.y <- max(y.axis)
#   yrange <- max.y-min.y
#   y.axis2 <- seq(min.y,max.y,length=5)[c(2:4)]
#   z.axis <- seq(min(predmat), max(predmat), length=5)[c(2:4)]
#   min.z <- min(predmat)
#   max.z <- max(predmat)
#   zrange <- max.z-min.z
#   xlab=paste("",var1_2,sep="")
#   ylab=paste("",var2_2,sep="")
#   zlab="Selection"
#   persp(x.axis,y.axis,predmat,
#         axes=F,
#         theta = 55, phi = 40, r = sqrt(10), d = 3, 
#         ticktype = "detailed", mgp = c(4, 1, 0)) -> res
#   # lines(trans3d(x.axis, min.y, min.z, res) , col="black",lwd=2)
#   # lines(trans3d(max.x, y.axis, min.z, res) , col="black",lwd=2)
#   # lines(trans3d(min.x, min.y, z.axis, res) , col="black",lwd=2)
#   
#   labels <- sprintf("%.1f",x.axis2)
#   label.pos <- trans3d(x.axis2, (min.y - yrange*0.3), min.z, res)
#   text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=30, cex=1)
#   
#   labels <- xlab
#   label.pos <- trans3d(x.axis2[1]-xrange*0.3, (min.y - yrange*0.4), min.z, res)
#   text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=30+100+180, cex=1)
#   
#   labels <- sprintf("%.1f",y.axis2)
#   label.pos <- trans3d((max.x + xrange*0.1), y.axis2, min.z, res)
#   text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), cex=1) 
#   
#   labels <- ylab
#   label.pos <- trans3d((max.x + xrange*0.3), y.axis2[2], min.z, res)
#   text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), cex=1,srt=25)
#   
#   labels <- sprintf("%.1f",z.axis)
#   label.pos <- trans3d(min.x, (min.y - 0.08), z.axis, res)
#   text(label.pos$x, label.pos$y, labels=labels, adj=c(1, NA), cex=1)
#   
#   labels <- zlab
#   label.pos <- trans3d(min.x, (min.y - yrange*0.6), z.axis[3]+zrange*0.2, res)
#   text(label.pos$x, label.pos$y, labels=labels, adj=c(1, NA), cex=1,srt=100)
#   
#   
#   var1_2 <- strsplit(var1,"\\.")[[1]][2]
#   var2_2 <- strsplit(var2,"\\.")[[1]][2]
#   subst <- data[sample(1:nrow(data),1000,replace = F),]
#   points(trans3d(subst[[var1_2]], subst[[var2_2]], max(predmat)-0.02, pmat = res), col = gray(0.8), pch = 19, cex=0.5)
#   
#   
#   if(svg) dev.off() 
# }

# VisualizeRelation_col <- function(data=allcollar[[season]],meancoefs=thisres,predvar="std.ndvi",varname="ndvi",allvars=allvar_c,season,svg=F){
#   len <- 100
#   
#   predvar2 <- strsplit(predvar,"\\.")[[1]][2]
#   
#   dataclasses <- sapply(data,class)
#   
#   dim <- data[,predvar]
#   range <- seq(min(dim),max(dim),length=len)
#   
#   realmean <- mean(data[[predvar2]])
#   realsd <- sd(data[[predvar2]])
#   
#   newdata <- data.frame(temp=range)
#   names(newdata) <- c(predvar)
#   
#   othervars <- allvars[!allvars%in%c(predvar,"Used")]
#   
#   
#   var = othervars[5]
#   for(var in othervars){
#     thisvar <- data[[var]]
#     if(is.factor(thisvar)){
#       tab <- table(thisvar)
#       vals <- names(tab)
#       levs <- levels(thisvar)
#       mostcom <- vals[which.max(tab)]
#       newvec <- factor(rep(mostcom,times=nrow(newdata)),levels=levs)
#       newdata[,var] <- newvec
#     }else{
#       newdata[,var] <- 0 #mean(thisvar)
#     }
#   }
#   
#   ndx <- match(names(newdata),meancoefs$names)
#   
#   int <- meancoefs$meancoefs[meancoefs$names=="(Intercept)"]
#   coef <- meancoefs$meancoefs[meancoefs$names==predvar]
#   
#   pred <- sapply(1:nrow(newdata), function(t)              # make prediction plots based on population-averaged model
#     plogis(int + as.matrix(newdata[t,]) %*% as.matrix(meancoefs$meancoefs[ndx]) ) )
#   
#   m <- c(int=int,coef=coef)
#   v <- diag(c(meancoefs$secoefs[meancoefs$names=="(Intercept)"],meancoefs$secoefs[meancoefs$names==predvar]))
#   
#   se_pred <- sapply(1:nrow(newdata), function(t)
#     as.numeric(car::deltaMethod(m,sprintf("1/(1+exp(-1*(int+%s*coef)))",newdata[t,predvar]),vcov.=v)["SE"]) )
#   
#   svg(sprintf("pdplots_col_%s_%s.svg",predvar2,season),width=4.5,height = 4.5)
#   
#   tcrit <- qnorm(0.975) #qt(0.95,nrow(meancoefs))
#   
#   par(mai=c(0.9,0.8,0.1,0.1))
#   plot(range,pred,xlab=predvar2,ylab="Selection propensity",type="l",lwd=2,xaxt="n",ylim=c(0,min(1,max(pred)*1.4)))
#   points(range,pred+tcrit*se_pred,type="l",lty=2)
#   points(range,pred-tcrit*se_pred,type="l",lty=2)
#   ats <- seq(min(range),max(range),length=6)
#   axis(1,ats,labels = round(realmean+ats*2*realsd,2))
#   rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
#   
#   
#   dev.off()
# }

