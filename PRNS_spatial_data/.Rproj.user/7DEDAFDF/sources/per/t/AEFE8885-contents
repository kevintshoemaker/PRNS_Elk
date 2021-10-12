###############
# Satellite-image analysis for PRNS elk!
###############


#########
# Questions for Kelley Lacey etc.
#########


#########
# Research questions
#########

# do elk tend to keep a distance from cattle?
# do cattle tend to keep a distance from elk?

# what is the probability of an elk-cattle interaction under different seasons, circumstances (stocking rates etc)

# do elk tend to spend more time in natural areas when there are more cattle in the pastures? 

#########
# Clear workspace
#########

rm(list=ls())


#########
# Load packages
#########

library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(spatstat)
library(maptools)
library(MASS)
library(spatialEco)
library(adehabitatHR)


#########
# Global vars
#########

GVars <- list()

GVars$projection <- CRS("+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
GVars$projectionll <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 ")

#########
# Load functions
#########




#########
# Directories
#########

thisdir <- getwd()
rasterdir <- sprintf("%s/R files",thisdir)
rasterdir
shapefiledir <- sprintf("%s/Shapefiles",thisdir)
shapefiledir

#########
# Load raster data
#########

files <- list.files(path = rasterdir, pattern = "\\.rds$", full.names = TRUE)
rastnames <- sapply(sapply(strsplit(files,split=c("/")),function(t) strsplit(t[[grep("rds",t)]],split="\\.") ),function(l) l[1] )
temp <- table(rastnames)
rastnames <- unlist(sapply(1:length(temp),function(t) paste(rep(names(temp[t]),times=temp[t]),1:temp[t],sep="" ) ))

r <- lapply(files, readRDS)

########
# only keep the raster data I want

keep <- c(
  "buildings1",
  #"cows1",
  #"cows2",
  #"cows3",
  #"cows4",
  "DEM1",
  "eastness1",
  #"fences1",
  "northness1",
  "ponds1",
  "ranch1",
  "ranch2",
  "roads1",
  "slope1",
  "streams1",
  "wetlands1"
)

names(r) <- rastnames
rasterCovs <- r[keep]


rasterCovs <- stack(rasterCovs)


rm(r,keep,rastnames)

GVars$master_extent <- extent(rasterCovs$ranch1)


#########
# Visualize raster data


for(i in 1:length(rasterCovs)) plot(rasterCovs[[i]],main=names(rasterCovs)[i])    # plot all rasters


#######
# Load vector data
#######

covarsdir <- sprintf("%s/Covariates_shp",shapefiledir)
covarsdir

files <- list.files(path = covarsdir, pattern = "\\.shp$", full.names = FALSE)
files2 <- gsub(".shp","",files)

v <- lapply(1:length(files), function(t) readOGR(dsn=covarsdir,layer=files2[t]))  # p4s=as.character(projection)
names(v) <- files2

moredir <- sprintf("%s/Other_shp",shapefiledir)
moredir

files <- list.files(path = moredir, pattern = "\\.shp$", full.names = FALSE)
files2 <- gsub(".shp","",files)

v2 <- lapply(1:length(files), function(t) readOGR(dsn=moredir,layer=files2[t]))  # p4s=as.character(projection)
names(v2) <- files2

v <- c(v,v2)

satdir <- sprintf("%s/Satellite_annotations_shp",shapefiledir)
satdir

files <- list.files(path = satdir, pattern = "\\.shp$", full.names = FALSE)
files2 <- gsub(".shp","",files)

v2 <- lapply(1:length(files), function(t) readOGR(dsn=satdir,layer=files2[t]))  # p4s=as.character(projection)
names(v2) <- files2

v <- c(v,v2)

#### clean up workspace

rm(temp,files,files2,covarsdir,moredir,satdir,v2)


#######
# crop all vectors to the proper extent and reproject to same coordinate system

v <- lapply(v,function(t) spTransform(t,GVars$projection))

v <- lapply(v,function(t) crop(t,GVars$master_extent) )


#######
# separate response and predictors/helpers

names(v)

covarnames <- c(
  "ranch_boundary",
  "trails_and_roads_clip2",
  "wilderness",
  "female.recursions.072018",
  "joint_home_range1a",
  "male.recursions.072018",
  "survey.pts"
)

satnames <- setdiff(names(v),covarnames)

satnames

length(satnames)   # 31 images!



vectorCovs <- v[covarnames]
satPoints <- v[satnames]

rm(satnames,covarnames,v)
rm(i)

vectorCovs$survey.pts$X


######
# Remove satellite points with NA or uncertain designation

#satPoints$`2015March18_cows`@data

satPoints <- lapply(satPoints,function(t){ keep <- !is.na(t@data$Name); t[keep,]  })
satPoints <- lapply(satPoints,function(t){ keep <- t@data$Name%in%c("Elk","Livestock","Possible livestock","Blank"); t[keep,]  })

######
# Ensure that only sat points in ranches are considered further

plot(rasterCovs$ranch1)

i=1
for(i in 1:length(satPoints)){
  keep <- extract(rasterCovs$ranch1,satPoints[[i]])==1
  satPoints[[i]] <- satPoints[[i]][keep,]
}
  
rm(i,keep)  

#######
# Explore sat files

allsats <- names(satPoints)

i=1
for(i in 1:length(allsats)){
  temp <- satPoints[[allsats[i]]]@data
  elkndx <- temp$Name=="Elk"
  cowndx <- temp$Name%in%c("Livestock","Possible livestock","Blank")
  plot(rasterCovs$ranch1,main=allsats[i])
  plot(vectorCovs$wilderness,add=T)
  plot(satPoints[[allsats[i]]][cowndx,],pch=20,col="red",cex=1,add=T)
  plot(satPoints[[allsats[i]]][elkndx,],add=T)
}

graphics.off()

rm(temp,allsats,i,elkndx,cowndx)


########
# Make lists of only elk and cow points

elkPoints <- lapply(satPoints,function(t) t[t@data$Name=="Elk",] )
cowPoints <- lapply(satPoints,function(t) t[t@data$Name%in%c("Livestock","Possible livestock","Blank"),] )


########
# remove any images that have no elk!

keep <- sapply(elkPoints,function(t) nrow(t@data))!=0


# sapply(cowPoints,function(t) nrow(t@data))!=0

elkPoints <- elkPoints[keep]
cowPoints <- cowPoints[keep]

length(elkPoints)     # 28 images with both elk and cattle


# total elk points

sum(sapply(elkPoints,function(t) nrow(t@data)))
min(sapply(elkPoints,function(t) nrow(t@data)))

########
# Make polygons of elk and cow groups: note area and number of elk in each polygon

elkGroups <- list()
cowGroups <- list()

# elkmcps <- list()
# cowmcps <- list()

covnames <- names(rasterCovs)

i=10
for(i in 1:length(elkPoints)){
  
      # ELK
  buf1 <- gBuffer(elkPoints[[i]], width=50, byid=TRUE)
  buf2 <- gUnaryUnion(buf1)
  buf <-  sp::disaggregate(gBuffer(buf2, width=-25))     # was -40
  elkPoints[[i]]@data$Group <- over(elkPoints[[i]], buf)
  df <- data.frame(id1 = 1:length(buf))
  elkGroups[[i]] <- SpatialPolygonsDataFrame(buf,df)
  elkGroups[[i]]@data$Elk <- table(elkPoints[[i]]@data$Group)
  elkGroups[[i]]@data$Area <- sapply(elkGroups[[i]]@polygons,function(t) t@area )
  elkGroups[[i]]@data$centerx <- coordinates(elkGroups[[i]])[,1]
  elkGroups[[i]]@data$centery <- coordinates(elkGroups[[i]])[,2]
  
  # temp <- SpatialPointsDataFrame(elkPoints[[i]]@coords,data.frame(id=elkPoints[[i]]@data$Group))
  # mcp1 <- mcp(temp,percent=100)
  # elkmcps[[i]] <- mcp1
  
      # COWS
  buf1 <- gBuffer(cowPoints[[i]], width=50, byid=TRUE)
  buf2 <- gUnaryUnion(buf1)
  buf <-  sp::disaggregate(gBuffer(buf2, width=-25))
  cowPoints[[i]]@data$Group <- over(cowPoints[[i]], buf)
  df <- data.frame(id1 = 1:length(buf))
  cowGroups[[i]] <- SpatialPolygonsDataFrame(buf,df)
  cowGroups[[i]]@data$Cows <- table(cowPoints[[i]]@data$Group)
  cowGroups[[i]]@data$Area <- sapply(cowGroups[[i]]@polygons,function(t) t@area )
  cowGroups[[i]]@data$centerx <- coordinates(cowGroups[[i]])[,1]
  cowGroups[[i]]@data$centery <- coordinates(cowGroups[[i]])[,2]
  
  # temp <- SpatialPointsDataFrame(cowPoints[[i]]@coords,data.frame(id=cowPoints[[i]]@data$Group))
  # mcp1 <- mcp(temp,percent=100)
  # cowmcps[[i]] <- mcp1
  # 

  # Extract covariates for each polygon
  for(j in 1:length(covnames)){
    elkGroups[[i]]@data[,covnames[j]] <- extract(rasterCovs[[j]],elkGroups[[i]],fun=median)
    cowGroups[[i]]@data[,covnames[j]] <- extract(rasterCovs[[j]],cowGroups[[i]],fun=median)
  }
}

names(elkGroups) <- names(elkPoints)
names(cowGroups) <- names(cowPoints)

# names(elkmcps) <- names(elkPoints)
# names(cowmcps) <- names(cowPoints)

rm(buf,buf1,buf2,i,j)

i=1
for(i in 1:length(elkPoints)){
  #plot(cowPoints[[i]],pch="x",cex=0.2)
  temp <- extent(cowGroups[[i]])
  temp <- raster(temp,resolution=10,vals=NA)
  plot(temp,main=names(cowGroups[i]))
  plot(cowGroups[[i]],lwd=0.1,add=T)
  plot(raster::buffer(elkGroups[[i]],50),col="green",add=T,lwd=3)
  plot(elkPoints[[i]],pch=20,cex=0.2,add=T,col="red")
}

rm(i,df,temp)

graphics.off()


#######
# Distance to wildland?

# apply(gDistance(spts, columbus,byid=TRUE),2,min)

# ?gDistance
# 
# plot(vectorCovs$wilderness,col="red")
# plot(elkPoints[[1]],add=T)

i=1
for(i in 1:length(elkPoints)){
  elkPoints[[i]]@data$distWild <- gDistance(elkPoints[[i]],vectorCovs$wilderness,byid = TRUE)[1,]
  elkGroups[[i]]@data$distWild <- gDistance(elkGroups[[i]],vectorCovs$wilderness,byid = TRUE)[1,]
  cowPoints[[i]]@data$distWild <- gDistance(cowPoints[[i]],vectorCovs$wilderness,byid = TRUE)[1,]
  cowGroups[[i]]@data$distWild <- gDistance(cowGroups[[i]],vectorCovs$wilderness,byid = TRUE)[1,]
}


#plot(elkmcps[[1]])


#######
# fit inhomogeneous clustered point process model for both elk and cattle
# use Matern point process model

plot(vectorCovs$ranch_boundary)
plot(rasterCovs$ranch1,add=T)
plot(cowPoints$`2015April25`,add=T)
plot(cowPoints$`2017June14`,add=T)

######
# Develop a mask for the ppm

temp <- unionSpatialPolygons(vectorCovs$ranch_boundary,rep(1,times=length(vectorCovs$ranch_boundary)))
temp <- sp::disaggregate(gBuffer(temp, width=10))
plot(temp)
GVars$vector_mask <- temp

GVars$raster_mask <- rasterCovs$ranch1
GVars$raster_template <- reclassify(GVars$raster_mask,rcl=c(-Inf,Inf,NA))

plot(GVars$raster_mask)
plot(GVars$vector_mask,add=T)

GVars$raster_mask2 <- rasterize(GVars$vector_mask,GVars$raster_template)
GVars$raster_mask2 <- reclassify(GVars$raster_mask2,rcl=c(-Inf,0.5,NA, 0.5,Inf,1))

window <- as.owin(temp)

plot(GVars$raster_mask2)
plot(GVars$raster_mask,add=T)
plot(GVars$vector_mask,add=T)


rm(temp)

######
# Set up the points for ppm

elkppp <- lapply(1:length(elkPoints),function(t)  ppp(elkPoints[[t]]@coords[,1],elkPoints[[t]]@coords[,2],window)   )
cowppp <- lapply(1:length(cowPoints),function(t)  ppp(cowPoints[[t]]@coords[,1],cowPoints[[t]]@coords[,2],window)   )

names(elkppp) <- names(elkPoints)
names(cowppp) <- names(cowPoints)

#######
# Develop a merged shapefile for all elk and cattle observations

allx <- unlist(lapply(1:length(elkppp),function(t) elkppp[[t]]$x ))
ally <- unlist(lapply(1:length(elkppp),function(t) elkppp[[t]]$y ))
allelkppp <- ppp(allx,ally,window)

allx <- unlist(lapply(1:length(cowppp),function(t) cowppp[[t]]$x ))
ally <- unlist(lapply(1:length(cowppp),function(t) cowppp[[t]]$y ))
allcowppp <- ppp(allx,ally,window)


plot(allcowppp,main="All cattle and elk points")
plot(allelkppp,add=T,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.2))
legend("bottomright",bty="n",pch=c(1,1),col=c("black","red"),legend=c("Cattle","Elk"),cex=1.2)


######
# Density contours for elk and cattle, showing areas of overlap?

?kde2d

elkkernel <- kde2d(x=allelkppp$x,y=allelkppp$y,n=500,lims=c(GVars$master_extent@xmin,GVars$master_extent@xmax,GVars$master_extent@ymin,GVars$master_extent@ymax))
cowkernel <- kde2d(x=allcowppp$x,y=allcowppp$y,n=500,lims=c(GVars$master_extent@xmin,GVars$master_extent@xmax,GVars$master_extent@ymin,GVars$master_extent@ymax))

elk_kde2 <- elkkernel$z/sum(elkkernel$z)
sum(elk_kde2)

cow_kde2 <- cowkernel$z/sum(cowkernel$z)
sum(cow_kde2)

#######
# find key contours of kernel density distribution

cover <- function(z,k=0.95){
  function(x){
    sum(z[z>x])-k
  }
}
z=elk_kde2
sum(z[z>0.00000000001])

contours <- c(0.25,0.5,0.75,0.95)
elkcontours <- c()
cowcontours <- c()

for(cont in contours){
  c1 <- cover(elk_kde2,cont)
  elkcontours <- c(elkcontours, uniroot(c1,interval=c(1e-9,max(elk_kde2)),tol=1e-9,trace=10)$root)
  c1 <- cover(cow_kde2,cont)
  cowcontours <- c(cowcontours, uniroot(c1,interval=c(1e-9,max(elk_kde2)),tol=1e-9,trace=10)$root)
}

names(elkcontours) <- contours
names(cowcontours) <- contours

rm(contours)


elkcontours <- elkcontours*sum(elkkernel$z)
cowcontours <- cowcontours*sum(cowkernel$z)



#1/sum(elkkernel$z)

svg("kernel overlap.svg",7,3)

layout(matrix(1:2,nrow=1))

par(mai=c(0.75,0.75,0,0))

plot(reclassify(rasterCovs$ranch1,rcl=c(-Inf,0.1,NA, 0.1,Inf,1)),legend=FALSE,col=gray(0.8))

contour(cowkernel,add=T,levels=cowcontours,lwd=c(1,1,1,2),labels="")
contour(elkkernel,add=T,col="red",levels=elkcontours,lwd=c(1,1,1,2))

plot(allcowppp,add=T,cex=0.1,col=rgb(red = .1, green = .1, blue = .1, alpha = 0.01),pch=19)
plot(allelkppp,add=T,cex=0.1,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.01),pch=19)


plot(reclassify(rasterCovs$ranch1,rcl=c(-Inf,0.1,NA, 0.1,Inf,1)),legend=FALSE,col=gray(0.8))
# plot(allcowppp,add=T,cex=0.1,col=rgb(red = .1, green = .1, blue = .1, alpha = 0.05),pch=19)
# plot(allelkppp,add=T,cex=0.1,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.05),pch=19)

plot(allcowppp,add=T,col=rgb(red = .1, green = .1, blue = .1, alpha = 0.02),pch=20,cex=0.5)
plot(allelkppp,add=T,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.05),pch=20,cex=0.5)
legend("bottomright",bty="n",pch=c(1,1),col=c("black","red"),legend=c("Cattle","Elk"),cex=1.2)

# contour(cowkernel,add=T,lwd=1.5,labels="")
# contour(elkkernel,add=T,col="red",nlevels=4,lwd=2)

# plot(allcowppp,main="")
# plot(allelkppp,col="red",add=T)

dev.off()

graphics.off()
#####
# Plot each image

svg("allpoints.svg",8,8)

layout(matrix(1:9,nrow=3))
par(mai=c(0,0,0.5,0))

i=3
for(i in 1:9){      # plot only the first 9 images
  plot(elkppp[[i]],main=names(elkppp)[i],pch="+",cex=0.1)
  contour(elkkernel,add=T,col="red",levels=elkcontours[4],lwd=2)
  plot(elkppp[[i]],add=T,col="red",pch="+",cex=1.4)
  plot(cowppp[[i]],add=T)
}

dev.off()


######
# real statistics: (1) mean distance to nearest cow, (2) minimum elk-cow distance, (3) percent overlap between elk groups and cow groups   

mean_mindist_obs <- sapply(c(1:length(elkPoints)),function(t) mean(apply(gDistance(elkPoints[[t]],cowPoints[[t]],byid = T),2,min)))
min_mindist_obs <- sapply(c(1:length(elkPoints)),function(t) min(apply(gDistance(elkPoints[[t]],cowPoints[[t]],byid = T),2,min)))

range(mean_mindist_obs)
range(min_mindist_obs)

mean_mindist2_obs <- mean(mean_mindist_obs)
min_mindist2_obs <- mean(min_mindist_obs)

#  real statistic : percent of elk groups that overlap with cow groups

real_overlap <- sapply(1:length(elkGroups), function(t){ 
  sum(apply(gIntersects(elkGroups[[t]],cowGroups[[t]],byid=TRUE),2,sum))/length(elkGroups[[t]]) } ) 

real_overlap2 <- mean(real_overlap)


######
# Random sampler- bootstrap!

nreps <- 200

mean_mindist <- list()
min_mindist <- list()

per_overlap <- list()

rep=1
for(rep in 1:nreps){
  null_centers <- expand.grid(elkkernel$x,elkkernel$y)
  names(null_centers) <- c("x","y")
  
  nullgrp <- list()
  mean_mindist[[rep]] <- numeric(length(cowPoints))
  min_mindist[[rep]] <- numeric(length(cowPoints))
  
  per_overlap[[rep]] <- numeric(length(cowPoints))
  
  thisimg <- 1
  for(thisimg in 1:length(cowPoints)){
    test <- FALSE
    while(!test){    # NOTE: maybe just sample centers from 
      centers <- sample(1:nrow(null_centers),size=length(elkGroups[[thisimg]]),replace=T,prob=as.numeric(elkkernel$z))
      test <- all(extract(GVars$raster_mask,SpatialPoints(null_centers[centers,]))==1)
    }
    xresids <- elkPoints[[thisimg]]@coords[,1] - elkGroups[[thisimg]]@data$centerx[elkPoints[[thisimg]]@data$Group]
    yresids <- elkPoints[[thisimg]]@coords[,2] - elkGroups[[thisimg]]@data$centery[elkPoints[[thisimg]]@data$Group]
    #hist(c(xresids,yresids))
    
    nullgrp[[thisimg]] <- data.frame(x=1:length(elkPoints[[thisimg]]@data$Group),y=NA)
    
    i=1
    for(i in 1:length(elkPoints[[thisimg]]@data$Group)){
      test=FALSE
      while(!test){
        propx <- null_centers[centers,]$x[elkPoints[[thisimg]]@data$Group][i]+c(sample(xresids,1))
        propy <- null_centers[centers,]$y[elkPoints[[thisimg]]@data$Group][i]+c(sample(yresids,1))
        prop <- data.frame(x=propx,y=propy)
        test <- extract(GVars$raster_mask,prop)==1
        test <- ifelse(is.na(test),FALSE,test)
      }
      nullgrp[[thisimg]][i,] <- prop[1,]
    }
    
    nullgrp[[thisimg]] <- SpatialPoints(nullgrp[[thisimg]],proj4string = GVars$projection)
    
    buf1 <- gBuffer(nullgrp[[thisimg]], width=50, byid=TRUE)
    buf2 <- gUnaryUnion(buf1)
    buf <-  sp::disaggregate(gBuffer(buf2, width=-40))
    
    per_overlap[[rep]][[thisimg]] <- 
      sum(apply(gIntersects(buf,cowGroups[[thisimg]],byid=TRUE),2,sum))/length(elkGroups[[thisimg]])
    
    #extract(GVars$raster_mask,nullgrp)
    
    mean_mindist[[rep]][thisimg] <- mean(apply(gDistance(nullgrp[[thisimg]],cowPoints[[thisimg]],byid = T),2,min))
    min_mindist[[rep]][thisimg] <- min(apply(gDistance(nullgrp[[thisimg]],cowPoints[[thisimg]],byid = T),2,min))
  }
  
}


names(nullgrp) <- names(elkPoints)

mean_mindist2 <- t(as.matrix(as.data.frame(mean_mindist)))
rownames(mean_mindist2) <- NULL

min_mindist2 <- t(as.matrix(as.data.frame(min_mindist)))
rownames(min_mindist2) <- NULL

mean_peroverlap2 <- t(as.matrix(as.data.frame(per_overlap)))
rownames(mean_peroverlap2) <- NULL

mean_mindist_sim <- apply(mean_mindist2,2,mean)
min_mindist_sim <- apply(min_mindist2,2,mean)

mean_peroverlap <- apply(mean_peroverlap2,2,mean)

graphics.off()
svg("grandmeans.svg",4,5)

layout(matrix(1:3,nrow=3))
par(mai=c(0.6,0.6,0,0.1))

meanz <- sapply(mean_mindist,mean)
p1 <- length(which(meanz>mean_mindist2_obs))/nreps
hist(meanz,xlab="Mean distance from elk to nearest cattle (m)",ylab="Proportion",freq=F,
                     ylim=c(0,0.011),xlim=c(0,850),main="")    # xlim=c(200,1500),
lines(density(meanz),col=gray(0.4),lwd=3,lty=2)
segments(mean_mindist2_obs,0,mean_mindist2_obs,0.008,lwd=4,col="blue")
#abline(v=mean_mindist2_obs,lwd=4,col="blue")
text(100,0.003,sprintf("p-value: %s",p1),cex=1.1)
text(10,0.009,"a)",cex=1.5)

meanz <- sapply(min_mindist,mean)
p1 <- length(which(meanz>min_mindist2_obs))/nreps
hist(meanz,xlab="Min distance from elk to nearest cattle (m)",ylab="Proportion",freq=F,
                    ylim=c(0,0.011),xlim=c(0,850),main="")    # xlim=c(0,800),
lines(density(meanz),col=gray(0.4),lwd=3,lty=2)
segments(min_mindist2_obs,0,min_mindist2_obs,0.008,lwd=4,col="blue")
#abline(v=min_mindist2_obs,lwd=4,col="blue")
text(600,0.003,sprintf("p-value: %s",p1),cex=1.1)
text(10,0.009,"b)",cex=1.5)

meanz <- sapply(per_overlap,mean)
p1 <- length(which(meanz<=real_overlap2))/nreps
hist(meanz,xlim=c(0,0.25),xlab="Percent elk-cattle group overlap",ylab="Proportion",freq=F,ylim=c(0,20),main="")
lines(density(meanz),col=gray(0.4),lwd=3,lty=2)
segments(real_overlap2,0,real_overlap2,18,lwd=4,col="blue")
#abline(v=real_overlap2,lwd=4,col="blue")
text(0.15,12,sprintf("p-value: %s",p1),cex=1.1)
text(0.003,18,"c)",cex=1.5)

dev.off()


graphics.off()
svg("imagemeans1.svg",8,8)

layout(matrix(1:9,nrow=3))

i=1
for(i in 1:9){
  meanz <- mean_mindist2[,i]
  hist(meanz,xlim=c(0,1500),xlab="Mean distance from elk to nearest cow",
       ylab="Proportion",freq=F,ylim=c(0,0.003),main=names(elkPoints)[i])
  lines(density(meanz),col=gray(0.4),lwd=3,lty=2)
  abline(v=mean_mindist_obs[i],lwd=4,col="blue")
  p1 <- length(which(meanz>mean_mindist_obs[i]))/nreps
  text(1000,0.002,sprintf("p-value: %s",p1))
  #text(200,0.006,"a)",cex=2)
}

dev.off()

graphics.off()
svg("imagemeans2.svg",8,8)

layout(matrix(1:9,nrow=3))

i=1
for(i in 10:18){
  meanz <- mean_mindist2[,i]
  hist(meanz,xlim=c(0,1500),xlab="Mean distance from elk to nearest cow",
       ylab="Proportion",freq=F,ylim=c(0,0.003),main=names(elkPoints)[i])
  lines(density(meanz),col=gray(0.4),lwd=3,lty=2)
  abline(v=mean_mindist_obs[i],lwd=4,col="blue")
  p1 <- length(which(meanz>mean_mindist_obs[i]))/nreps
  text(1000,0.002,sprintf("p-value: %s",p1))
  #text(200,0.006,"a)",cex=2)
}

dev.off()

graphics.off()
svg("imagemeans3.svg",8,8)

layout(matrix(1:9,nrow=3))

i=1
for(i in 19:27){
  meanz <- mean_mindist2[,i]
  hist(meanz,xlim=c(0,1500),xlab="Mean distance from elk to nearest cow",
       ylab="Proportion",freq=F,ylim=c(0,0.003),main=names(elkPoints)[i])
  lines(density(meanz),col=gray(0.4),lwd=3,lty=2)
  abline(v=mean_mindist_obs[i],lwd=4,col="blue")
  p1 <- length(which(meanz>mean_mindist_obs[i]))/nreps
  text(1000,0.002,sprintf("p-value: %s",p1))
  #text(200,0.006,"a)",cex=2)
}

dev.off()


svg("imageoverlaps.svg",8,8)

layout(matrix(1:9,nrow=3))

i=1
for(i in 1:9){
  meanz <- mean_peroverlap2[,i]
  hist(meanz,xlim=c(0,1),xlab="Mean group elk-cattle overlap",
       ylab="Proportion",freq=F,ylim=c(0,45),main=names(elkPoints)[i])
  lines(density(meanz),col=gray(0.4),lwd=3,lty=2)
  abline(v=real_overlap[i],lwd=4,col="blue")
  p1 <- length(which(meanz==real_overlap[i]))/nreps
  text(0.6,20,sprintf("p-value: %s",p1))
  #text(200,0.006,"a)",cex=2)
}

dev.off()

min(min_mindist_obs)    # minimum elk-cattle distance ever observed was 44 meters. 



plot(GVars$vector_mask)
plot(nullgrp[[thisimg]],add=T)
plot(elkPoints[[thisimg]],add=T,col="red")
plot(cowPoints[[thisimg]],add=T,col="green")

#grp <- elkPoints[[1]][elkPoints[[1]]@data$Group==1,]
#ovl <- GVars$raster_mask*








rm(grp,temp,prop,i)

######
# Run the ppm

?kppm



######
# new idea- sample (with clustering) from the kernel density distribution, use custom sampling algorithm... 


?rMatClust

elkppm <- kppm(elkppp[[2]],trend=~1,clusters="MatClust")  #test

temp <- rMatClust(kappa = elkppm$par[1],scale = elkppm$par[2],mu=elkppm$mu, nsim=1,win=window)
plot(elkppp[[2]])
plot(temp,add=T,col="red")


rm(covnames,i)







#######
# establish random polygons

#######
# visualize shapefiles




#######
# Elk use of pasture vs cow abundance



    ### Note this assumes constant elk population across all images


elkuse <- sapply(elkPoints,length)

cowuse <- sapply(cowPoints,length)

mod1 <- lm(elkuse~cowuse)

mod2 <- lm(elkuse~cowuse+I(cowuse^2))

summary(mod1)
summary(mod2)    # non-significant

newdata <- data.frame(
  cowuse=seq(100,1500,length=100)
)

pred <- predict(mod1,newdata=newdata,interval="confidence",level=0.95)

svg("elkuse_vs_cattle.svg",4,4)
plot(elkuse~cowuse, pch=19,cex=1.2,xlab="# Cattle in pasture",ylab="# Elk in pasture" )

polygon(x=c(newdata$cowuse,rev(newdata$cowuse)),y=c(pred[,2],rev(pred[,3])),col=gray(0.7),lty=0)

points(cowuse, elkuse,pch=19,cex=1.2)
points(newdata$cowuse,pred[,1],type="l",lwd=2,lty=2)
dev.off()



##############
# look at only cows within the elk home range area?


a <- contourLines(elkkernel,levels=elkcontours["0.95"])  # "0.95"   # do again for 50... 
a2 <- lapply(a,function(t) as.matrix(as.data.frame(t[c(2,3)])))
names(a2) <- c("p1","p2","p3")
hr1 <- lapply(a2,function(t) list(Polygon(t,hole=F))  )      #Polygon(a3,hole=F)
hr2 <- lapply(names(a2),function(t) Polygons(hr1[[t]],t) )
hr3 <- SpatialPolygons(hr2,proj4string = GVars$projection)

plot(hr3)


plot(cowPoints[[1]])
plot(hr3,add=T)

sapply(cowPoints,length)

CowsinHR95 <- sapply(cowPoints,function(t)  length(which(!is.na(t%over%hr3) )) )
ElkinHR95 <- sapply(elkPoints,function(t)  length(which(!is.na(t%over%hr3) )) )

mod1 <- lm(ElkinHR95~CowsinHR95)
mod2 <- lm(ElkinHR95~CowsinHR95+I(CowsinHR95^2))
summary(mod1)
summary(mod2)    # non-significant

newdata95 <- data.frame(
  CowsinHR95=seq(0,300,length=100)
)

pred95 <- predict(mod1,newdata=newdata95,interval="confidence",level=0.95)

svg("elkuse_vs_cattle95.svg",4,4)
plot(ElkinHR95~CowsinHR95, pch=19,cex=1.2,xlab="# Cattle in Elk HR",ylab="# Elk in HR" )

polygon(x=c(newdata95$CowsinHR95,rev(newdata95$CowsinHR95)),y=c(pred95[,2],rev(pred95[,3])),col=gray(0.7),lty=0)

points(CowsinHR95, ElkinHR95,pch=19,cex=1.2)
points(newdata95$CowsinHR95,pred[,1],type="l",lwd=2,lty=2)
dev.off()

a <- contourLines(elkkernel,levels=elkcontours["0.5"])  # "0.95"   # do again for 50... 
a2 <- lapply(a,function(t) as.matrix(as.data.frame(t[c(2,3)])))
names(a2) <- c("p1","p2","p3")
hr1 <- lapply(a2,function(t) list(Polygon(t,hole=F))  )      #Polygon(a3,hole=F)
hr2 <- lapply(names(a2),function(t) Polygons(hr1[[t]],t) )
hr3 <- SpatialPolygons(hr2,proj4string = GVars$projection)

plot(hr3)


plot(cowPoints[[1]])
plot(hr3,add=T)

sapply(cowPoints,length)

CowsinHR50 <- sapply(cowPoints,function(t)  length(which(!is.na(t%over%hr3) )) )
ElkinHR50 <- sapply(elkPoints,function(t)  length(which(!is.na(t%over%hr3) )) )

mod1 <- lm(ElkinHR50~CowsinHR50)
mod2 <- lm(ElkinHR50~CowsinHR50+I(CowsinHR50^2))
summary(mod1)
summary(mod2)    # non-significant

newdata <- data.frame(
  CowsinHR50=seq(0,300,length=100)
)

pred <- predict(mod1,newdata=newdata,interval="confidence",level=0.95)

svg("elkuse_vs_cattle50.svg",4,4)
plot(ElkinHR50~CowsinHR50, pch=19,cex=1.2,xlab="# Cattle in Elk HR",ylab="# Elk in HR" )

polygon(x=c(newdata$CowsinHR50,rev(newdata$CowsinHR50)),y=c(pred[,2],rev(pred[,3])),col=gray(0.7),lty=0)

points(CowsinHR50, ElkinHR50,pch=19,cex=1.2)
points(newdata$CowsinHR50,pred[,1],type="l",lwd=2,lty=2)
dev.off()


svg("elkuse_vs_cattle95_50.svg",6,4)
  layout(matrix(c(1:2),nrow=1))
  plot(ElkinHR95~CowsinHR95, pch=19,cex=1.2,xlab="# Cattle in Elk HR",ylab="# Elk in HR",main="Elk Home Range" )
  polygon(x=c(newdata95$CowsinHR95,rev(newdata95$CowsinHR95)),y=c(pred95[,2],rev(pred95[,3])),col=gray(0.7),lty=0)
  points(CowsinHR95, ElkinHR95,pch=19,cex=1.2)
  points(newdata95$CowsinHR95,pred95[,1],type="l",lwd=2,lty=2)
  
  plot(ElkinHR50~CowsinHR50, pch=19,cex=1.2,xlab="# Cattle in Elk CA",ylab="# Elk in CA",main="Elk Core Area" )
  polygon(x=c(newdata$CowsinHR50,rev(newdata$CowsinHR50)),y=c(pred[,2],rev(pred[,3])),col=gray(0.7),lty=0)
  points(CowsinHR50, ElkinHR50,pch=19,cex=1.2)
  points(newdata$CowsinHR50,pred[,1],type="l",lwd=2,lty=2)
dev.off()


plot(cowPoints[[4]])
plot(elkPoints[[4]],add=T,col="red")
plot(elkGroups[[4]],add=T,col="red")



##########
# VISUALIZE ELK CORE AREAS VS CATTLE LOCATIONS

a <- contourLines(elkkernel,levels=elkcontours["0.5"])  # "0.95"   # do again for 50... 
a2 <- lapply(a,function(t) as.matrix(as.data.frame(t[c(2,3)])))
names(a2) <- c("p1","p2","p3")
hr1 <- lapply(a2,function(t) list(Polygon(t,hole=F))  )      #Polygon(a3,hole=F)
hr2 <- lapply(names(a2),function(t) Polygons(hr1[[t]],t) )
hr3 <- SpatialPolygons(hr2,proj4string = GVars$projection)

#plot(hr3)

for(i in 1:length(cowPoints)){
  plot(cowPoints[[i]],main=names(cowPoints)[i])
  plot(hr3,add=T)
}


forfig <- c("2017Apr22b","2016Mar25","2015April25")
names(forfig) <- c("High","Med","Low")
svg("elkuse_vs_cattle_lowmedhigh.svg",6.5,3.5)
  layout(matrix(c(1:3),nrow=1))
  par(mai=c(0.1,0.1,0.7,0.1))
  for(f in 1:3){
    plot(cowPoints[[forfig[f]]],main=names(forfig)[f])
    plot(hr3,add=T,col="lightblue")
    plot(cowPoints[[forfig[f]]],add=T)
  }
dev.off()


names(cowPoints)



###########
# Visualize elk and cattle locations with ranch boundaries in background

graphics.off()

#vectorCovs$ranch_boundary
i=1
for(i in 1:length(cowPoints)){
  plot(vectorCovs$ranch_boundary,main=names(cowPoints)[i])
  plot(cowPoints[[i]],add=T,pch=20,col="red",cex=0.4)
  plot(elkPoints[[i]],add=T,pch=20,col="blue",cex=0.4)
  
}


##########
# TODO: test for overlap between elk MCPs and cattle MCPs







##########
# Make csv file for Kelley- for blossom software (MRPP analysis)

cowPoints$`2015Apr25`@coords

unique(substr(names(unlist(sapply(cowPoints,function(t) t@coords[,1]))),1,9))



df_cow <- data.frame(
  x = as.numeric(unlist(sapply(cowPoints,function(t) t@coords[,1]))),
  y = as.numeric(unlist(sapply(cowPoints,function(t) t@coords[,2]))),
  species = "COW",
  img = rep(names(cowPoints),times=sapply(cowPoints,length))
)

df_elk <- data.frame(
  x = as.numeric(unlist(sapply(elkPoints,function(t) t@coords[,1]))),
  y = as.numeric(unlist(sapply(elkPoints,function(t) t@coords[,2]))),
  species = "ELK",
  img = rep(names(elkPoints),times=sapply(elkPoints,length))
)

df <- rbind(df_cow,df_elk)

write.csv(df,file="allpoints.csv",row.names = FALSE)









##############
# Shapefiles for cow and elk contours
##############

a <- lapply(names(elkcontours), function(t) contourLines(elkkernel,levels=elkcontours[t]))
names(a) <- names(elkcontours)
a2 <- lapply(a,function(v) lapply(v,function(t) as.matrix(as.data.frame(t[c(2,3)]))))
temp <- lapply(1:length(a2),function(t) names(a2[[t]]) <<- paste("p",1:length(a2[[t]]),sep="") )
a2

hr1 <- lapply(a2,function(v) lapply(v,function(t) list(Polygon(t,hole=F))  ))      #Polygon(a3,hole=F)
hr2 <- lapply(hr1,function(v) lapply(names(v),function(t) Polygons(v[[t]],t) ))
hr3 <- lapply(hr2,function(v) SpatialPolygons(v,proj4string = GVars$projection))

pids <- lapply(hr3,function(v) sapply(slot(v, "polygons"), function(x) slot(x, "ID")) )
pdfs <- lapply(names(hr3), function(t) data.frame( ID=1:length(hr3[[t]]), row.names = pids[[t]]) ) 
names(pdfs) <- names(hr3)
hr4 <- lapply(names(hr3), function(t) SpatialPolygonsDataFrame(hr3[[t]], pdfs[[t]]))
names(hr4) <- names(hr3)

plot(hr4[[4]])
plot(hr4[[3]],add=T)
plot(hr4[[2]],add=T)
plot(hr4[[1]],add=T)

#?writeOGR
#setwd(shapefiledir)
temp <- lapply(names(hr4),function(t) 
  writeOGR(obj=hr4[[t]],dsn=shapefiledir,layer=paste("elkcontour",t,sep=""), driver="ESRI Shapefile")
) 

### and cow contours

a <- lapply(names(elkcontours), function(t) contourLines(elkkernel,levels=elkcontours[t]))
names(a) <- names(elkcontours)
a2 <- lapply(a,function(v) lapply(v,function(t) as.matrix(as.data.frame(t[c(2,3)]))))
temp <- lapply(1:length(a2),function(t) names(a2[[t]]) <<- paste("p",1:length(a2[[t]]),sep="") )
a2

hr1 <- lapply(a2,function(v) lapply(v,function(t) list(Polygon(t,hole=F))  ))      #Polygon(a3,hole=F)
hr2 <- lapply(hr1,function(v) lapply(names(v),function(t) Polygons(v[[t]],t) ))
hr3 <- lapply(hr2,function(v) SpatialPolygons(v,proj4string = GVars$projection))

pids <- lapply(hr3,function(v) sapply(slot(v, "polygons"), function(x) slot(x, "ID")) )
pdfs <- lapply(names(hr3), function(t) data.frame( ID=1:length(hr3[[t]]), row.names = pids[[t]]) ) 
names(pdfs) <- names(hr3)
hr4 <- lapply(names(hr3), function(t) SpatialPolygonsDataFrame(hr3[[t]], pdfs[[t]]))
names(hr4) <- names(hr3)

plot(hr4[[4]])
plot(hr4[[3]],add=T)
plot(hr4[[2]],add=T)
plot(hr4[[1]],add=T)

#?writeOGR
#setwd(shapefiledir)
temp <- lapply(names(hr4),function(t) 
  writeOGR(obj=hr4[[t]],dsn=shapefiledir,layer=paste("elkcontour",t,sep=""), driver="ESRI Shapefile")
) 


#  




