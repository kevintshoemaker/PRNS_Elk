###Script for estimating sensitivity of RSF coefficients to the size of "available" sample points. Created by Jared Stabach (stabachj@si.edu). Please do not distribute without permission.

#rm(list=ls())

#Set WD
#setwd("/Users/Lacey/Box Sync/PhD_Project/R.Analyses/Elk/Raw data")


###For collar data, after subsetting by individual ID, I went into this script to test sensitivity of various sample sizes. ###


#load dataset and re-name as deer1 to run through script below.
deer1 <- m1a

Total.Avail.Samples <- c(0.25, 0.5, 1, 2)# I started with 10X number of used points, so working backwards to see how taking smaller samples would affect sensitivity

NumberCoefficients <- 9 #change this to be 1 + number of coefficients being run

#Subset "all" into used and available
deer.use <- subset(deer1, Used == 1)#swap males 2 for "deer1" in original
deer.use <- deer.use[c(19:27)]#subset only for coefficients that will be included plus "Used" column

deer.avail <- subset(deer1, Used == 0)#swap all.std for "deer1" in original
deer.avail <- deer.avail[c(19:27)]#subset only for coefficients that will be included plus "Used" column

Sensitivity.RSF <- function(USE,AVAIL,NumberSamples,NumberCoefficients){
  
  # Number of Iterations
  n.iter=100
  
  # Create matrix to hold everything 
  # *** Edit number of columns to create…..based on the number of covariates (below) ***
  beta.1 <- matrix(NA, nrow = n.iter, ncol = NumberCoefficients)
  
  # Use simulation to test the sensitivity of the number of samples.
  for(i in 1:n.iter){
    # Draw a sample from the random points file to test the sensitivity
    s.index <- sample(nrow(AVAIL),nrow(USE) * NumberSamples,replace=TRUE)# NOTE: I replaced "FALSE" with "TRUE" here to make it run- need to understand what's happening to do this before trusting these results. Also need to confirm final glm.
    x.a <- AVAIL[s.index,]
    
    # Bind together
    x.glm <-rbind(USE, x.a)
    
    # Run GLM # *******Edit covariates ********
    out <- glm (Used ~ std.pasture + std.stocking + std.ponds + std.elevation + 
                  std.insolation + std.scrub + std.slope +std.ndvi,
                data= x.glm, family="binomial" (link= "logit"))
    
    
    # Output coefficients # *** Remove beta.t values as appropriate….should match columns in matrix
    beta.1[i,1] <-out$coeff[1]
    beta.1[i,2] <-out$coeff[2]
    beta.1[i,3] <-out$coeff[3]
    beta.1[i,4] <-out$coeff[4]
    beta.1[i,5] <-out$coeff[5]
    beta.1[i,6] <-out$coeff[6]
    beta.1[i,7] <-out$coeff[7]
    beta.1[i,8] <-out$coeff[8]
    #beta.1[i,9] <-out$coeff[9]
    #beta.1[i,10] <-out$coeff[10]
    #beta.1[i,11] <-out$coeff[11]
    #beta.1[i,12] <-out$coeff[12]
    #beta.1[i,13] <-out$coeff[13]
    #beta.1[i,14] <-out$coeff[14]
    #beta.1[i,15] <-out$coeff[15]
    #beta.1[i,16] <-out$coeff[16]
  }
  # Return the results
  return(beta.1)
}

#======ANALYSIS=================
Sens.List <- vector("list")

# Run loop to create samples availability and model for use in Sensitivity analysis
# Output results to a list (Sens.List)
for (i in 1:length(Total.Avail.Samples)){
  Sens.1Sample <- Sensitivity.RSF(deer.use, deer.avail, Total.Avail.Samples[i], NumberCoefficients)
  Sens.List[[i]] <- Sens.1Sample
}

#=============VISUALIZATION=====================
# Set up blank matrices to store result summaries
beta.1.up=matrix(0, length(Total.Avail.Samples),NumberCoefficients)
beta.1.low=matrix(0, length(Total.Avail.Samples),NumberCoefficients)
beta.1.mean=matrix(0, length(Total.Avail.Samples),NumberCoefficients)

# Variables being monitored
Variables <- c("Beta", "Pasture", "Stocking rate", "Ponds", "Elevation", "Insolation", "NDVI", "Scrub", "Slope")

for (i in 1:length(Total.Avail.Samples)){
  beta.1.up[i,] <-apply(Sens.List[[i]],2,quantile,prob=0.975, na.rm=TRUE)
  beta.1.low[i,] <-apply(Sens.List[[i]],2,quantile,prob=0.025, na.rm=TRUE)
  beta.1.mean[i,] <-apply(Sens.List[[i]],2,mean,na.rm=TRUE)
}

# Transpose to plot
beta.1.up <-t(beta.1.up)
beta.1.low <- t(beta.1.low)
beta.1.mean <- t(beta.1.mean)

#remove NAs- if NA's are present in "available samples"- which they should be due to the dune mask for NDVI.
beta.1.up <- na.omit(beta.1.up)
beta.1.low <- na.omit(beta.1.low)
beta.1.mean <- na.omit(beta.1.mean)

# Plot the results
for (i in 1:nrow(beta.1.up)){
  plot(Total.Avail.Samples, beta.1.up[i,], typ="l", xaxt="n", ylim=c(min(beta.1.low[i,]), max(beta.1.up[i,])),lty=2, main=paste0(Variables[i]," Sensitivity"), ylab="Coefficient",xlab="Availability Sample Size")
  axis(side=1, at=Total.Avail.Samples, label=c("0.25x", "0.5x", "1x", "2x"))#XLab
  lines(Total.Avail.Samples, beta.1.low[i,], lty=2)
  lines(Total.Avail.Samples, beta.1.mean[i,])
}

