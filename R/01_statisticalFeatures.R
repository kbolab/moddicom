#######################################################################################################################
################################ STATISTICAL FEATURES #################################################################
#######################################################################################################################

statisticalFeatures <- function(imgObj){ 

  nVoxel <- dim(imgObj)[1]*dim(imgObj)[2]*dim(imgObj)[3] - sum(is.na(imgObj))
  n <- numeric()
  # for (i in seq(1:length(table(imgObj)))){
  #   n[i] <- table(imgObj)[i]
  # }
  n <- table(imgObj)
  p <- n/nVoxel

  #Initialise data table for storing Statistical features
  featNames <- c("F_stat.mean","F_stat.var","F_stat.skew","F_stat.kurt","F_stat.median","F_stat.min","F_stat.10thpercentile",
                 "F_stat.90thpercentile","F_stat.max","F_stat.iqr","F_stat.range","F_stat.mad",
                 "F_stat.rmad","F_stat.energy","F_stat.rms","F_stat.entropy","F_stat.uniformity",
                 "F_stat.Nic.entropy", "F_stat.Nic.kurt", "F_stat.Nic.Skew")
  F_stat <- data.frame(matrix(NA, ncol=length(featNames)))
  colnames(F_stat) <- featNames

  F_stat$F_stat.mean <- sum(imgObj, na.rm = T)/nVoxel
  F_stat$F_stat.var <- sum((imgObj -   F_stat$F_stat.mean)^2, na.rm = T) / nVoxel
  F_stat$F_stat.skew <- (sum((imgObj -   F_stat$F_stat.mean)^3, na.rm = T) / nVoxel) / (sum((imgObj -   F_stat$F_stat.mean)^2, na.rm = T) / nVoxel)^(3/2)
  F_stat$F_stat.kurt <- ((sum((imgObj -   F_stat$F_stat.mean)^4, na.rm = T) / nVoxel) / (sum((imgObj -   F_stat$F_stat.mean)^2, na.rm = T) / nVoxel)^2 ) - 3
  F_stat$F_stat.median <- median(imgObj, na.rm = T)
  F_stat$F_stat.min <- min(imgObj, na.rm = T)
  F_stat$F_stat.10thpercentile <- as.numeric(quantile(imgObj,0.10,na.rm = T))
  F_stat$F_stat.90thpercentile <- as.numeric(quantile(imgObj,0.90,na.rm = T))
  F_stat$F_stat.max <- max(imgObj, na.rm = T)
  F_stat$F_stat.iqr <- as.numeric(quantile(imgObj,0.75,na.rm = T)) - as.numeric(quantile(imgObj,0.25,na.rm = T))
  F_stat$F_stat.range <-   F_stat$F_stat.max -   F_stat$F_stat.min
  F_stat$F_stat.mad <- sum(abs(imgObj -   F_stat$F_stat.mean), na.rm=T) / nVoxel

#### STUFF for RMAD (ROBUST MEAN ABSOLUTE DEVIATION)
nVoxel_P10_P90 <- length(which(as.numeric(imgObj) >=   F_stat$F_stat.10thpercentile & as.numeric(imgObj) <=   F_stat$F_stat.90thpercentile))
rmadVoxels <- as.numeric(imgObj)[which(as.numeric(imgObj) >=   F_stat$F_stat.10thpercentile & as.numeric(imgObj) <=   F_stat$F_stat.90thpercentile)]
rmadVoxels_mean <- mean(rmadVoxels)
####

F_stat$F_stat.rmad <- sum(abs(rmadVoxels - rmadVoxels_mean)) / nVoxel_P10_P90
F_stat$F_stat.energy <- sum(imgObj^2, na.rm=T)
F_stat$F_stat.rms <- sqrt(sum(imgObj^2, na.rm=T)/nVoxel)
F_stat$F_stat.entropy <- -sum(p * log2(p))
F_stat$F_stat.uniformity <- sum(p^2)

F_stat$F_stat.Nic.entropy <- entropy(y = discretize(x = imgObj[which(!is.na(imgObj))], numBins = 25))
F_stat$F_stat.Nic.skew <- skewness(imgObj, na.rm = T)
F_stat$F_stat.Nic.kurt <- kurtosis( imgObj ,na.rm = T)

return(F_stat)
}


statisticalFeatures.int <- function(imgObj){

  Nv <- length(imgObj[which(!is.na(imgObj))])
  p_i<- numeric()
  for(j in 1:length(table(imgObj))){
    p_i[j] <- table(imgObj)[j]/Nv
  }

  F_stat.int$F_stat.entropy <- - sum(p_i * log2(p_i))
  F_stat.int$F_stat.uniformity <- sum(p_i^2)

  return(F_stat.int)
}
