
################################################################################
# Function for calculating custom new data for raw or scaled predictors
# Heather Kenny-Duddela
# March 16, 2025
################################################################################

custom.newdata <- function(old.data, data.type, scaled, xvar, xunscale) {
  
  # Pull columns from female table or male table
  if (data.type=="female") {
    list <- colnames(old.data)[c(8:11,21:24,30:31,36,40:50)]  
  } else {
    list <- colnames(old.data)[c(8:11,17:20,31:32,36,40:50)]
  }
  
  # vector for range of xvar
  seq.xvar <- seq(range(old.data[ ,which(colnames(old.data)==xvar)], 
                        na.rm=T)[1], 
                      range(old.data[ ,which(colnames(old.data)==xvar)], 
                            na.rm=T)[2], 0.001)
  xvar.n <- length(seq.xvar)  
  
  ## Different variable names for females and males
  if (data.type=="female") {
    
    # Different new data for scaled or unscaled variables
    if (scaled==T) {
      
    ## For scaled variables
      # generate new data
      newdata <- 
        data.frame(
          tail_scaled = rep(mean(old.data$tail_scaled, na.rm=T), xvar.n),
          site.size = rep("medium", xvar.n),
          socM_tail_scaled = rep(mean(old.data$socM_tail_scaled, na.rm=T), xvar.n),
          ci_1_julian_scaled = rep(mean(old.data$ci_1_julian_scaled, na.rm=T), xvar.n),
          throat_avg_bright_scaled = rep(mean(old.data$throat_avg_bright_scaled, na.rm=T), xvar.n),
          socM_t.avg.bright_scaled = rep(mean(old.data$socM_t.avg.bright_scaled, na.rm=T), xvar.n),
          breast_avg_bright_scaled = rep(mean(old.data$breast_avg_bright_scaled, na.rm=T), xvar.n),
          socM_r.avg.bright_scaled = rep(mean(old.data$socM_r.avg.bright_scaled, na.rm=T), xvar.n))
      
      # add sequence for xvar
      newdata[, which(colnames(newdata)==xvar)] <- seq.xvar
    
      # un-scale the predictor
      newdata$unscale <- (newdata[, which(colnames(newdata)==xvar)] * 
                            sd(old.data[ ,which(colnames(old.data)==xunscale)], 
                               na.rm=T) +
                            mean(old.data[ ,which(colnames(old.data)==xunscale)], 
                                 na.rm=T))
      # informative colname
      colnames(newdata)[9] <- xunscale
      
    } else {
      
    ## For unscaled variables
      
      # generate new data
      newdata <- data.frame(
        tail = rep(mean(old.data$tail, na.rm=T), xvar.n),
        socM_tail = rep(mean(old.data$socM_tail, na.rm=T), xvar.n),
        site.size = rep ("medium", xvar.n), 
        ci_1_julian = rep(mean(old.data$ci_1_julian, na.rm=T), xvar.n),
        throat_avg_bright = rep(mean(old.data$throat_avg_bright, na.rm=T), xvar.n),
        socM_t.avg.bright = rep(mean(old.data$socM_t.avg.bright, na.rm=T), xvar.n),
        breast_avg_bright = rep(mean(old.data$breast_avg_bright, na.rm=T), xvar.n),
        socM_r.avg.bright = rep(mean(old.data$socM_r.avg.bright, na.rm=T), xvar.n))
      
      # add sequence for xvar
      newdata[, which(colnames(newdata)==xvar)] <- seq.xvar
    }
  
  # For male data table    
  } else {
    
    # Different new data for scaled or unscaled variables
    if (scaled==T) {
      
      ## For scaled variables
      # generate new data
      newdata <- 
        data.frame(
          tail_scaled = rep(mean(old.data$tail_scaled, na.rm=T), xvar.n),
          site.size = rep("medium", xvar.n),
          socF_tail_scaled = rep(mean(old.data$socF_tail_scaled, na.rm=T), xvar.n),
          ci_julian_scaled = rep(mean(old.data$ci_julian_scaled, na.rm=T), xvar.n),
          throat_avg_bright_scaled = rep(mean(old.data$throat_avg_bright_scaled, na.rm=T), xvar.n),
          socF_t.avg.bright_scaled = rep(mean(old.data$socF_t.avg.bright_scaled, na.rm=T), xvar.n),
          breast_avg_bright_scaled = rep(mean(old.data$breast_avg_bright_scaled, na.rm=T), xvar.n),
          socF_r.avg.bright_scaled = rep(mean(old.data$socF_r.avg.bright_scaled, na.rm=T), xvar.n))
      
      # add sequence for xvar
      newdata[, which(colnames(newdata)==xvar)] <- seq.xvar
      
      # un-scale the predictor
      newdata$unscale <- (newdata[, which(colnames(newdata)==xvar)] * 
                            sd(old.data[ ,which(colnames(old.data)==xunscale)], na.rm=T) + 
                            mean(old.data[ ,which(colnames(old.data)==xunscale)], na.rm=T))
      # informative colname
      colnames(newdata)[9] <- xunscale

    } else {
      
      ## For unscaled variables
      
      # generate new data
      newdata <- data.frame(
        tail = rep(mean(old.data$tail, na.rm=T), xvar.n),
        socF_tail = rep(mean(old.data$socF_tail, na.rm=T), xvar.n),
        site.size = rep ("medium", xvar.n), 
        ci_julian = rep(mean(old.data$ci_julian, na.rm=T), xvar.n),
        throat_avg_bright = rep(mean(old.data$throat_avg_bright, na.rm=T), xvar.n),
        socF_t.avg.bright = rep(mean(old.data$socF_t.avg.bright, na.rm=T), xvar.n),
        breast_avg_bright = rep(mean(old.data$breast_avg_bright, na.rm=T), xvar.n),
        socF_r.avg.bright = rep(mean(old.data$socF_r.avg.bright, na.rm=T), xvar.n))
      
      # add sequence for xvar
      newdata[, which(colnames(newdata)==xvar)] <- seq.xvar
      
    }
  } 
  return(newdata)
} 
