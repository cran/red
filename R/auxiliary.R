###############################################################################
##############################AUX FUNCTIONS####################################
###############################################################################

longlat2utm <- function(longlat){
  ### Added exception if point data.
  if (is(longlat, "matrix") || is(longlat, "data.frame")){
    if(is(longlat, "data.frame")){
      longlat = as.matrix(longlat)
    }
    minlong = min(longlat[,1])
    zone = floor((minlong + 180) / 6) + 1
    res = terra::project(longlat,
                         from = "+proj=longlat +datum=WGS84",
                         to = paste0("+proj=utm +zone=", zone," ellps=WGS84"))
    return(res)
  }
  
  if (is(longlat, "SpatRaster")){
    minlong = terra::xmin(longlat)
    zone = floor((minlong + 180) / 6) + 1
    res <- terra::project(longlat, paste0("+proj=utm +zone=", zone," ellps=WGS84"))
    res <- data.frame(x = terra::crds(res)[,1], y = terra::crds(res)[,2])
    return(res)
  }
  
  longlat = as.matrix(longlat)
  minlong = min(longlat[,1])
  zone = floor((minlong + 180) / 6) + 1
  res = terra::project(longlat, paste0("+proj=utm +zone=", zone," ellps=WGS84"))
  
  res <- data.frame(x = terra::crds(res)[,1],
                    y = terra::crds(res)[,2])
  return(res)
}

utm2longlat <- function(utm, zone){
  if(is(utm, "SpatRaster")){
    if(!is.null(zone))
      terra::crs(utm) <- paste("+proj=utm +zone=", zone, sep="") ### raster::crs()
      
    res <- terra::project(utm, y = "+proj=longlat +datum=WGS84")
  } else {
    utm <- SpatialPoints(utm, CRS(paste("+proj=utm +zone=", zone,sep="")))
    res <- as.data.frame(spTransform(utm,CRS(paste("+proj=longlat"))))
  }
  return(res)
}

##detect which layers are categorical by checking if all values are integers and if the max is less than 50 (may fail, just an attempt)
find.categorical <- function(layers){
  categorical = c()
  for(l in 1:(dim(layers)[3])){
    lay <- terra::as.matrix(layers[[l]])
    lay <- as.vector(lay)
    lay <- lay[!is.na(lay)]
    if(sum(floor(lay)) == sum(lay) && length(unique(lay)) < 50)
      categorical = c(categorical, l)
  }
  return(categorical)
}

##basic function to calculate the rli of any group of species
rli.calc <- function(spData, tree = NULL, boot = FALSE, dd = FALSE, runs = 1000){
  if(all(is.na(spData)))
    return(NA)
  spData <- rli.convert(spData)                ##call function to convert spData to a 0-1 scale
  
  if(is.null(tree)){                           ##if not weighted by PD or FD
    if(!boot){                                 ##if no bootstrap to be made
      return (mean(spData, na.rm = TRUE))
    } else {
      run <- rep(NA, runs)
      if(!dd){
        for(i in 1:runs){
          rnd <- sample(spData, replace = TRUE) ##bootstrap with all species
          run[i] <- mean(rnd, na.rm = TRUE)
        }
      } else {                                       ##bootstrap with only DD species
        nDD = sum(is.na(spData))                     ##number of DD species
        rliBase = sum(spData, na.rm = TRUE)
        for(i in 1:runs){
          rnd <- sample(spData[!is.na(spData)], nDD, replace = TRUE)
          run[i] <- (rliBase + sum(rnd)) / length(spData)
        }
      }
      res <- matrix(quantile(run, c(0.025, 0.5, 0.975)), nrow = 1)
      colnames(res) <- c("LowCL", "Median", "UpCL")
      return(res)
    }
  } else {                                     ##if weighted by PD or FD, still to work, not available at the moment!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    comm <- matrix(1, nrow = 2, ncol = length(spData))
    contrib <- BAT::contribution(comm, tree, relative = TRUE)[1,]
    contrib <- contrib/sum(contrib[!is.na(spData)]) #needed to standardize the contribution by the total contribution of species living in the community
    if(!boot){                                 ##if no bootstrap to be made
      return(sum(spData * contrib, na.rm = TRUE))
    } else {
      run <- rep(NA, runs)
      for(i in 1:runs){
        rndSpp <- sample(length(spData), replace = TRUE)
        rndComm <- spData[rndSpp]
        rndContrib <- contrib[rndSpp]/sum(contrib[rndSpp])
        run[i] <- sum(rndComm * rndContrib, na.rm = TRUE)
      }
      res <- matrix(quantile(run, c(0.025, 0.5, 0.975)), nrow = 1)
      colnames(res) <- c("LowCL", "Median", "UpCL")
      return(res)
    }
  }
}

##function to convert strings to numbers in the RLI
rli.convert <- function(spData){
  if(!is.numeric(spData)){                                ##if letters are given, convert to [0,1]
    spData <- replace(spData, which(spData == "EX" ), 0)
    spData <- replace(spData, which(spData == "EW" ), 0)
    spData <- replace(spData, which(spData == "RE" ), 0)
    spData <- replace(spData, which(spData == "CR" ), 0.2)
    spData <- replace(spData, which(spData == "CR(PE)" ), 0.2)
    spData <- replace(spData, which(spData == "EN" ), 0.4)
    spData <- replace(spData, which(spData == "VU" ), 0.6)
    spData <- replace(spData, which(spData == "NT" ), 0.8)
    spData <- replace(spData, which(spData == "LC" ), 1)
    spData <- replace(spData, which(spData == "DD" ), NA)
    spData <- as.numeric(spData)
  } else if (all(spData == floor(spData))){  #if all integers, a scale [0,5] is given, convert to [0,1]
    spData <- 1 - spData/5
  }
  return(spData)
}

#required for kml
spline.poly <- function(xy, vertices, k=3, ...) {
  # Assert: xy is an n by 2 matrix with n >= k.
  
  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  
  # Spline the x and y coordinates.
  data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y
  
  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}
