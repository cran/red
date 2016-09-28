#####RED - IUCN Redlisting Tools
#####Version 0.1.1 (2016-09-28)
#####By Pedro Cardoso
#####Maintainer: pedro.cardoso@helsinki.fi
#####Reference: Cardoso, P.(in prep.) An R package to facilitate species red list assessments.
#####Changed from v0.1.0:
#####Notation of mcp in function map.draw
#####clarified a number of points in documentation

#####data origins:
#####climate -> Hijmans, R.J., Cameron, S.E, Parra, J.L., Jones, P.G. & Jarvis A. (2005) Very high resolution interpolated climate surfaces for global land areas. International Journal of Climatology, 25: 1965-1978.
#####altitude -> Farr, T. G., et al. (2007), The Shuttle Radar Topography Mission, Rev. Geophys., 45, RG2004
#####landcover -> Tuanmu, M.-N. & Jetz, W. (2014) A global 1-km consensus land-cover product for biodiversity and ecosystem modeling. Global Ecology and Biogeography, 23: 1031-1045.

#####RED Stats:
#####library("cranlogs")
#####day <- cran_downloads(package = "red", from = "2016-08-19", to = "2016-08-019")
#####group <- matrix(day$count, 140, byrow = TRUE)
#####plot(rowSums(group), type = "n")
#####lines(rowSums(group))

#####required packages
library("BAT")
library("dismo")
library("geosphere")
library("graphics")
library("grDevices")
library("igraph")
library("maptools")
library("raster")
library("rgdal")
library("rgeos")
library("rJava")
library("sp")
library("stats")
library("utils")
#' @importFrom geosphere areaPolygon
#' @import graphics
#' @importFrom grDevices chull dev.copy dev.off pdf
#' @import maptools
#' @importFrom raster area cellStats clump crop extent extract getValues layerStats mask raster rasterize rasterToPoints reclassify sampleRandom scalebar terrain trim writeRaster
#' @import rgdal
#' @import rgeos
#' @import sp
#' @import stats
#' @import utils

###############################################################################
##############################AUX FUNCTIONS####################################
###############################################################################

##detect which layers are categorical by checking if all values are integers and if the max is less than 100 (may fail, just an attempt)
find.categorical <- function(layers){
  categorical = c()
  for(l in 1:(dim(layers)[3])){
    lay <- raster::as.matrix(layers[[l]])
    lay[is.na(lay)] <- 0
    if(sum(as.integer(lay)) == sum(lay) && max(lay < 100))
      categorical = c(categorical, l)
  }
  return(categorical)
}

##basic function to calculate the rli of any group of species
rli.calc <- function(spData, tree = NULL, boot = FALSE, runs = 1000){
  if(max(spData) > 5){                                ##if letters are given, convert to [0,1]
    spData <- replace(spData, which(spData == "EX" ), 0)
    spData <- replace(spData, which(spData == "EW" ), 0)
    spData <- replace(spData, which(spData == "RE" ), 0)
    spData <- replace(spData, which(spData == "CR" ), 0.2)
    spData <- replace(spData, which(spData == "EN" ), 0.4)
    spData <- replace(spData, which(spData == "VU" ), 0.6)
    spData <- replace(spData, which(spData == "NT" ), 0.8)
    spData <- replace(spData, which(spData == "LC" ), 1)
    spData <- replace(spData, which(spData == "DD" ), 2)
    spData <- as.numeric(spData)
    spData <- subset(spData, spData < 2)
  } else if (max(spData) > 1){                       ##if a scale [0,5] is given, convert to [0,1]
    spData <- 1 - spData / 5
  }
  if(is.null(tree)){                           ##if not weighted by PD or FD
    if(!boot){                                 ##if no bootstrap to be made
      return (mean(spData))
    } else {
      run <- rep(NA, runs)
      for(i in 1:runs){
        rnd <- sample(spData, length(spData), replace = TRUE) ##bootstrap
        run[i] <- mean(rnd)
      }
      res <- matrix(quantile(run, c(0.025, 0.5, 0.975)), nrow = 1)
      colnames(res) <- c("LowCL", "Median", "UpCL")
      return(res)
    }
  } else {                                     ##if weighted by PD or FD
    comm <- matrix(1, nrow = 2, ncol = length(spData))
    contrib <- BAT::contribution(comm, tree, relative = TRUE)[1,]
    if(!boot){                                 ##if no bootstrap to be made
      return (sum(contrib * spData))
    } else {
      for(i in 1:runs){
        rnd <- sample(spData, length(spData), replace = TRUE) ##bootstrap
        for(i in 1:length(rnd))
          run[i] <- run[i] + contrib[rnd[i]] * spData[rnd[i]]
      }
      return(quantile(run, c(0.025, 0.5, 0.975)))
    }
  }
}

##################################################################################
##################################MAIN FUNCTIONS##################################
##################################################################################

#' Red package setup.
#' @description Setup red to work with species distribution modelling and layers available online.
#' @details Please check that you have at least 20Gb free in your disk (and a fast internet connection) to download all files. In the end of the process "only" 3Gb will be left though. This function will:
#' 1. Check if maxent.jar is available in the dismo package directory.
#' 2. Create a folder named "gis" in the working directory.
#' 3. Download global bioclim and elevation files (20) from http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/.
#' 4. Download landcover files (12) from http://data.earthenv.org/consensus_landcover/without_DISCover/.
#' 5. Unzip all files and delete the originals.
#' 6. Create a new layer (1) with the dominant land cover at each cell.
#' 7. Resample all files (33) to 2x2km (as required by IUCN) and 10x10km (for use with widespread species) grid cells.
#' 8. Delete the 1x1km files.
#' Sit back and enjoy, this should take a while.
#' @export
red.setup <- function(){

  ##test if maxent.jar is in the right directory
  if(!file.exists(paste(.libPaths()[[1]], "/dismo/java/maxent.jar", sep=""))){
    return(warning("RED could not find maxent.jar.
1. Download the latest version of maxent from:
https://www.cs.princeton.edu/~schapire/maxent/
2. Move the file maxent.jar to the java directory inside dismo package
(there should be a file named dismo.jar already there)
3. Install the latest version of java runtime environment (JRE) with the same architecture (32 or 64 bits) as your version of R:
http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html"))
  }

  ##basic setup
  dir.create("gis")
  pb <- txtProgressBar(min = 0, max = 33, style = 3)

  ##download and process bioclim 1-9
  download.file("http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio1-9_30s_bil.zip", "gis/bio1-9_30s_bil.zip")
  unzip(zipfile = "gis/bio1-9_30s_bil.zip", exdir = "gis")
  file.remove("gis/bio1-9_30s_bil.zip")
  for(i in 1:9){
    setTxtProgressBar(pb, i)
    rast <- raster(paste("gis/bio_", i, ".bil", sep=""))
    rast <- crop(rast, c(-180, 180, -56, 90))
    rast2 <- aggregate(rast, 2)
    writeRaster(rast2, paste("gis/red_2km_", i, ".tif", sep=""))
    rast10 <- aggregate(rast, 10)
    writeRaster(rast10, paste("gis/red_10km_", i, ".tif", sep=""))
    file.remove(paste("gis/bio_", i, ".bil", sep=""))
    file.remove(paste("gis/bio_", i, ".hdr", sep=""))
    gc()
  }

  ##download and process bioclim 10-19
  download.file("http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio10-19_30s_bil.zip", "gis/bio10-19_30s_bil.zip")
  unzip(zipfile = "gis/bio10-19_30s_bil.zip", exdir = "gis")
  file.remove("gis/bio10-19_30s_bil.zip")
  for(i in 10:19){
    setTxtProgressBar(pb, i)
    rast <- raster(paste("gis/bio_", i, ".bil", sep=""))
    rast <- crop(rast, c(-180, 180, -56, 90))
    rast2 <- aggregate(rast, 2)
    writeRaster(rast2, paste("gis/red_2km_", i, ".tif", sep=""))
    rast10 <- aggregate(rast, 10)
    writeRaster(rast10, paste("gis/red_10km_", i, ".tif", sep=""))
    file.remove(paste("gis/bio_", i, ".bil", sep=""))
    file.remove(paste("gis/bio_", i, ".hdr", sep=""))
    gc()
  }

  ##download and process altitude
  setTxtProgressBar(pb, 20)
  download.file("http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/alt_30s_bil.zip", "gis/alt_30s_bil.zip")
  unzip(zipfile = "gis/alt_30s_bil.zip", exdir = "gis")
  file.remove("gis/alt_30s_bil.zip")
  rast <- raster("gis/alt.bil")
  rast <- crop(rast, c(-180, 180, -56, 90))
  rast2 <- aggregate(rast, 2)
  writeRaster(rast2, "gis/red_2km_20.tif")
  rast10 <- aggregate(rast, 10)
  writeRaster(rast10, "gis/red_10km_20.tif")
  file.remove("gis/alt.bil")
  file.remove("gis/alt.hdr")
  gc()

  ##download and process land cover
  max2 <- raster::stack()
  max10 <- raster::stack()
  for(i in 1:12){
    setTxtProgressBar(pb, (i+20))
    download.file(paste("http://data.earthenv.org/consensus_landcover/without_DISCover/Consensus_reduced_class_", i, ".tif", sep=""), destfile = paste("gis/Consensus_reduced_class_", i, ".tif", sep=""), mode = "wb")
    rast <- raster(paste("gis/Consensus_reduced_class_", i, ".tif", sep=""))
    rast2 <- aggregate(rast, 2)
    writeRaster(rast2, paste("gis/red_2km_", (i+20), ".tif", sep=""))
    rast10 <- aggregate(rast, 10)
    writeRaster(rast10, paste("gis/red_10km_", (i+20), ".tif", sep=""))
    file.remove(paste("gis/Consensus_reduced_class_", i, ".tif", sep=""))
    ##write stack for later use of maxStack
    max2 <- raster::stack(max2, rast2)
    max10 <- raster::stack(max10, rast10)
    gc()
  }
  remove(rast, rast2, rast10)

  ##create new rasters with most common landcover at each cell
  setTxtProgressBar(pb, 33)
  max2 <- which.max(max2)
  writeRaster(max2, "gis/red_2km_33")
  max10 <- which.max(max10)
  writeRaster(max10, "gis/red_10km_33")
  remove(max2, max10)
  gc()

  ##Now the files should be named as:
  ##red_2km_1.tif
  ##...
  ##red_10km_33.tif
  ##Where 1 to 19 are the corresponding bioclim variables, 20 is altitude, 21 to 32 are landcover proportion and 33 is most common landcover per cell
}

#' Visual detection of outliers.
#' @description Draws plots of sites in geographical (latlong) and environmental (2-axis PCA) space.
#' @param longlat Matrix of longitude and latitude (two columns) of species occurrence records.
#' @param layers Raster* object as defined by package raster. It can be any set of environmental layers thought to allow the identification of environmental outliers.
#' @details Erroneous data sources or errors in transcriptions may introduce outliers that can be easily detected by looking at simple graphs of geographical or environmental space.
#' @return No values are returned, two plots are drawn for visual inspection.
# @examples data(data.records)
# data(data.layers)
# outliers(data.records, data.layers[[1:3]])
#' @export
outliers <- function(longlat, layers){
  pca <- raster.reduce(layers, n = 2)
  ##extract pca values from longlat
  pca <- extract(pca, longlat)[,1:2]
  par(mfrow = c(1,2))
  graphics::plot(longlat, main = "Geographical")
  raster::plot(pca, main = "Environmental")
}

#' Spatial thinning of occurrence records.
#' @description Thinning of records with minimum distances either absolute or relative to the species range.
#' @param longlat Matrix of longitude and latitude (two columns) of species occurrence records.
#' @param distance Distance either in relative terms (proportion of maximum distance between any two records) or in raster units.
#' @param relative If TRUE, represents the proportion of maximum distance between any two records. If FALSE, is in raster units.
#' @param runs Number of runs
#' @details Clumped distribution records due to ease of accessibility of sites, emphasis of sampling on certain areas in the past, etc. may bias species distribution models.
#' The algorithm used here eliminates records closer than a given thres to any other record. The choice of records to eliminate is random, so a number of runs are made and the one keeping more of the original records is chosen.
#' @return A matrix of longitude and latitude (two columns) of species occurrence records separated by at least the given distance.
#' @examples records <- matrix(sample(100), ncol = 2)
#' par(mfrow=c(1,2))
#' graphics::plot(records)
#' records <- thin(records, 0.1)
#' graphics::plot(records)
#' @export
thin <- function(longlat, distance = 0.01, relative = TRUE, runs = 100){
  longlat = longlat[!duplicated(longlat),]                #first, remove duplicate rows
  nSites = nrow(longlat)
  if(nSites < 4)
    return(longlat)

  ##if relative, calculate maxDist between any two points
  if(relative){
    maxDist = 0
    for(x in 1:(nSites-1)){
      for(y in (x+1):nSites){
        maxDist = max(maxDist,((longlat[x,1]-longlat[y,1])^2+(longlat[x,2]-longlat[y,2])^2)^.5)
      }
    }
    distance = maxDist*distance
  }

  listSites = matrix(longlat[1,], ncol=2, byrow = TRUE)
  for (r in 1:runs){
    longlat = longlat[sample(nSites),]       ##shuffle rows (sites)
    rndSites = longlat[1,]                   ##start with first random site
    for(newSite in 2:nSites){
      for(oldSite in 1:(newSite-1)){
        addSite = TRUE
        dist = ((longlat[newSite,1]-longlat[oldSite,1])^2+(longlat[newSite,2]-longlat[oldSite,2])^2)^.5
        if(dist < distance){
          addSite = FALSE
          break
        }
      }
      if(addSite)
        rndSites = rbind(rndSites, longlat[newSite,])
    }
    if(nrow(rndSites) > nrow(listSites))
      listSites = rndSites
  }
  return(as.matrix(listSites))
}

#' Read and buffer raster layers.
#' @description Read raster layers of environmental or other variables and crop them to a given extent around the known occurrences.
#' @param longlat Matrix of longitude and latitude (two columns) of species occurrence records.
#' @param layers Raster* object as defined by package raster.
#' @param ext Buffer around the known records used to crop layers. It is relative to the maximum distance between any two records.
#' @details If layers are not given, the function will read either 1 arc-minute (approx. 2km) or 5 arc-minutes (approx. 10km) resolution rasters from worldclim (Hijmans et al. 2005) and landcover (Tuanmu & Jetz 2014) if these are provided in the data folder of the package (see read.me file on how to do it).
#' @return A RasterStack object (Variables 1-19 = bioclim, 20 = elevation, 21-32 = proportion landcover, 33 = most common landcover).
#' @examples data(data.layers)
#' data(data.records)
#' par(mfrow=c(1,2))
#' raster::plot(data.layers[[1]])
#' points(data.records)
#' croppedLayers <- raster.read(data.records, data.layers, 0.1)
#' raster::plot(croppedLayers[[1]])
#' points(data.records)
#' @export
raster.read <- function(longlat, layers = NULL, ext = 1){

  xmin = min(longlat[,1])
  xmax = max(longlat[,1])
  xlen = xmax - xmin
  ymin = min(longlat[,2])
  ymax = max(longlat[,2])
  ylen = ymax - ymin

  if(is.null(layers)){          ##if no layers are provided read the ones available
    ##calculate species range and buffer around it
    if((xlen * ylen) < 10){
      layers <- raster::stack(raster("gis/red_2km_1.tif"))
      for(i in 2:33)
        layers <- raster::stack(layers, raster(paste("gis/red_2km_", i, ".tif", sep = "")))
    } else {
      layers <- raster::stack(raster("gis/red_10km_1.tif"))
      for(i in 2:33)
        layers <- raster::stack(layers, raster(paste("gis/red_10km_", i, ".tif", sep = "")))
    }
    ##determine longitude limits of species to check if crop and paste are needed around longitude 180 for Pacific species
    if(xmin < -90 && xmax > 90 && sum(longlat[longlat[,1] < 90 && longlat[,1] > -90,]) != 0){
      ##crop and merge layers
      rightHalf = crop(layers, c(0,180,raster::extent(layers)@ymin,raster::extent(layers)@ymax))
      raster::extent(rightHalf) <- c(-180,0,raster::extent(layers)@ymin,raster::extent(layers)@ymax)
      leftHalf = crop(layers, c(-180,0,raster::extent(layers)@ymin,raster::extent(layers)@ymax))
      raster::extent(leftHalf) <- c(0,180,raster::extent(layers)@ymin,raster::extent(layers)@ymax)
      layers <- merge(rightHalf, leftHalf)
      ##modify longlat
      for(i in 1:nrow(longlat))
        if(longlat[i,1] > 0)
          longlat[i,1] = longlat[i,1] - 180
      else
        longlat[i,1] = longlat[i,1] + 180
    }
  }

  if(length(ext) == 4)                          ##if absolute extent is given crop and return, else calculate buffer
    return(crop(layers, ext))

  if(xlen == 0)      ##in case some dimensions are inexistent consider equal to extent
    xlen = ext
  if(ylen == 0)
    ylen = ext

  ##calculate new extent of layers and crop
  ext = max(1, ((xlen + ylen) * ext))
  xmin <- max(raster::extent(layers)@xmin, xmin-ext)
  xmax <- min(raster::extent(layers)@xmax, xmax+ext)
  ymin <- max(raster::extent(layers)@ymin, ymin-ext)
  ymax <- min(raster::extent(layers)@ymax, ymax+ext)
  return(crop(layers, c(xmin,xmax,ymin,ymax)))
}

#' Uniformize raster layers.
#' @description Crop raster layers to minimum size and uniformize NA values accross layers.
#' @param layers Raster* object as defined by package raster.
#' @details Excludes all marginal rows and columns with only NA values and change values to NA if they are NA in any of the layers.
#' @return A Raster* object, same class as layers.
#' @examples data(data.layers)
#' raster::plot(raster.clean(data.layers))
#' @export
raster.clean <- function(layers){

  ##apply mask to have NAs everywhere where any layer has NAs
  maskLayer <- sum(layers)
  maskLayer[!is.na(maskLayer)] <- 1
  layers <- mask(layers, maskLayer)

  ##crop by excluding external rows and columns with NAs only
  layers <- trim(layers)

  return(layers)
}

#' Reduce dimensionality of raster layers.
#' @description Reduce the number of layers by either performing a PCA on them or by eliminating highly correlated ones.
#' @param layers Raster* object as defined by package raster.
#' @param method Either Principal Components Analysis (pca, default) or Pearson's correlation (cor).
#' @param n Number of layers to reduce to.
#' @param thres Value for pairwise Pearson's correlation above which one of the layers (randomly selected) is eliminated.
#' @details Using a large number of explanatory variables in models with few records may lead to overfitting. This function allows to avoid it as much as possible.
#' If both n and thres are given, n has priority. If method is not recognized and layers come from raster.read function, only landcover is reduced by using only the dominating landuse of each cell.
#' @return A RasterStack object.
#' @export
raster.reduce <- function(layers, method = "pca", n = NULL, thres = NULL){
  ##method = "pca, cor", if unrecognized method only reduce landcover but not climate
  out <- raster::stack()

  if(dim(layers)[3] == 33){          ##check if layers are obtained with raster.read
    out <- raster::stack(layers[[33]])
    layers = layers[[1:19]]
  }
  if(method == "cor"){                       ##if correlation
    if(is.null(n)){
      if(is.null(thres))
        thres = 0.7
      for(i in 1:dim(layers)[3]){                  ##delete layers until none are correlated above threshold
        cor = as.matrix(as.dist(layerStats(layers, 'pearson', na.rm = TRUE)[[1]]))
        if(max(cor) < thres)
          break
        corLayer = sample(which(cor == max(cor), arr.ind = TRUE)[,1],1)
        layers = layers[[-corLayer]]
      }
    } else {
      while (dim(layers)[3] > n){                   ##delete layers until reaching n layers
        cor = abs(as.matrix(as.dist(layerStats(layers, 'pearson', na.rm = TRUE)[[1]])))
        corLayer = sample(which(cor == max(cor), arr.ind = TRUE)[,1],1)
        layers = layers[[-corLayer]]
      }
    }
  } else if(method == "pca"){                                  ##if pca
    if(is.null(n))
      n = 3
    if(sum(!is.na(getValues(layers[[1]]))) > 2000)
      sr <- sampleRandom(layers, 1000)
    else
      sr <- sampleRandom(layers, as.integer(sum(!is.na(getValues(layers[[1]])))/2))
    pca <- prcomp(sr)
    layers <- raster::predict(layers, pca, index = 1:n)
    for(i in 1:n)
      names(layers[[i]]) <- paste("pca",i)
  }
  out <- raster::stack(out, layers)
  return(out)
}

#' Create longitude layer.
#' @description Create a layer depicting longitude based on any other.
#' @param layer RasterLayer object as defined by package raster.
#' @details Using longitude (and latitude) in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return A RasterLayer object.
#' @examples data(data.layers)
#' raster::plot(raster.long(data.layers[[1]]))
#' @export
raster.long <- function(layer){
  x <- rasterToPoints(layer)[,1:2]
  long <- rasterize(x, layer, x[,1])
  long <- mask(long, layer)
  names(long) <- "longitude"
  return(long)
}

#' Create latitude layer.
#' @description Create a layer depicting latitude based on any other.
#' @param layer RasterLayer object as defined by package raster.
#' @details Using latitude (and longitude) in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return A RasterLayer object.
#' @examples data(data.layers)
#' raster::plot(raster.lat(data.layers[[1]]))
#' @export
raster.lat <- function(layer){
  x <- rasterToPoints(layer)[,1:2]
  lat <- rasterize(x, layer, x[,2])
  lat <- mask(lat, layer)
  names(lat) <- "latitude"
  return(lat)
}

#' Create northness layer.
#' @description Create a layer depicting northness based on an elevation layer.
#' @param layer RasterLayer object of elevation as defined by package raster.
#' @details Using elevation, aspect can be calculated. Yet, it is a circular variable (0 = 360) and has to be converted to northness and eastness to be useful for modelling.
#' @return A RasterLayer object.
#' @examples data(data.sp)
#' raster::plot(raster.north(data.sp))
#' @export
raster.north <- function(layer){
  asp <- terrain(layer, opt = "aspect")
  return(cos(asp))
}

#' Create eastness layer.
#' @description Create a layer depicting eastness based on an elevation layer.
#' @param layer RasterLayer object of elevation as defined by package raster.
#' @details Using elevation, aspect can be calculated. Yet, it is a circular variable (0 = 360) and has to be converted to northness and eastness to be useful for modelling.
#' @return A RasterLayer object.
#' @examples data(data.sp)
#' raster::plot(raster.east(data.sp))
#' @export
raster.east <- function(layer){
  asp <- terrain(layer, opt = "aspect")
  return(sin(asp))
}

#' Predict species distribution.
#' @description Prediction of potential species distributions using maximum entropy (maxent).
#' @param longlat Matrix of longitude and latitude (two columns) of each occurrence record.
#' @param layers Raster* object as defined by package raster.
#' @param bg Background data as a matrix of longitude and latitude (two columns). If not defined 1000 points will be randomly selected.
#' @param categorical Vector of layer indices of categorical (as opposed to quantitative) data. If NULL the package will try to find them automatically based on the data itself.
#' @param thres Threshold of logistic output used for conversion of probabilistic to binary (presence/absence) maps. If 1 this will be the value that maximizes the sum of sensitivity and specificity.
#' @param polygon Used for a precautionary approach. If TRUE, all areas predicted as present but outside the minimum convex hull polygon encompassing all occurrence records are converted to absence. Only cells connected to other areas inside the polygon are kept.
#' @param eval Build a matrix with AUC, Kappa, TSS, EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).
#' @param jack If 0 no jackknife is performed. If 1, a jackknife with number of runs equivalent to the number of records is made. If > 1 these are the number of runs. For each run one random record is left out at a time and a new set of 1000 background points is chosen.
#' @details Builds maxent (maximum entropy) species distribution models (Phillips et al. 2004, 2006; Elith et al. 2011) using function maxent from R package dismo (Hijmans et al. 2016). Dismo requires the MaxEnt species distribution model software, a java program that can be downloaded from https://www.cs.princeton.edu/~schapire/maxent/. Put the file 'maxent.jar' in the 'java' folder of the dismo package. That is the folder returned by system.file("java", package="dismo"). You need MaxEnt version 3.3.3b or higher. Please note that this program (maxent.jar) cannot be redistributed or used for commercial or for-profit purposes.
#' @return Either one or two raster files (depending if jackknifes are performed, in which case the second is a probabilistic map from all the runs) and possibly a matrix with AUC, Kappa, TSS, EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model). Aggregate values are taken from maps after transformation of probabilities to incidence, with presence predicted for cells with presence for 50% or more of the runs.
#' @references Hijmans, R.J., Phillips, S., Leathwick, J., Elith, J. (2016) dismo: Species Distribution Modeling. R package version 1.0-15. https://CRAN.R-project.org/package=dismo
#' @references Phillips, S.J., Dudik, M., Schapire, R.E. (2004) A maximum entropy approach to species distribution modeling. Proceedings of the Twenty-First International Conference on Machine Learning. p. 655-662.
#' @references Phillips, S.J., Anderson, R.P., Schapire, R.E. (2006) Maximum entropy modeling of species geographic distributions. Ecological Modelling 190:231-259.
#' @references Elith, J., Phillips, S.J., Hastie, T., Dudik, M., Chee, Y.E., Yates, C.J. (2011) A statistical explanation of MaxEnt for ecologists. Diversity and Distributions 17:43-57.
#' @export
map.sdm <- function(longlat, layers, bg = NULL, categorical = NULL, thres = 1, polygon = TRUE, eval = TRUE, jack = 0){

  ##if jackknife is to be done
  if(jack > 0){
    if(jack == 1)
      jack = nrow(longlat)
    if(eval)
      jackEval = matrix(NA, nrow = 1, ncol = 7)
    jackMap <- rasterize(longlat, layers[[1]], field = 0, background = 0)
    pb <- txtProgressBar(min = 0, max = jack, style = 3)
    for(i in 1:jack){
      jackData <- longlat[-(sample.int(nrow(longlat),1)),]
      jackRun <- map.sdm(jackData, layers, bg, categorical, thres, polygon, eval, jack = 0)
      if(eval){
        jackEval <- rbind(jackEval, jackRun[[2]])
        jackRun <- jackRun[[1]]
      }
      jackMap <- jackMap + jackRun / jack
      setTxtProgressBar(pb, i)
    }
    jackMap01 <- reclassify(jackMap, matrix(c(0,0.499,0,0.499,1,1), ncol = 3, byrow = TRUE))
    if(eval){
      jackEval <- jackEval[-1,]
      clEval <- matrix(NA, nrow = 4, ncol = 7)
      colnames(clEval) <- c("AUC", "Kappa", "TSS", "EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)")
      rownames(clEval) <- c("Aggregate", "LowCL", "Median", "UpCL")
      clEval[1,1:3] <- colMeans(jackEval[,1:3])
      clEval[2,] <- apply(jackEval, 2,  quantile, probs= 0.025, na.rm = TRUE)
      clEval[3,] <- apply(jackEval, 2,  quantile, probs= 0.5, na.rm = TRUE)
      clEval[4,] <- apply(jackEval, 2,  quantile, probs= 0.975, na.rm = TRUE)
      clEval[1,4] <- eoo(longlat)
      clEval[1,5] <- eoo(jackMap01)
      clEval[1,6] <- aoo(jackMap01, longlat)
      clEval[1,7] <- aoo(jackMap01)
      return(list(jackMap01, jackMap, clEval))
    } else {
      return (jackMap01)
    }
  }

  if(is.null(bg))
    bg <- dismo::randomPoints(layers, 1000)                                ##extract background points (to use as absence)
  ##if no categorical variables are given try to figure out which
  if(is.null(categorical))
    categorical <- find.categorical(layers)

  model <- dismo::maxent(layers, longlat, a = bg, factors = categorical) ##build model
  p <- raster::predict(model, layers)                                     ##do prediction
  e <- dismo::evaluate(longlat, bg, model, layers)                       ##do evaluation of model
  if(thres >= 1)
    thres <- dismo::threshold(e)$spec_sens                                   ##extract threshold from evaluation
  p <- reclassify(p, matrix(c(0,thres,0,thres,1,1), nrow=2, byrow = TRUE))  ##convert to presence/absence

  if(polygon){
    vertices <- chull(longlat)
    vertices <- c(vertices, vertices[1])
    vertices <- longlat[vertices,]
    poly = Polygon(vertices)
    poly = Polygons(list(poly),1)
    poly = SpatialPolygons(list(poly))    ##original EOO, before modelling
    patches <- clump(p, gaps=FALSE)       ##individual patches, numbered
    selPatches <- unique(extract(patches, poly, df = TRUE, weights = TRUE)$clumps) ##which patches are inside original EOO
    #p <- replace(p, !(patches %in% selPatches), 0)
    p[!(as.vector(patches) %in% as.vector(selPatches))] <- 0
  }

  if(eval){
    e <- dismo::evaluate(longlat, bg, model, layers, thres)                  ##do evaluation of model with threshold
    auc <- e@auc
    kappa <- e@kappa
    sensitivity <- as.numeric(e@TPR/(e@TPR+e@FNR))
    specificity <- as.numeric(e@TNR/(e@TNR+e@FPR))
    tss <- sensitivity + specificity - 1
    eooRaw <- eoo(longlat)
    eooModel <- eoo(p)
    aooRaw <- aoo(p, longlat)
    aooModel <- aoo(p)
    txtEval <- matrix(c(auc, kappa, tss, eooRaw, eooModel, aooRaw, aooModel), nrow = 1)
    colnames(txtEval) <- c("AUC", "Kappa", "TSS", "EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)")
    return(list(p, txtEval))
  } else {
    return(p)
  }
}

#' Map species distribution of habitat specialist.
#' @description Mapping of all habitat areas where the species is known to occur.
#' @param longlat Matrix of longitude and latitude (two columns) of each occurrence record.
#' @param layer RasterLayer object representing the presence/absence (1/0) of a single habitat type.
#' @param polygon If TRUE, all habitat patches inside the minimum convex hull polygon encompassing all occurrence records are converted to presence.
#' @param eval Build a matrix with EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).
#' @details In many cases a species has a very restricted habitat and we generally know where it occurs. In such cases using the distribution of the habitat patches may be enough to map the species.
#' @return One raster file and possibly a matrix with EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).
#' @export
# @examples
# data(data.records)
# data(data.sp)
# raster::plot(map.habitat(data.records, data.sp, eval = FALSE))
# points(data.records)
map.habitat <- function(longlat, layer, polygon = TRUE, eval = TRUE){
  if(polygon){
    vertices <- chull(longlat)
    vertices <- c(vertices, vertices[1])
    vertices <- longlat[vertices,]
    poly = Polygon(vertices)
    poly = Polygons(list(poly),1)
    poly = SpatialPolygons(list(poly))    ##minimum convex polygon
    patches <- clump(layer, gaps=FALSE)       ##individual patches, numbered
    selPatches <- unique(extract(patches, poly, df = TRUE, weights = TRUE)$clumps) ##which patches are inside polygon
  } else {
    patches <- clump(layer, gaps=FALSE)       ##individual patches, numbered
    selPatches <- unique(extract(patches, longlat, df = TRUE, weights = TRUE)$clumps) ##which patches have the species
  }
  layer <- mask(patches%in%selPatches, layer) * layer
  if(eval){
    eooRaw <- eoo(longlat)
    eooModel <- eoo(layer)
    aooRaw <- aoo(layer, longlat)
    aooModel <- aoo(layer)
    txtEval <- matrix(c(eooRaw, eooModel, aooRaw, aooModel), nrow = 1)
    colnames(txtEval) <- c("EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)")
    return(list(layer, txtEval))
  } else {
    return(layer)
  }
}

#' Map recorded species distribution of species.
#' @description Mapping of all cells where the species is known to occur.
#' @param longlat Matrix of longitude and latitude (two columns) of each occurrence record.
#' @param layers Raster* object as defined by package raster. Any raster with the relevant extent and cell size can be used.
#' @param eval Build a matrix with EOO and AOO calculated from occurrence records only.
#' @details To be used if either information on the species is very scarce (and it is not possible to model the species distribution) or, on the contrary, complete (and there is no need to model the distribution).
#' @return One raster file and possibly a matrix with EOO and AOO.
#' @export
# @examples
# data(data.records)
# data(data.layers)
# raster::plot(map.points(data.records, data.layers, eval = FALSE))
# points(data.records)
map.points <- function(longlat, layers, eval = TRUE){
  p <- rasterize(longlat, layers[[1]], field = 1, background = 0)
  maskLayer <- sum(layers)
  maskLayer[!is.na(maskLayer)] <- 1
  p <- mask(p, maskLayer)
  if(eval){
    eooRaw <- eoo(longlat)
    aooRaw <- aoo(p, longlat)
    txtEval <- matrix(c(eooRaw, aooRaw), nrow = 1)
    colnames(txtEval) <- c("EOO", "AOO")
    return(list(p, txtEval))
  } else {
    return(p)
  }
}

## #' Predict species distribution change in time.
## #' @description Prediction and projection in time of potential species distribution using maximum entropy (maxent).
## #' @param longlat Matrix of longitude and latitude (two columns) of each occurrence record.
## #' @param layers Raster* object as defined by package raster, used for modelling current conditions.
## #' @param projectLayers Raster* object as defined by package raster, used for projection into the future.
## #' @param bg Background data as a matrix of longitude and latitude (two columns). If not defined 1000 points will be randomly selected.
## #' @param categorical Vector of layer indices of categorical (as opposed to quantitative) data. If NULL the package will try to find them automatically based on the data itself.
## #' @param thres Threshold of logistic output used for conversion of probabilistic to binary (presence/absence) maps. If 1 this will be the value that maximizes the sum of sensitivity and specificity.
## #' @param polygon Used for a precautionary approach. If TRUE, all areas predicted as present but outside the minimum convex hull polygon encompassing all occurrence records are converted to absence. Only cells connected to other areas inside the polygon are kept for both present and future projections.
## #' @param jack If 0 no jackknife is performed. If 1, a jackknife with number of runs equivalent to the number of records is made. If > 1 these are the number of runs. For each run one random record is left out at a time and a new set of 1000 background points is chosen.
## #' @details Builds maxent (maximum entropy) species distribution models (Phillips et al. 2004, 2006; Elith et al. 2011) using function maxent from R package dismo (Hijmans et al. 2016) for both the present and future.
## #' @return Either three or six RasterLayer (depending if jackknifes are performed, in which case the second trio are probabilistic maps from all the runs) and a matrix with Present AOO, Future AOO, Gain, Keep and Loss.
## #' @references Hijmans, R.J., Phillips, S., Leathwick, J., Elith, J. (2016) dismo: Species Distribution Modeling. R package version 1.0-15. https://CRAN.R-project.org/package=dismo
## #' @references Phillips, S.J., Dudik, M., Schapire, R.E. (2004) A maximum entropy approach to species distribution modeling. Proceedings of the Twenty-First International Conference on Machine Learning. p. 655-662.
## #' @references Phillips, S.J., Anderson, R.P., Schapire, R.E. (2006) Maximum entropy modeling of species geographic distributions. Ecological Modelling 190:231-259.
## #' @references Elith, J., Phillips, S.J., Hastie, T., Dudik, M., Chee, Y.E., Yates, C.J. (2011) A statistical explanation of MaxEnt for ecologists. Diversity and Distributions 17:43-57.
## #' @export
## map.change <- function(longlat, layers, projectLayers, bg = NULL, categorical = NULL, thres = 1, polygon = FALSE, jack = 0){
##   options(warn=-1)
#   ##if jackknife is to be done
#   if(jack > 0){
#     if(jack == 1)
#       jack = nrow(longlat)
#     jackEval = matrix(NA, nrow = 1, ncol = 5)
#     jackMap <- rasterize(longlat, layers[[1]], field = 0, background = 0)
#     jackMap <- raster::stack(jackMap, jackMap, jackMap)
#     jackMap01 <- jackMap
#     pb <- txtProgressBar(min = 0, max = jack, style = 3)
#     for(i in 1:jack){
#       jackData <- longlat[-(sample.int(nrow(longlat),1)),]
#       jackRun <- map.change(jackData, layers, projectLayers, bg, categorical, thres, polygon, jack = 0)
#       jackEval <- rbind(jackEval, jackRun[[4]])
#       jackRun <- list(jackRun[[1]], jackRun[[2]], jackRun[[3]])
#       for(j in 1:3)
#         jackMap[[j]] <- jackMap[[j]] + jackRun[[j]] / jack
#       setTxtProgressBar(pb, i)
#     }
#     for(i in 1:2)
#       jackMap01[[i]] <- reclassify(jackMap[[i]], matrix(c(0,0.5,0,0.5,1,1), ncol = 3, byrow = TRUE))
#     jackMap01[[3]] <- jackMap01[[2]] * 2 - jackMap01[[1]] ##gain = 2, kept = 1, loss = -1, never exists = 0
#     jackEval <- jackEval[-1,]
#     clEval <- matrix(NA, nrow = 4, ncol = 5)
#     colnames(clEval) <- c("Present AOO", "Future AOO", "Gain", "Keep", "Loss")
#     rownames(clEval) <- c("Aggregate", "LowCL", "Median", "UpCL")
#     clEval[1,1] <- aoo(jackMap01[[1]])
#     clEval[1,2] <- aoo(jackMap01[[2]])
#     clEval[1,3] <- cellStats((raster::area(jackMap01[[3]]) * subs(jackMap01[[3]], as.data.frame(matrix(c(0,0,1,0,-1,0,2,1), ncol = 2, byrow = TRUE)))),sum)
#     clEval[1,4] <- cellStats((raster::area(jackMap01[[3]]) * subs(jackMap01[[3]], as.data.frame(matrix(c(0,0,1,1,-1,0,2,0), ncol = 2, byrow = TRUE)))),sum)
#     clEval[1,5] <- cellStats((raster::area(jackMap01[[3]]) * subs(jackMap01[[3]], as.data.frame(matrix(c(0,0,1,0,-1,1,2,0), ncol = 2, byrow = TRUE)))),sum)
#     clEval[2,] <- apply(jackEval, 2,  quantile, probs= 0.025, na.rm = TRUE)
#     clEval[3,] <- apply(jackEval, 2,  quantile, probs= 0.5, na.rm = TRUE)
#     clEval[4,] <- apply(jackEval, 2,  quantile, probs= 0.975, na.rm = TRUE)
#     options(warn=0)
#     return(list(jackMap01, jackMap, clEval))
#   }
#
#   ##if no background points are given randomly sample them
#   if(is.null(bg))
#     bg <- randomPoints(layers, 1000)                                ##extract background points (to use as absence)
#   ##if no categorical variables are given try to figure out which
#   if(is.null(categorical))
#     categorical = find.categorical(layers)
#
#   model <- dismo::maxent(layers, longlat, a = bg, factors = categorical) ##build model
#   present <- raster::predict(model, layers)                             ##do prediction for the present
#   future <- raster::predict(model, projectLayers)
#   e <- dismo::evaluate(longlat, bg, model, layers)                       ##do evaluation of model
#   if(thres >= 1)
#     thres <- threshold(e)$spec_sens                                   ##extract threshold from evaluation
#   present <- reclassify(present, matrix(c(0,thres,0,thres,1,1), nrow=2, byrow = TRUE))  ##convert to presence/absence
#   future <- reclassify(future, matrix(c(0,thres,0,thres,1,1), nrow=2, byrow = TRUE))  ##convert to presence/absence
#
#   if(polygon){                          ##if species is limited in dispersal ability
#     vertices <- chull(longlat)
#     vertices <- c(vertices, vertices[1])
#     vertices <- longlat[vertices,]
#     poly = Polygon(vertices)
#     poly = Polygons(list(poly),1)
#     poly = SpatialPolygons(list(poly))    ##original EOO, before modelling
#     ##present
#     patches <- clump(present, gaps=FALSE)       ##individual patches, numbered, present
#     selPatches <- unique(extract(patches, poly, df = TRUE, weights = TRUE)$clumps) ##which patches are inside original EOO
#     present <- replace(present, !(patches %in% selPatches), 0)
#     ##future
#     patches <- clump(future, gaps=FALSE)       ##individual patches, numbered, future
#     selPatches <- unique(extract(patches, poly, df = TRUE, weights = TRUE)$clumps) ##which patches are inside original EOO
#     future <- replace(future, !(patches %in% selPatches), 0)
#   }
#
#   spDiff <- future * 2 - present ##gain = 2, kept = 1, loss = -1, never exists = 0
#   txtChange <- rep(NA,5)
#   names(txtChange) <- c("Present AOO", "Future AOO", "Gain", "Keep", "Loss")
#   txtChange[1] <- cellStats((raster::area(present) * present), sum)
#   txtChange[2] <- cellStats((raster::area(future) * future), sum)
#   txtChange[3] <- cellStats((raster::area(spDiff) * subs(spDiff, as.data.frame(matrix(c(0,0,1,0,-1,0,2,1), ncol = 2, byrow = TRUE)))),sum)
#   txtChange[4] <- cellStats((raster::area(spDiff) * subs(spDiff, as.data.frame(matrix(c(0,0,1,1,-1,0,2,0), ncol = 2, byrow = TRUE)))),sum)
#   txtChange[5] <- cellStats((raster::area(spDiff) * subs(spDiff, as.data.frame(matrix(c(0,0,1,0,-1,1,2,0), ncol = 2, byrow = TRUE)))),sum)
#   return(list(present, future, spDiff, txtChange))
# }

#' Predict species distribution made simple.
#' @description Prediction of species distributions, multiple species simultaneously, and output of maps, klms and relevant data to files. All using a single step.
#' @param longlat data.frame of species names, longitude and latitude (three columns) of each occurrence record.
#' @param layers Raster* object as defined by package raster. If NULL they are read from data files.
#' @param file Name of output csv file with all results. If NULL it is named "Results_All.csv".
#' @param minimum Minimum number of occurrence records to perform a maxent model. If these are lower than the minimum, the function will return a map with presences points only.
#' @param jack If 0 no jackknife is performed. If 1, a jackknife with number of runs equivalent to the number of records is made. If > 1 these are the number of runs. For each run one random record is left out at a time and a new set of 1000 background points is chosen.
#' @details Builds maxent species distribution models and outputs maps in both pdf and kml format, plus a file with EOO, AOO and a list of countries where the species is predicted to be present.
#' @return Writes maps, kmls and all information to a file.
#' @export
map.easy <- function(longlat, layers = NULL, file = NULL, minimum = 3, jack = 0){
  spNames <- unique(longlat[,1])
  nSp <- length(spNames)
  if (jack > 0) {
    res <- matrix(NA, nrow = nSp, ncol = 11)
    colnames(res) <- c("EOO (raw)", "EOO (aggregate)", "EOO (LowCL)", "EOO (Median)", "EOO (UpCL)", "AOO (raw)", "AOO (aggregate)", "AOO (LowCL)", "AOO (Median)", "AOO (UpCL)", "Countries")
  } else {
    res <- matrix(NA, nrow = nSp, ncol = 5)
    colnames(res) <- c("EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)", "Countries")
  }
  rownames(res) <- spNames
  if(is.null(layers))
    newLayers <- TRUE
  else
    newLayers <- FALSE
  for(s in 1:nSp){
    spData <- longlat[longlat[,1] == spNames[s], -1]
    cat("\nModelling species", s, "of", nSp, "-", toString(spNames[s]), "\n")
    spData <- thin(spData)
    if(newLayers){
      ##check if worldclim and landcover layers are available, if not ask to run red.setup
      if(!file.exists("gis/red_2km_1.tif")){
        return(warning("Please check if worldclim and land cover layers are available in the gis directory under the working directory of r. If not, run function red.setup."))
      }
      if(!file.exists(paste(.libPaths()[[1]], "/dismo/java/maxent.jar", sep=""))){
        return(warning("Please check if maxent.jar is available in the dismo package directory named java. If not, run function red.setup."))
      }
      layers <- raster.read(spData)
      layers <- raster.reduce(layers)
      layers <- raster.clean(layers)
      layers <- raster::stack(layers, raster.long(layers[[1]]), raster.lat(layers[[1]]))
    }
    if(nrow(spData) >= minimum){
      p <- map.sdm(spData, layers, jack = jack)
    } else {
      p <- map.points(spData, layers)
    }
    writeRaster(p[[1]], paste(toString(spNames[s]), ".asc", sep=""), overwrite = TRUE)
    map.draw(spData, p[[1]], spNames[s], print = TRUE)
    if(nrow(spData) >= minimum)
      kml(p[[1]], paste(toString(spNames[s]), ".kml", sep=""))
    else
      kml(spData, paste(toString(spNames[s]), ".kml", sep=""), minimum)
    countryList <- countr(p[[1]])
    if(nrow(spData) >= minimum && jack > 0){
      writeRaster(p[[2]], paste(toString(spNames[s]), "_prob.asc", sep=""), overwrite = TRUE)
      map.draw(spData, p[[2]], paste(toString(spNames[s]), "_prob", sep = ""), legend = TRUE, print = TRUE)
      res[s,] <- c(p[[3]][1,4], p[[3]][1:4,5], p[[3]][1,6], p[[3]][1:4,7], toString(countryList))
    } else if (nrow(spData) >= minimum && jack <= 0){
      res[s,] <- c(p[[2]][1,4:7], toString(countryList))
    } else if (jack > 0){
      res[s,] <- c(p[[2]][1,c(1,1,1,1,1,2,2,2,2,2)], toString(countryList))
    } else if (jack <= 0){
      res[s,] <- c(p[[2]][1,c(1,1,2,2)], toString(countryList))
    }
      write.csv(res[s,], paste(toString(spNames[s]), ".csv", sep = ""))
  }
  if(is.null(file))
    write.csv(res, "Results_All.csv")
  else
    write.csv(res, toString(file))
  return(res)
}

#' Map creation.
#' @description Creates maps ready to print in pdf or other formats.
#' @param longlat Matrix of longitude and latitude (two columns) of each occurrence record.
#' @param layer RasterLayer object representing the presence/absence map for the species.
#' @param spName String of species name.
#' @param countries If TRUE country borders are drawn.
#' @param scale If TRUE a distance scale in km is drawn.
#' @param legend If TRUE the legend for the map is drawn.
#' @param sites If TRUE the record locations are drawn.
#' @param mcp If TRUE the minimum convex polygon representing the Extent of Occurrence is drawn.
#' @param print If TRUE a pdf is saved instead of the output to the console.
#' @examples data(data.records)
#' data(data.sp)
#' par(mfrow = c(1,2))
#' map.draw(data.records, layer = data.sp, mcp = TRUE)
#' @export
map.draw <- function(longlat = NULL, layer, spName,  countries = FALSE, scale = TRUE, legend = FALSE, sites = TRUE, mcp = FALSE, print = FALSE){
  data.worldborders <- NULL
  data(data.worldborders, envir = environment())
  if (countries){
    layer[layer == 0] <- NA
    raster::plot(layer, main = toString(spName), legend = legend, xlab = "longitude", ylab = "latitude", col = "forestgreen")
    lines(data.worldborders)
  } else {
    raster::plot(layer, main = spName, legend = legend, colNA = "lightblue", xlab = "longitude", ylab = "latitude")
  }
  if (scale)
    scalebar(type="bar", divs = 2, below = "km")
  if (sites && !is.null(longlat))
    points(longlat)
  if (mcp){
    e <- rasterToPoints(layer, fun = function(dat){dat == 1})   ##convert raster to points
    vertices <- chull(e[,1], e[,2])
    vertices <- c(vertices, vertices[1])
    vertices <- e[vertices,c(1,2)]
    poly <- SpatialPolygons(list(Polygons(list(Polygon(vertices)),1)))
    raster::plot(poly, add = TRUE)
  }
  if(print){
    dev.copy(device = pdf, file = paste(toString(spName), ".pdf", sep=""))
    dev.off()
  }
}

#' Extent of Occurrence (EOO).
#' @description Calculates the Extent of Occurrence of a species based on either records or predicted distribution.
#' @param spData Either a matrix of longitude and latitude (two columns) of each occurrence record or a presence/absence (1/0) map as a RasterLayer object.
#' @details EOO is calculated as the minimum convex polygon covering all known or predicted sites for the species.
#' @return A single value in km2.
#' @examples data(data.records)
#' data(data.sp)
#' eoo(data.records)
#' eoo(data.sp)
#' @export
eoo <- function(spData){
  if(class(spData) == "RasterLayer"){
    e <- rasterToPoints(spData, fun = function(dat){dat == 1})   ##convert raster to points
    vertices <- chull(e[,1], e[,2])
    vertices <- c(vertices, vertices[1])
    vertices <- e[vertices,c(1,2)]
  } else {
    vertices <- chull(spData)
    vertices <- c(vertices, vertices[1])
    vertices <- spData[vertices,]
  }
  return(areaPolygon(vertices)/1000000)
}

#' Area of Occupancy (AOO).
#' @description Calculates the Area of Occupancy of a species based on either known records or predicted distribution.
#' @param layer RasterLayer object. If AOO is to be calculated based on known records, any raster with the relevant extent and cell size can be used. If AOO is to be calculated based on predicted distribution, a raster of presence/absence (1/0) is needed.
#' @param longlat Matrix of longitude and latitude (two columns) of each occurrence record. Only needed if AOO is to be calculated based on known records.
#' @details AOO is calculated as the area of all known or predicted cells (usually 2x2 km) for the species. Known cells are used if long is given, predicted cells
#' @return A single value in km2.
#' @examples data(data.records)
#' data(data.sp)
#' aoo(data.sp, data.records)
#' aoo(data.sp)
#' @export
aoo <- function(layer, longlat = NULL){
  if(!is.null(longlat)){
    layer = rasterize(longlat, layer, field = 1, background = 0)
  }
  return(cellStats((raster::area(layer) * layer), sum))
}

#' Countries of occurrence.
#' @description Extracts the names or ISO codes of countries of occurrence of a species based on either records or predicted distribution.
#' @param spData Either a matrix of longitude and latitude (two columns) of each occurrence record or a presence/absence map as a RasterLayer object.
#' @param ISO Outputs either country names (FALSE) or ISO codes (TRUE).
#' @details Country boundaries and designations are based on data(data.worldborders) from package maptools.
#' @return A vector with country names or codes.
#' @examples data(data.records)
#' data(data.sp)
#' countr(data.records)
#' countr(data.sp, ISO = TRUE)
#' @export
countr <- function(spData, ISO = FALSE){
  data.worldborders <- NULL
  data(data.worldborders, envir = environment())
  if(class(spData) == "RasterLayer")
    spData <- rasterToPoints(spData, fun = function(dat){dat == 1})   ##convert raster to points
  spCountries <- sp::over(sp::SpatialPoints(spData), sp::SpatialPolygons(data.worldborders@polygons))
  if(ISO)
    spCountries <- unique(data.worldborders@data[spCountries,])$ISO2
  else
    spCountries <- unique(data.worldborders@data[spCountries,])$NAME
  spCountries <- sort(as.vector(spCountries[!is.na(spCountries)]))
  return(spCountries)
}

#' Output kml files.
#' @description Creates kml files for Google Maps as required by IUCN guidelines.
#' @param spData Either a matrix of longitude and latitude (two columns) of each occurrence record or a presence/absence map as a RasterLayer object.
#' @param filename The name of file to save, should end with .kml.
#' @param minimum If the number of records is lower than the minimum no polygon in drawn, only circles of about 10km radius around each point.
#' @export
kml <- function(spData, filename, minimum = 3){
  if(nrow(spData) < minimum){
    poly = list()
    for(i in 1:nrow(spData)){
      rad = 0.1 # radius
      pts = seq(0, 2 * pi, length.out = 100)
      xy = cbind(spData[i, 1] + rad * sin(pts), spData[i, 2] + rad * cos(pts))
      poly[[i]] = Polygon(xy)
    }
    poly = Polygons(poly,1)
  } else {
    if (class(spData) == "RasterLayer"){
      e <- rasterToPoints(spData, fun = function(dat){dat == 1})   ##convert raster to points
      vertices <- chull(e[,1], e[,2])
      vertices <- c(vertices, vertices[1])
      vertices <- e[vertices,c(1,2)]
    } else {
      vertices <- chull(spData)
      vertices <- c(vertices, vertices[1])
      vertices <- spData[vertices,]
    }
    poly = Polygon(vertices)
    poly = Polygons(list(poly),1)
  }
  kmlPolygon(poly, filename)
}

#' Red List Index.
#' @description Calculates the Red List Index (RLI) for a group of species.
#' @param spData Either a vector with species assessment categories for a single point in time or a matrix with two points in time in different columns (species x date). Values can be textual (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).
#' @param boot If TRUE bootstrapping for statistical significance is performed on both values per date and the trend between dates.
#' @param runs Number of runs for bootstrapping
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to Extinct/Extinct in the Wild (5).
#' Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
#' Importantly, the RLI is based on true improvements or deteriorations in the status of species, i.e. genuine changes. It excludes category changes resulting from, e.g., new knowledge (Butchart et al. 2007).
#' The RLI approach helps to develop a better understanding of which taxa, regions or ecosystems are declining or improving.
#' Juslen et al. (2016a, b) suggested the use of bootstrapping to search for statistical significance when comparing taxa or for trends in time of the index
#' @return Either a vector (if no two dates are given) or a matrix with the RLI values and, if bootstrap is performed, their confidence limits and significance.
#' @references Butchart, S.H.M., Stattersfield, A.J., Bennun, L.A., Shutes, S.M., Akcakaya, H.R., Baillie, J.E.M., Stuart, S.N., Hilton-Taylor, C. & Mace, G.M. (2004) Measuring global trends in the status of biodiversity: Red List Indices for birds. PloS Biology,
#' @references Butchart, S.H.M., Akcakaya, H.R., Chanson, J., Baillie, J.E.M., Collen, B., Quader, S., Turner, W.R., Amin, R., Stuart, S.N. & Hilton-Taylor, C. (2007) Improvements to the Red List index. PloS One,
#' @references Juslen, A., Cardoso, P., Kullberg, J., Saari, S. & Kaila, L. (2016a) Trends of extinction risk for Lepidoptera in Finland: the first national Red List Index of butterflies and moths. Insect Conservation and Diversity, 9: 118-123.
#' @references Juslen, A., Pykala, J., Kuusela, S., Kaila, L., Kullberg, J., Mattila, J., Muona, J., Saari, S. & Cardoso, P. (2016b) Application of the Red List Index as an indicator of habitat change. Biodiversity and Conservation, 25: 569-585.
#' @examples rliData <- matrix(c("LC","LC","EN","EN","EX","EX","LC","CR","CR","EX"), ncol = 2, byrow = TRUE)
#' colnames(rliData) <- c("2000", "2010")
#' rli(rliData[,1])
#' rli(rliData[,1], boot = TRUE)
#' rli(rliData)
#' rli(rliData, boot = TRUE)
#' @export
rli <- function (spData, boot = FALSE, runs = 1000){
  ##RLI with phylogenetic or functional data to be implemented soon
  ##rli <- function (spData, tree = NULL, boot = FALSE, runs = 1000){
  tree = NULL   ##to add soon

  ##if only one point in time is given
  if(is.null(dim(spData)))
    return(rli.calc(spData, tree, boot, runs))  ##return either 1 or 3 values

  ##if two points in time are given
  ts <- apply(spData, 2, function(x) rli.calc(x, tree, boot = FALSE))
  sl <- ts[2] - ts[1]
  if(!boot){
    res <- matrix(c(ts, sl), nrow = 1)
    colnames(res) <- c(colnames(spData), "Change")
    rownames(res) <- c("Raw")
    return(res)
  } else {
    tr <- apply(spData, 2, function(x) rli.calc(x, tree, boot, runs))
    p = 0
    rndSl = rep(NA, runs)
    for(r in 1:runs){
      rndSl[r] <- rli.calc(spData[,2], tree, boot, 1)[2] - rli.calc(spData[,1], tree, boot, 1)[2]
      if(sign(sl) < sign(rndSl[r]) || sign(sl) > sign(rndSl[r]))
        p = p + 1
    }
    p = p / runs
    rndSl = quantile(rndSl, c(0.025, 0.5, 0.975))
    res <- matrix(c(ts[1], tr[,1], ts[2], tr[,2], sl, rndSl), nrow = 4, ncol = 3)
    colnames(res) <- c(colnames(spData), "Change")
    rownames(res) <- c("Raw", "LowCL", "Median", "UpCL")
    return(list("Values" = res, "P_change" = p))
  }
}

#' Red List Index for multiple groups.
#' @description Calculates the Red List Index (RLI) for multiple groups of species.
#' @param spData A matrix with group names (first column) and species assessment categories for one or two points in time (remaining columns). Values can be textual (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).
#' @param boot If TRUE bootstrapping for statistical significance is performed on both values per date and the trend between dates.
#' @param runs Number of runs for bootstrapping
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to Extinct/Extinct in the Wild (5).
#' Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
#' Importantly, the RLI is based on true improvements or deteriorations in the status of species, i.e. genuine changes. It excludes category changes resulting from, e.g., new knowledge (Butchart et al. 2007).
#' The RLI approach helps to develop a better understanding of which taxa, regions or ecosystems are declining or improving.
#' Juslen et al. (2016a, b) suggested the use of bootstrapping to search for statistical significance when comparing taxa or for trends in time of the index
#' @return A matrix with the RLI values and, if bootstrap is performed, their confidence limits and significance.
#' @references Butchart, S.H.M., Stattersfield, A.J., Bennun, L.A., Shutes, S.M., Akcakaya, H.R., Baillie, J.E.M., Stuart, S.N., Hilton-Taylor, C. & Mace, G.M. (2004) Measuring global trends in the status of biodiversity: Red List Indices for birds. PloS Biology,
#' @references Butchart, S.H.M., Akcakaya, H.R., Chanson, J., Baillie, J.E.M., Collen, B., Quader, S., Turner, W.R., Amin, R., Stuart, S.N. & Hilton-Taylor, C. (2007) Improvements to the Red List index. PloS One,
#' @references Juslen, A., Cardoso, P., Kullberg, J., Saari, S. & Kaila, L. (2016a) Trends of extinction risk for Lepidoptera in Finland: the first national Red List Index of butterflies and moths. Insect Conservation and Diversity, 9: 118-123.
#' @references Juslen, A., Pykala, J., Kuusela, S., Kaila, L., Kullberg, J., Mattila, J., Muona, J., Saari, S. & Cardoso, P. (2016b) Application of the Red List Index as an indicator of habitat change. Biodiversity and Conservation, 25: 569-585.
#' @examples rliData <- matrix(c("LC","LC","EN","EN","EX","EX","LC","CR","CR","EX"), ncol = 2, byrow = TRUE)
#' colnames(rliData) <- c("2000", "2010")
#' rliData <- cbind(c("Arthropods","Arthropods","Birds","Birds","Birds"), rliData)
#' rli.multi(rliData[,1:2])
#' rli.multi(rliData[,1:2], boot = TRUE)
#' rli.multi(rliData)
#' rli.multi(rliData, boot = TRUE)
#' @export
rli.multi <- function (spData, boot = FALSE, runs = 1000){
  ##RLI with phylogenetic or functional data to be implemented soon
  ##rli.multi <- function (spData, tree = NULL, boot = FALSE, runs = 1000){
  tree = NULL ##to add soon

  groups <- unique(spData[,1])
  nGroups <- length(groups)
  if(ncol(spData) == 2 && !boot){
    res <- matrix(NA, nrow = nGroups, ncol = 1)
  } else if((ncol(spData) == 2 && boot) || (ncol(spData) == 3 && !boot)){
    res <- matrix(NA, nrow = nGroups, ncol = 3)
  } else {
    res <- matrix(NA, nrow = nGroups, ncol = 13)
    colnames(res) <- c(paste(colnames(spData)[2], "(raw)"), paste(colnames(spData)[2], "(lowCL)"), paste(colnames(spData)[2], "(median)"), paste(colnames(spData)[2], "(upCL)"), paste(colnames(spData)[3], "(raw)"), paste(colnames(spData)[3], "(lowCL)"), paste(colnames(spData)[3], "(median)"), paste(colnames(spData)[3], "(upCL)"), "Change (raw)", "Change (lowCL)", "Change (median)", "Change (upCL)", "p (change)")
  }
  row.names(res) <- groups
  for(g in 1:nGroups){
    if(is.null(tree)){
      v <- rli(spData[spData[,1] == groups[g],-1], boot = boot, runs = runs)
      if(ncol(res) < 13){
        res[g,] <- v
        colnames(res) <- colnames(v)
      } else {
        res[g,1:4] <- v$Values[,1]
        res[g,5:8] <- v$Values[,2]
        res[g,9:12] <- v$Values[,3]
        res[g,13] <- v$P_change
      }
    } else {
      i = 1
    }
  }
  return(res)
}

#' Occurrence records for a fictional species in Finland.
#'
#' Occurrence records for a fictional species in Finland.
#'
#' @docType data
#' @keywords datasets
#' @name data.records
#' @usage data(data.records)
#' @format Matrix of longitude and latitude (two columns) of species occurrence records.
NULL

#' Geographic range for a fictional species in Finland.
#'
#' Geographic range for a fictional species in Finland.
#'
#' @docType data
#' @keywords datasets
#' @name data.sp
#' @usage data(data.sp)
#' @format RasterLayer object as defined by package raster.
NULL

#' Environmental layers for Finland.
#'
#' Average annual temperature, total annual precipitation and landcover for Finland (Hijmans et al. 2005, Tuanmu & Jetz 2014).
#'
#' @docType data
#' @keywords datasets
#' @name data.layers
#' @usage data(data.layers)
#' @format RasterStack object as defined by package raster.
#' @references Hijmans, R.J., Cameron, S.E, Parra, J.L., Jones, P.G. & Jarvis A. (2005) Very high resolution interpolated climate surfaces for global land areas. International Journal of Climatology, 25: 1965-1978.
#' @references Tuanmu, M.-N. & Jetz, W. (2014) A global 1-km consensus land-cover product for biodiversity and ecosystem modeling. Global Ecology and Biogeography, 23: 1031-1045.
NULL
#'
#'
#' World country borders.
#'
#' World country borders.
#'
#' @docType data
#' @keywords datasets
#' @name data.worldborders
#' @usage data(data.worldborders)
#' @format SpatialPolygonsDataFrame.
NULL
