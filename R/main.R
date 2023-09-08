#####RED - IUCN Redlisting Tools
#####Version 1.6.1 (2023-08-01)
#####By Pedro Cardoso & Vasco Branco
#####Maintainer: vasco.branco@helsinki.fi
#####Reference: Cardoso, P.(2017) An R package to facilitate species red list assessments according to the IUCN criteria. Biodiversity Data Journal 5: e20530 doi: 10.3897/BDJ.5.e20530
#####Changed from v1.6.0:
#####rgeos import fix, terra update fixes

################################################################################
##################################MAIN FUNCTIONS################################
################################################################################

#' Download taxon records from GBIF.
#' @description Downloads species or higher taxon data from GBIF and outputs non-duplicate records with geographical coordinates.
#' @param taxon Taxon name.
#' @details As always when using data from multiple sources the user should be careful and check if records "make sense". This can be done by either ploting them in a map (e.g. using red::map.draw()) or using red::outliers().
#' @return A data.frame with longitude and latitude, plus species names if taxon is above species.
#' @examples
#' rec = records("Nephila senegalensis")
#' plot(rec)
#' @export
records <- function(taxon){
  taxon = unlist(strsplit(taxon, split = " ")[[1]])
  dat <- dismo::gbif(taxon[1], paste(taxon[2], "*", sep = ""))
  dat <- dat[c("species","lon","lat")] #filter columns
  dat <- dat[!(is.na(dat$lon) | is.na(dat$lat)),] #filter rows
  dat <- unique(dat)       #delete duplicate rows
  colnames(dat) <- c("Species", "long", "lat")
  if (length(taxon) == 1){      #if genus
    dat[which(is.na(dat[,1])),1] <- paste(taxon, "sp.")
  } else {                #if species
    dat <- dat[,-1]
  }
  return(dat)
}

#' Move records to closest non-NA cell.
#' @description Identifies and moves presence records to cells with environmental values.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers SpatRaster* object as defined by package raster.
#' @param buffer Maximum distance in map units that a record will move. If 0 all NA records will be changed.
#' @details Often records are in coastal or other areas for which no environmental data is available. This function moves such records to the closest cells with data so that no information is lost during modelling.
#' @return A matrix with new coordinate values.
#' @examples
#' rast <- terra::rast(matrix(c(rep(NA,100), rep(1,100), rep(NA,100)), ncol = 15))
#' pts <- cbind(runif(100, 0, 5), runif(100, 0, 15))
#' terra::plot(rast)
#' points(pts)
#' 
#' pts <- move(pts, rast)
#' terra::plot(rast)
#' points(pts)
#' @export
move <- function(longlat, layers, buffer = 0) {
  if (dim(layers)[3] > 1) {
    layers <- layers[[1]]
  }

  if (is(longlat, "matrix")) {
    longlat <- as.data.frame(longlat)
  }

  if (terra::crs(layers) == "") {
    terra::crs(layers) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
  }

  values <- terra::extract(layers, longlat, ID = FALSE) # get values of each record

  suppressWarnings(
    for (i in which(is.na(values))) { # if a value is NA, move it
      # Distance does not work when the crs is "". needs an exception
      distRaster <- terra::distance(
        layers,
        terra::vect(as.data.frame(longlat)[i, ], # remove if enforced at start
          geom = colnames(longlat),
          crs = terra::crs(layers)
        )
      )
      distRaster <- terra::mask(distRaster, layers)
      vmin <- terra::where.min(distRaster)

      if (buffer <= 0 || buffer > vmin) {
        # vmin = terra::as.points(distRaster, function(x) x == vmin)
        longlat[i, ] <- terra::xyFromCell(distRaster, vmin[2]) # vmin[1,1:2]
      }
    }
  )
  return(longlat)
}

#' Visual detection of outliers.
#' @description Draws plots of sites in geographical (longlat) and environmental (2-axis PCA) space.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers Raster* object as defined by package raster. It can be any set of environmental layers thought to allow the identification of environmental outliers.
#' @details Erroneous data sources or errors in transcriptions may introduce outliers that can be easily detected by looking at simple graphs of geographical or environmental space.
#' @return A data.frame with coordinate values and distance to centroid in pca is returned. Two plots are drawn for visual inspection. The environmental plot includes row numbers for easy identification of possible outliers.
#' @examples
#' records = red.examples("red.records")
#' layers = red.examples("red.layers")
#' outliers(records, layers[[1:3]])
#' @export
outliers <- function(longlat, layers){
  userpar <- par(no.readonly = TRUE) 
  on.exit(par(userpar))
  if(dim(layers)[3] == 33)      #if layers come from read
    pca <- raster.reduce(layers[[1:19]], n = 2)
  else
    pca <- raster.reduce(layers, n = 2)
  
  ##extract pca values from longlat
  pca <- as.data.frame(terra::extract(pca, longlat))
  goodRows <-  which(!is.na(pca[,1]))
  pca <- pca[goodRows,]
  longlat <- longlat[goodRows,]
  par(mfrow = c(1,2))
  map.draw(longlat, layers[[1]], spName = "Geographical")
  plot(pca, main = "Environmental")
  centroid = colMeans(pca)
  text(centroid[1], centroid[2], label = "X")
  for(i in 1:nrow(pca)){
    text(pca[i,1], pca[i,2], label = row.names(longlat)[i])
  }
  
  ##build new matrix ordered by distance to centroid
  dist2centroid = apply(pca, 1, function(x) dist(rbind(x, centroid)))
  out = as.data.frame(cbind(longlat, dist2centroid))
  out = out[order(-dist2centroid),]
  return(out)
}

#' Spatial thinning of occurrence records.
#' @description Thinning of records with minimum distances either absolute or relative to the species range.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param distance Distance either in relative terms (proportion of maximum distance between any two records) or in raster units.
#' @param relative If TRUE, represents the proportion of maximum distance between any two records. If FALSE, is in raster units.
#' @param runs Number of runs
#' @details Clumped distribution records due to ease of accessibility of sites, emphasis of sampling on certain areas in the past, etc. may bias species distribution models.
#' The algorithm used here eliminates records closer than a given distance to any other record. The choice of records to eliminate is random, so a number of runs are made and the one keeping more of the original records is chosen.
#' @return A matrix of species occurrence records separated by at least the given distance.
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
    if(nSites < 40){ #if limited number of sites use all data
      maxDist = 0
      for(x in 1:(nSites-1)){
        for(y in (x+1):nSites){
          maxDist = max(maxDist,((longlat[x,1]-longlat[y,1])^2+(longlat[x,2]-longlat[y,2])^2)^.5)
        }
      }
    } else { #if many sites use hypothenusa of square encompassing all of them
      horiDist = max(longlat[,1]) - min(longlat[,1])
      vertDist = max(longlat[,2]) - min(longlat[,2])
      maxDist = (horiDist^2 + vertDist^2)^0.5
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
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers Raster* object as defined by package raster.
#' @param ext Either extent of map or buffer around the known records used to crop layers. If buffer, it is relative to the maximum distance between any two records.
#' @details If layers are not given, the function will read either 30 arc-second (approx. 1km) or 5 arc-minutes (approx. 10km) resolution rasters from worldclim (Fick & Hijmans 2017) and landcover (Tuanmu & Jetz 2014) if red.setup() is run previously.
#' @return A RasterStack object (If no layers are given: Variables 1-19 = bioclim, 20 = elevation, 21-32 = proportion landcover, 33 = most common landcover).
#' @references Fick, S.E. & Hijmans, R.J. (2017) Worldclim 2: new 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology, in press.
#' @references Tuanmu, M.-N. & Jetz, W. (2014) A global 1-km consensus land-cover product for biodiversity and ecosystem modeling. Global Ecology and Biogeography, 23: 1031-1045.
#' @examples
#' layers = red.examples("red.layers")
#' records = red.examples("red.records")
#' par(mfrow=c(1,2))
#' 
#' terra::plot(layers[[1]])
#' points(records)
#' 
#' croppedLayers <- raster.read(records, layers, 0.1)
#' terra::plot(croppedLayers[[1]])
#' points(records)
#' @export
raster.read <- function(longlat, layers = NULL, ext = 1){

  xmin = min(longlat[,1])
  xmax = max(longlat[,1])
  xlen = xmax - xmin
  ymin = min(longlat[,2])
  ymax = max(longlat[,2])
  ylen = ymax - ymin

  if(is.null(layers)){          ##if no layers are provided read the ones available
    gisdir = red.getDir()
    layers = terra::rast()
    ##calculate species range and buffer around it
    if (eoo(longlat) < 200000) {
      layers <- c(terra::rast(paste(gisdir, "red_1km_1.tif", sep = "")))
      ### These iterations over the downloaded layers need to be checked.
      for (i in 2:33) {
        layers <- c(layers, terra::rast(paste(gisdir, "red_1km_", i, ".tif", sep = "")))
      }
    } else {
      layers <- c(terra::rast(paste(gisdir, "red_10km_1.tif", sep = "")))
      for (i in 2:33) {
        layers <- c(layers, terra::rast(paste(gisdir, "red_10km_", i, ".tif", sep = "")))
      }
    }
    ##determine longitude limits of species to check if crop and paste are needed around longitude 180 for Pacific species
    if(xmin < -90 && xmax > 90 && sum(longlat[longlat[,1] < 90 && longlat[,1] > -90,]) != 0){
      ##crop and merge layers
      rightHalf <- terra::crop(layers, c(
        0, 180,
        terra::ymin(layers),
        terra::ymax(layers)
      ))

      terra::ext(rightHalf) <- c(
        -180, 0,
        terra::ymin(layers),
        terra::ymax(layers)
      )

      leftHalf <- terra::crop(layers, c(
        -180, 0,
        terra::ymin(layers),
        terra::ymax(layers)
      ))
      terra::ext(leftHalf) <- c(
        0, 180,
        terra::ymin(layers),
        terra::ymax(layers)
      )

      layers <- terra::merge(rightHalf, leftHalf)
      ##modify longlat
      for(i in 1:nrow(longlat))
        if(longlat[i,1] > 0)
          longlat[i,1] = longlat[i,1] - 180
      else
        longlat[i,1] = longlat[i,1] + 180
    }
  }
  
  ##if absolute extent is given crop and return, else calculate buffer
  if(length(ext) == 4) return(crop(layers, ext))
  
  ##in case some dimensions are inexistent consider equal to extent
  if(xlen == 0) xlen = ext
  if(ylen == 0) ylen = ext

  ##calculate new extent of layers and crop
  ext = max(1, ((xlen + ylen) * ext))
  xmin <- max(terra::xmin(layers), xmin-ext)
  xmax <- min(terra::xmax(layers), xmax+ext)
  ymin <- max(terra::ymin(layers), ymin-ext)
  ymax <- min(terra::ymax(layers), ymax+ext)
  layers <- crop(layers, c(xmin,xmax,ymin,ymax))
  return(layers)
}

#' Uniformize raster layers.
#' @description Crop raster layers to minimum size possible and uniformize NA values across layers.
#' @param layers Raster* object as defined by package raster.
#' @details Excludes all marginal rows and columns with only NA values and change values to NA if they are NA in any of the layers.
#' @return A Raster* object, same class as layers.
#' @examples
#' layers = red.examples("red.layers")
#' terra::plot(raster.clean(layers))
#' @export
raster.clean <- function(layers){

  ##apply mask to have NAs everywhere where any layer has NAs
  maskLayer <- sum(layers)
  maskLayer[!is.na(maskLayer)] <- 1
  layers <- terra::mask(layers, maskLayer)

  ##crop by excluding external rows and columns with NAs only
  layers <- terra::trim(layers)

  return(layers)
}

#' Reduce dimensionality of raster layers.
#' @description Reduce the number of layers by either performing a PCA on them or by eliminating highly correlated ones.
#' @param layers Raster* object as defined by package raster.
#' @param method Either Principal Components Analysis ("pca", default) or Pearson's correlation ("cor").
#' @param n Number of layers to reduce to.
#' @param thres Value for pairwise Pearson's correlation above which one of the layers (randomly selected) is eliminated.
#' @details Using a large number of explanatory variables in models with few records may lead to overfitting. This function allows to avoid it as much as possible.
#' If both n and thres are given, n has priority. If method is not recognized and layers come from raster.read function, only landcover is reduced by using only the dominating landuse of each cell.
#' @return A RasterStack object.
#' @export
raster.reduce <- function(layers, method = "pca", n = NULL, thres = NULL){
  ##method = "pca, cor", if unrecognized method only reduce landcover but not climate

  if(dim(layers)[3] == 33){          ##check if layers are obtained with read
    out <- c(layers[[33]])
    layers = layers[[1:19]]
  }
  
  # HANDLE COR
  if(method == "cor"){                       ##if correlation
    if(is.null(n)){
      if(is.null(thres))
        thres = 0.7
      for(i in 1:dim(layers)[3]){                  ##delete layers until none are correlated above threshold
        cor = as.matrix(
          as.dist(terra::layerCor(layers, 'pearson', na.rm = TRUE)[[1]])
        )
        
        if(max(cor) < thres)
          break
        corLayer = sample(which(cor == max(cor), arr.ind = TRUE)[,1],1)
        layers = layers[[-corLayer]]
      }
    } else {
      while (dim(layers)[3] > n){                   ##delete layers until reaching n layers
        cor = abs(as.matrix(as.dist(
          terra::layerCor(layers, 'pearson', na.rm = TRUE)[[1]])))
        
        corLayer = sample(which(cor == max(cor), arr.ind = TRUE)[,1],1)
        layers = layers[[-corLayer]]
      }
    }
  } else if(method == "pca"){                                  ##if pca
    if(is.null(n))
      n = 3
    if(sum(!is.na(terra::values(layers[[1]], mat = FALSE))) > 2000)
      sr <- terra::spatSample(layers, 1000)
    else
      sr <- terra::spatSample(layers, as.integer(sum(!is.na(terra::values(layers[[1]], mat = FALSE)))/2), na.rm = TRUE) # added na.rm
    pca <- prcomp(sr)
    layers <- terra::predict(layers, pca, index = 1:n)
    for(i in 1:n){
      names(layers[[i]]) <- paste("pca",i)
    }
  }
  if(dim(layers)[3] == 33){
    out <- c(layers, out)
  } else {
    out <- layers
  }
  return(out)
}

#' Create distance layer.
#' @description Creates a layer depicting distances to records using the minimum, average, distance to the minimum convex polygon or distance taking into account a cost surface.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers Raster* object as defined by package raster to serve as model to create distance layer. Cost surface in case of param ="cost".
#' @param type text string indicating whether the output should be the "minimum", "average", or "mcp" distance to all records. "mcp" means the distance to the minimum convex polygon encompassing all records.
#' @details Using distance to records in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return A RasterLayer object.
#' @examples
#' layers = red.examples("red.layers")
#' alt = layers[[3]]
#' records = red.examples("red.records")
#' par(mfrow=c(3,2))
#' terra::plot(alt)
#' points(records)
#' 
#' terra::plot(raster.distance(records, alt))
#' terra::plot(raster.distance(records, alt, type = "average"))
#' terra::plot(raster.distance(records, alt, type = "mcp"))
#' @export
raster.distance <- function(longlat, layers, type = "minimum"){
  ### terra::distance() works substantially different from how it did in raster.
  ### see documentation and red::move()
  
  if(dim(layers)[3] > 1)
    layers <- layers[[1]]
  
  layers_template = terra::classify(!is.na(layers), c(TRUE, 0))
  
  if(type == "average"){
    for(d in 1:nrow(longlat)){
      layers <- c(
        layers,
        terra::distance(
          layers_template,
          terra::vect(longlat[d, ],
                      geom = colnames(longlat),
                      crs = terra::crs(layers)
          )
        )
      )
    }
    layers = layers[[2:dim(layers)[3]]]
    layers = terra::mean(layers)
    layers <- terra::mask(layers, layers_template) ###
    names(layers) <- "average distance"
  } else if (type == "mcp"){
    vertices <- grDevices::chull(longlat)
    vertices <- c(vertices, vertices[1])
    vertices <- longlat[vertices,]
    poly = sp::Polygon(vertices)
    poly = sp::Polygons(list(poly),1)
    poly = sp::SpatialPolygons(list(poly))    ##minimum convex polygon
    longlat = terra::as.points(terra::rasterize(terra::vect(poly), layers))
    layers <- terra::mask(terra::distance(layers, longlat), layers)
    names(layers) <- "mcp distance"
  } else if (type == "cost"){
    ### These do not work with terra objects! an alternative is needed.
    layers <- gdistance::transition(layers, function(x) 1/mean(x), 8) ###
    layers <- gdistance::geoCorrection(layers) 
    layers <- gdistance::accCost(layers, as.matrix(longlat))
    names(layers) <- "cost distance"
  } else {
    distRaster <- terra::distance(
      layers,
      terra::vect(longlat,
                  geom = colnames(longlat),
                  crs = terra::crs(layers)
      )
    )
    layers <- terra::mask(distRaster, layers)
    names(layers) <- "minimum distance"
  }
  return(layers)
}

#' Create longitude layer.
#' @description Create a layer depicting longitude based on any other.
#' @param layers Raster* object as defined by package raster.
#' @details Using longitude (and latitude) in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return A RasterLayer object.
#' @examples
#' layers = red.examples("red.layers")
#' terra::plot(raster.long(layers))
#' @export
raster.long <- function(layers){
  if(dim(layers)[3] > 1) {
    layers <- layers[[1]]
  }
  x <- terra::as.points(layers)[,1:2]
  long <- terra::rasterize(terra::crds(x), y = layers, values = terra::crds(x)[,2] )
  long <- terra::mask(long, layers)
  names(long) <- "longitude"
  return(long)
}

#' Create latitude layer.
#' @description Create a layer depicting latitude based on any other.
#' @param layers Raster* object as defined by package raster.
#' @details Using latitude (and longitude) in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return A RasterLayer object.
#' @examples
#' layers = red.examples("red.layers")
#' terra::plot(raster.lat(layers[[1]]))
#' @export
raster.lat <- function(layers){
  if(dim(layers)[3] > 1){
    layers <- layers[[1]]
  }
  x <- terra::as.points(layers)[,1:2]
  lat <- terra::rasterize(terra::crds(x), y = layers, values = terra::crds(x)[,1] )
  lat <- terra::mask(lat, layers)
  names(lat) <- "latitude"
  return(lat)
}

#' Create eastness layer.
#' @description Create a layer depicting eastness based on an elevation layer.
#' @param dem RasterLayer object of elevation (a digital elevation model - DEM) as defined by package raster.
#' @details Using elevation, aspect can be calculated. Yet, it is a circular variable (0 = 360) and has to be converted to northness and eastness to be useful for modelling.
#' @return A RasterLayer object.
#' @examples layers = red.examples("red.layers")
#' terra::plot(raster.east(layers[[3]]))
#' @export
raster.east <- function(dem){
  asp <- terra::terrain(dem, v = "aspect")
  return(sin(asp))
}

#' Create northness layer.
#' @description Create a layer depicting northness based on an elevation layer.
#' @param dem RasterLayer object of elevation (a digital elevation model - DEM) as defined by package raster.
#' @details Using elevation, aspect can be calculated. Yet, it is a circular variable (0 = 360) and has to be converted to northness and eastness to be useful for modelling.
#' @return A RasterLayer object.
#' @examples
#' layers = red.examples("red.layers")
#' terra::plot(raster.north(layers[[3]]))
#' @export
raster.north <- function(dem){
  asp <- terra::terrain(dem, v = "aspect")
  return(cos(asp))
}

#' Extent of Occurrence (EOO).
#' @description Calculates the Extent of Occurrence of a species based on either records or predicted distribution.
#' @param spData spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (either 0/1 or probabilistic values).
#' @details EOO is calculated as the minimum convex polygon covering all known or predicted sites for the species.
#' @return A single value in km2 or a vector with lower confidence limit, consensus and upper confidence limit (probabilities 0.975, 0.5 and 0.025 respectively).
#' @examples
#' records = red.examples("red.records")
#' range = red.examples("red.range")
#' eoo(records)
#' eoo(range)
#' @export
eoo <- function(spData){
  if(is(spData, "SpatRaster")){
    if(!all(terra::as.matrix(spData) == floor(terra::as.matrix(spData)), na.rm = TRUE)){ #if probabilistic map
      upMap <- terra::classify(
        spData,
        matrix(c(0, 0.025, 0, 0.025, 1, 1),
               ncol = 3, byrow = TRUE
        )
      )
      consensusMap <- terra::classify(
        spData,
        matrix(c(0, 0.499, 0, 0.499, 1, 1),
               ncol = 3, byrow = TRUE
        )
      )
      downMap <- terra::classify(
        spData,
        matrix(c(0, 0.975, 0, 0.975, 1, 1),
               ncol = 3, byrow = TRUE
        )
      )
      
      area <- c(eoo(downMap), eoo(consensusMap), eoo(upMap))
    } else {
      if (terra::xmax(spData) <= 180) {  #if longlat data
        e <- as.data.frame(terra::as.points(spData), geom = "XY")
        e <- e[ e[,1] == 1 , ]
        
        
        vertices <- chull(e[,2], e[,3]) # from 1, 2
        if(length(vertices) < 3) return(0)
        vertices <- c(vertices, vertices[1])
        vertices <- e[vertices, c(2, 3)] # from 1, 2. seems redundant
        area = geosphere::areaPolygon(vertices)/1000000
      } else {
        spData[spData < 1] <- NA
        spData <- as.data.frame(terra::as.points(spData))
        vertices <- chull(spData)
        if(length(vertices) < 3) return(0)
        vertices <- c(vertices, vertices[1])
        vertices <- spData[vertices,]
        area = 0
        for(i in 1:(nrow(vertices)-1))
          area = area + (as.numeric(vertices[i,1])*as.numeric(vertices[(i+1),2]) - as.numeric(vertices[i,2])*as.numeric(vertices[(i+1),1]))
        area = abs(area/2000000)
      }
    }
  } else if (ncol(spData) == 2){
    vertices <- chull(spData)
    if(length(vertices) < 3) return(0)
    vertices <- c(vertices, vertices[1])
    vertices <- spData[vertices,]
    if(max(spData) <= 180) {  #if longlat data
      area = geosphere::areaPolygon(vertices)/1000000
    } else { #if square data in meters
      area = 0
      for(i in 1:(nrow(vertices)-1))
        area = area + (as.numeric(vertices[i,1])*as.numeric(vertices[(i+1),2]) - as.numeric(vertices[i,2])*as.numeric(vertices[(i+1),1]))
      area = abs(area/2000000)
    }
  } else {
    return(warning("Data format not recognized"))
  }
  return(round(area))
}

#' Area of Occupancy (AOO).
#' @description Calculates the Area of Occupancy of a species based on either known records or predicted distribution.
#' @param spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (either 0/1 or probabilistic values).
#' @details AOO is calculated as the area of all known or predicted cells for the species. The resolution will be 2x2km as required by IUCN.
#' @return A single value in km2 or a vector with lower confidence limit, consensus and upper confidence limit (probabilities 0.975, 0.5 and 0.025 respectively).
#' @examples
#' range = red.examples("red.range")
#' aoo(range)
#' @export
aoo <- function(spData){
  if (is(spData, "SpatRaster")){ #if rasterlayer
    if(sum(terra::minmax(spData)) == 0){  #if no data (empty raster)
      area = 0
    } else if(!all(terra::as.matrix(spData) == floor(terra::as.matrix(spData)), na.rm = TRUE)){ #if probabilistic map
      upMap <- terra::classify(
        spData,
        matrix(c(0, 0.025, 0, 0.025, 1, 1),
               ncol = 3, byrow = TRUE
        )
      )
      consensusMap <- terra::classify(
        spData,
        matrix(c(0, 0.499, 0, 0.499, 1, 1),
               ncol = 3, byrow = TRUE
        )
      )
      downMap <- terra::classify(
        spData,
        matrix(c(0, 0.975, 0, 0.975, 1, 1),
               ncol = 3, byrow = TRUE
        )
      )
      area <- c(eoo(downMap), eoo(consensusMap), eoo(upMap))
    } else {
      if (terra::xmax(spData) <= 180) {  #if longlat data
        if(terra::res(spData)[1] > 0.05){ #if resolution is > 1km use area of cells rounded to nearest 4km
          area = round(terra::global(
            (terra::cellSize(spData, unit = "km") * spData), fun = "sum", na.rm = T)/4)*4  
        } else {
          spData[spData < 1] <- NA
          spData <- terra::as.points(spData)
          spData <- data.frame(x = terra::crds(spData)[,1],
                               y = terra::crds(spData)[,2])

          if(nrow(unique(spData)) == 1){
            area = 4
          } else {
            spData <- longlat2utm(spData)
            spData = floor(spData/2000)
            ncells = nrow(unique(spData))
            area = ncells * 4
          }
        }
      } else { #if square data in meters
        spData[spData < 1] <- NA
        spData <- terra::as.points(spData) # rasterToPoints(spData)
        spData <- data.frame(x = terra::crds(spData)[,1],
                             y = terra::crds(spData)[,2])
        spData = floor(spData/2000)
        ncells = nrow(unique(spData))
        area = ncells * 4
      }
    }
  } else if (ncol(spData) == 2){
    if (max(spData) <= 180) {  #if longlat data
      spData <- longlat2utm(spData)
      spData = floor(spData/2000)
      ncells = nrow(unique(spData))
      area = ncells * 4
    } else { #if square data in meters
      spData = floor(spData/2000)
      ncells = nrow(unique(spData))
      area = ncells * 4
    }
  } else {
    return(warning("Data format not recognized!"))
  }
  return(round(area))
}

#' Elevation limits.
#' @description Calculates the elevation (or depth) limits (range) of a species based on either known records or predicted distribution.
#' @param spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (0/1 values).
#' @param dem RasterLayer object. Should be a digital elevation model (DEM) of the relevant area. If not given the function will try to read it from base data, only works with longlat data.
#' @details Maximum and minimum elevation are calculated based on the DEM.
#' @return A vector with two values (min and max) in meters above (or below) sea level.
#' @examples
#' records = red.examples("red.records")
#' range = red.examples("red.range")
#' layers = red.examples("red.layers")
#' dem = layers[[3]]
#' elevation(records, dem)
#' elevation(range, dem)
#' @export
elevation <- function(spData, dem = NULL){
  if(is(spData, "data.frame")){
    spData = as.matrix(spData)
  }
  if(!is(spData, "SpatRaster")){ #if no rasterlayer is given but just a matrix of longlat.
    if(is.null(dem) && max(spData) <= 180){
      gisdir = red.getDir()
      dem <- terra::rast(paste(gisdir, "red_1km_20.tif", sep =""))
      dem <- terra::crop(dem, c(min(spData[,1])-0.1, max(spData[,1]+0.1), min(spData[,2])-0.1, max(spData[,2])+0.1))
    }
    spData = terra::rasterize(spData, dem, background = NA) #create a layer of presence based on the dem
  } else if (is.null(dem)){
    gisdir = red.getDir()
    dem <- terra::rast(paste(gisdir, "red_1km_20.tif", sep = ""))
    dem <- terra::crop(dem, spData)
  }
  spData <- terra::lapp(c(spData, dem), fun = function(x,y){(x*y)} )
  out <- c(terra::minmax(spData)[1], terra::minmax(spData)[2])
  names(out) <- c("Min", "Max")
  return(round(out))
}

#' Countries of occurrence.
#' @description Extracts the names or ISO codes of countries of occurrence of a species based on either records or predicted distribution.
#' @param spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (0/1 values).
#' @param zone UTM zone if data is in metric units.
#' @param ISO Outputs either country names (FALSE) or ISO codes (TRUE).
#' @details Country boundaries and designations are based on data(worldborders) from package maptools.
#' @return A vector with country names or codes.
#' @examples
#' records = red.examples("red.records")
#' range = red.examples("red.range")
#' countries(records)
#' countries(range, ISO = TRUE)
#' @export
countries <- function(spData, zone = NULL, ISO = FALSE){
  if ((is(spData, "SpatRaster") && terra::xmax(spData) > 180) || (!is(spData, "SpatRaster") && max(spData) > 180))   ##if need to project to longlat
    spData <- utm2longlat(spData, zone)

  worldborders <- red.examples("worldborders")
  
  if(is(spData, "SpatRaster")){
    spData <- terra::as.points(spData)   ##convert raster to points
    spData <- terra::as.data.frame(spData, geom = "XY")
    spData <- spData[ spData[,1] == 1 , 2:3]
  }
  
  if(ISO){
    countryList <- unique(terra::extract(worldborders, spData))$ISO2
  } else {
    countryList <- unique(terra::extract(worldborders, spData))$NAME
  }
  countryList <- unique(sort(as.vector(countryList[!is.na(countryList)])))
  return(countryList)
}

#' Output kml files.
#' @description Creates kml files for Google Maps as required by IUCN guidelines.
#' @param spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (0/1 values).
#' @param zone UTM zone if data is in metric units.
#' @param filename The name of file to save, should end with .kml.
#' @param mapoption Type of representation, any of "points", "eoo" or "aoo".
#' @param smooth Smooths the kml lines as per IUCN guidelines. Higher values represent smoother polygons.
#' @param rad radius of circles in degrees if mapoption is "points". It can be the same value for all points or a vector with length equal to number of records in spData representing associated error. The default is about 10km (0.1 degrees) as per IUCN guidelines.
#' @return A kml with polygon or circles around records.
#' @export
kml <- function(spData, zone = NULL, filename, mapoption = "aoo", smooth = 0, rad = 0.1){
  
  worldborders = red.examples("worldborders")
  if ((is(spData, "SpatRaster") && terra::xmax(spData) > 180) || (!is(spData, "SpatRaster") && max(spData) > 180))   ##if need to project to longlat
    spData <- utm2longlat(spData, zone)
  
  if(mapoption == "aoo" && is(spData, "SpatRaster")){
    spData[spData != 1] <- NA
    spData <- terra::as.polygons(spData, dissolve = TRUE)

    #simplify
    if(smooth > 0){
      trytol <- c(seq(0.001,0.01,0.001),
                  seq(0.02,0.1,0.01),
                  seq(0.2,1,0.1),
                  2:10,
                  seq(20,100,10),
                  seq(200,1000,100),
                  seq(2000,10000,1000),
                  seq(20000,100000,10000),
                  seq(200000,1000000,100000))
      for (i in trytol){
        if(!is(try(terra::simplifyGeom(spData, tolerance = (1 / i)), silent = TRUE), "try-error")){ ### rgeos::gSimplify
          spData <- terra::simplifyGeom(spData, tolerance = (smooth / (i*10)))
          break
        }
      }

      #cut to coast
      spData = terra::intersect(spData, worldborders) ### is this the intended effect?

      #round
      smooth = smooth * 100
      polys = methods::slot(spData@polygons[[1]], "Polygons")
      
      terra::plot(spData[,2])
      
      spData = as(spData, "Spatial")
      polys = methods::slot(spData@polygons[[1]], "Polygons")
      
      spData <- sp::SpatialPolygons(
        Srl = lapply(1:length(polys),
                     function(x){
                       p <- polys[[x]]

                       #applying spline.poly function for smoothing polygon edges
                       px <- methods::slot(polys[[x]], "coords")[,1]
                       py <- methods::slot(polys[[x]], "coords")[,2]
                       bz <- spline.poly(methods::slot(polys[[x]], "coords"),smooth, k=3)
                       bz <- rbind(bz, bz[1,])
                       methods::slot(p, "coords") <- bz

                       # create Polygons object
                       poly <- sp::Polygons(list(p), ID = x)
                     }
        )
      )
      spData <- sp::SpatialPolygonsDataFrame(spData, data = data.frame(ID = 1:length(spData)))
      
      terra::writeVector(spData, filename, filetype = "KML")
    } else {
      terra::writeVector(spData, filename, filetype = "KML")
    }
    
  } else if(mapoption == "points" || ( is(spData, "SpatRaster") && aoo(spData) <= 8) || nrow(spData) < 3){
    poly = list()
    for(i in 1:nrow(spData)){
      pts = seq(0, 2 * pi, length.out = 100)
      if(length(rad) == 1)
        xy = cbind(spData[i, 1] + rad * sin(pts), spData[i, 2] + rad * cos(pts))
      else
        xy = cbind(spData[i, 1] + rad[i] * sin(pts), spData[i, 2] + rad[i] * cos(pts))
      poly[[i]] = Polygon(xy)
    }
    poly = sp::Polygons(poly,1)
    terra::writeVector(poly, filename, filetype = "KML")
  } else {
    if (is(spData, "SpatRaster")){
      e <- terra::as.points(spData, fun = function(dat){dat == 1})   ##convert raster to points
      vertices <- chull(e[,1], e[,2])
      vertices <- c(vertices, vertices[1])
      vertices <- e[vertices,c(1,2)]
    } else {
      vertices <- chull(spData)
      vertices <- c(vertices, vertices[1])
      vertices <- spData[vertices,]
    }
    poly = sp::Polygon(vertices)
    poly = sp::Polygons(list(poly),1)
    terra::writeVector(poly, filename, filetype = "KML")
  }
}