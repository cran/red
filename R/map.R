#' Predict species distributions.
#' @description Prediction of potential species distributions using maximum entropy (maxent).
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of each occurrence record.
#' @param layers Predictor variables, a Raster* object as defined by package raster.
#' @param error Vector of spatial error in longlat (one element per row of longlat) in the same unit as longlat. Used to move any point randomly within the error radius.
#' @param year Vector of sampling years in longlat (one element per row of longlat). Used to exclude old records with a given probability proportional to time passed since sampling (never excluded only for current year).
#' @param idconf Vector of identification confidence in longlat (one element per row of longlat). Used to exclude uncertain records with a given probability. Can be on any scale where max values are certain (e.g. from 1 - very uncertain to 10 - holotype).
#' @param categorical Vector of layer indices of categorical (as opposed to quantitative) data. If NULL the package will try to find them automatically based on the data.
#' @param thres Threshold of logistic output used for conversion of probabilistic to binary (presence/absence) maps. If 0 this will be the value that maximizes the sum of sensitivity and specificity.
#' @param testpercentage Percentage of records used for testing only. If 0 all records will be used for both training and testing.
#' @param mcp Used for a precautionary approach. If TRUE, all areas predicted as present but outside the minimum convex hull polygon encompassing all occurrence records are converted to absence. Exceptions are cells connected to other areas inside the polygon.
#' @param points If TRUE, force map to include cells with presence records even if suitable habitat was not identified.
#' @param eval If TRUE, build a matrix with AUC, Kappa, TSS, EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).
#' @param runs If <= 0 no ensemble modelling is performed. If > 0, ensemble modelling with n runs is made. For each run, a new random sample of occurrence records (if testpercentage > 0), background points and predictive variables (if subset > 0) are chosen. In the ensemble model, each run is weighted as max(0, (runAUC - 0.5)) ^ 2.
#' @param subset Number of predictive variables to be randomly selected from layers for each run if runs > 0. If <= 0 all layers are used on all runs. Using a small number of layers is usually better than using many variables for rare species, with few occurrence records (Lomba et al. 2010, Breiner et al. 2015).
#' @details Builds maxent (maximum entropy) species distribution models (Phillips et al. 2004, 2006; Elith et al. 2011) using function maxent from R package dismo (Hijmans et al. 2017). Dismo requires the MaxEnt species distribution model software, a java program that can be downloaded from http://biodiversityinformatics.amnh.org/open_source/maxent. Copy the file 'maxent.jar' into the 'java' folder of the dismo package. That is the folder returned by system.file("java", package="dismo"). You need MaxEnt version 3.3.3b or higher. Please note that this program (maxent.jar) cannot be redistributed or used for commercial or for-profit purposes.
#' @return List with either one or two raster objects (depending if ensemble modelling is performed, in which case the second is a probabilistic map from all the runs) and, if eval = TRUE, a matrix with AUC, Kappa, TSS, EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model). Aggregate values are taken from maps after transformation of probabilities to incidence, with presence predicted for cells with ensemble values > 0.5.
#' @references Breiner, F.T., Guisan, A., Bergamini, A., Nobis, M.P. (2015) Overcoming limitations of modelling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6: 1210-1218.
#' @references Hijmans, R.J., Phillips, S., Leathwick, J., Elith, J. (2017) dismo: Species Distribution Modeling. R package version 1.1-4. https://CRAN.R-project.org/package=dismo
#' @references Lomba, A., Pellissier, L., Randin, C.F., Vicente, J., Moreira, F., Honrado, J., Guisan, A. (2010) Overcoming the rare species modelling paradox: a novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation, 143: 2647-2657.
#' @references Phillips, S.J., Dudik, M., Schapire, R.E. (2004) A maximum entropy approach to species distribution modeling. Proceedings of the Twenty-First International Conference on Machine Learning. p. 655-662.
#' @references Phillips, S.J., Anderson, R.P., Schapire, R.E. (2006) Maximum entropy modeling of species geographic distributions. Ecological Modelling, 190: 231-259.
#' @references Elith, J., Phillips, S.J., Hastie, T., Dudik, M., Chee, Y.E., Yates, C.J. (2011) A statistical explanation of MaxEnt for ecologists. Diversity and Distributions, 17: 43-57.
#' @export
map.sdm <- function(longlat, layers, error = NULL, year = NULL, idconf = NULL, categorical = NULL, thres = 0, testpercentage = 0, mcp = TRUE, points = FALSE, eval = TRUE, runs = 0, subset = 0){
  
  terra::terraOptions(memmax = 2e+09) # raster::rasterOptions(maxmemory = 2e+09)
  origLonglat = longlat
  
  ##if ensemble is to be done
  if(runs > 0){
    longlat = origLonglat
    
    #if there is spatial error randomly move points within its radius
    if(!is.null(error)){
      for(i in 1:nrow(longlat)){
        #move up to given error (angular movement converted to x and y)
        rndAngle = sample(1:360, 1)
        rndDist = runif(1, 0, error[i])
        longlat[i,1] = longlat[i,1] + rndDist * cos(rndAngle)
        longlat[i,2] = longlat[i,2] + rndDist * sin(rndAngle)
      }
    }
    
    #if there is year
    if(!is.null(year)){
      for(i in 1:nrow(longlat)){
        if(year[i] < sample(min(year):as.integer(substr(Sys.Date(), 1, 4)), 1))
          longlat = longlat[-i,]
      }
    }
    
    #if there is idconf
    if(!is.null(idconf)){
      for(i in 1:nrow(longlat)){
        if(idconf[i] < sample(1:max(idconf), 1))
          longlat = longlat[-i,]
      }
    }
    
    if(eval)
      runEval = matrix(NA, nrow = 1, ncol = 7)
    runMap <- terra::rasterize(as.matrix(longlat), layers[[1]], background = 0) ### removed field = 0,
    pb <- txtProgressBar(min = 0, max = runs, style = 3)
    totalAUC = 0
    for(i in 1:runs){
      if(subset > 0 && subset < dim(layers)[3]){
        runLayers <- layers[[sample.int(dim(layers)[3], subset)]]
        thisRun <- map.sdm(longlat, runLayers, error = NULL, year = NULL, idconf = NULL, categorical, thres, testpercentage, mcp, points, eval, runs = 0, subset = 0)
      } else {
        thisRun <- map.sdm(longlat, layers, error = NULL, year = NULL, idconf = NULL, categorical, thres, testpercentage, mcp, points, eval, runs = 0, subset = 0)
      }
      runAUC = 1
      if(eval){
        runAUC <- thisRun[[2]][1]
        runAUC <- max(0, (runAUC - 0.5)) ^ 2      #weight the map by its AUC above 0.5 to the square
        runEval <- rbind(runEval, thisRun[[2]])
        thisRun <- thisRun[[1]]
      }
      totalAUC = totalAUC + runAUC
      runMap <- runMap + (thisRun * runAUC)
      setTxtProgressBar(pb, i)
    }
    runMap <- terra::app(runMap, function(x) {x/totalAUC}) # raster::calc
    
    ### This should be an aux_function.
    upMap <- terra::classify(
      runMap,
      matrix(c(0, 0.025, 0, 0.025, 1, 1),
             ncol = 3, byrow = TRUE
      )
    )
    consensusMap <- terra::classify(
      runMap,
      matrix(c(0, 0.499, 0, 0.499, 1, 1),
             ncol = 3, byrow = TRUE
      )
    )
    downMap <- terra::classify(
      runMap,
      matrix(c(0, 0.975, 0, 0.975, 1, 1),
             ncol = 3, byrow = TRUE
      )
    )
    
    if(mcp && aoo(consensusMap) >= 4)
      consensusMap <- map.habitat(longlat, consensusMap, mcp = TRUE, eval = FALSE)
    
    if(eval){
      runEval <- runEval[-1,]
      clEval <- matrix(NA, nrow = 3, ncol = 7)
      colnames(clEval) <- c("AUC", "Kappa", "TSS", "EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)")
      rownames(clEval) <- c("UpCL", "Consensus", "LowCL")
      clEval[1,] <- apply(runEval, 2,  quantile, probs= 0.975, na.rm = TRUE)
      clEval[2,] <- apply(runEval, 2,  quantile, probs= 0.5, na.rm = TRUE)
      clEval[3,] <- apply(runEval, 2,  quantile, probs= 0.025, na.rm = TRUE)
      clEval[1:3,4] <- eoo(longlat)
      clEval[1:3,6] <- aoo(longlat)
      clEval[1,5] <- eoo(upMap)
      clEval[1,7] <- aoo(upMap)
      clEval[2,5] <- eoo(consensusMap)
      clEval[2,7] <- aoo(consensusMap)
      clEval[3,5] <- eoo(downMap)
      clEval[3,7] <- aoo(downMap)
      return(list(consensusMap, runMap, clEval))
    } else {
      return (consensusMap)
    }
  }
  
  longlat <- move(longlat, layers)  #move all records falling on NAs
  
  nPoints = min(1000, sum(!is.na(as.vector(layers[[1]])), na.rm = TRUE)/4)
  
  ### dismo has not been updated to allow extraction from terra objects! It is to
  ### be replaced with the package predicts, so you should use predicts::backgroundSample
  ### However, this function is not working as intended. Checking the function it
  ### just calls SpatSample until the number required is achieved and samples leftovers.
  ### bg <- dismo::randomPoints(layers, nPoints)
  ### bg <- predicts::backgroundSample(mask = terra::rast(red.range), n = nPoints)
  
  bg = data.frame(x = c(), y = c())
  while (nrow(bg) < nPoints){
    bg = rbind(bg, terra::spatSample(layers[[1]], nPoints, xy = TRUE, values = FALSE) )  
    vals = terra::extract(layers[[1]], bg, ID = FALSE)
    bg = bg[-c(which(is.na(vals))),]
    
    if (nrow(bg) > nPoints){
      bg = bg[-c(sample(1:nrow(bg), size = (nrow(bg) - nPoints))),]
    }
  }
  
  ##if no categorical variables are given try to figure out which are
  if(is.null(categorical))
    categorical <- find.categorical(layers)
  
  llTrain <- longlat
  llTest <- longlat
  if(testpercentage > 0){
    testRecords <- sample(1:nrow(longlat), ceiling(nrow(longlat)*testpercentage/100))
    llTrain <- longlat[-testRecords,]
    llTest <- longlat[testRecords,]
  }
  mod <- predicts::MaxEnt(layers, llTrain, a = bg, factors = categorical) ##build model
  
  ### p <- raster::predict(mod, layers)                                     
  p <- terra::predict(mod, layers)  ##do prediction
  
  e <- predicts::pa_evaluate(p = llTrain, a = bg, model = mod, x = layers) ##do evaluation of model
  
  ### There are two spec_sens now, max_spenc_sens and equal_sens_spec, check if correct call
  if(thres == 0)
    thres <- predicts::threshold(e)$max_spec_sens ##extract threshold from evaluation
  
  p <- terra::classify( ##convert to presence/absence
    p,
    matrix(c(0, thres, 0, thres, 1, 1),
           ncol = 3, byrow = TRUE
    )
  )
  if(mcp && aoo(p) >= 4)
    p <- map.habitat(longlat, p, mcp = TRUE, eval = FALSE)
  
  if(points)
    p <- max(p, map.points(longlat, p, eval = FALSE))
  
  if(eval){
    e <- predicts::pa_evaluate(p = llTest, a = bg, model = mod, x = layers, tr = thres) ##do evaluation of model with threshold
    auc <- e@stats$auc ### e@auc
    kappa <- e@tr_stats$kappa ### e@kappa
    
    sensitivity <- as.numeric(e@tr_stats$TPR/(e@tr_stats$TPR+e@tr_stats$FNR)) ### e@TPR, e@FNR
    specificity <- as.numeric(e@tr_stats$TNR/(e@tr_stats$TNR+e@tr_stats$FPR)) ### e@TNR, e@FPR
    tss <- sensitivity + specificity - 1
    eooRaw <- eoo(longlat)
    aooRaw <- aoo(longlat)
    aooModel <- aoo(p)
    if(aooModel > 8)
      eooModel <- eoo(p)
    else
      eooModel = aooModel
    txtEval <- matrix(c(auc, kappa, tss, eooRaw, eooModel, aooRaw, aooModel), nrow = 1)
    colnames(txtEval) <- c("AUC", "Kappa", "TSS", "EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)")
    return(list(p, txtEval))
  } else {
    return(p)
  }
}

#' Map species distribution of habitat specialist.
#' @description Mapping of all habitat patches where the species is known to occur.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of each occurrence record.
#' @param layer RasterLayer object representing the presence/absence (1/0) of a single habitat type.
#' @param move If TRUE, identifies and moves presence records to closest cells with suitable habitat. Use when spatial error might put records outside the correct patch.
#' @param mcp If TRUE, all habitat patches inside the minimum convex hull polygon encompassing all occurrence records are converted to presence.
#' @param points If TRUE, force map to include cells with presence records even if suitable habitat was not identified.
#' @param eval If TRUE, build a matrix with EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).
#' @details In many cases a species has a very restricted habitat and we generally know where it occurs. In such cases using the distribution of the known habitat patches may be enough to map the species.
#' @return One raster object and, if eval = TRUE, a matrix with EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).
#' @export
map.habitat <- function(longlat, layer, move = TRUE, mcp = FALSE, points = FALSE, eval = TRUE){
  
  if(points)
    layer <- max(layer, map.points(longlat, layer, eval = FALSE))
  
  if(move){
    moveLayer <- layer
    moveLayer[moveLayer == 0] <- NA
    longlat <- move(longlat, moveLayer)
    remove(moveLayer)
  }
  
  if(mcp){
    vertices <- grDevices::chull(longlat)
    vertices <- c(vertices, vertices[1])
    vertices <- longlat[vertices,]
    poly = sp::Polygon(vertices)
    poly = sp::Polygons(list(poly),1)
    poly = sp::SpatialPolygons(list(poly))    ##minimum convex polygon
    patches <- terra::patches(layer, allowGaps=FALSE)       ##individual patches, numbered
    selPatches <- terra::unique(extract(patches, poly, df = TRUE, weights = TRUE)$patches) ##which patches are inside polygon
  } else {
    patches <- terra::patches(layer, allowGaps=FALSE)       ##individual patches, numbered
    selPatches <- terra::unique(extract(patches, longlat, df = TRUE, weights = TRUE)$patches) ##which patches have the species
  }
  
  selPatches <- selPatches[!is.na(selPatches)]
  allPatches <- terra::unique(patches)
  allPatches <- as.data.frame(cbind(allPatches, rep(0, length(allPatches))))
  colnames(allPatches) <- c("patches", "selected")
  allPatches[selPatches, 2] <- 1
  patches <- terra::classify(patches, allPatches)
  
  layer <- terra::mask(layer, patches, maskvalue = 0, updatevalue = 0)
  if(eval){
    eooRaw <- eoo(longlat)
    eooModel <- eoo(layer)
    aooRaw <- aoo(longlat)
    aooModel <- aoo(layer)
    txtEval <- matrix(c(eooRaw, eooModel, aooRaw, aooModel), nrow = 1)
    colnames(txtEval) <- c("EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)")
    return(list(layer, txtEval))
  } else {
    return(layer)
  }
}

#' Map recorded distribution of species.
#' @description Mapping of all cells where the species is known to occur.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of each occurrence record.
#' @param layers Raster* object as defined by package raster. Any raster with the relevant extent and cell size can be used.
#' @param eval If TRUE, build a matrix with EOO and AOO calculated from occurrence records only.
#' @details To be used if either information on the species is very scarce (and it is not possible to model the species distribution) or, on the contrary, complete (and there is no need to model the distribution).
#' @return One raster object and, if EVAL = TRUE, a matrix with EOO and AOO.
#' @examples
#' records = red.examples("red.records")
#' layers = red.examples("red.layers")
#' terra::plot(map.points(records, layers, eval = FALSE))
#' points(records)
#' @export
map.points <- function(longlat, layers, eval = TRUE){
  if(is(longlat, "data.frame")){ ### Added exception for data.frame
    longlat = as.matrix(longlat)
  }
  p <- terra::rasterize(longlat, layers[[1]], background = 0)
  maskLayer <- sum(layers)
  maskLayer[!is.na(maskLayer)] <- 1
  p <- terra::mask(p, maskLayer)
  if(eval){
    eooRaw <- eoo(longlat)
    aooRaw <- aoo(longlat)
    txtEval <- matrix(c(eooRaw, aooRaw), nrow = 1)
    colnames(txtEval) <- c("EOO", "AOO")
    return(list(p, txtEval))
  } else {
    return(p)
  }
}

#' Species distributions made easy (multiple species).
#' @description Single step for prediction of multiple species distributions. Output of maps (in pdf format), klms (for Google Earth) and relevant data (in csv format).
#' @param longlat data.frame of taxon names, longitude and latitude or eastness and northness (three columns in this order) of each occurrence record.
#' @param layers If NULL analyses are done with environmental layers read from data files of red.setup(). If a Raster* object as defined by package raster, analyses use these.
#' @param habitat Raster* object as defined by package raster. Habitat extent layer (0/1) used instead of layers if any species is an habitat specialist.
#' @param zone UTM zone if data is in metric units. Used only for correct placement of kmls and countries.
#' @param thin boolean defining if species data should be thinned before modeling (only for SDMs).
#' @param error Vector of spatial error in longlat (one element per row of longlat) in the same unit as longlat. Used to move any point randomly within the error radius.
#' @param move If TRUE, identifies and moves presence records to closest cells with environmental data. Use when spatial error might put records outside such data.
#' @param dem RasterLayer object. It should be a digital elevation model for calculation of elevation limits of the species. If NULL, dem from red.setup() is used if possible, otherwise it will be 0.
#' @param pca Number of pca axes for environmental data reduction. If 0 (default) no pca is made.
#' @param filename Name of output csv file with all results. If NULL it is named "Results_All.csv".
#' @param mapoption Vector of values within options: points, habitat and sdm; each value corresponding to the function to be used for each species (map.points, map.habitat, map.sdm). If a single value, all species will be modelled according to it. If NULL, the function will perform analyses using map.points. Species values must be in same order as latlong.
#' @param testpercentage Percentage of records used for testing only. If 0 all records will be used for both training and testing.
#' @param mintest Minimim number of total occurrence records of any species to set aside a test set. Only used if testpercentage > 0.
#' @param points If TRUE, force map to include cells with presence records even if suitable habitat was not identified.
#' @param runs If <= 0 no ensemble modelling is performed. If > 0, ensemble modelling with n runs is made. For each run, a new random sample of occurrence records (if testpercentage > 0), background points and predictive variables (if subset > 0) are chosen. In the ensemble model, each run is weighted as max(0, (runAUC - 0.5)) ^ 2.
#' @param subset Number of predictive variables to be randomly selected from layers for each run if runs > 0. If <= 0 all layers are used on all runs. Using a small number of layers is usually better than using many variables for rare species, with few occurrence records (Lomba et al. 2010, Breiner et al. 2015).
#' @return Outputs maps in asc, pdf and kml format, plus a file with EOO, AOO and a list of countries where the species is predicted to be present if possible to extract.
#' @references Breiner, F.T., Guisan, A., Bergamini, A., Nobis, M.P. (2015) Overcoming limitations of modelling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6: 1210-1218.
#' @references Lomba, A., Pellissier, L., Randin, C.F., Vicente, J., Moreira, F., Honrado, J., Guisan, A. (2010) Overcoming the rare species modelling paradox: a novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation, 143: 2647-2657.
#' @export
map.easy <- function(longlat, layers = NULL, habitat = NULL, zone = NULL, thin = TRUE, error = NULL, move = TRUE, dem = NULL, pca = 0, filename = NULL, mapoption = NULL, testpercentage = 0, mintest = 20, points = FALSE, runs = 0, subset = 0){
  
  try(dev.off(), silent = TRUE)
  spNames <- unique(longlat[,1])
  nSp <- length(spNames)
  
  if(is.null(mapoption))
    mapoption = rep("points", nSp)
  else if(length(mapoption) == 1)
    mapoption = rep(mapoption, nSp)
  else if(length(mapoption) != nSp)
    return(warning("Number of species different from length of mapoption"))
  
  if (all(mapoption == rep("points", nSp))){
    res <- matrix(NA, nrow = nSp, ncol = 5)
    colnames(res) <- c("EOO", "AOO", "Min elevation", "Max elevation", "Countries")
  } else if (("sdm" %in% mapoption) && runs > 0) {
    res <- matrix(NA, nrow = nSp, ncol = 11)
    colnames(res) <- c("EOO (raw)", "EOO (LowCL)", "EOO (Consensus)", "EOO (UpCL)", "AOO (raw)", "AOO (LowCL)", "AOO (Consensus)", "AOO (UpCL)", "Min elevation", "Max elevation", "Countries")
  } else {
    res <- matrix(NA, nrow = nSp, ncol = 7)
    colnames(res) <- c("EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)", "Min elevation", "Max elevation", "Countries")
  }
  rownames(res) <- spNames
  
  if(is.null(layers))
    newLayers <- TRUE
  else
    newLayers <- FALSE
  if(is.null(dem))
    newDem <- TRUE
  else
    newDem <- FALSE
  
  rad = 0.1
  for(s in 1:nSp){
    cat("\nSpecies", s, "of", nSp, "-", toString(spNames[s]),"\n")
    spData <- longlat[longlat[,1] == spNames[s], -1]
    if(!is.null(error)){
      spError <- error[longlat[,1] == spNames[s]]
      if(max(spError) > 1)
        rad <- spError/100000
      else
        rad <- spError
    } else {
      spError <- NULL
    }
    if(newLayers){
      layers <- raster.read(spData)
      if(newDem)
        dem <- layers[[20]]
      if(pca > 0)
        layers <- raster.reduce(layers, n = pca)
    }
    
    if(mapoption[s]  == "sdm" && aoo(move(spData, layers)) > 8){
      if(move)
        spData <- move(spData, layers)
      if(thin)
        spData <- thin(spData)
      if(testpercentage > 0)
        p <- map.sdm(spData, layers, spError, testpercentage = testpercentage, mcp = TRUE, points = points, runs = runs, subset = subset)
      else
        p <- map.sdm(spData, layers, spError, testpercentage = 0, mcp = TRUE, points = points, runs = runs, subset = subset)
    } else if (mapoption[s] == "habitat"){
      p <- map.habitat(spData, habitat, move, points = points)
    } else {
      mapoption[s] = "points"
      p <- map.points(spData, layers)
    }
    terra::writeRaster(p[[1]], paste(toString(spNames[s]), ".asc", sep=""), overwrite = TRUE)
    map.draw(spData, p[[1]], spNames[s], sites = FALSE, print = TRUE)
    if(mapoption[s] != "points"){
      kml(p[[1]], zone = zone, paste(toString(spNames[s]), ".kml", sep=""), mapoption = "aoo")
      countryList <- countries(p[[1]], zone = zone)
      if(is.null(dem))
        elev <- c(0, 0)
      else
        elev <- elevation(p[[1]], dem)
    } else {
      kml(spData, zone = zone, paste(toString(spNames[s]), ".kml", sep=""), mapoption = "points", rad = rad)
      countryList <- countries(spData, zone = zone)
      if(is.null(dem))
        elev <- c(0, 0)
      else
        elev <- elevation(spData, dem)
    }
    
    if(mapoption[s]  == "sdm" && aoo(spData) > 8 && runs > 0){
      terra::writeRaster(p[[2]], paste(toString(spNames[s]), "_prob.asc", sep=""), overwrite = TRUE)
      map.draw(spData, p[[2]], paste(toString(spNames[s]), "_prob", sep = ""), legend = TRUE, print = TRUE)
    }
    
    ##write output values to csv
    spRes = p[[length(p)]]
    if(ncol(res) == 5){      #colnames(res) <- c("EOO", "AOO", "Min elevation", "Max elevation", "Countries")
      res[s,] <- c(spRes, elev, toString(countryList))
    }
    if(ncol(res) == 7){      #colnames(res) <- c("EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)", "Min elevation", "Max elevation", "Countries")
      if(length(spRes) == 7)
        res[s,] <- c(spRes[4:7], elev, toString(countryList))
      else  #if length(spRes) < 7
        res[s,] <- c(spRes[c(1,1,2,2)], elev, toString(countryList))
    }
    if(ncol(res) == 11){     #colnames(res) <- c("EOO (raw)", "EOO (LowCL)", "EOO (Consensus)", "EOO (UpCL)", "AOO (raw)", "AOO (LowCL)", "AOO (Consensus)", "AOO (UpCL)", "Min elevation", "Max elevation", "Countries")
      if(length(spRes) == 2)
        res[s,] <- c(spRes[c(1,1,1,1,2,2,2,2)], elev, toString(countryList))
      else if(length(spRes) == 4)
        res[s,] <- c(spRes[c(1,2,2,2,3,4,4,4)], elev, toString(countryList))
      else if(is.null(dim(spRes)))
        res[s,] <- c(spRes[4:7], elev, toString(countryList))
      else   #if matrix
        res[s,] <- c(spRes[2,4], spRes[3:1,5], spRes[2,6], spRes[3:1,7], elev, toString(countryList))
    }
    write.csv(res[s,], paste(toString(spNames[s]), ".csv", sep = ""))
    if(mapoption[s]  == "sdm" && aoo(spData) > 8){
      if(runs > 0)
        write.csv(p[[3]], paste(toString(spNames[s]), "_detail.csv", sep = ""))
      else
        write.csv(p[[2]], paste(toString(spNames[s]), "_detail.csv", sep = ""))
    }
    
    
  }
  if(is.null(filename))
    write.csv(res, "Results_All.csv")
  else
    write.csv(res, toString(filename))
  return(as.data.frame(res))
}

#' Map creation.
#' @description Creates maps ready to print in pdf or other formats.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of each occurrence record.
#' @param layer RasterLayer object representing the presence/absence map for the species.
#' @param spName String of species name.
#' @param borders If TRUE country borders are drawn.
#' @param scale If TRUE a distance scale in km is drawn.
#' @param legend If TRUE the legend for the map is drawn.
#' @param sites If TRUE the record locations are drawn.
#' @param mcp If TRUE the minimum convex polygon representing the Extent of Occurrence is drawn.
#' @param print If TRUE a pdf is saved instead of the output to the console.
#' @examples
#' records = red.examples("red.records")
#' range = red.examples("red.range")
#' par(mfrow = c(1,2))
#' map.draw(records, layer = range, mcp = TRUE)
#' @export
map.draw <- function(longlat = NULL, layer, spName, borders = FALSE, scale = TRUE, legend = FALSE, sites = TRUE, mcp = FALSE, print = FALSE){
  worldborders = red.examples("worldborders")
  if (borders){
    layer[layer == 0] <- NA
    terra::plot(layer, main = spName, legend = legend, xlab = "longitude", ylab = "latitude", col = "forestgreen")
    lines(worldborders)
  } else {
    terra::plot(layer, main = spName, legend = legend, colNA = "lightblue", xlab = "longitude", ylab = "latitude")
  }
  if (scale){
    width = (terra::xmax(layer) - terra::xmin(layer))
    d = round(width/10^(nchar(width)-1))*10^(nchar(width)-2)
    terra::sbar(d = d, type="bar", divs = 2)
  }
  if (sites && !is.null(longlat))
    points(longlat, pch = 19)
  if (mcp){
    e <- as.data.frame(terra::as.points(layer), geom = "XY")
    e <- e[ e[,1] == 1 , ]
    
    
    vertices <- chull(e[,1], e[,2])
    vertices <- c(vertices, vertices[1])
    vertices <- e[vertices,c(1,2)]
    poly <- sp::SpatialPolygons(list(Polygons(list(Polygon(vertices)),1)))
    terra::plot(poly, add = TRUE)
  }
  if(print){
    dev.copy(device = pdf, file = paste(toString(spName), ".pdf", sep=""))
    dev.off()
  }
}
