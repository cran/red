#' Red List Index.
#' @description Calculates the Red List Index (RLI) for a group of species.
#' @param spData Either a vector with species assessment categories for a single point in time or a matrix with two points in time in different columns (species x date). Values can be text (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).
#' @param tree An hclust or phylo object (used when species are weighted by their unique contribution to phylogenetic or functional diversity).
#' @param boot If TRUE bootstrapping for statistical significance is performed on both values per date and the trend between dates.
#' @param dd bootstrap among all species (FALSE) or Data Deficient species only (TRUE).
#' @param runs Number of runs for bootstrapping
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to Extinct/Extinct in the Wild (5).
#' Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
#' Each species weight can further be influenced by how much it uniquely contributes to the phylogenetic or functional diversity of the group (Cardoso et al. in prep.).
#' To incorporate Importantly, the RLI is based on true improvements or deteriorations in the status of species, i.e. genuine changes. It excludes category changes resulting from, e.g., new knowledge (Butchart et al. 2007).
#' The RLI approach helps to develop a better understanding of which taxa, regions or ecosystems are declining or improving.
#' Juslen et al. (2016a, b) suggested the use of bootstrapping to search for statistical significance when comparing taxa or for trends in time of the index and this approach is here implemented.
#' @return Either a vector (if no two dates are given) or a matrix with the RLI values and, if bootstrap is performed, their confidence limits and significance.
#' @references Butchart, S.H.M., Stattersfield, A.J., Bennun, L.A., Shutes, S.M., Akcakaya, H.R., Baillie, J.E.M., Stuart, S.N., Hilton-Taylor, C. & Mace, G.M. (2004) Measuring global trends in the status of biodiversity: Red List Indices for birds. PloS Biology, 2: 2294-2304.
#' @references Butchart, S.H.M., Akcakaya, H.R., Chanson, J., Baillie, J.E.M., Collen, B., Quader, S., Turner, W.R., Amin, R., Stuart, S.N. & Hilton-Taylor, C. (2007) Improvements to the Red List index. PloS One, 2: e140.
#' @references Juslen, A., Cardoso, P., Kullberg, J., Saari, S. & Kaila, L. (2016a) Trends of extinction risk for Lepidoptera in Finland: the first national Red List Index of butterflies and moths. Insect Conservation and Diversity, 9: 118-123.
#' @references Juslen, A., Pykala, J., Kuusela, S., Kaila, L., Kullberg, J., Mattila, J., Muona, J., Saari, S. & Cardoso, P. (2016b) Application of the Red List Index as an indicator of habitat change. Biodiversity and Conservation, 25: 569-585.
#' @examples
#' rliData <- matrix(c("LC","LC","EN","EN","EX","EX","LC","CR","DD","DD"), ncol = 2, byrow = TRUE)
#' colnames(rliData) <- c("2000", "2010")
#' rli(rliData[,1])
#' rli(rliData[,1], boot = TRUE)
#' rli(rliData)
#' rli(rliData, boot = TRUE, dd = TRUE)
#' @export
rli <- function (spData, tree = NULL, boot = FALSE, dd = FALSE, runs = 1000){
  
  ##if only one point in time is given
  if(is.null(dim(spData)))
    return(rli.calc(spData, tree, boot, dd, runs))  ##return either 1 or 3 values
  
  ##if two points in time are given
  ts <- apply(spData, 2, function(x) rli.calc(x, tree, boot = FALSE))
  sl <- (ts[2] - ts[1]) / (as.numeric(colnames(spData))[2] - as.numeric(colnames(spData))[1])
  if(!boot){
    res <- matrix(c(ts, sl), nrow = 1)
    colnames(res) <- c(colnames(spData), "Change/year")
    rownames(res) <- c("Raw")
    return(res)
  } else {
    tr <- apply(spData, 2, function(x) rli.calc(x, tree, boot, dd, runs))
    p = 0
    rndSl = rep(NA, runs)
    for(r in 1:runs){
      rndSl[r] <- rli.calc(spData[,2], tree, boot, dd, runs = 1)[2] - rli.calc(spData[,1], tree, boot, dd, runs = 1)[2]
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
#' @param spData A matrix with group names (first column) and species assessment categories for one or two points in time (remaining columns). Values can be text (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).
#' @param tree A list of hclust or phylo objects, each corresponding to a tree per group (used when species are weighted by their unique contribution to phylogenetic or functional diversity).
#' @param boot If TRUE bootstrapping for statistical significance is performed on both values per date and the trend between dates.
#' @param dd bootstrap among all species (FALSE) or Data Deficient species only (TRUE).
#' @param runs Number of runs for bootstrapping
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to 5 (Extinct/Extinct in the Wild).
#' Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
#' Each species weight can further be influenced by how much it uniquely contributes to the phylogenetic or functional diversity of the group (Cardoso et al. in prep.).
#' Importantly, the RLI is based on true improvements or deteriorations in the status of species, i.e. genuine changes. It excludes category changes resulting from, e.g., new knowledge (Butchart et al. 2007).
#' The RLI approach helps to develop a better understanding of which taxa, regions or ecosystems are declining or improving.
#' Juslen et al. (2016a, b) suggested the use of bootstrapping to search for statistical significance when comparing taxa or for trends in time of the index and this approach is here implemented.
#' @return A matrix with the RLI values and, if bootstrap is performed, their confidence limits and significance.
#' @references Butchart, S.H.M., Stattersfield, A.J., Bennun, L.A., Shutes, S.M., Akcakaya, H.R., Baillie, J.E.M., Stuart, S.N., Hilton-Taylor, C. & Mace, G.M. (2004) Measuring global trends in the status of biodiversity: Red List Indices for birds. PloS Biology, 2: 2294-2304.
#' @references Butchart, S.H.M., Akcakaya, H.R., Chanson, J., Baillie, J.E.M., Collen, B., Quader, S., Turner, W.R., Amin, R., Stuart, S.N. & Hilton-Taylor, C. (2007) Improvements to the Red List index. PloS One, 2: e140.
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
rli.multi <- function (spData, tree = NULL, boot = FALSE, dd = FALSE, runs = 1000){
  
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
    if(is.null(tree))
      v <- rli(spData[spData[,1] == groups[g],-1], tree = NULL, boot = boot, dd = dd, runs = runs)
    else
      v <- rli(spData[spData[,1] == groups[g],-1], tree[[g]], boot = boot, dd = dd, runs = runs)
    if(ncol(res) < 13){
      res[g,] <- v
      colnames(res) <- colnames(v)
    } else {
      res[g,1:4] <- v$Values[,1]
      res[g,5:8] <- v$Values[,2]
      res[g,9:12] <- v$Values[,3]
      res[g,13] <- v$P_change
    }
  }
  return(res)
}

#' Prediction of Red List Index.
#' @description Linearly interpolates and extrapolates RLI values to any years.
#' @param rliValue Should be a vector with RLI values and names as the corresponding year numbers.
#' @param from Starting year of the sequence to predict.
#' @param to Ending year of the sequence to predict.
#' @param rliPlot Plots the result
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' @return A matrix with the RLI values and confidence limits.
#' @examples rliValue <- c(4.5, 4.3, 4.4, 4.2, 4.0)
#' names(rliValue) <- c(2000, 2004, 2008, 2011, 2017)
#' rli.predict(rliValue, 1990, 2020)
#' @export
rli.predict <- function(rliValue, from = NA, to = NA, rliPlot = FALSE){
  year = as.numeric(c(names(rliValue)))
  rliTable = data.frame(rliValue, year)
  if(is.na(from))
    from = min(year)
  if(is.na(to))
    to = max(year)
  newYear = data.frame(year = seq(from = from, to = to, by = 1))
  lmOut = predict(lm(rliValue ~ year, data = rliTable), newYear, interval = "confidence", level = 0.95)
  res = lmOut[,c(2,1,3)]
  colnames(res) = c("LowCL", "Fitted RLI", "UpCL")
  rownames(res) = newYear$year
  
  if(rliPlot){
    plot(year, rliValue, xlab="Year", ylab="Fitted RLI", xlim = c(from, to), ylim = c(0,5))
    abline(lm(rliValue ~ year, data = rliTable), col = "red")
    matlines(newYear, lmOut[,2:3], col = "blue", lty = 2)
  }
  
  return(res)
}

#' Sampled Red List Index.
#' @description Calculates accumulation curve of confidence limits in sampled RLI.
#' @param spData A vector with species assessment categories for a single point in time. Values can be text (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).
#' @param tree An hclust or phylo object (used when species are weighted by their unique contribution to phylogenetic or functional diversity).
#' @param p p-value of confidence limits (in a two-tailed test).
#' @param runs Number of runs for smoothing accumulation curves.
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to Extinct/Extinct in the Wild (5).
#' Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
#' Yet, in many groups, it is not possible to assess all species due to huge diversity and/or lack of resources. In such case, the RLI is estimated from a randomly selected sample of species - the Sampled Red List Index (SRLI; Stuart et al. 2010).
#' This function allows to calculate how many species are needed to reach a given maximum error of the SRLI around the true value of the RLI (with all species included) for future assessments of the group.
#' @return A vector with the accumulation of the error of the SRLI around the true value of the RLI (with all species included).
#' @references Butchart, S.H.M., Stattersfield, A.J., Bennun, L.A., Shutes, S.M., Akcakaya, H.R., Baillie, J.E.M., Stuart, S.N., Hilton-Taylor, C. & Mace, G.M. (2004) Measuring global trends in the status of biodiversity: Red List Indices for birds. PLoS Biology, 2: 2294-2304.
#' @references Butchart, S.H.M., Akcakaya, H.R., Chanson, J., Baillie, J.E.M., Collen, B., Quader, S., Turner, W.R., Amin, R., Stuart, S.N. & Hilton-Taylor, C. (2007) Improvements to the Red List index. PLoS One, 2: e140.
#' @references Stuart, S.N., Wilson, E.O., McNeely, J.A., Mittermeier, R.A. & Rodriguez, J.P. (2010) The barometer of Life. Science 328, 117.
#' @examples rliData <- c("LC","LC","EN","EN","EX","EX","LC","CR","CR","EX")
#' rli.sampled(rliData)
#' @export
rli.sampled <- function (spData, tree = NULL, p = 0.05, runs = 1000){
  
  nSpp <- length(spData)
  accum <- rep(NA, nSpp)
  for(n in 1:nSpp){               #test with n species from the entire set
    diff = rep(NA, runs)          #try runs times each species
    for(r in 1:runs){             #do r runs for each n species
      rndComm = rep(NA, nSpp)
      rndSpp = sample(nSpp, n)
      rndComm[rndSpp] = spData[rndSpp]
      diff[r] = abs(rli.calc(spData, tree, FALSE, FALSE, runs = 1) - rli.calc(rndComm, tree, FALSE, FALSE, runs = 1))  #calculate absolute difference between true and sampled rli for each run
    }
    accum[n] = quantile(diff, (1-p))
  }
  return(accum)   #returns the accumulation curve of confidence limit of sampled RLI
}

#' Mapping the Red List Index.
#' @description Creates a map for the red list index according to species distribution and threat status.
#' @param spData Either a vector with species assessment categories for a single point in time or a matrix with two points in time in different columns (species x date). Values can be text (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).
#' @param layers Species distributions (0/1), a Raster* object as defined by package raster.
#' @param layers2 Species distributions (0/1) on the second point in time, a Raster* object as defined by package raster. If there are two dates but no layers2, the distributions are assumed to be kept constant in time.
#' @param tree An hclust or phylo object (used when species are weighted by their unique contribution to phylogenetic or functional diversity).
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to Extinct/Extinct in the Wild (5).
#' Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
#' Each species weight can further be influenced by how much it uniquely contributes to the phylogenetic or functional diversity of the group (Cardoso et al. in prep.).
#' @return A RasterLayer with point values  (if a single date is given) or change per cell (if two dates are given).
#' @references Butchart, S.H.M., Stattersfield, A.J., Bennun, L.A., Shutes, S.M., Akcakaya, H.R., Baillie, J.E.M., Stuart, S.N., Hilton-Taylor, C. & Mace, G.M. (2004) Measuring global trends in the status of biodiversity: Red List Indices for birds. PloS Biology, 2: 2294-2304.
#' @references Butchart, S.H.M., Akcakaya, H.R., Chanson, J., Baillie, J.E.M., Collen, B., Quader, S., Turner, W.R., Amin, R., Stuart, S.N. & Hilton-Taylor, C. (2007) Improvements to the Red List index. PloS One, 2: e140.
#' @examples sp1 <- terra::rast(matrix(c(1,1,1,0,0,0,0,0,NA), ncol = 3))
#' sp2 <- terra::rast(matrix(c(1,0,0,1,0,0,1,0,NA), ncol = 3))
#' sp3 <- terra::rast(matrix(c(1,0,0,0,0,0,0,0,NA), ncol = 3))
#' sp4 <- terra::rast(matrix(c(0,1,1,1,1,1,1,1,NA), ncol = 3))
#' layers <- c(sp1, sp2, sp3, sp4)
#' spData <- c("CR","EN","VU","LC")
#' terra::plot(rli.map(spData, layers))
#' @export
rli.map <- function (spData, layers, layers2 = NULL, tree = NULL){
  
  if(!is.null(dim(spData))){             #if to calculate change call this same function twice
    if(is.null(layers2)){
      layers2 <- layers
    }
    map1 <- rli.map(spData[,1], layers = layers, tree = tree)
    map2 <- rli.map(spData[,2], layers = layers2, tree = tree)
    return(map2 - map1)
  }
  
  #convert rasters to array
  layers = terra::as.array(layers)
  
  #get data for each cell (row by row)
  cells = matrix(NA, (nrow(layers) * ncol(layers)), dim(layers)[3])
  i = 0
  for (r in 1:nrow(layers)){
    for(c in 1:ncol(layers)){
      i = i+1
      cells[i,] = layers[r,c,]
    }
  }
  
  #RLI of each cell
  rliCells = rep(NA, nrow(cells))
  for (i in 1:nrow(cells)){
    rliNA <- ifelse(cells[i,] == 1, spData, NA) #only consider species present in each cell
    rliCells[i] = rli.calc(rliNA, tree = tree)
  }
  
  #create RLI map
  rliMap = terra::rast(matrix(rliCells, nrow = nrow(layers), byrow = T))
  return(rliMap)
}