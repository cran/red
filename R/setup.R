#' Setup GIS directory.
#' @description Setup directory where GIS files are stored.
#' @param gisPath Path to the directory where the gis files are stored.
#' @details Writes a txt file in the red directory allowing the package to always access the world GIS files directory.
#' @export
red.setDir <- function(gisPath = NULL){
  if (is.null(gisPath)) {
    gisPath <- readline("Input directory for storing world gis layers:")
  }
  gisPath <- paste(gisPath, "/", sep = "")
  redFile <- paste(find.package("red"), "/red.txt", sep = "")
  dput(gisPath, redFile)
}

#' Read GIS directory.
#' @description Read directory where GIS files are stored.
#' @details Reads a txt file pointing to where the world GIS files are stored.
#' @export
red.getDir <- function(){
  redFile <- paste(find.package("red"), "/red.txt", sep = "")
  if (file.exists(redFile)) { # if there is already a file read from it
    dir <- dget(redFile)
  } else {
    warning(paste(redFile, "not found, please run red.setDir()"))
    return()
  }
  return(dir)
}

#' Download and setup GIS files.
#' @description Setup red to work with species distribution modelling and layers available online.
#' @details Please check that you have at least 50Gb free in your disk (and a fast internet connection) to download all files. In the end of the process "only" 17.4Gb will be left though. This function will:
#' 1. Check if maxent.jar is available in the dismo package directory.
#' 2. Ask user input for GIS directory.
#' 3. Download global bioclim and elevation files (20) from http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_30s_bio.zip.
#' 4. Download landcover files (12) from http://data.earthenv.org/consensus_landcover/without_DISCover/.
#' 5. Unzip all files and delete the originals.
#' 6. Create a new layer (1) with the dominant land cover at each cell.
#' 7. Resample all files (33) to approximately 10x10km (for use with widespread species) grid cells.
#' Sit back and enjoy, this should take a while.
#' @export
red.setup <- function(){
  worldclim_refs = list(url = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip",
                        zip_name = "bioclim2.zip",
                        layer_prefix = "wc2.1_30s_bio_",
                        altitude_url = "http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/alt_30s_bil.zip",
                        altitude_zip = "alt_30s_bil.zip"
  ) 
  
  wopt = list()
  
  oldwd = getwd()
  on.exit(expr = setwd(oldwd))
  gisdir = red.setDir()
  setwd(gisdir)
  
  ##basic setup
  pb <- txtProgressBar(min = 0, max = 33, style = 3)
  
  ##download and process bioclim
  download.file(worldclim_refs$url, worldclim_refs$zip_name)
  ### This is prone to timeout. We should maybe have a comment saying you can choose
  ### to download the file yourself in your browser and place it in the user's red
  ### directory
  
  unzip(zipfile = worldclim_refs$zip_name) # was bioclim.zip
  file.remove(worldclim_refs$zip_name)  # bioclim.zip
  
  for (i in 1:19) {
    ### Bioclim no longer uses double digit notation, e.g.: 01, 02 so code for
    ### that isn't needed. In the future, instead of making if() print() we
    ### can just use this: paste0(rep("0", c(2 - nchar(i))), i).
    setTxtProgressBar(pb, i)
    rast_name <- paste(worldclim_refs$layer_prefix,
                       i,
                       ".tif",
                       sep = ""
    )
    
    rast <- terra::rast(rast_name)
    rast <- terra::crop(rast, c(-180, 180, -56, 90),
                        filename = paste0("red_1km_", i, ".tif"),
                        datatype = "FLT4S", filetype = "GTiff",
                        gdal = c("COMPRESS=LZW"), overwrite = FALSE
    )
    rast <- terra::aggregate(rast, 10,
                             filename = paste("red_10km_", i, ".tif", sep = ""),
                             datatype = "FLT4S", filetype = "GTiff",
                             gdal = c("COMPRESS=LZW"), overwrite = FALSE
    )
    file.remove(rast_name)
  }
  
  ## download and process altitude
  setTxtProgressBar(pb, 20)
  download.file(worldclim_refs$altitude_url, worldclim_refs$altitude_zip)
  unzip(zipfile = worldclim_refs$altitude_zip)
  file.remove(worldclim_refs$altitude_zip)
  
  rast <- terra::rast("alt.bil")
  rast <- terra::crop(rast, c(-180, 180, -56, 90), filename = "red_1km_20.tif")
  rast <- terra::aggregate(rast, 10, filename = "red_10km_20.tif")
  
  file.remove("alt.bil")
  file.remove("alt.hdr")
  gc()
  
  ##download and process land cover
  altmask1 = terra::rast("red_1km_20.tif")
  altmask10 =  terra::rast("red_10km_20.tif")
  for (i in 1:12) {
    setTxtProgressBar(pb, (i + 20))
    ### should probably have bigger timeout period
    download.file(
      paste0(
        "http://data.earthenv.org/consensus_landcover/",
        "without_DISCover/Consensus_reduced_class_", i, ".tif"
      ),
      destfile = paste0("Consensus_reduced_class_", i, ".tif"),
      mode = "wb"
    )
    rast <- terra::rast(paste("Consensus_reduced_class_", i, ".tif", sep = ""))
    rast <- terra::mask(rast, altmask1,
                        filename = paste("red_1km_", (i + 20), ".tif", sep = "")
    )
    rast <- terra::aggregate(rast, 10) ### maskLayer <- sum(altmask, rast)
    rast <- terra::mask(rast, altmask10, filename = paste("red_10km_", (i + 20),
                                                          ".tif",
                                                          sep = ""
    ))
    file.remove(paste("Consensus_reduced_class_", i, ".tif", sep = ""))
    gc()
  }
  remove(rast)
  
  ##create new rasters with most common landcover at each cell
  setTxtProgressBar(pb, 33)
  max1 <- terra::rast()
  max10 <- terra::rast()
  for(i in 21:32){
    rast <- terra::rast(paste("red_1km_", i, ".tif", sep=""))
    max1 <- c(max1, rast)
    rast <- terra::rast(paste("red_10km_", i, ".tif", sep=""))
    max10 <- c(max10, rast)
  }
  
  max1 <- which.max(max1)
  terra::writeRaster(max1, "red_1km_33.tif")
  max10 <- which.max(max10)
  terra::writeRaster(max10, "red_10km_33.tif")
  remove(max1, max10)
  gc()
  setwd(oldwd)
}