#####required packages
library("BAT")
library("dismo")
library("gdistance")
library("geosphere")
library("graphics")
library("grDevices")
library("jsonlite")
library("methods")
library("predicts")
library("sp")
library("stats")
library("terra")
library("utils")
#' @import gdistance
#' @import graphics
#' @import jsonlite
#' @import sp
#' @import stats
#' @import utils
#' @importFrom BAT contribution
#' @importFrom dismo gbif
#' @importFrom geosphere areaPolygon
#' @importFrom grDevices chull dev.copy dev.off pdf
#' @importFrom methods slot as is
#' @importFrom predicts MaxEnt pa_evaluate threshold
#' @importFrom terra cellSize global patches crop ext extract values layerCor mask rast rasterize as.points as.polygons as.array classify res spatSample sbar terrain trim writeRaster xmax xmin ymax ymin distance nrow terraOptions minmax project simplifyGeom intersect writeVector predict merge
NULL
#> NULL 

############################################################################
##################################DATASETS##################################
############################################################################

#' Example data packaged with *red*
#' @description Load data included in the package. This includes *red.records*,
#' a matrix of longitude and latitude (two columns) occurrence records for
#' Hogna maderiana (Walckenaer, 1837); *red.range*, a SpatRaster object, as
#' defined by package terra, of the geographic range of Hogna maderiana
#' (Walckenaer, 1837); *red.layers*, a SpatRaster object with layers 
#' representing the average annual temperature, total annual precipitation,
#' altitude and landcover for Madeira Island
#' (Fick & Hijmans 2017, Tuanmu & Jetz 2014); and *worldborders* is a small vector
#' of global country borders.
#' @param data Name of data in quotes. E.g.: `"red.records"`
#' If `NULL`, the example files will be listed.
#' @examples
#' red.examples()
#' red.examples("red.range")
#' @source This function is inspired by `palmerpanguins::path_to_file()`
#' which in turn is based on `readxl::readxl_example()`.
#' @export
red.examples <- function(data = NULL) {
  # worldborders is to possibly be replaced by RedList_countries (v. 2022.1)
  # by Victor Cazalis in a later version (https://github.com/victorcazalis/RedList_countries).
  if (is.null(data)) {
    print(
      c(
        "red.records", "red.range", "red.layers", "worldborders" 
      )
    )
    return(NULL)
  } else {
    if(data == "red.records"){
      path = system.file(paste0("extdata/red.records.csv"), package = "red")
      out = read.csv(path)
    } else if (data == "red.range") {
      path = system.file(paste0("extdata/red.range.tif"), package = "red")
      out = terra::rast(x = path)
    } else if (data == "red.layers") {
      path = system.file(paste0("extdata/red.layers.", c(1:4), ".tif"), package = "red")
      out = terra::rast(x = path)
    } else if (data == "worldborders") {
      path = system.file(paste0("extdata/worldborders"), package = "red")
      out = terra::vect(x = path)
    }
  }
  return(out)
}
