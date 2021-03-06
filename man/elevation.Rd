% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/red.R
\name{elevation}
\alias{elevation}
\title{Elevation limits.}
\usage{
elevation(spData, dem = NULL)
}
\arguments{
\item{spData}{One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (0/1 values).}

\item{dem}{RasterLayer object. Should be a digital elevation model (DEM) of the relevant area. If not given the function will try to read it from base data, only works with longlat data.}
}
\value{
A vector with two values (min and max) in meters above (or below) sea level.
}
\description{
Calculates the elevation (or depth) limits (range) of a species based on either known records or predicted distribution.
}
\details{
Maximum and minimum elevation are calculated based on the DEM.
}
\examples{
data(red.records)
data(red.range)
data(red.layers)
dem = red.layers[[3]]
elevation(red.records, dem)
elevation(red.range, dem)
}
