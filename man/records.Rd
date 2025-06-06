% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{records}
\alias{records}
\title{Download taxon records from GBIF.}
\usage{
records(taxon)
}
\arguments{
\item{taxon}{Taxon name.}
}
\value{
A data.frame with longitude and latitude, plus species names if taxon is above species.
}
\description{
Downloads species or higher taxon data from GBIF and outputs non-duplicate records with geographical coordinates.
}
\details{
As always when using data from multiple sources the user should be careful and check if records "make sense". This can be done by either ploting them in a map (e.g. using red::map.draw()) or using red::outliers().
}
\examples{
\dontrun{
rec = records("Nephila senegalensis")
plot(rec)
}
}
