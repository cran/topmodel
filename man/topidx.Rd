\name{topidx}
\alias{topidx}
\title{Calculation of the topographic index from a DEM raster}
\description{Calculation of the topographic index from a DEM raster}
\usage{
topidx(DEM, resolution)
}
\arguments{
  \item{DEM}{A matrix representing a DEM with equally sized pixels and equal NS and EW resolution}
  \item{resolution}{Resolution of the DEM (m) }
}
\details{The function returns a raster with ln(a/tanb) values for each grid cell}
%\value{}
%\references{}
\author{Wouter Buytaert, University of Bristol}
%\note{}

\seealso{\code{\link{topmodel}}}
\examples{}
\keyword{models}
