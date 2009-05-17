\name{topidx}
\alias{topidx}
\title{Calculation of the topographic index from a DEM raster}
\description{Calculation of the topographic index from a DEM raster}
\usage{
topidx(DEM, resolution, river = NA)
}
\arguments{
  \item{DEM}{A matrix representing a DEM with equally sized pixels and equal NS and EW resolution}
  \item{resolution}{Resolution of the DEM (m)}
  \item{river}{A matrix representing a river map (optional). All cells with values higher than 0 will be treated as river cells.}
}
\details{The river map should not contain negative values}
\value{The function returns a list, containing 2 rasters. Atb are the topographic index values of each grid cell. Area contains the contributing area}
%\references{}
\author{Wouter Buytaert, Imperial College London, based on routines developed by the hydrology group of Lancaster University}
%\note{}

\seealso{\code{\link{topmodel}}}
\examples{
  data(huagrahuma.dem)
  topidx <- topidx(huagrahuma.dem, resolution= 25)$atb
  image(topidx)
}
\keyword{models}
