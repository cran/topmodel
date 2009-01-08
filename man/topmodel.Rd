\name{topmodel}
\alias{topmodel}
\title{implementation of the hydrological model TOPMODEL}
\description{
  This is an implementation of the 1995 Fortran version of TOPMODEL by Keith Beven.
}
\usage{
topmodel(parameters, topidx, delay, rain,ET0,verbose = F, Qobs=NA)
}
\arguments{
  \item{parameters}{A vector or a matrix containing the input parameters (see below for the exact structure)}
  \item{topidx}{A 2 column matrix with respectively the topographic index classes and values (see below for the exact structure)}
  \item{delay}{Delay function for overland flow (see below)}
  \item{rain}{A vector of rain data (m per time unit dt)}
  \item{ET0}{A vector of ET0 data (m per time unit dt)}
  \item{verbose}{If set to TRUE, returns besides predicted discharge also overland flow, base flow and storage. The default is FALSE}
  \item{Qobs}{Shortcut: If Qobs is given, normal output is suppressed and only the Nash and Sutcliffe efficiency is returned (m per time unit dt)}
}
\details{
topmodel() automatically implements a Monte Carlo simulation. If the parameter argument is a single vector, the model is run once. If the parameter argument is a matrix, each row should represent a parameter set. In that case, the model is run with each parameter set (see the examples below).

A single parameter set consists of: \code{c(qs0,lnTe,m,Sr0,SrMax,td,vch,vr,k0,psi,dtheta,dt)}, with:

\tabular{ll}{
qs0 \tab Initial subsurface flow per unit area [m]\cr
lnTe 	\tab log of the areal average of T0 [m2/h]\cr
m	\tab Model parameter controlling the rate of decline of transmissivity in the soil profile, see Beven, 1984\cr
Sr0 	\tab Initial root zone storage deficit [m]\cr
Srmax 	\tab Maximum root zone storage deficit [m]\cr
td 	\tab Unsaturated zone time delay per unit storage deficit\cr
vch 	\tab channel flow outside the catchment catchment [m/h]\cr
vr 	\tab channel flow inside catchment [m/h]\cr
k0 	\tab Surface hydraulic conductivity [m/h]\cr
CD 	\tab capillary drive, see Morel-Seytoux and Khanji (1974)\cr
dt 	\tab The timestep [h]\cr
}

The topidx dataframe can be derived conveniently with \code{make.classes()}. It should contain 2 columns. The first column should give the lower boundary of each topographic index class, and the second column should give the respective area fraction. The second column must sum to 1.

Parameters k0 is only used in the unsaturated zone time delay funtion if td<0. psi and dtheta are currently not used, because the infiltration excess module is not yet functional.

Overland flow is routed through a delay function which represents the time spent in the channel system. The parameter \code{delay} is used for this. Delay is a matrix with 2 columns. The first column gives the cumulative relative area. The second column gives the average distance towards the outlet (m).

}
\value{
The function returns an array of observed discharges. If more than one parameter set is given, a matrix is returned, with each column representing a discharge set coinciding with the parameter sets. If Qobs is given, the function returns an array of Nash-Sutcliffe efficiencies, 1 for each parameter sets.

If verbose output is requested, a list is returned, with the modelled discharge (Q), overland flow (qo), subsurface flow (qs), storage (S), infiltration excess overland flow (fex), and actual evapotranspiration (Ea) for each time step.

Be aware that invoking \code{topmodel()} without Q for a large number of runs, may require a large amount of memory.
}
\references{Beven, K. J., Kirkby, M. J., 1979. A physically based variable contributing area model of basin hydrology. Hydrol. Sci. Bull. 24, 43-69.

Beven K, Lamb R, Quinn P, Romanowicz R, Freer J, 1995. TOPMODEL. In: Sing VP (Ed), Computer Models of Watershed Hydrology. Water Resources Publications, Colorado. pp. 627-668.

Morel-Seytoux, H.J., Khanji, J., 1974. Derivation of an Equation of Infiltration. Water Resources Research, 10, 795-800.
Beven, K., 1984. Infiltration into a Class of Vertically Non-Uniform Soils. Hydrological Sciences Journal 29, 425-434.

See also \url{http://source.ggy.bris.ac.uk/wiki/Hydrology_in_R} for a more examples on how to run topmodel in R.
}
\author{Wouter Buytaert, University of Bristol}
\note{}

% \seealso{\code{\link{print}}}
\examples{
data(huagrahuma)
attach(huagrahuma)

## returns the simulated runoff (Qobs not given)
Qsim <- topmodel(parameters,topidx,delay,rain,ET0)

## returns a list of simulated runoff (Qobs), overland flow (qo), subsurface flow (qs) and storage (S):
topmodel(parameters,topidx,delay,rain,ET0,verbose = TRUE)

## plot observed and simulated discharge:
plot(Qobs)
points(Qsim, col="red", type="l")

## For a Monte carlo run with random sampling, we construct a parameter matrix:

runs<-10

qs0   <- runif(runs)*4e-5
lnTe  <- runif(runs)*3-2
m     <- runif(runs)*0.2
Sr0   <- runif(runs)*0.02
Srmax <- runif(runs)*2
td    <- runif(runs)*3-3
vch   <- 1000
vr    <- 100+runif(runs)*2400
k0    <- runif(runs)*0.01
CD    <- runif(runs)*5
dt    <- 0.25

parameters<-cbind(qs0,lnTe,m,Sr0,Srmax,td,vch,vr,k0,CD,dt)

## returns an array of 10 Nash Sutcliffe efficiencies; one for each parameter set:
result<-topmodel(parameters,topidx,delay,rain,ET0,Qobs = Qobs)

}
\keyword{models}
