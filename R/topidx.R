topidx <- function(DEM, resolution, river=NA) {

  ## data preparation

  nrow <- dim(DEM)[1]
  ncol <- dim(DEM)[2]

  ## data checking

  if(!is(DEM, "matrix"))
    stop("The DEM should be a matrix")
  if(min(as.vector(DEM)) < -9000)
     stop("DEM contains unrealistic values (< -9000)")
  if(!is.na(river) && !is(DEM, "matrix"))
     stop("The river should be a matrix")
  
  if(is(river, "matrix") && (min(river) < 0))
     stop("Error: the object 'river' should only contain positive values")
  else river = rep(0,nrow*ncol)

  ## calling the function

  result <- .C("topidx",
               PACKAGE = "topmodel",
               as.double(DEM),
               as.integer(river),
               as.integer(nrow),
               as.integer(ncol),
               as.double(resolution),
               as.double(resolution),
               result = double(length(DEM)*2))$result

  ## formatting of the results

  atb  <- matrix(result[1:(nrow*ncol)],nrow=nrow)
  area <- matrix(result[(nrow*ncol+1):(nrow*ncol*2)],nrow=nrow)

  return(list(atb = atb,area = area))

}
