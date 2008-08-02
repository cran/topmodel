topidx <- function(DEM, resolution) {

if(!is(DEM, "matrix"))
	print("DEM should be a matrix")

else {

#   ew_res <- DEM@grid@cellsize[1]
#   ns_res <- DEM@grid@cellsize[2]
   cols <- dim(DEM)[2]
   rows <- dim(DEM)[1]

# check for values below -9000

   if(length(DEM[DEM[!is.na(DEM)]< -8999]) != 0)
	print("Check map for irrealistic values")

# replace NA's

   DEM[is.na(DEM)] = -9999

   topidxmap <- .C("topidx",
		PACKAGE = "topmodel",
		as.double(t(DEM)),
		as.integer(rows),
		as.integer(cols),
		as.double(resolution),
		as.double(resolution),
		topidxmap = double(rows*cols))$topidxmap

   topidxmap[topidxmap == -9999]<-NA
   
   topidxmap <- matrix(topidxmap, nrow = rows, ncol = cols, byrow=T)

   return(topidxmap)

   }
}
