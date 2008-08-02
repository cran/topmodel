topidx2 <- function(map, resolution) {

if(!is(map, "matrix"))
	print("map should be a matrix")

else {

#   ew_res <- map@grid@cellsize[1]
#   ns_res <- map@grid@cellsize[2]
   cols <- dim(map)[2]
   rows <- dim(map)[1]

# check for values below -9000

   if(length(map[map[!is.na(map)]< -8999]) != 0)
	print("Check map for irrealistic values")

# replace NA's

   map[is.na(map)] = -9999

# NOTE: as.double fills the vector per matrix COLUMN
# we want rows, so we use t(matrix)

   topidxmap <- .C("topidx2",
		PACKAGE = "topmodel",
		as.double(t(map)),
		as.integer(rows),
		as.integer(cols),
		as.double(resolution),
		as.double(resolution),
		topidxmap = double(rows*cols))$topidxmap

   topidxmap[topidxmap == -9999]<-NA
   
   topidxmap <- matrix(t(topidxmap), nrow = rows, ncol = cols, byrow = T)

   return(topidxmap)

   }
}
