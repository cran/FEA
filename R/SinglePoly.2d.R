#' @title SinglePoly.2d
#'
#' @description Generates a mesh for polygon with a single continuous geometry
#'
#' @usage SinglePoly.2d(x, y, ptDS, ptDL)
#'
#' @param x X-coordinates for geometry.
#' @param y Y-coordinates for geometry.
#' @param ptDS Density of points desired within the geometry.
#' @param ptDL Density of points desired at the perimeter of the geometry.
#'
#' @return Coordinate points of nodes distributed within and on the line of a given geometry.
#' \item{AllCoords}{all coordinate points distributed across the geometry.}
#' \item{Within}{all coordinate points within the geometry ONLY.}
#' \item{Line}{all coordinate points that lay on the perimeter of the geometry ONLY.}
#'
#' @examples
#' data(Cart)
#'
#' x = Cart[,1]
#' y= Cart[,2]
#' ptDS = 30
#' ptDL = 20
#'
#' polyshape = SinglePoly.2d(x, y, ptDS, ptDL)
#'
#'@export

SinglePoly.2d = function(x, y, ptDS, ptDL){
  poly <- sp::Polygon(matrix(rbind(x, y),  nrow = NROW(x), ncol =2, byrow = T))
  pts <- as.data.frame(sp::spsample(poly, n= ptDS, "regular"), pch = 3) #change n to reflect desired point density
  names(pts)[1] <- "x"
  names(pts)[2] <- "y"
  line1 = as.data.frame(sp::spsample(sp::SpatialLines(list(sp::Lines(sp::Line(sp::SpatialPoints(poly)), ID="a"))), n= ptDL, offset =c(0,1), "regular"), pch = 3)
  names(line1)[1] <- "x"
  names(line1)[2] <- "y"
  coords= rbind(pts,line1)
  plot(coords, pch = 46, type = "p") #polygon coordinates

  Rlist = list("AllCoords" = coords, "Within" = pts, "Line" = line1)}
