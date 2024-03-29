#' @title AutoAdjust.2d
#'
#' @description Allows for automatic refinement of the triangular mesh generated based on given parameters. Will remove elements that are outside the margin of the geometry.
#'
#' @usage AutoAdjust.2d(meshP, meshT, edge, centroid, AspectR, AR)
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param edge Coordinate points of the initial geometry.
#' @param centroid Matrix (2 x n) of triangle elements.
#' @param AspectR Aspect ratio of each triangle element.
#' @param AR maximum desired aspect ratio, numeric value.
#'
#' @importFrom graphics text
#' @importFrom stats na.omit
#'
#' @return Generates new mesh and centroid tables
#' \item{Meshpts}{Includes both new mesh coordinate points and triangulation of points.}
#' \item{Centroids}{Centroid positions for each triangle element.}
#'
#' @examples
#' data(triMesh)
#' data(polyshape)
#' data(dime)
#'
#' meshP = triMesh$MeshPts$p
#' meshT = triMesh$MeshPts$T
#' edge =  polyshape$Line
#' centroid = triMesh$Centroids
#' AspectR = dime$AspectRatio
#' AR = 10
#'
#' auto = AutoAdjust.2d(meshP, meshT, edge, centroid, AspectR, AR)
#'
#' @export

AutoAdjust.2d= function(meshP, meshT, edge, centroid, AspectR, AR){
  #Is the point outside the edge? (in the space and need to be removed if = 1)
  edge = as.matrix(edge)
  outO = ptinpoly::pip2d(edge, centroid)
  v = NROW(centroid)
  total1 = numeric(v)
  total2 = numeric(v)
  for (m in 1:v){
    if (outO[m] <= 0) {total1[m] = 1}
    else {total1[m] = NA}}

  #clean up high aspect ratios
  for (m in 1:v)
  {if (AspectR[m] <= AR) {total2[m] = 1}
    else {total2[m] = NA}}

  total = total1*total2
  pd1 = meshT*total
  pd2 = na.omit(pd1)   #Remove triangles outside shape
  Meshpt = geometry::trimesh(pd2, meshP)     #Trimesh with inner regions removed
  Pmid = matrix(c(((meshP[pd2[,1],1]+ meshP[pd2[,2],1] + meshP[pd2[,3],1])/3),
                  (meshP[pd2[,1],2]+ meshP[pd2[,2],2] + meshP[pd2[,3],2])/3), ncol = 2) #triangle centroids
  data1 = c(1:NROW(Meshpt$p))
  data2 = c(1:NROW(Meshpt$T))
  text(Meshpt$p[,1],Meshpt$p[,2], data1, col = "red", cex= 0.6, font = 2)
  text(Pmid[,1],Pmid[,2], data2, col = "blue", cex= 0.6, font = 1)

  Rlist = list("MeshPts" = Meshpt, "Centroids" = Pmid)}
