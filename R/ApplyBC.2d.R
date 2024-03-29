#' @title ApplyBC.2d
#'
#' @description Boundary constraint for element centroids based on coordinate points. For the x & y direction per centroid create matrix with boundary 1(unfixed) or 0(fixed).
#'
#' @usage ApplyBC.2d(meshP, BoundConx, BoundCony)
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh
#' @param BoundConx Boundary constraint for nodes in the x-direction
#' @param BoundCony Boundary constraint for nodes in the y-direction
#'
#' @return A data frame with constraint parameters applied to each node in the x and y directions. Formatted for use in reduced element matrix.
#' \item{NodeKnownL}{Constraint parameters}
#'
#' @examples
#' data(triMesh)
#'
#' meshP = triMesh$MeshPts$p
#' BoundConx = BoundCony = numeric(NROW(meshP))
#' BoundConx[1:NROW(meshP)] = BoundCony[1:NROW(meshP)] = 1
#' BoundConx[c(10, 11, 12)] = BoundCony[c(10, 11, 12)] = 0
#'
#' bound = ApplyBC.2d(meshP, BoundConx, BoundCony)
#'
#' @export

ApplyBC.2d = function(meshP, BoundConx, BoundCony){
  m= o=  n = 2*NROW(meshP) # m=col, n=row, z=element#

  BoundConxy = matrix(rbind(BoundConx, BoundCony), byrow = TRUE, ncol = 1)
  Apply_BC = function(BoundConxy, n){
    for (m in 1:n){if (BoundConxy[m] < 0){BoundConxy[m] = 0} else {next}}
    NodeKnownL = (1:n)*BoundConxy} #ONLY if loads are already listed as in the example

  NodeKnownL = Apply_BC(BoundConxy, n)}
