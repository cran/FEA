#' @title beamApplyBC
#'
#' @description Boundary constraint for element centroids based on coordinate points. For the x & y direction per centroid create matrix with boundary 1(unfixed) or 0(fixed).
#'
#' @usage beamApplyBC(beamP, BCtran, BCrot)
#'
#' @param beamP Matrix (2 x n) of beam coordinates.
#' @param BCtran Boundary constraint for nodes to translate in x or y directions. Input as a non-matrix column.
#' @param BCrot Boundary constraint for nodes to rotate. Input as a non-matrix column.
#'
#' @return A data frame with constraint parameters applied to each node for directional translation and rotation. Formatted for use in reduced element matrix.
#' \item{NodeKnownL}{Matrix (1 x n) of constraint parameters}
#'
#' @examples
#' data(beamGeo)
#'
#' beamBC = beamApplyBC(beamGeo$beamP, beamGeo$BCtran, beamGeo$BCrot)
#'
#' @export


beamApplyBC = function(beamP, BCtran, BCrot){
n = 2*NROW(beamP) # m=col, n=row, z=element#

  ApplyBC = function(BCtran, BCrot){
  BoundContr = matrix(rbind(BCtran, BCrot), byrow = TRUE, ncol = 1)
  Apply_BC = (1:n) * BoundContr }

NodeKnownL = ApplyBC(BCtran, BCrot)

return(NodeKnownL)}
