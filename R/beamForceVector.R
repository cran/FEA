#' @title beamForceVector
#'
#' @description Creates a matrix of loads for beams in the x & y direction for each load unconstrained node.
#'
#' @usage beamForceVector(beamP, fx, fy, NodeKnownL)
#'
#' @param beamP Matrix (2 x n) of beam coordinates.
#' @param fx Load vector (newtons) in the x-direction.
#' @param fy Load vector (newtons) in the y-direction.
#' @param NodeKnownL Data frame with constraint parameters applied to each node in the x and y directions. Formatted for use in reduced element matrix. Generated from ApplyBC function.
#'
#' @return Produces a matrix with loading parameters for each node.
#' \item{ReducedFV}{Reduced force vector matrix containing the model load parameters.}
#'
#' @examples
#' data(beamGeo)
#' data(beamUDL)
#'
#' NodeKnownL = beamBC
#' FV = beamForceVector(beamGeo$beamP, beamGeo$fx, beamGeo$fy, NodeKnownL)
#'
#' @export

beamForceVector = function(beamP, fx, fy, NodeKnownL){
  z= 2*NROW(beamP) # m=col, n=row, z=# of beam elements
  f = matrix(cbind(fx, fy), ncol = 2)
  Vec = matrix((rbind(f[,1], f[,2])), ncol = 1)
  ReducedFV = matrix(Vec[c(NodeKnownL)], ncol = 1)}
