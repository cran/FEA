#' @title ForceVector.2d
#'
#' @description Creates a matrix of loads in the x & y direction for each load unconstrained node.
#'
#' @usage ForceVector.2d(Fx, Fy, RSF, meshP, NodeKnownL)
#'
#' @param Fx Load vector for the x-direction
#' @param Fy Load vector for the y-direction
#' @param RSF If surface traction is present assign value as the ReducedSF matrix; if there is no surface traction set RSF = 0
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param NodeKnownL data frame with constraint parameters applied to each node in the x and y directions. Formatted for use in reduced element matrix. Generated from ApplyBC function.
#'
#' @return Produces a matrix with loading parameters for each node.
#' \item{ReducedFV}{Reduced force vector matrix containing the model load parameters.}
#'
#' @examples
#' data(triMesh)
#' data(reduc_SF)
#' data(bound)
#'
#' meshP = triMesh$MeshPts$p
#' RSF = reduc_SF
#' Fx = 10
#' Fy = 10
#' NodeKnownL = bound
#'
#' load = ForceVector.2d(Fx, Fy, RSF, meshP, NodeKnownL)
#'
#' @export

ForceVector.2d = function(Fx, Fy, RSF, meshP, NodeKnownL){
  z= NROW(meshP) # m=col, n=row, z=element#
  Vec = matrix(cbind(Fx, Fy), nrow = z, ncol =2, byrow = FALSE)
  F.vector = matrix((rbind(Vec[,1], Vec[,2])), ncol = 1, byrow = TRUE) + RSF
  ReducedFV = matrix(F.vector[c(NodeKnownL)], ncol = 1)}
