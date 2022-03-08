#' @title NodeDis
#'
#' @description Calculates global nodal displacements
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param REM Reduced element matrix, returned from function ReducedEM.
#' @param ForceV Reduced force vector matrix containing the model load parameters. Returned from function ForceVector.
#' @param NodeKnownL data frame with constraint parameters applied to each node in the x and y directions. Formatted for use in reduced element matrix. Generated from ApplyBC function.
#'
#' @return Produces tables with new node coordinates that are produced by the geometry under an applied load.
#' \item{NodeDis}{Nodal displacement}
#' \item{GlobalND}{Nodal displacement in the global environment}
#'
#' @examples
#' data(triMesh)
#' data(load)
#' data(reduc_EM)
#' data(bound)
#'
#' meshP = triMesh$MeshPts$p
#' REM = reduc_EM
#' ForceV = load
#' NodeKnownL = bound
#'
#' displacN = NodeDis(meshP, REM, ForceV, NodeKnownL)
#'
#' @import MASS
#' @export

#Global Nodal Displacement
NodeDis = function(meshP, REM, ForceV, NodeKnownL){
  m= n= z= o= NROW(meshP)
  TDOF = 2*z

  UKNodeDisplace = MASS::ginv(REM) %*% ForceV #unknown node displacements where forces are known

  GlobalND = matrix(rep(0,TDOF),byrow=T)
  GlobalND[NodeKnownL] = UKNodeDisplace

  GlobalND[is.na(GlobalND)] <- 0

  Rlist = list("NodalDisplacement" = UKNodeDisplace, "GlobalND" = GlobalND)}
