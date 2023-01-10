#' @title beamNodeDis
#'
#' @description Calculates global nodal displacements of beam.
#'
#' @usage beamNodeDis(beamP, BCtran, BCrot, REM, NodeKnownL, ForceV)
#'
#' @param beamP Matrix (2 x n) of beam coordinates.
#' @param BCtran Boundary constraint for nodes to translate in x or y directions.
#' @param BCrot Boundary constraint for nodes to rotate.
#' @param REM Reduced element matrix, returned from function ReducedEM.
#' @param NodeKnownL data frame with constraint parameters applied to each node in the x and y directions. Formatted for use in reduced element matrix. Generated from ApplyBC function.
#' @param ForceV Reduced force vector matrix containing the model load parameters. Returned from function ForceVector.
#'
#' @return Produces tables with new node coordinates that are produced by the geometry under an applied load.
#' \item{NodeDis}{Nodal displacement}
#' \item{GlobalND}{Nodal displacement in the global environment}
#' \item{GlobalNDMatrix}{Nodal displacement in the global environment as a reduced matrix}
#'
#' @examples
#' data(beamGeo)
#' data(beamFV)
#' data(beamREM)
#' data(beamBC)
#'
#' ForceV = beamFV
#' REM = beamREM
#' NodeKnownL = beamBC
#'
#' beamND = beamNodeDis(beamGeo$beamP, beamGeo$BCtran, beamGeo$BCrot, REM, NodeKnownL, ForceV)
#'
#' @export

beamNodeDis = function(beamP, BCtran, BCrot, REM, NodeKnownL, ForceV){
  n = NROW(beamP) * 2

  BC = matrix(rbind(BCtran, BCrot), byrow = TRUE, ncol = 1)

  UKNodeDisplace = MASS::ginv(REM) %*% ForceV

  GlobalND = matrix(rep(0,n),byrow=T)
  GlobalND[NodeKnownL] = UKNodeDisplace
  GlobalND[is.na(GlobalND)] <- 0

  GNDmat = matrix(GlobalND, ncol = 2, byrow= TRUE)

  Rlist = list("BEAMDisplacement" = UKNodeDisplace, "GlobalND" = GlobalND, "GlobalNDMatrix" = GNDmat)}
