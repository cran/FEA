#' @title beamElementMat
#'
#' @description Generates element stiffness matrix for beams.
#'
#' @usage beamElementMat(beamP, beamT, Y, Length, MoI)
#'
#' @param beamP Matrix (2 x n) of beam coordinates.
#' @param beamT Matrix (2 x n) containing the number of the coordinate point as shown in beamP that connect to form a given beam (Discretization table).
#' @param Y Elastic modulus value for material.
#' @param Length Length of beams.
#' @param MoI Moment of inertia for each beam segment.
#'
#' @return Generates initial element matrix needed for the finite element model.
#' \item{beamEmat}{An element matrix of the beam}
#'
#' @examples
#' data(beamGeo)
#' data(beamDime)
#'
#' Length = beamDime$Length
#' MoI = beamDime$MomentofInertia
#'
#' beamEmat = beamElementMat(beamGeo$beamP, beamGeo$beamT, beamGeo$Y, Length, MoI)
#'
#' @export

beamElementMat = function(beamP, beamT, Y, Length, MoI){
  n = NROW(beamP)
  DOF = n

  beamEMfx = function(Y, Length, MoI, DOF){
  Stiffness = (Y[m]*MoI[m]/Length[m]^3)
  coefficients=c(12, 6*Length[m], -12, 6*Length[m],
                 6*Length[m], 4*Length[m]^2, -6*Length[m], 2*Length[m]^2,
                 -12, -6*Length[m], 12, -6*Length[m],
                 6*Length[m], 2*Length[m]^2, -6*Length[m], 4*Length[m]^2)
  ematrix= matrix(coefficients, nrow=DOF, byrow=TRUE) * Stiffness; return(ematrix)}

n = NROW(beamT)
ElementMatrix = list()
  for (m in 1:n){ ElementMatrix[[m]] = beamEMfx(Y, Length, MoI, DOF)}

return("beamEmat" = ElementMatrix) }
