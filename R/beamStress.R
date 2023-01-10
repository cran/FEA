#' @title beamStress
#'
#' @description Calculates local stress and strain for beam elements
#'
#' @usage beamStress(beamP, beamT, Y, Length, MoI, RotAng, BND)
#'
#' @param beamP Matrix (2 x n) of beam coordinates.
#' @param beamT Matrix (2 x n) containing the number of the coordinate point as shown in beamP that connect to form a given beam (Discretization table).
#' @param Y Value of Young's (Elastic) modulus
#' @param Length Length of beam
#' @param MoI Moment of Inertia
#' @param RotAng Angle of rotation
#' @param BND Global nodal displacement matrix, return from function beamNodeDis
#'
#' @return Completes FEM by calculating values of stress and strain, produces three (3) [3 x n] matrix.
#' \item{BendingStress}{Bending Stress}
#'
#' @examples
#' data(beamGeo)
#' data(beamGLforce)
#'
#' Length = beamDime$Length
#' MoI = beamDime$MomentofInertia
#' RotAng = beamDime$Angle
#' BND = beamND
#'
#' beamBendStress = beamStress(beamGeo$beamP, beamGeo$beamT, beamGeo$Y, Length, MoI, RotAng, BND)
#'
#' @export

beamStress = function(beamP, beamT, Y, Length, MoI, RotAng, BND){
  DOF = NROW(beamP)

  BendingStress = function(beamT, Y, Length, RotAng, GNDmat) {
    GND = matrix(c(GNDmat[beamT[m,1],], GNDmat[beamT[m,2],]), ncol = 1, nrow = DOF)
    Tmatrix= matrix(c(-cos(RotAng[m]), -sin(RotAng[m]), cos(RotAng[m]), sin(RotAng[m])), nrow=1, byrow=T);
    local_stress= (Y[m]/Length[m])*Tmatrix %*% GND; return(local_stress)}

  n = NROW(beamT)
  GNDmat = BND$GlobalNDMatrix
  Bend = numeric(n)
  for (m in 1:n){
    Bend[m] = BendingStress(beamT, Y, Length, RotAng, GNDmat)}
  return(Bend)}
