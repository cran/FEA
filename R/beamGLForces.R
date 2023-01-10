#' @title beamGLForces
#'
#' @description Uses nodal displacements to determine global and local forces at each node
#'
#' @usage beamGLForces(beamP, beamT, Y, MoI, Length, GMat, BUDL, BND)
#'
#' @param beamP Matrix (2 x n) of beam coordinates.
#' @param beamT Matrix (2 x n) containing the number of the coordinate point as shown in beamP that connect to form a given beam (Discretization table).
#' @param Y Elastic Modulus of material
#' @param MoI Moment of Inertia
#' @param Length Length of beam
#' @param GMat Global stiffness matrix
#' @param BUDL Column matrix for beam distributed load
#' @param BND beam nodal displacement, output from function "beamNodeDis"
#'
#' @return Matrices of global and local forces. (Results in kN)
#' \item{GForce}{Large global force matrix.}
#' \item{Lforce}{Large local force matrix.}
#'
#' @examples
#' data(beamGeo)
#' data(beamDime)
#' data(beamsUDL)
#' data(beamND)
#' data(beamGloMat)
#'
#' Length = beamDime$Length
#' MoI = beamDime$MomentofInertia
#' BUDL = beamsUDL
#' BND = beamND
#' GMat = beamGloMat
#'
#' GLforce = beamGLForces(beamGeo$beamP, beamGeo$beamT, beamGeo$Y, MoI, Length, GMat, BUDL, BND)
#'
#' @export

beamGLForces = function(beamP, beamT, Y, MoI, Length, GMat, BUDL, BND){
  DOF = NROW(beamP)
  TDOF = 2 * NROW(beamP)

# Global force moment
  GND = BND$GlobalND
  GFMoment = function(GMat, GND){
    GForce = GMat %*% GND }
  GForce = GFMoment(GMat, GND)

# Local force
  LFMoment = function(beamT, Y, MoI, Length, GNDmat, LoadDis, DOF, TDOF){
    GND = matrix(c(GNDmat[beamT[m,1],], GNDmat[beamT[m,2],]), ncol = 1, nrow = NROW(beamP))
    Bendingstiff = Y[m]*MoI[m]/Length[m]^3
    coefficients=c(12, 6*Length[m], -12, 6*Length[m],
                   6*Length[m], 4*Length[m]^2, -6*Length[m], 2*Length[m]^2,
                   -12, -6*Length[m], 12, -6*Length[m],
                   6*Length[m], 2*Length[m]^2, -6*Length[m], 4*Length[m]^2)
    ematrix= matrix(coefficients, nrow=DOF, byrow=T) * Bendingstiff
    localL = (ematrix %*% GND) - LoadDis[[m]]
    return(localL)}

  n = NROW(beamT)
  GNDmat = BND$GlobalNDMatrix
  LoadDis = BUDL$DLMatrix
  LocalL = list()
  for (m in 1:n){
    LocalL[[m]] = LFMoment(beamT, Y, MoI, Length, GNDmat, LoadDis, DOF, TDOF)}

  Rlist = list( "Gforce" = GForce, "Lforce" = LocalL)}
