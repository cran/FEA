#' @title beamDimensions
#'
#' @description Calculates input dimensions needed for beam finite element.
#'
#' @usage beamDimensions(Y, G, Nu, beamP, beamT, thick,  fx, fy)
#'
#' @param Y Elastic modulus value for material (Pa).
#' @param G Shear modulus value for material (Pa). If using Euler-Bernoulli model, G = 0.
#' @param Nu Poisson's ratio value for material.
#' @param beamP Matrix (2 x n) of beam coordinates.
#' @param beamT Matrix (2 x n) containing the number of the coordinate point as shown in beamP that connect to form a given beam (Discretization table).
#' @param thick Thickness of the beam
#' @param fx Load value (newtons) in the x direction.
#' @param fy Load value (newtons) in the y direction.
#'
#' @return Calculates values needed for both Timoshenko-Ehrenfest and Euler-Bernoulli beam theories.
#' \item{k}{Timoshenko-Ehrenfest correction}
#' \item{Length}{Beam length}
#' \item{Angle}{Beam angle within the plane}
#' \item{MomentofInertia}{Moment of Inertia for each beam}
#' \item{Displacement}{Displacement under Timoshenko-Ehernfest beam theory}
#' \item{RotationAngle}{Angle of rotation}
#' \item{StiffnessAngle}{Stiffness angle}
#'
#' @examples
#' data(beamGeo)
#'
#' DOF = 4
#' n = NROW(beamGeo$beamT)
#' thick = matrix(c(0.039149, 0.03, 0.0246625), ncol = 1, nrow = n) #height(thickness) of beam
#'
#' beamDime = beamDimensions(beamGeo$Y, beamGeo$G, beamGeo$Nu, beamGeo$beamP, beamGeo$beamT,
#'                           beamGeo$thick, beamGeo$fx, beamGeo$fy)
#'
#' @export

beamDimensions = function(Y, G, Nu, beamP, beamT, thick,  fx, fy){
  k = function(v){
    TBk = (10*(1+Nu[m]))/(12+ (11*Nu[m]))}

  L <- function(beamT, beamP){
    Length = sqrt(((beamP[beamT[m,2], 2] - beamP[beamT[m,1], 2])^2) + ((beamP[beamT[m,2], 1] - beamP[beamT[m,1], 1])^2))}

  matH <- function(beamT, beamP){
    matH = ((180/pi) * atan((beamP[beamT[m,2], 2] - beamP[beamT[m,1], 2]) / (beamP[beamT[m,2], 1] - beamP[beamT[m,1], 1])))}

  MomentI = function(Length, thick){
    MoI = (Length[m]*(thick[m]^3)) * (1/12)}

  Displace = function(Y, G, MoI, Length, thick, f){
    dis = (1+((48*Y[m]*MoI[m])/(5*(thick[m]*Length[m])*G[m]*Length[m]^2)))* ((5*f[m]*(Length[m]^4))/(38*Y[m]*MoI[m]))}

  RotAng = function(Y, G, MoI, Length, thick){
    phi = (12*Y[m]*MoI[m])/(G[m]*(thick[m]*Length[m])*Length[m]^2)}

  StifAng = function(Y, G, MoI, Length, thick, TBk){
    ta = ((12*Y[m]*MoI[m])/(TBk[m]*G[m]*(thick[m]*Length[m])*(Length[m]^2)))}

  f = rowSums(fx, fy)
  n = NROW(beamT)
  TBk = numeric(n)        #TB correction
  Length = numeric(n)     #Beam length
  mH = numeric(n)         #Beam initial angle
  MoI = numeric(n)        #moment of inertia for individual rectangular beams
  w = numeric(n)          #displacement
  phi = numeric(n)        #Rotation
  theta = numeric(n)      #stiffness angle
  for (m in 1:n) {
    TBk[m] = k(Nu)
    Length[m] = L(beamT, beamP)
    mH[m] = matH(beamT, beamP)
    MoI[m] = MomentI(Length, thick)
    w[m] = Displace(Y, G, MoI, Length, thick, f)
    phi[m] = RotAng(Y, G, MoI, Length, thick)
    theta[m] = StifAng(Y, G, MoI, Length, thick, TBk)}

  Rlist = list("k (TBcorrection)" = TBk, "Length" = Length, "Angle" = mH, "MomentofInertia" = MoI, "Displacement" = w, "RotationAngle" = phi, "StiffnessAngle" = theta)
}
