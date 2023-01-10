#' @title EulerBeamFEA
#'
#' @description Calculates stress in beam structures using the Euler-Bernoulli beam theory.
#'
#' @usage EulerBeamFEA(Y, beamP, beamT, fx, fy, BCtran, BCrot, Length, MoI, RotAng)
#'
#' @param Y Elastic modulus value for material (Pa).
#' @param beamP Matrix (2 x n) of beam coordinates.
#' @param beamT Matrix (2 x n) containing the number of the coordinate point as shown in beamP that connect to form a given beam (Discretization table).
#' @param fx Load value (newtons) in the x direction.
#' @param fy Load value (newtons) in the y direction.
#' @param BCtran Boundary constraint for nodes to translate in x or y directions. Input as a non-matrix column.
#' @param BCrot Boundary constraint for nodes to rotate. Input as a non-matrix column.
#' @param Length Length of beam.
#' @param MoI Moment of inertia for each beam segment.
#' @param RotAng Angle of rotation
#'
#' @return Calculates local forces and stresses, as well as bending stress for beams following the Euler-Bernoulli beam theory.
#' \item{Stress}{Local stress at node}
#' \item{LocalLoad}{Local load at node}
#' \item{BendingStress}{Bending Stress}
#'
#' @examples
#' data(beamGeo)
#' data(beamDime)
#'
#' Length = beamDime$Length
#' MoI = beamDime$MomentofInertia
#' RotAng = beamDime$Angle
#'
#' beamFEA = EulerBeamFEA(beamGeo$Y, beamGeo$beamP, beamGeo$beamT, beamGeo$fx, beamGeo$fy,
#'                        beamGeo$BCtran, beamGeo$BCrot, Length, MoI, RotAng)
#'
#' @export

EulerBeamFEA = function(Y, beamP, beamT, fx, fy, BCtran, BCrot, Length, MoI, RotAng){

  #Element Matrix
  p2 = beamElementMat(beamP, beamT, Y, Length, MoI)

  #Expanded Element Matrix
  ElementMat = p2
  p3 = beamExpandEM(beamP, beamT, ElementMat)

  #Global element matrix (reducedEM)
  beamExEM = p3
  p4 = beamGlobalEM(beamExEM)

  #Apply boundary conditions
  p5 = beamApplyBC(beamP, BCtran, BCrot)

  #Reduced element matrix
  GMat = p4
  NodeKnownL = p5

  p6 = beamReducedEM(GMat, NodeKnownL)

  #Load Distribution
  p7 = beamUDL(beamP, beamT, Length, fx, fy)

  #Force vector
  p8 = beamForceVector(beamP, fx, fy, NodeKnownL)

  #Node displacement
  ForceV = p8
  REM = p6

  p9 = beamNodeDis(beamP, BCtran, BCrot, REM, NodeKnownL, ForceV)

  #Global Local forces
  BUDL = p7
  BND = p9
  p10 = beamGLForces(beamP, beamT, Y, MoI, Length, GMat, BUDL, BND)

  #Stress
  p11 = beamStress(beamP, beamT, Y, Length, MoI, RotAng, BND)

Rlist = list("GlobalForces"= p10$Gforce, "LocalForces" = p10$Lforce, "BendingStress" = p11)}
