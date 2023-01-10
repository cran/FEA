#' @title FEMStress.2d
#'
#' @description Creates a complete finite element model using stress for a given 2D mesh under specified boundary conditions (constrain and load).
#'
#' @usage FEMStress.2d(meshP, meshT, centroid, BoundConx, BoundCony, SFShear,
#' SFTensile, Length, area, Fx, Fy,  Y, Nu, Thick)
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param centroid Matrix (2 x n) containing coordinate points of the centroid of each triangular element.
#' @param BoundConx Boundary constraint for nodes in the x-direction
#' @param BoundCony Boundary constraint for nodes in the y-direction
#' @param SFTensile Magnitude of tensile surface traction; if there is no surface traction then SFTensile = 0
#' @param SFShear Magnitude of positive shear traction; if there is no surface traction then SFShear = 0
#' @param Length Truss length
#' @param Thick Triangle element thickness
#' @param area Triangle element area
#' @param Fx Load vector for the x-direction
#' @param Fy Load vector for the y-direction
#' @param Y Value of Young's (Elastic) modulus
#' @param Nu Value of Poisson's ratio
#' @param Thick Value of the thickness of the mesh, a value must be given.
#'
#' @return Completes the FEM to generate values of stress and strain and nodal displacement.
#' \item{NodeDisplacement}{Node displacement on each axis}
#' \item{LocalStress}{Stress as calucated from stress, strain, and stress from strain. Three (3) [3 x n] matrices where [x, y, tau]}
#'
#' @examples
#' \donttest{
#' data(triMesh)
#' data(dime)
#'
#' meshP = triMesh$MeshPts$p
#' meshT = triMesh$MeshPts$T
#' centroid = triMesh$Centroids
#' Y = matrix(20e9, nrow = NROW(meshT))
#' Nu = matrix(0.45, nrow = NROW(meshT))
#' Thick = 0.001
#' DOF = 6
#' BoundConx = BoundCony = numeric(NROW(meshP))
#' BoundConx[1:NROW(meshP)] = BoundCony[1:NROW(meshP)] = 1
#' BoundConx[c(10, 11, 12)] = BoundCony[c(10, 11, 12)] = 0
#' SFShear = 0
#' SFTensile = 0
#' Length = dime$TrussLength
#' area = dime$Area
#' Fx = 10
#' Fy = 10
#'
#' fea_stress = FEMStress.2d(meshP, meshT, centroid, BoundConx, BoundCony, SFShear, SFTensile,
#'              Length, area, Fx, Fy, Y, Nu, Thick)
#' }
#'
#' @export

#Complete FEM
FEMStress.2d = function(meshP, meshT, centroid, BoundConx, BoundCony, SFShear,
                             SFTensile, Length, area, Fx, Fy, Y, Nu, Thick){
  FEAt1 = ElementMat.2d(meshP, meshT, Nu, Y, Thick)

  #Expanded element matrix represents the contribution of individual finite elements towards the global structural matrix
  EMatrixlist = FEAt1$EMPStress
  FEAt2 = ExpandEM.2d(meshP, meshT, centroid, EMatrixlist)

  #Global stiffness matrix - once established expanded must be combined to create the global structural stiffness matrix (by adding the expanded matrices)
  ExEM = FEAt2 #expanded element matrix
  FEAt3 = GlobalMat.2d(meshP, meshT, ExEM)

  #Boundary conditions for element centroids based on coordinate points. For the x & y direction per centroid create matrix with boundary 1(unfixed) or 0(fixed)
  FEAt4 = ApplyBC.2d(meshP, BoundConx, BoundCony)

  #Reduced stiffness matrix - use boundary condition to reduce matrix to smaller form by removing systems that are bound.
  GMat = FEAt3
  NodeKnownL = FEAt4
  FEAt5 = ReducedEM.2d(GMat, NodeKnownL)

  #Element Surface Traction - generates the column matric for the uniformly distributed, SFTensile and SFShear must be given values. load.
  FEAt6 = SurfaceTraction.2d(meshP,SFTensile, SFShear, Length, Thick, area)

  #Expanded Surface Force
  SurfTrac = FEAt6 #surface traction
  FEAt7 = ExpandSFT.2d(meshP, meshT, SurfTrac)

  #Reduced Surface force
  ExSurf = FEAt7
  FEAt8 = ReducedSF.2d(meshP, ExSurf)

  #Reduced Load vector - values must be known or given in the problem statement and applied to the model.
  RSF = FEAt8 #reduced surface traction. If none is present RSF = 0
  NodeKnownL = FEAt4
  FEAt9 = ForceVector.2d(Fx, Fy, RSF, meshP, NodeKnownL)

  #Global Nodal Displacement
  ForceV = FEAt9 #reduced force vector
  REM = FEAt5 #reduced element martix
  NodeKnownL = FEAt4
  FEAt10 = NodeDis.2d(meshP, REM, ForceV, NodeKnownL)

  #Global and Local Forces
  GMat = FEAt3
  GlobalND = FEAt10$GlobalND
  EMatrixlist = FEAt1$EMPStress
  FEAt11 = GLForces.2d(meshP, meshT, GMat, GlobalND, EMatrixlist)

  #Stresses
  GlobalND = FEAt10$GlobalND
  FEAt12 = LocalStress.2d(meshP, meshT, Y, Nu, GlobalND)

  Rlist = list("NodeDisplacement" = FEAt10, "Global/Local_Force" = FEAt11, "LocalStress" = FEAt12)
}
