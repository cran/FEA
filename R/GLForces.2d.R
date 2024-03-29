#' @title GLForces.2d
#'
#' @description Uses nodal displacements to determine global and local forces at each node
#'
#' @usage GLForces.2d(meshP, meshT, GMat, GlobalND, EMatrixlist)
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param GMat Global matrix
#' @param GlobalND Global nodal displacement
#' @param EMatrixlist Element matrix list
#'
#' @return Matrices of global and local forces
#' \item{GForce}{Large global force matrix.}
#' \item{Lforce}{Large local force matrix.}
#'
#' @examples
#' data(triMesh)
#' data(gloMat)
#' data(displacN)
#' data(fea_EM)
#'
#' meshP = triMesh$MeshPts$p
#' meshT = triMesh$MeshPts$T
#' GMat = gloMat
#' GlobalND = displacN$GlobalND
#' EMatrixlist = fea_EM$EMPStress
#'
#' glfor = GLForces.2d(meshP, meshT, GMat, GlobalND, EMatrixlist)
#'
#' @export

#Global and local forces
GLForces.2d = function(meshP, meshT, GMat, GlobalND, EMatrixlist){
  z= NROW(meshT) # m=col, n=row, z=element#

  G_Force = function(GMat, GlobalND){
    GMat %*% GlobalND}
  GlobalForce = G_Force(GMat, GlobalND)

  L_Force = function(EMatrixlist, meshT, GlobalND){
    r1=2*meshT[m,1]-1
    r2=2*meshT[m,1]
    r3=2*meshT[m,2]-1
    r4=2*meshT[m,2]
    r5=2*meshT[m,3]-1
    r6=2*meshT[m,3]
    EMatrixlist[[m]] * GlobalND[c(r1, r2, r3, r4, r5, r6)] }
  LocalF = list()
  for (m in 1:z){
    LocalF[[m]] = L_Force(EMatrixlist, meshT, GlobalND)}

  Rlist = list("GForce" = GlobalForce, "Lforce" = LocalF)}
