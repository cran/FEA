#' @title GlobalMat.2d
#'
#' @description Generates global stiffness matrix - once established, the expanded element matrix must be combined to create the global structural stiffness matrix by adding the expanded matrices.
#'
#' @usage GlobalMat.2d(meshP, meshT, ExEM)
#'
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param ExEM Expanded element matrix
#'
#' @return Produces large (n x n) global matrix
#' \item{GlobalMat}{Global matrix}
#'
#' @examples
#' data(triMesh)
#' data(fea_ExEM)
#'
#' meshP = triMesh$MeshPts$p
#' meshT = triMesh$MeshPts$T
#' ExEM = fea_ExEM
#'
#' gloMat = GlobalMat.2d(meshP, meshT, ExEM)
#'
#' @export

GlobalMat.2d = function(meshP, meshT, ExEM){
  m= n= o= 2*NROW(meshP)
  z= NROW(meshT) # m=col, n=row, z=element#
  GlobalMat = matrix(0, nrow = o, ncol = o)
  for (m in 1:o){for (n in 1:o){
    cell = numeric(z)
    for (z in 1:z){
      cell[z] = ExEM[[z]][m, n]}
    GlobalMat[m,n] = sum(cell)}}
  return(GlobalMat)}
