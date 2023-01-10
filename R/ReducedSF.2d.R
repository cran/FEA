#' @title ReducedSF.2d
#'
#' @description Reduced matrix of surface forces
#'
#' @usage ReducedSF.2d(meshP, ExSurf)
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param ExSurf Expanded surface matrix, output from ExpandSFT
#'
#' @return Produces a large matrix.
#' \item{RSF}{Produces a large, reduced surface force matrix}
#'
#' @examples
#' data(triMesh)
#' data(expSurf)
#' meshP = triMesh$MeshPts$p
#' ExSurf = expSurf
#' reduc_SF = ReducedSF.2d(meshP, ExSurf)
#'
#' @export

#Reduced Surface force
ReducedSF.2d = function(meshP, ExSurf){
  TDOF = 2*NROW(meshP)
  z= NROW(meshP)

  RSF = matrix(0, nrow = TDOF)
  for (n in 1:TDOF){
    cell = numeric(TDOF)
    for (m in 1:z){
      cell[m] = ExSurf[[m]][n,]}
    RSF[n,] = sum(cell) }
  return(RSF)}
