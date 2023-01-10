#' @title beamExpandEM
#'
#' @description Expanded element matrix for beam.
#'
#' @usage beamExpandEM(beamP, beamT, ElementMat)
#'
#' @param beamP Matrix (2 x n) of beam coordinates.
#' @param beamT Matrix (2 x n) containing the number of the coordinate point as shown in beamP that connect to form a given beam (Discretization table).
#' @param ElementMat Element stiffness matrix list.
#'
#' @return produces large (n x n) element matrix from initial element matrix.
#' \item{beamExMat}{The expanded element matrix}
#'
#' @examples
#' data(beamGeo)
#' data(beamEmat)
#'
#' ElementMat = beamEmat
#' beamExMat =  beamExpandEM(beamGeo$beamP, beamGeo$beamT, ElementMat)
#'
#' @export

beamExpandEM = function(beamP, beamT, ElementMat){
  n = NROW(beamP)
  TDOF = 2 * n
  EM = ElementMat

  beamExEm = function(TDOF, EM, ik, jk){
    r1 = (ik-1)+ik; r2 = (ik-1)+(ik+1); r3 = (jk-2)+(jk+1); r4 = (jk-2)+(jk+2);
    bigmatrix = matrix(vector(length = TDOF*TDOF), nrow = TDOF, byrow = TRUE);
    bigmatrix[c(r1, r2, r3, r4), c(r1, r2, r3, r4)] = EM[[m]]; return(bigmatrix)}

ExpandedMat = list()
n = NROW(beamT)
for (m in 1:n){
  ik = beamT[m,1]
  jk = beamT[m,2]
  ExpandedMat[[m]] = beamExEm(TDOF, EM, ik, jk)}

return("beamExMat" = ExpandedMat)}
