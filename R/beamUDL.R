#' @title beamUDL
#'
#' @description Uniformly distributes load over the length of the beam.
#'
#' @usage beamUDL(beamP, beamT, Length, fx, fy)
#'
#' @param beamP Matrix (2 x n) of beam coordinates.
#' @param beamT Matrix (2 x n) containing the number of the coordinate point as shown in beamP that connect to form a given beam (Discretization table).
#' @param Length Length of beam.
#' @param fx Load value (newtons) in the x direction.
#' @param fy Load value (newtons) in the y direction.
#'
#' @return Produces matrix representing uniformly distributed load on beam
#' \item{DLMatrix}{Column matrix for beam distributed load}
#' \item{ExpandedDLMatrix}{Expanded beam distribution load}
#' \item{ReductedDLMatrix}{Reduced beam distribution load}
#'
#' @examples
#' data(beamGeo)
#' data(beamDime)
#'
#' Length = beamDime$Length
#' beamUDL = beamUDL(beamGeo$beamP, beamGeo$beamT, Length, beamGeo$fx, beamGeo$fy)
#'
#' @export

beamUDL = function(beamP, beamT, Length, fx, fy) {
  DOF = NROW(beamP)
  TDOF = 2*NROW(beamP)
  f = rowSums(fx, fy)

#Column matrix for beam distributed load
  BeamUDL = function(DOF, f, Length){
    LoadReduc = matrix(f, ncol = 1)
    DisLMat = matrix(c(-Length[m]/2, -(Length[m]^2)/12, -Length[m]/2, (Length[m]^2)/12), nrow = DOF, byrow = TRUE)
    DisLoad = LoadReduc[m]*DisLMat
    return(DisLoad)}

#Expanded beam distributed load
  BeamUDL_EE = function(TDOF, BUDL, ik, jk)
    {LoadReduc = matrix(f, ncol = 1)
    r1 = (ik-1)+ik; r2 = (ik-1)+(ik+1); r3 = (jk-2)+(jk+1); r4 = (jk-2)+(jk+2);
    bigmatrix = matrix(vector(length = TDOF*TDOF), nrow = TDOF, byrow = TRUE);
    bigmatrix[c(r1, r2, r3, r4), c(r1, r2, r3, r4)] = BUDL[[m]]; return(bigmatrix)}

  n = NROW(beamT)
  BUDL = list()
    for (m in 1:n){BUDL[[m]] = BeamUDL(DOF, f, Length)}
  BK = list()
    for (m in 1:n){
      ik = beamT[m,1]
      jk = beamT[m,2]
      BK[[m]] = BeamUDL_EE(TDOF, BUDL, ik, jk)}

#Reduced beam distributed load
DistributeL = Reduce('+', BK)

Rlist = list("DLMatrix" = BUDL, "ExpandedDLMatrix" = BK, "ReducedDLMatrix" = DistributeL)
}
