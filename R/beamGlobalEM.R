#' @title beamGlobalEM
#'
#' @description Generates global stiffness matrix for beams.
#'
#' @usage beamGlobalEM(beamExEM)
#'
#' @param beamExEM Expanded Element Matrix
#'
#' @return Produces large (n x n) global matrix
#' \item{GlobalMat}{Global matrix}
#'
#' @examples
#' data(beamExMat)
#'
#' beamExEM = beamExMat
#' GMat = beamGlobalEM(beamExEM)
#'
#' @export

beamGlobalEM = function(beamExEM){
  bGloMat = Reduce('+', beamExEM)

  return(bGloMat) }
