#' @title beamReducedEM
#'
#' @description Reduced stiffness matrix - use boundary condition to reduce matrix to smaller form by removing systems that are bound.
#'
#' @usage beamReducedEM(GMat, NodeKnownL)
#'
#' @param GMat Global stiffness matrix
#' @param NodeKnownL data frame with constraint parameters applied to each node in the x and y directions. Formatted for use in reduced element matrix. Generated from ApplyBC function.
#'
#' @return Produces a large matrix.
#' \item{ReducedEM}{Reduced element matrix.}
#'
#' @examples
#' data(beamBC)
#' data(beamGloMat)
#'
#' NodeKnownL = beamBC
#' GMat = beamGloMat
#' beamREM = beamReducedEM(GMat, NodeKnownL)
#'
#' @export

beamReducedEM = function(GMat, NodeKnownL){GMat[c(NodeKnownL), c(NodeKnownL)]}
