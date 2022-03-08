#' @title ReducedEM
#'
#' @description Reduced stiffness matrix - use boundary condition to reduce matrix to smaller form by removing systems that are bound.
#'
#' @param GMat Global stiffness matrix
#' @param NodeKnownL data frame with constraint parameters applied to each node in the x and y directions. Formatted for use in reduced element matrix. Generated from ApplyBC function.
#'
#' @return Produces a large matrix.
#' \item{ReducedEM}{Reduced element matrix.}
#'
#' @examples
#' data(gloMat)
#' data(bound)
#' GMat = gloMat
#' NodeKnownL = bound
#' reduc_EM = ReducedEM(GMat, NodeKnownL)
#'
#' @export

ReducedEM = function(GMat, NodeKnownL){GMat[c(NodeKnownL), c(NodeKnownL)]}
