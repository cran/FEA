% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReducedSF.R
\name{ReducedSF}
\alias{ReducedSF}
\title{ReducedSF}
\usage{
ReducedSF(meshP, ExSurf)
}
\arguments{
\item{meshP}{Matrix (2 x n) containing coordinate points of the mesh nodes.}

\item{ExSurf}{Expanded surface matrix, output from ExpandSFT}
}
\value{
Produces a large matrix.
\item{RSF}{Produces a large, reduced surface force matrix}
}
\description{
Reduced matrix of surface forces
}
\examples{
data(triMesh)
data(expSurf)
meshP = triMesh$MeshPts$p
ExSurf = expSurf
reduc_SF = ReducedSF(meshP, ExSurf)

}
