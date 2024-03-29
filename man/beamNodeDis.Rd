% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beamNodeDis.R
\name{beamNodeDis}
\alias{beamNodeDis}
\title{beamNodeDis}
\usage{
beamNodeDis(beamP, BCtran, BCrot, REM, NodeKnownL, ForceV)
}
\arguments{
\item{beamP}{Matrix (2 x n) of beam coordinates.}

\item{BCtran}{Boundary constraint for nodes to translate in x or y directions.}

\item{BCrot}{Boundary constraint for nodes to rotate.}

\item{REM}{Reduced element matrix, returned from function ReducedEM.}

\item{NodeKnownL}{data frame with constraint parameters applied to each node in the x and y directions. Formatted for use in reduced element matrix. Generated from ApplyBC function.}

\item{ForceV}{Reduced force vector matrix containing the model load parameters. Returned from function ForceVector.}
}
\value{
Produces tables with new node coordinates that are produced by the geometry under an applied load.
\item{NodeDis}{Nodal displacement}
\item{GlobalND}{Nodal displacement in the global environment}
\item{GlobalNDMatrix}{Nodal displacement in the global environment as a reduced matrix}
}
\description{
Calculates global nodal displacements of beam.
}
\examples{
data(beamGeo)
data(beamFV)
data(beamREM)
data(beamBC)

ForceV = beamFV
REM = beamREM
NodeKnownL = beamBC

beamND = beamNodeDis(beamGeo$beamP, beamGeo$BCtran, beamGeo$BCrot, REM, NodeKnownL, ForceV)

}
