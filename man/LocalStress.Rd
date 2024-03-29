% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LocalStress.R
\name{LocalStress}
\alias{LocalStress}
\title{LocalStress}
\usage{
LocalStress(meshP, meshT, Y, Nu, GlobalND)
}
\arguments{
\item{meshP}{Matrix (2 x n) containing coordinate points of the mesh nodes.}

\item{meshT}{Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.}

\item{Y}{Value of Young's (Elastic) modulus}

\item{Nu}{Value of Poisson's ratio}

\item{GlobalND}{Global nodal displacement, return from function NodeDis}
}
\value{
Completes FEM by calculating values of stress and strain, produces three (3) [3 x n] matrix.
\item{Strain}{Calculated strain. [x, y, tau]}
\item{Stress}{Calculated stress in pascals. [x, y, tau]}
\item{StressStrain}{Stress as calucated from strain. [x, y, tau]}
}
\description{
Calculates local stress and strain for triangular elements of the mesh
}
\examples{
\donttest{
data(triMesh)
data(displacN)

meshP = triMesh$MeshPts$p
meshT = triMesh$MeshPts$T
Y = matrix(20e9, nrow = NROW(meshT))
Nu = matrix(0.45, nrow = NROW(meshT))
GlobalND = displacN$GlobalND

fea_result = LocalStress(meshP, meshT, Y, Nu, GlobalND)
}

}
