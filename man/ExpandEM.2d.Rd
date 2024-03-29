% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExpandEM.2d.R
\name{ExpandEM.2d}
\alias{ExpandEM.2d}
\title{ExpandEM.2d}
\usage{
ExpandEM.2d(meshP, meshT, centroid, EMatrixlist)
}
\arguments{
\item{meshP}{Matrix (2 x n) containing coordinate points of the mesh nodes.}

\item{meshT}{Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.}

\item{centroid}{Matrix (2 x n) containing coordinate points of the centroid of each triangular element.}

\item{EMatrixlist}{EMPStress or EMPStrain generated from ElementMat function. List of element matrices.}
}
\value{
Produces large (n x n) matrix.
\item{ExpandedMat}{The expanded element matrix}
}
\description{
Generates the expanded element matrix, which represents the contribution of individual finite elements towards the global structural matrix
}
\examples{
data(triMesh)
data(fea_EM)

meshP = triMesh$MeshPts$p
meshT = triMesh$MeshPts$T
centroid = triMesh$Centroids
EMatrixlist = fea_EM$EMPStress

fea_ExEM = ExpandEM.2d(meshP, meshT, centroid, EMatrixlist)

}
