% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GlobalMat.2d.R
\name{GlobalMat.2d}
\alias{GlobalMat.2d}
\title{GlobalMat.2d}
\usage{
GlobalMat.2d(meshP, meshT, ExEM)
}
\arguments{
\item{meshP}{Matrix (2 x n) containing coordinate points of the mesh nodes.}

\item{meshT}{Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.}

\item{ExEM}{Expanded element matrix}
}
\value{
Produces large (n x n) global matrix
\item{GlobalMat}{Global matrix}
}
\description{
Generates global stiffness matrix - once established, the expanded element matrix must be combined to create the global structural stiffness matrix by adding the expanded matrices.
}
\examples{
data(triMesh)
data(fea_ExEM)

meshP = triMesh$MeshPts$p
meshT = triMesh$MeshPts$T
ExEM = fea_ExEM

gloMat = GlobalMat.2d(meshP, meshT, ExEM)

}
