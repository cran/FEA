% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ManualAdjust.2d.R
\name{ManualAdjust.2d}
\alias{ManualAdjust.2d}
\title{ManualAdjust.2d}
\usage{
ManualAdjust.2d(meshP, meshT, edge, centroid, loc)
}
\arguments{
\item{meshP}{Matrix (2 x n) containing coordinate points of the mesh nodes.}

\item{meshT}{Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.}

\item{edge}{Coordinate points of the initial geometry.}

\item{centroid}{Matrix (2 x n) of triangle elements.}

\item{loc}{String containing the number of the meshT matrix row of the triangle chosen to be removed.}
}
\value{
Generates new mesh and centroid tables
\item{Meshpts}{Includes both new mesh coordinate points and triangulation of points.}
\item{Centroids}{Centroid positions for each triangle element.}
}
\description{
Allows for manual refinement of the triangular mesh generated based on given parameters. Will remove triangle elements that are identified in the input (loc).
}
\examples{
data(triMesh)
data(polyshape)

meshP = triMesh$MeshPts$p
meshT = triMesh$MeshPts$T
edge =  polyshape$Line
centroid = triMesh$Centroids
loc = c(7, 35, 17)

ManualAdjust.2d(meshP, meshT, edge, centroid, loc)

}
