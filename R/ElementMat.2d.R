#' @title ElementMat.2d
#'
#' @description Generates an element stiffness matrix
#'
#' @usage ElementMat.2d(meshP, meshT, Nu, Y, Thick)
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param Nu Value of Poisson's ratio for each element
#' @param Y Value of Young's (Elastic) modulus for each element
#' @param Thick Value of the thickness of the mesh, a positive value must be given.
#'
#' @return Generates initial element matrix needed for the finite element model.
#' \item{EMPStress}{An element matrix of the geometry under stress.}
#' \item{EMPStrain}{An element matrix of the geometry under strain.}
#'
#' @examples
#' data(triMesh)
#'
#' meshP = triMesh$MeshPts$p
#' meshT = triMesh$MeshPts$T
#' Y = matrix(20e9, nrow = NROW(meshT))
#' Nu = matrix(0.45, nrow = NROW(meshT))
#' Thick = 0.001
#' DOF = 6
#'
#' fea_EM = ElementMat.2d(meshP, meshT, Nu, Y, Thick)
#'
#' @export

ElementMat.2d = function(meshP, meshT, Nu, Y, Thick){
  CST_EM = function(meshP, meshT, Nu, Y, Thick){
    x1 = meshP[meshT[m,1],1]
    y1 = meshP[meshT[m,1],2]
    x2 = meshP[meshT[m,2],1]
    y2 = meshP[meshT[m,2],2]
    x3 = meshP[meshT[m,3],1]
    y3 = meshP[meshT[m,3],2]

    A2 = x3*(y1-y2) + x2*(y3-y1) + x1*(y2-y3)
    A = A2/2

    B = matrix(c(((y2-y3)/A2), 0, ((y3-y1)/A2), 0, ((y1-y2)/A2), 0,
                 0, ((x3-x2)/A2), 0, ((x1-x3)/A2), 0, ((x2-x1)/A2),
                 ((x3-x2)/A2), ((y2-y3)/A2), ((x1-x3)/A2),((y3-y1)/A2), ((x2-x1)/A2),((y1-y2)/A2)), nrow = 3, byrow = TRUE)

    #for plane stress
    d1 = Y[m]/(1-Nu[m]^2)
    d2 = matrix(c(1, Nu[m], 0,
                  Nu[m], 1, 0,
                  0, 0, ((1-Nu[m])/2)), nrow = 3, byrow = TRUE)
    D1 = d1*d2

    #for plane strain
    d3 = (Y[m])/((1+Nu[m])*(1-2*Nu[m]))
    d4 = matrix(c((1-Nu[m]), Nu[m], 0,
                  Nu[m], (1-Nu[m]), 0,
                  0, 0, ((1-2*Nu[m])/2)), nrow = 3, byrow = TRUE)
    D2 = d3 * d4

    #Element matrix
    Emat1 = (Thick*A)* t(B) %*% D1 %*% B #Element matrix for plane stress
    Emat2 = (Thick*A)* t(B) %*% D2 %*% B #Element matrix for plane strain

    Rlist = list("Pstress" = Emat1, "Pstrain" = Emat2)
    return(Rlist) }

  #Run for acquiring individual element matrix
  m= n= z= o= NROW(meshT) # m=col, n=row, z=element#
  EMPStress = list()
  EMPStrain = list()
  for (m in 1:z){
    test6 = CST_EM(meshP, meshT, Nu, Y, Thick)
    EMPStress[[m]] = test6$Pstress
    EMPStrain[[m]] = test6$Pstrain}

  Rlist = list("EMPStress" = EMPStress, "EMPStrain" = EMPStrain)}
