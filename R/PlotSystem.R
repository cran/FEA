#' @title PlotSystem
#'
#' @description Generates heat map for given stress or strain on the geometry. Threshold values for the color must be assigned.
#' @usage PlotSystem(meshP, meshT, PlotVal, a, b, c, d, e, f, g, h, i, j,
#'                   Oc, ac, bc, cc, dc, ec, fc, gc, hc, ic, jc)
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param PlotVal Value to be plotted, either stress or strain, return from function LocalStress function.
#' @param a Threshold 1
#' @param b Threshold 2
#' @param c Threshold 3
#' @param d Threshold 4
#' @param e Threshold 5
#' @param f Threshold 6
#' @param g Threshold 7
#' @param h Threshold 8
#' @param i Threshold 9
#' @param j Threshold 10
#'
#' @param Oc Color for all zero values
#' @param ac Color 1
#' @param bc Color 2
#' @param cc Color 3
#' @param dc Color 4
#' @param ec Color 5
#' @param fc Color 6
#' @param gc Color 7
#' @param hc Color 8
#' @param ic Color 9
#' @param jc Color 10
#'
#' @return Plot of colored polygon with mesh colored based on the plot value
#'
#' @examples
#' \donttest{
#' data(triMesh)
#' data(fea_result)
#'
#' meshP = triMesh$MeshPts$p
#' meshT = triMesh$MeshPts$T
#' PlotVal = abs(fea_result$Stress[,1])
#' Oc = "slateblue"; ac = "steelblue2"; bc = "cyan2"; cc = "palegreen2";
#' dc = "darkolivegreen1"; ec = "lemonchiffon"; fc = "lightgoldenrod1"; gc = "gold";
#' hc= "lightsalmon"; ic= "tomato"; jc= "firebrick3"
#' a = 1e5;  b = 5e5;  c = 1e6;  d = 5e6;  e = 1e7;  f = 5e7;  g = 1e8;  h = 5e8; i = 1e9; j =5e9
#'
#' PlotSystem(meshP, meshT, PlotVal, a, b, c, d, e, f, g, h, i, j,
#'            Oc, ac, bc, cc, dc, ec, fc, gc, hc, ic, jc)
#'}
#'
#' @export

#Plot colored elements
PlotSystem = function(meshP, meshT, PlotVal, a, b, c, d, e, f, g, h, i, j,
                                Oc, ac, bc, cc, dc, ec, fc, gc, hc, ic, jc){
  m= n= o= NROW(meshP)*2
  z= NROW(meshT) # m=col, n=row, z=element#
  ColorMap = matrix(nrow=z)

  for (m in 1:z){
      if (PlotVal[m] == 0) {ColorMap[m] = Oc}
      else if (PlotVal[m] < a) {ColorMap[m] = ac}
      else if (PlotVal[m] < b) {ColorMap[m] = bc}
      else if (PlotVal[m] < c) {ColorMap[m] = cc}
      else if (PlotVal[m] < d) {ColorMap[m] = dc}
      else if (PlotVal[m] < e) {ColorMap[m] = ec}
      else if (PlotVal[m] < f) {ColorMap[m] = fc}
      else if (PlotVal[m] < g) {ColorMap[m] = gc}
      else if (PlotVal[m] < h) {ColorMap[m] = hc}
      else if (PlotVal[m] < i) {ColorMap[m] = ic}
      else if (PlotVal[m] < j) {ColorMap[m] = jc}
      else {ColorMap[m] = "black"}
    }

  PolyCo = list()
  PolyColor= function(meshP, meshT){
    x1 = meshP[meshT[m,1],1]
    y1 = meshP[meshT[m,1],2]
    x2 = meshP[meshT[m,2],1]
    y2 = meshP[meshT[m,2],2]
    x3 = meshP[meshT[m,3],1]
    y3 = meshP[meshT[m,3],2]
    sp::Polygons(list(sp::Polygon(cbind(c(x1, x2, x3, x1), c(y1, y2, y3, y1)))), m)}
  for (m in 1:z){PolyCo[[m]]=PolyColor(meshP, meshT)}
  PolyCo2 = sp::SpatialPolygons(PolyCo)
  sp::plot(PolyCo2, col = ColorMap) }
