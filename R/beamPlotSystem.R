#' @title beamPlotSystem
#'
#' @description Generates heat map for given stress or strain on the beam geometry. Threshold values for the color must be assigned.
#'
#' @usage beamPlotSystem(beamP, beamT, PlotVal, a, b, c, d, e, f, g, h, i, j,
#'                   Oc, ac, bc, cc, dc, ec, fc, gc, hc, ic, jc, LWD)
#'
#' @param beamP Matrix (2 x n) of beam coordinates.
#' @param beamT Matrix (2 x n) containing the number of the coordinate point as shown in beamP that connect to form a given beam (Discretization table).
#' @param PlotVal Value to be plotted, either stress or strain, return from function beamLocalStress function.
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
#' @param LWD Line (beam) width
#'
#' @return Plot of colored beam based on the plot value
#'
#' @examples
#' data(beamGeo)
#' data(beamStressResult)
#'
#' PlotVal = beamStressResult
#'
#' Oc = "slateblue"; ac = "steelblue2"; bc = "cyan2"; cc = "palegreen2";
#' dc = "darkolivegreen1"; ec = "lemonchiffon"; fc = "lightgoldenrod1";
#' gc = "gold"; hc= "lightsalmon"; ic= "tomato"; jc= "firebrick3"
#'
#' a = 1e5;  b = 5e5;  c = 1e6;  d = 5e6;  e = 1e7;  f = 5e7;  g = 1e8;  h = 5e8; i = 1e9; j =5e9
#' beamPlotSystem(beamGeo$beamP, beamGeo$beamT, PlotVal, a, b, c, d, e, f, g, h, i, j, Oc,
#' ac, bc, cc, dc, ec, fc, gc, hc, ic, jc, LWD = 4)
#'
#' @export

beamPlotSystem = function(beamP, beamT, PlotVal, a, b, c, d, e, f, g, h, i, j,
                           Oc, ac, bc, cc, dc, ec, fc, gc, hc, ic, jc, LWD){
  m= z= NROW(beamT)
  ColorMap = matrix(nrow=m)

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
    else {ColorMap[m] = "black"}}

  for (m in 1:z){
    Beam = c(beamP[beamT[m,1],], beamP[beamT[m,2],])
    graphics::segments(Beam[1],Beam[2], Beam[3],Beam[4], col = ColorMap[m], lwd= LWD)}}
