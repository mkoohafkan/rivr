#' This package is designed as an educational tool for students and 
#' instructors of undergraduate courses in open channel hydraulics. 
#' Functions are provided for computing normal and critical depths, 
#' steady (e.g. backwater curves) and unsteady (flood wave routing) 
#' flow computations for prismatic trapezoidal channels of arbitrary 
#' geometry.
#' @name rivr-package
#' @docType package
#' @useDynLib rivr
#' @importFrom Rcpp evalCpp
NULL

#' @title Gradually-varied flow profiles
#' @description Compute the gradually-varied flow profile of a prismatic channel.
#' @param So Channel slope [\eqn{L L^{-1}}].
#' @param n Manning's roughness coefficient.
#' @param Q Flow rate [\eqn{L^3 T^{-1}}].
#' @param y0 The water depth at the control section [\eqn{L}].
#' @param Cm Unit conversion coefficient for Manning's equation. For SI units, Cm = 1.
#' @param g Gravitational acceleration [\eqn{L T^{-2}}].
#' @param B Channel bottom width [\eqn{L}].
#' @param SS Channel sideslope [\eqn{L L^{-1}}].
#' @param z Elevation reference datum at control section [\eqn{L}]. Default is 0.
#' @param x Distance reference at control section [\eqn{L}]. Default is 0.
#' @param stepdist The spatial interval used in the Standard step method [\eqn{L}].
#' @param totaldist The total distance upstream (or downstream) to compute the profile [\eqn{L}].
#' @return data.frame with columns:
#'   \item{x}{Along-channel distance.}
#'   \item{z}{Elevation.}
#'   \item{y}{Flow depth.}
#'   \item{v}{Flow velocity.}
#'   \item{A}{Flow area.}
#'   \item{Sf}{Friction slope.}
#'   \item{E}{Total energy.}
#'   \item{Fr}{Froude Number.}
#' @details Computes the longitudinal water surface profile of a prismatic 
#'   channel using the standard step method by solving the non-linear ODE 
#'   \eqn{\frac{dy}{dx} = \frac{S_0 - S_f}{1 - Fr^2}}. The standard-step 
#'   method operates by stepping along the channel by a constant distance 
#'   interval, starting from a cross-section where the flow depth is known 
#'   (the control section). The flow depth is computed at the adjacent 
#'   cross-section (target section). The computed value at the target is then 
#'   used as the basis for computing flow depth at the next cross-section, i.e. 
#'   the previous target section becomes the new control section for each step. 
#'   A Newton-Raphson scheme is used each step to compute the flow depth and 
#'   friction slope. Technically, the average friction slope of the control and
#'   target section is used to compute the flow depth at the target section.
#' @examples
#' # example M1 profile
#' compute_profile(0.001, 0.045, 250, 2.7, 1.486, 32.2, 100, 0, stepdist=50, totaldist=3000)
#' # example M2 profile
#' compute_profile(0.001, 0.045, 250, 0.64, 1.486, 32.2, 100, 0, stepdist=50, totaldist=3000)
#' @export
compute_profile = function(So, n, Q, y0, Cm, g, B, SS, z=0, x=0, stepdist, totaldist){
  stepsize = stepdist
  yc = critical_depth(Q, y0, g, B, SS)
  if(yc < y0){
    message('flow at control section is subcritical. ', 
	  'Upstream profile will be computed.')
	stepsize = -abs(stepdist)
  } else if(yc > y0) {
    message('flow at control section is supercritical. ',
	  'Downstream profile will be computed.')
	stepsize = abs(stepdist)
  }
  else {
    stop('flow at control section is critical. Starting flow depth must ',
      'be greater than or less than critical depth.')
  }
  res = as.data.frame(loop_step(So, n, Q, Cm, g, y0, B, SS, z, x, stepsize, totaldist))
  names(res) = c("x", "z", "y", "v", "A", "Sf", "E", "Fr")
  return(res)
}
