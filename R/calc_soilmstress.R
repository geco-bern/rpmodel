#' Calculates an empirical soil moisture stress factor
#'
#' Calculates an empirical soil moisture stress factor as a function of relative 
#' soil moisture (fraction of field capacity).
#'
#' @param soilm Relative soil moisture as a fraction 
#' of field capacity (unitless). Defaults to 1.0 (no soil moisture stress).
#' @param meanalpha Local annual mean ratio of 
#' actual over potential evapotranspiration, measure for average aridity. Defaults to 1.0. 
#' @param apar_soilm (Optional, used only if \code{do_calc_soilmstress==TRUE}) Parameter determining the 
#' sensitivity of the empirical soil moisture stress function. Defaults to 0.0, the empirically fitted value
#' as presented in Stocker et al. (2019) Geosci. Model Dev. for model setup 'FULL' (corresponding to a setup
#' with \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_calc_soilmstress=TRUE}).
#' @param bpar_soilm (Optional, used only if \code{do_calc_soilmstress==TRUE}) Parameter determining the 
#' sensitivity of the empirical soil moisture stress function. Defaults to 0.685, the empirically fitted value
#' as presented in Stocker et al. (2019) Geosci. Model Dev. for model setup 'FULL' (corresponding to a setup
#' with \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_calc_soilmstress=TRUE}).
#'
#' @details The soil moisture stress factor is calculated using a quadratic function that
#' is 1 above \code{soilm} = 0.6 and has a sensitivity, given by the y-axis cutoff, 
#' (zero soil moisture), determined  by average aridity (argument \code{meanalpha}) as:
#' \deqn{
#'     \beta = q(\theta - \theta*)^2 + 1
#' }
#' for \eqn{\theta < \theta*} and \eqn{\beta = 1.0} otherwise. \eqn{\theta*} is fixed at 0.6.
#' \eqn{q} is the sensitivity parameter and is calculated as a linear function of average aridity, 
#' quantified by the local annual mean ratio of actual over potential evapotranspiration, termed 
#' \eqn{\alpha}:
#' \deqn{
#'      q=(\beta0-1)/(\theta*-\theta0)^2
#' }
#' \eqn{\theta0} is 0.0, and 
#' \deqn{
#'      \beta0 = a + b \alpha
#' }
#' \eqn{a} is given by argument \code{apar}, \eqn{b} is given by argument \code{bpar}.
#' 
#' @return A numeric value for \eqn{\beta}
#'
#' @examples 
#' ## Relative reduction (%) in GPP due to soil moisture stress at 
#' ## relative soil water content ('soilm') of 0.2:
#' print((calc_soilmstress(0.2)-1)*100 )
#' 
#' @references  Stocker, B. et al. Geoscientific Model Development Discussions (in prep.)
#'
#' @export
calc_soilmstress <- function( soilm, meanalpha = 1.0, apar_soilm = 0.0, bpar_soilm = 0.685 ){

  # Fixed parameters
  x0 <- 0.0
  x1 <- 0.6

  if (soilm > x1) {

    outstress <- 1.0

  } else {

    y0 <- (apar_soilm + bpar_soilm * meanalpha)

    beta <- (1.0 - y0) / (x0 - x1)^2
    outstress <- 1.0 - beta * ( soilm - x1 )^2
    outstress <- max( 0.0, min( 1.0, outstress ) )

  }

  return(outstress)
}
