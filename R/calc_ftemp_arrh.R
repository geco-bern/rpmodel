#' Arrhenius temperature response
#'
#' Given a kinetic rate at a reference temperature (argument \code{tkref})
#' this function calculates its temperature-scaling factor following Arrhenius kinetics.
#'
#' @param tk Temperature (Kelvin)
#' @param dha Activation energy (J mol-1)
#' @param tkref Reference temperature (Kelvin)
#'
#' @details To correct for effects by temperature following Arrhenius kinetics,
#' and given a reference temperature \eqn{T0}, \eqn{f} calculates the temperature
#' scaling. Arrhenius kinetics are described by an equation of form
#' \eqn{x(T)= exp(c - \Delta Ha / (T R))}. The temperature-correction
#' function \eqn{f(T, \Delta Ha)} is thus given by \eqn{f=x(T)/x(T0)} which is:
#' \deqn{
#'      f = exp( \Delta Ha (T - T0) / (T0 R T_K) )
#' }
#' \eqn{\Delta Ha} is given by argument \code{dha}. \eqn{T} is given by argument
#' \code{tk} and has to be provided in Kelvin. \eqn{R} is the universal gas constant
#' and is 8.3145 J mol-1 K-1. Note that this is equivalent to 
#' \deqn{
#'      f = exp( (\Delta Ha/R) (1/T0 - 1/T) )  
#' }
#'
#' @return A numeric value for \eqn{f}
#'
#' @examples ftemp <- calc_ftemp_arrh(283.15, 100000, tkref = 298.15)
#'
#' @export
#'
calc_ftemp_arrh <- function( tk, dha, tkref = 298.15 ){

  # Note that the following forms are equivalent:
  # ftemp = exp( dha * (tk - 298.15) / (298.15 * kR * tk) )
  # ftemp = exp( dha * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) )
  # ftemp = exp( (dha/kR) * (1/298.15 - 1/tk) )
  #-----------------------------------------------------------------------
  kR   <- 8.3145     # Universal gas constant, J/mol/K
  ftemp <- exp( dha * (tk - tkref) / (tkref * kR * tk) )

  return(ftemp)
}
