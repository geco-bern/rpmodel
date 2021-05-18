#' CO2 partial pressure
#'
#' Calculates CO2 partial pressure from concentration in ppm.
#'
#' @param co2 Atmospheric CO2 concentration (ppm)
#' @param patm Atmospheric pressure (Pa).
#'
#'
#' @return CO2 partial pressure in Pa.
#'
#' @export
#'
co2_to_ca <- function( co2, patm ){
  # Pa, atms. CO2
  ca   <- ( 1.0e-6 ) * co2 * patm
  return( ca )
}
