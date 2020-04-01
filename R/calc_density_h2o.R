#' Density of water
#'
#' Calculates the density of water as a function of temperature and atmospheric
#' pressure, using the Tumlirz Equation.
#'
#' @param tc numeric, air temperature (tc), degrees C
#' @param p numeric, atmospheric pressure (p), Pa
#'
#' @return numeric, density of water, kg/m^3
#'
#' @examples print("Density of water at 20 degrees C and standard atmospheric pressure:")
#' print(calc_density_h2o(20, 101325))
#'
#' @references  F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of
#' pure water and sea water, Tech. Rept., Marine Physical
#' Laboratory, San Diego, CA.
#'
#' @export
#'
calc_density_h2o <- function( tc, p ){

  # Calculate lambda, (bar cm^3)/g:
  my_lambda <- 1788.316 +
    21.55053*tc +
    -0.4695911*tc*tc +
    (3.096363e-3)*tc*tc*tc +
    -(7.341182e-6)*tc*tc*tc*tc

  # Calculate po, bar
  po <- 5918.499 +
    58.05267*tc +
    -1.1253317*tc*tc +
    (6.6123869e-3)*tc*tc*tc +
    -(1.4661625e-5)*tc*tc*tc*tc

  # Calculate vinf, cm^3/g
  vinf <- 0.6980547 +
    -(7.435626e-4)*tc +
    (3.704258e-5)*tc*tc +
    -(6.315724e-7)*tc*tc*tc +
    (9.829576e-9)*tc*tc*tc*tc +
    -(1.197269e-10)*tc*tc*tc*tc*tc +
    (1.005461e-12)*tc*tc*tc*tc*tc*tc +
    -(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc +
    (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc +
    -(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc

  # Convert pressure to bars (1 bar <- 100000 Pa)
  pbar <- (1e-5)*p

  # Calculate the specific volume (cm^3 g^-1):
  v <- vinf + my_lambda/(po + pbar)

  # Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
  rho <- (1e3/v)

  return(rho)
}


