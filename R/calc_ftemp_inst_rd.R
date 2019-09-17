#' Calculates the temperature response of dark respiration
#'
#' Given the dark respiration at the reference temperature 25 degress Celsius,
#' this function calculates its temperature-scaling factor following Heskel et al. 2016.
#'
#' @param tc Temperature (degrees Celsius)
#'
#' @details To correct for effects by temperature Heskel et al. 2016,
#' and given the reference temperature \eqn{T0 =} 25 deg C, this calculates the temperature
#' scaling factor to calculate dark respiration at temperature \eqn{T} (argument \code{tc}) as:
#' \deqn{
#'      fr = exp( 0.1012 (T0-T) - 0.0005(T0^2 - T^2) ) 
#' }
#' where \eqn{T} is given in degrees Celsius.
#' 
#' @return A numeric value for \eqn{fr}
#' 
#' @references  Heskel,  M.,  O’Sullivan,  O.,  Reich,  P.,  Tjoelker,  M.,  Weerasinghe,  L.,  Penillard,  A.,Egerton, J., 
#'              Creek, D., Bloomfield, K., Xiang, J., Sinca, F., Stangl, Z., Martinez-De La Torre, A., Griffin, K., 
#'              Huntingford, C., Hurry, V., Meir, P., Turnbull, M.,and Atkin, O.:  Convergence in the temperature response 
#'              of leaf respiration across biomes and plant functional types, Proceedings of the National Academy of Sciences,
#'              113,  3832–3837,  doi:10.1073/pnas.1520282113,2016.
#'
#' @examples 
#' ## Relative change in Rd going (instantaneously, i.e. not 
#' ## acclimatedly) from 10 to 25 degrees (percent change):
#' print( (calc_ftemp_inst_rd(25)/calc_ftemp_inst_rd(10)-1)*100 )
#'
#' @export
#'
calc_ftemp_inst_rd <- function( tc ){

  # loal parameters
  apar <- 0.1012
  bpar <- 0.0005

  fr <- exp( apar * (tc - 25.0) - bpar * (tc^2 - 25.0^2) )

  return(fr)
}
