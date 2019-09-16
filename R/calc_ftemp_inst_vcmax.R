#' Calculates the instantaneous temperature response of Vcmax
#'
#' Given Vcmax at a reference temperature (argument \code{tcref})
#' this function calculates its temperature-scaling factor following modified Arrhenius
#' kinetics based on Kattge & Knorr (2007). Calculates \eqn{f} for the conversion
#' \deqn{
#'    V = f Vref
#' }
#'
#' @param tcleaf Leaf temperature, or in general the temperature relevant for photosynthesis
#' (degrees Celsius)
#' @param tcgrowth (Optional) Growth temperature, in the P-model, taken to be equal to \code{tcleaf}
#' (in degrees Celsius). Defaults to \code{tcgrowth = tcleaf}.
#' @param tcref Reference temperature (in degrees Celsius)
#'
#' @details The function is given by Kattge & Knorr (2007) as
#' \deqn{
#' 		fv = f(T, \Delta Hv) A/B
#' }
#' where \eqn{f(T, \Delta Hv)} is a regular Arrhenius-type temperature response function (see
#' \link{calc_ftemp_arrh}) with \eqn{Hv=71513} J mol-1,
#' \deqn{
#' 		A = 1 + exp( (T0 \Delta S - Hd) / (T0 R) )
#' }
#' and
#' \deqn{
#' 	    B = 1 + exp( (T \Delta S - Hd) / (TK R) )
#' }
#' Here, \eqn{T} is in Kelvin, \eqn{T0=293.15} K, \eqn{Hd = 200000} J mol-1 is the deactivation
#' energy and \eqn{R} is the universal gas constant and is 8.3145 J mol-1 K-1, and
#' \deqn{
#' 		\Delta S = aS - bS T
#' }
#' with \eqn{aS = 668.39} J mol-1 K-1, and \eqn{bS = 1.07} J mol-1 K-2, and \eqn{T} given in
#' degrees Celsius (!)
#'
#' @references Kattge, J. and Knorr, W.:  Temperature acclimation in a biochemical model of
#' photosynthesis: a reanalysis of data from 36 species, Plant, Cell and Environment, 30,1176–1190, 2007.
#'
#' @return A numeric value for \eqn{fv}
#'
#' @examples 
#' ## Relative change in Vcmax going (instantaneously, i.e. 
#' ## not acclimatedly) from 10 to 25 degrees (percent change):
#' print((calc_ftemp_inst_vcmax(25)/calc_ftemp_inst_vcmax(10)-1)*100 )
#'
#' @export
#'
calc_ftemp_inst_vcmax <- function( tcleaf, tcgrowth = tcleaf, tcref = 25.0 ){
  
  # loal parameters
  Ha    <- 71513  # activation energy (J/mol)
  Hd    <- 200000 # deactivation energy (J/mol)
  Rgas  <- 8.3145 # universal gas constant (J/mol/K)
  a_ent <- 668.39 # offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
  b_ent <- 1.07   # slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)

  tkref <- tcref + 273.15  # to Kelvin

  # conversion of temperature to Kelvin, tcleaf is the instantaneous leaf temperature in degrees C.
  tkleaf <- tcleaf + 273.15

  # calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin !!!
  dent <- a_ent - b_ent * tcgrowth   # 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
  fva <- calc_ftemp_arrh( tkleaf, Ha, tkref = tkref )
  fvb <- (1 + exp( (tkref * dent - Hd)/(Rgas * tkref) ) ) / (1 + exp( (tkleaf * dent - Hd)/(Rgas * tkleaf) ) )
  fv  <- fva * fvb

  return( fv )
}
