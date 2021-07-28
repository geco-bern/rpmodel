#' Calculates the instantaneous temperature response of Jmax
#'
#' Given Jmax at a reference temperature (argument \code{tcref})
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
#' \link{calc_ftemp_arrh}) with \eqn{Hv=49884} J mol-1,
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
#' with \eqn{aS = 659.70} J mol-1 K-1, and \eqn{bS = 0.75} J mol-1 K-2, and \eqn{T} given in
#' degrees Celsius (!)
#'
#' @references Kattge, J. and Knorr, W.:  Temperature acclimation in a biochemical model of
#' photosynthesis: a reanalysis of data from 36 species, Plant, Cell and Environment, 30,1176–1190, 2007.
#'
#' @return A numeric value for \eqn{fv}
#'
#' @examples
#' ## Relative change in Jmax going (instantaneously, i.e.
#' ## not acclimatedly) from 10 to 25 degrees (percent change):
#' print((calc_ftemp_inst_jmax(25)/calc_ftemp_inst_jmax(10)-1)*100 )
#'
#' @export
#'
calc_ftemp_inst_jmax <- function( tcleaf, tcgrowth = tcleaf, tchome = NA, tcref = 25.0, kuma_par = T ){

  # loal parameters
  Hd    <- 200000 # deactivation energy (J/mol)
  tkref <- tcref + 273.15  # to Kelvin
  tkleaf <- tcleaf + 273.15  # conversion of temperature to Kelvin, tcleaf is the instantaneous leaf temperature in degrees C.

  # Kattge2007 Parametrization
  Ha    <- 49884  # activation energy (J/mol)
  Rgas  <- 8.3145 # universal gas constant (J/mol/K)
  a_ent <- 659.70 # offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
  b_ent <- 0.75   # slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)

  # calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin !!!
  dent <- a_ent - b_ent * tcgrowth   # 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's

  if(kuma_par){
    #-------------------------
    # Kumarathunge2019 Implementation:
    # local parameters
    Ha    = 40710  # activation energy (J/mol)
    a_ent = 658.77 # offset of entropy vs. temperature relationship (J/mol/K)
    b_ent = 0.84   # slope of entropy vs. temperature relationship (J/mol/K^2)
    c_ent = 0.52   # 2nd slope of entropy vs. temperature (J/mol/K^2)

    # Entropy calculation, equations given in Celsius, not in Kelvin
    dent = a_ent - (b_ent * tchome) - c_ent * (tcgrowth - tchome)
    # print*,'Kumarathunge jmax, dent:', dent
    #-------------------------
  }

  fva <- calc_ftemp_arrh( tkleaf, Ha, tkref = tkref )
  fvb <- (1 + exp( (tkref * dent - Hd)/(Rgas * tkref) ) ) / (1 + exp( (tkleaf * dent - Hd)/(Rgas * tkleaf) ) )
  fv  <- fva * fvb

  return( fv )
}
