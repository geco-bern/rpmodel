#' Dampen inputs of rpmodel
#'
#' Applies an exponential dampening input time series with specified time scale.
#'
#' @param vec A numeric vector for the time series of a daily meteorological
#'  variable used as input for \code{rpmodel} (temperature, 
#'  vapour pressure deficit, CO2, or atmospheric pressure). The length of
#'  \code{x} must be at least 365, i.e., corresponding to one year.
#' @param tau The time scale of dampening (e-folding time scale of a 
#'  perturbation). Must be smaller or equal to 365 d.
#' @return A numeric vector of equal length as \code{x} with damped variation.
#'  The dampening is calculated as:
#'  \deqn{
#'   S(t+1) - S(t) = (X(t+1) - S(t)) / \tau
#'  }
#'  Where \eqn{X} is the daily varying time series given by argument \code{x},
#'  \eqn{S} is the dampened time returned by this function, and \eqn{\tau} is 
#'  the decay time scale of a perturbation, given by argument \code{tau}.
#'
#' @export
#'
#' @examples 
#' \dontrun{
#' dampen_vec(
#'  vec = 20 * (sin(doy*pi/(365)))^2 + rnorm(365, mean = 0, sd = 5),
#'  tau = 40 
#'  )
#' }

dampen_vec <- function( vec, tau ){

  if (length(vec)<365){

    stop("dampen_vec(): Aborting.
         Length of argument 'vec' must be at least 365.")

  } else {

    ## Add one year's data to the head of 'vec' to avoid "boundary effects"
    vec <- c(vec[1:365], vec)

    ## apply dampening
    vec_damped <- rep(NA, length(vec))
    vec_damped[1] <- vec[1] # this is the unwanted "boundary effect"
    for (idx in 2:length(vec)){
      dvar <- (1.0/tau) * (vec[idx] - vec_damped[idx - 1])
      vec_damped[idx] <- vec_damped[idx - 1] + dvar
    }

    ## remove the first year again, that was added artifically before
    vec_damped <- vec_damped[366:length(vec_damped)]

  }

  return( vec_damped )
}