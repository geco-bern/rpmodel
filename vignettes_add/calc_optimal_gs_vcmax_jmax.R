calc_optimal_gs_vcmax_jmax <- function( kmm, gammastar, ns_star, ca, vpd, ppfd, fapar, kphio, beta, c_cost, vcmax_start, gs_start, jmax_start ){

  optimise_this_gs_vcmax_jmax <- function( par, args, iabs, kphio, beta, c_cost, maximize=FALSE, return_all=FALSE ){
    
    kmm       <- args[1]
    gammastar <- args[2]
    ns_star   <- args[3]
    ca        <- args[4]
    vpd       <- args[5]
    
    vcmax <- par[1]
    gs    <- par[2]
    jmax  <- par[3]
    
    ## Electron transport is limiting
    ## Solve Eq. system
    ## A = gs (ca- ci)
    ## A = kphio * iabs (ci-gammastar)/ci+2*gammastar) * L
    ## L = 1 / sqrt(1 + ((4 * kphio * iabs)/jmax)^2)
    
    ## This leads to a quadratic equation:
    ## A * ci^2 + B * ci + C  = 0
    ## 0 = a + b*x + c*x^2
    
    ## with
    L <- 1.0 / sqrt(1.0 + ((4.0 * kphio * iabs)/jmax)^2)
    A <- -gs
    B <- gs * ca - 2 * gammastar * gs - L * kphio * iabs
    C <- 2 * gammastar * gs * ca + L * kphio * iabs * gammastar
    
    ci_j <- QUADM(A, B, C)
    a_j  <- kphio * iabs * (ci_j - gammastar)/(ci_j + 2 * gammastar) * L
    
    
    ## Rubisco is limiting
    ## Solve Eq. system
    ## A = gs (ca- ci)
    ## A = Vcmax * (ci - gammastar)/(ci + Kmm)
    
    ## This leads to a quadratic equation:
    ## A * ci^2 + B * ci + C  = 0
    ## 0 = a + b*x + c*x^2
    
    ## with
    A <- -1.0 * gs
    B <- gs * ca - gs * kmm - vcmax
    C <- gs * ca * kmm + vcmax * gammastar
    
    ci_c <- QUADM(A, B, C)
    a_c <- vcmax * (ci_c - gammastar) / (ci_c + kmm)
    
    ## Take minimum of the two assimilation rates and maximum of the two ci
    assim <- min( a_j, a_c )
    ci <- max(ci_c, ci_j)
    
    ## only cost ratio is defined. for this here we need absolute values. Set randomly
    cost_transp <- 1.6 * ns_star * gs * vpd
    cost_vcmax  <- beta * vcmax
    cost_jmax   <- c_cost * jmax
    
    ## Option B: This is equivalent to the P-model with its optimization of ci:ca.
    if (assim<=0) {
      net_assim <- -(999999999.9)
    } else {
      net_assim <- -(cost_transp + cost_vcmax + cost_jmax) / assim
    }
    
    if (maximize) net_assim <- -net_assim

    # print(par)
    # print(net_assim)
    # 
    if (return_all){
      return( tibble( vcmax_mine = vcmax, jmax_mine = jmax, gs_mine = gs, ci_mine = ci, chi_mine = ci/ca, a_c_mine = a_c, a_j_mine = a_j, assim = assim, ci_c_mine = ci_c, ci_j_mine = ci_j, cost_transp = cost_transp, cost_vcmax = cost_vcmax, cost_jmax = cost_jmax, net_assim = net_assim  ) )
    } else {
      return( net_assim )
    }
  }
  
  out_optim <- optimr::optimr(
    par        = c( vcmax_start, gs_start, jmax_start ), # starting values
    lower      = c( vcmax_start*0.001, gs_start*0.001, jmax_start*0.001 ),
    upper      = c( vcmax_start*1000, gs_start*1000, jmax_start*1000 ),
    fn         = optimise_this_gs_vcmax_jmax,
    args       = c(kmm, gammastar, ns_star, ca, vpd),
    iabs       = (ppfd * fapar),
    kphio      = kphio,
    beta       = beta,
    c_cost     = c_cost/4,
    method     = "L-BFGS-B",
    maximize   = TRUE,
    control    = list( maxit = 100000 )
  )
  
  varlist <- optimise_this_gs_vcmax_jmax( par=out_optim$par, args=c(kmm, gammastar, ns_star, ca, vpd), iabs=(fapar*ppfd), kphio, beta, c_cost/4, maximize=FALSE, return_all=TRUE )

  return(varlist)
}
