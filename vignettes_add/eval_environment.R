eval_environment <- function(fn, args, kphio, nsteps = 50, method_jmaxlim = "none", show = "chi"){

  calc_aj <- function(kphio, ppfd, ca, chi, gammastar){
    aj <- kphio * ppfd * (ca * chi  - gammastar) / (ca * chi + 2 * gammastar)
    return(aj)
  }

  df_tc <- tibble(
    tc   = seq(from = 0, to = 40, length.out = nsteps),
    vpd  = rep(vpd, nsteps),
    co2  = rep(co2, nsteps),
    ppfd = rep(ppfd, nsteps),
    patm = rep(rpmodel::calc_patm(0), nsteps)
    ) %>%
    mutate( out_pmodel = purrr::pmap(., rpmodel,
        fapar          = 1.0,
        kphio          = 0.05,
        beta           = 146,
        method_optci   = "prentice14",
        method_jmaxlim = method_jmaxlim,
        do_ftemp_kphio = FALSE
        ) ) %>%
    mutate( out_pmodel = purrr::map(out_pmodel, ~as_tibble(.))) %>%
    unnest( out_pmodel ) %>%
    mutate( out = ifelse(fn == "calc_optimal_chi_num", 
                    purrr::pmap(dplyr::select(., kmm, gammastar, ns_star, ca, vpd), calc_optimal_chi_num, beta = 146),
                    ifelse(fn == "calc_optimal_jmax_num",
                      purrr::pmap(dplyr::select(., kmm, gammastar, ns_star, ca, vpd, chi, ppfd), calc_optimal_jmax_num, kphio, beta = 146, c_cost ),
                      NA ))) %>%
    unnest( out ) %>%
    mutate( aj_num = purrr::pmap(dplyr::select(., ppfd, ca, chi = chi_num, gammastar), calc_aj, kphio = kphio)) %>%
    unnest( aj_num ) %>%
    mutate( a_ = gpp / 12.0107 )

  df_vpd <- tibble(
    tc   = rep(tc, nsteps),
    vpd  = seq(from = 0, to = 5000, length.out = nsteps),
    co2  = rep(co2, nsteps),
    ppfd = rep(ppfd, nsteps),
    patm = rep(rpmodel::calc_patm(0), nsteps)
    ) %>%
    mutate( out_pmodel = purrr::pmap(., rpmodel,
        fapar          = 1.0,
        kphio          = 0.05,
        beta           = 146,
        method_optci   = "prentice14",
        method_jmaxlim = method_jmaxlim,
        do_ftemp_kphio = FALSE
        ) ) %>%
    mutate( out_pmodel = purrr::map(out_pmodel, ~as_tibble(.))) %>%
    unnest( out_pmodel ) %>%
    mutate( out = ifelse(fn == "calc_optimal_chi_num", 
                    purrr::pmap(dplyr::select(., kmm, gammastar, ns_star, ca, vpd), calc_optimal_chi_num, beta = 146),
                    ifelse(fn == "calc_optimal_jmax_num",
                      purrr::pmap(dplyr::select(., kmm, gammastar, ns_star, ca, vpd, chi, ppfd), calc_optimal_jmax_num, kphio, beta = 146, c_cost ),
                      NA ))) %>%
    unnest( out ) %>%
    mutate( aj_num = purrr::pmap(dplyr::select(., ppfd, ca, chi = chi_num, gammastar), calc_aj, kphio = kphio)) %>%
    unnest( aj_num ) %>%
    mutate( a_ = gpp / 12.0107 )

  df_co2 <- tibble(
    tc   = rep(tc, nsteps),
    vpd  = rep(vpd, nsteps),
    co2  = seq(from = 150, to = 1000, length.out = nsteps),
    ppfd = rep(ppfd, nsteps),
    patm = rep(rpmodel::calc_patm(0), nsteps)
    ) %>%
    mutate( out_pmodel = purrr::pmap(., rpmodel,
        fapar          = 1.0,
        kphio          = 0.05,
        beta           = 146,
        method_optci   = "prentice14",
        method_jmaxlim = method_jmaxlim,
        do_ftemp_kphio = FALSE
        ) ) %>%
    mutate( out_pmodel = purrr::map(out_pmodel, ~as_tibble(.))) %>%
    unnest( out_pmodel ) %>%
    mutate( out = ifelse(fn == "calc_optimal_chi_num", 
                    purrr::pmap(dplyr::select(., kmm, gammastar, ns_star, ca, vpd), calc_optimal_chi_num, beta = 146),
                    ifelse(fn == "calc_optimal_jmax_num",
                      purrr::pmap(dplyr::select(., kmm, gammastar, ns_star, ca, vpd, chi, ppfd), calc_optimal_jmax_num, kphio, beta = 146, c_cost ),
                      NA ))) %>%
    unnest( out ) %>%
    mutate( aj_num = purrr::pmap(dplyr::select(., ppfd, ca, chi = chi_num, gammastar), calc_aj, kphio = kphio)) %>%
    unnest( aj_num ) %>%
    mutate( a_ = gpp / 12.0107 )

  df_ppfd <- tibble(
    tc   = rep(tc, nsteps),
    vpd  = rep(vpd, nsteps),
    co2  = rep(co2, nsteps),
    ppfd = seq(from = 0, to = 3000, length.out = nsteps),
    patm = rep(rpmodel::calc_patm(0), nsteps)
    ) %>%
    mutate( out_pmodel = purrr::pmap(., rpmodel,
        fapar          = 1.0,
        kphio          = 0.05,
        beta           = 146,
        method_optci   = "prentice14",
        method_jmaxlim = method_jmaxlim,
        do_ftemp_kphio = FALSE
        ) ) %>%
    mutate( out_pmodel = purrr::map(out_pmodel, ~as_tibble(.))) %>%
    unnest( out_pmodel ) %>%
    mutate( out = ifelse(fn == "calc_optimal_chi_num", 
                    purrr::pmap(dplyr::select(., kmm, gammastar, ns_star, ca, vpd), calc_optimal_chi_num, beta = 146),
                    ifelse(fn == "calc_optimal_jmax_num",
                      purrr::pmap(dplyr::select(., kmm, gammastar, ns_star, ca, vpd, chi, ppfd), calc_optimal_jmax_num, kphio, beta = 146, c_cost ),
                      NA ))) %>%
    unnest( out ) %>%
    mutate( aj_num = purrr::pmap(dplyr::select(., ppfd, ca, chi = chi_num, gammastar), calc_aj, kphio = kphio)) %>%
    unnest( aj_num ) %>%
    mutate( a_ = gpp / 12.0107 )

  df_patm <- tibble(
    tc   = rep(tc, nsteps),
    vpd  = rep(vpd, nsteps),
    co2  = rep(co2, nsteps),
    ppfd = rep(ppfd, nsteps),
    patm = seq(from = rpmodel::calc_patm(0), to = rpmodel::calc_patm(4000), length.out = nsteps)
    ) %>%
    mutate( out_pmodel = purrr::pmap(., rpmodel,
        fapar          = 1.0,
        kphio          = 0.05,
        beta           = 146,
        method_optci   = "prentice14",
        method_jmaxlim = method_jmaxlim,
        do_ftemp_kphio = FALSE
        ) ) %>%
    mutate( out_pmodel = purrr::map(out_pmodel, ~as_tibble(.))) %>%
    unnest( out_pmodel ) %>%
    mutate( out = ifelse(fn == "calc_optimal_chi_num", 
                    purrr::pmap(dplyr::select(., kmm, gammastar, ns_star, ca, vpd), calc_optimal_chi_num, beta = 146),
                    ifelse(fn == "calc_optimal_jmax_num",
                      purrr::pmap(dplyr::select(., kmm, gammastar, ns_star, ca, vpd, chi, ppfd), calc_optimal_jmax_num, kphio, beta = 146, c_cost ),
                      NA ))) %>%
    unnest( out ) %>%
    mutate( aj_num = purrr::pmap(dplyr::select(., ppfd, ca, chi = chi_num, gammastar), calc_aj, kphio = kphio)) %>%
    unnest( aj_num ) %>%
    mutate( a_ = gpp / 12.0107 )


  if (show=="chi"){

    gg_tc <- df_tc %>%
      rename(analytical = chi, numerical = chi_num) %>%
      pivot_longer(cols = c(analytical, numerical), names_to = "approach", values_to = "chi") %>%
      ggplot(aes(x = tc, y = chi, color = approach, linetype = approach)) +
      geom_line(size = 1.0) +
      labs(title = "Temperature sensitivity", subtitle = paste0("VPD = ", df_tc$vpd[1], " Pa;  CO2 = ", df_tc$co2[1], " ppm;  PPFD = ", df_tc$ppfd[1], " mol m-2 d-1;  p = ", df_tc$patm[1], " Pa"))

    gg_vpd <- df_vpd %>%
      rename(analytical = chi, numerical = chi_num) %>%
      pivot_longer(cols = c(analytical, numerical), names_to = "approach", values_to = "chi") %>%
      ggplot(aes(x = vpd, y = chi, color = approach, linetype = approach)) +
      geom_line(size = 1.0) +
      labs(title = "VPD sensitivity", subtitle = paste0("T = ", df_tc$tc[1], " deg C;  CO2 = ", df_tc$co2[1], " ppm;  PPFD = ", df_tc$ppfd[1], " mol m-2 d-1;  p = ", df_tc$patm[1], " Pa"))

    gg_co2 <- df_co2 %>%
      rename(analytical = chi, numerical = chi_num) %>%
      pivot_longer(cols = c(analytical, numerical), names_to = "approach", values_to = "chi") %>%
      ggplot(aes(x = co2, y = chi, color = approach, linetype = approach)) +
      geom_line(size = 1.0) +
      labs(title = "CO2 sensitivity", subtitle = paste0("VPD = ", df_tc$vpd[1], " Pa;  T = ", df_tc$tc[1], " deg C;  PPFD = ", df_tc$ppfd[1], " mol m-2 d-1;  p = ", df_tc$patm[1], " Pa"))

    gg_ppfd <- df_ppfd %>%
      rename(analytical = chi, numerical = chi_num) %>%
      pivot_longer(cols = c(analytical, numerical), names_to = "approach", values_to = "chi") %>%
      ggplot(aes(x = ppfd, y = chi, color = approach, linetype = approach)) +
      geom_line(size = 1.0) +
      labs(title = "PPFD sensitivity", subtitle = paste0("VPD = ", df_tc$vpd[1], " Pa;  CO2 = ", df_tc$co2[1], " ppm;  T = ", df_tc$tc[1], " deg C;  p = ", df_tc$patm[1], " Pa"))

    gg_patm <- df_patm %>%
      rename(analytical = chi, numerical = chi_num) %>%
      pivot_longer(cols = c(analytical, numerical), names_to = "approach", values_to = "chi") %>%
      ggplot(aes(x = patm, y = chi, color = approach, linetype = approach)) +
      geom_line(size = 1.0) +
      labs(title = "Pressure sensitivity", subtitle = paste0("VPD = ", df_tc$vpd[1], " Pa;  CO2 = ", df_tc$co2[1], " ppm;  PPFD = ", df_tc$ppfd[1], " mol m-2 d-1;  T = ", df_tc$tc[1], " deg C"))

  } else if (show=="assim"){

    gg_tc <- df_tc %>%
      rename(analytical = a_, numerical = aj_num) %>%
      pivot_longer(cols = c(analytical, numerical), names_to = "approach", values_to = "A") %>%
      ggplot(aes(x = tc, y = A, color = approach, linetype = approach)) +
      geom_line(size = 1.0) +
      labs(title = "Temperature sensitivity", subtitle = paste0("VPD = ", df_tc$vpd[1], " Pa;  CO2 = ", df_tc$co2[1], " ppm;  PPFD = ", df_tc$ppfd[1], " mol m-2 d-1;  p = ", df_tc$patm[1], " Pa"))

    gg_vpd <- df_vpd %>%
      rename(analytical = a_, numerical = aj_num) %>%
      pivot_longer(cols = c(analytical, numerical), names_to = "approach", values_to = "A") %>%
      ggplot(aes(x = vpd, y = A, color = approach, linetype = approach)) +
      geom_line(size = 1.0) +
      labs(title = "VPD sensitivity", subtitle = paste0("T = ", df_tc$tc[1], " deg C;  CO2 = ", df_tc$co2[1], " ppm;  PPFD = ", df_tc$ppfd[1], " mol m-2 d-1;  p = ", df_tc$patm[1], " Pa"))

    gg_co2 <- df_co2 %>%
      rename(analytical = a_, numerical = aj_num) %>%
      pivot_longer(cols = c(analytical, numerical), names_to = "approach", values_to = "A") %>%
      ggplot(aes(x = co2, y = A, color = approach, linetype = approach)) +
      geom_line(size = 1.0) +
      labs(title = "CO2 sensitivity", subtitle = paste0("VPD = ", df_tc$vpd[1], " Pa;  T = ", df_tc$tc[1], " deg C;  PPFD = ", df_tc$ppfd[1], " mol m-2 d-1;  p = ", df_tc$patm[1], " Pa"))

    gg_ppfd <- df_ppfd %>%
      rename(analytical = a_, numerical = aj_num) %>%
      pivot_longer(cols = c(analytical, numerical), names_to = "approach", values_to = "A") %>%
      ggplot(aes(x = ppfd, y = A, color = approach, linetype = approach)) +
      geom_line(size = 1.0) +
      labs(title = "PPFD sensitivity", subtitle = paste0("VPD = ", df_tc$vpd[1], " Pa;  CO2 = ", df_tc$co2[1], " ppm;  T = ", df_tc$tc[1], " deg C;  p = ", df_tc$patm[1], " Pa"))

    gg_patm <- df_patm %>%
      rename(analytical = a_, numerical = aj_num) %>%
      pivot_longer(cols = c(analytical, numerical), names_to = "approach", values_to = "A") %>%
      ggplot(aes(x = patm, y = A, color = approach, linetype = approach)) +
      geom_line(size = 1.0) +
      labs(title = "Pressure sensitivity", subtitle = paste0("VPD = ", df_tc$vpd[1], " Pa;  CO2 = ", df_tc$co2[1], " ppm;  PPFD = ", df_tc$ppfd[1], " mol m-2 d-1;  T = ", df_tc$tc[1], " deg C"))

  }



  gg_out <- list(tc = gg_tc, vpd = gg_vpd, co2 = gg_co2, ppfd = gg_ppfd, patm = gg_patm)
  return(gg_out)
}
