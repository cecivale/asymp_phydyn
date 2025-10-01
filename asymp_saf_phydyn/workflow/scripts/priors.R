# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee 
#        V\ Y /V    Get priors 95% interval
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

# Load libraries and functions -------------------------------------------------
library(tidyverse)

get_prior_table <- function(param, distr, meanl, sdl, upper, lower, mean, value, 
                            n_round, n_samples = 25000) {
  if (distr == "LogNormal") {
    t <- tibble(parameter = param,
                distribution = paste0(distr, "(", meanl, ", ", sdl, ")"),
                median = round(qlnorm(0.5, meanl, sdl), n_round),
                density95 = paste0("(", round(qlnorm(0.025, meanl, sdl), n_round)," - ", 
                                   round(qlnorm(0.975, meanl, sdl), n_round), ")"),
                samples =  list(rlnorm(n_samples, meanlog = meanl, sdlog = sdl))
                )
    } else if (distr == "Uniform") {
      t <- tibble(parameter = param,
                  distribution = paste0(distr, "(", lower, ", ", upper, ")"),
                  median = round(qunif(0.5, lower, upper), n_round),
                  density95 = paste0("(", round(qunif(0.025, lower, upper), n_round)," - ", 
                                     round(qunif(0.975, lower, upper), n_round), ")"),
                  samples = list(runif(n_samples, min = lower, max = upper))
                  )
    } else if (distr == "OneOnX") {
      t <- tibble(parameter = param,
                  distribution = paste0(distr, "(", lower, ", ", upper, ")"),
                  #median = round(qunif(0.5, lower, upper), n_round),
                  #density95 = paste0("(", round(qunif(0.025, lower, upper), n_round)," - ", 
                  #                   round(qunif(0.975, lower, upper), n_round), ")")
                  )
    } else if (distr == "Exponential") {
      t <- tibble(parameter = param,
                  distribution = paste0(distr, "(", mean, ")"),
                  median = round(qexp(0.5, mean), n_round),
                  density95 = paste0("(", round(qexp(0.025, mean), n_round)," - ", 
                                     round(qexp(0.975, mean), n_round), ")"),
                  samples = list(rexp(n_samples, rate = mean))
                  )
    } else if (distr == "Fixed") {
      t <- tibble(parameter = param,
                  distribution = paste0(distr, " to ", value))
    } else {
      t <- tibble(parameter = param,
                  distribution = distr)
    }
    return(t)
  }


# Priors model3p_asymp_symp_samplingunif.xml --------------------------

phydyn_asymp_symp_3p <- bind_rows(
  get_prior_table("Origin", "LogNormal", meanl = -0.6, sdl =  0.3, n_round = 3), #origin
  get_prior_table("Origin type", "0.5 - 0.5"), #origin
  get_prior_table("Death rate", "Fixed", value = 36.5), #death
  get_prior_table("Transmission rate", "LogNormal", meanl = log(36.5) - (0.5 * 0.8 * 0.8), sdl =  0.8, n_round = 3), #transmission
  get_prior_table("fsymp", "LogNormal", meanl = 0, sdl =  1.0, n_round = 3), #fsymp
  get_prior_table("pasymp", "Uniform", lower = 0.0, upper =  1, n_round = 3), #pasymp
  get_prior_table("Sampling rate", "LogNormal", meanl = log(0.35579) - (0.5 * 0.2 * 0.2), sdl =  0.8, n_round = 3), #sampling rate
  get_prior_table("Rho sampling", "Uniform", lower = 4.2e-4, upper =  0.0023, n_round = 6), #rho
  get_prior_table("Gamma", "Exponential", mean = 0.5, n_round = 3), #gamma
  get_prior_table("Kappa", "LogNormal", meanl = 1, sdl =  1.25, n_round = 3), #kappa
  get_prior_table("Frequencies", "Empirical")) #frequencies

print(xtable::xtable(phydyn_asymp_symp_3p %>% select(-samples), type = "latex"), file = "results/report/priors_unif.tex")

