# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee Phylodynamics
#        V\ Y /V    Take medians value from trace file and create json params
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

library(tidyverse)
library(yaml)
library(jsonlite)
source("./workflow/scripts/talking_to_beast.R")

trace_to_config <- function(trace_file, burnin, median_value = TRUE) {
  trace <- read_trace(trace_file, burnin)
  if (median_value) {
    trace_value <- trace %>% summarise_all(median) 
  } else {
    trace_value <- trace %>% tail(1)
  }
  trace_value <- trace_value %>% 
      select(sim_time = originBDMMPrime,
           f_type2 = f_symp,
           p_type1 = p_asymp,
           birth_type1_0 = birthRateSVi0_asymp,
           birth_type1_1 = birthRateSVi1_asymp,
           birth_type1_2 = birthRateSVi2_asymp,
           birth_type1_3 = birthRateSVi3_asymp,
           birth_changetime1 = birthRateSVi0_endtime,
           birth_changetime2 = birthRateSVi1_endtime,
           birth_changetime3 = birthRateSVi2_endtime,
           sampling_type2 = samplingRateSVi1_symp,
           sampling_changetime = samplingRateSVi0_endtime,
           rho_type1_0 = rhoSamplingSVe0asymp,
           rho_type1_1 = rhoSamplingSVe1asymp,
           rho_type1_2 = rhoSamplingSVe2asymp,
           rho_time1 = rhoSamplingSVe0_time,
           rho_time2 = rhoSamplingSVe1_time,
           rho_time3 = rhoSamplingSVe2_time,
           death_type1 = deathRateSVasymp,
           death_type2 = deathRateSVsymp
    )
}


trace <- snakemake@params[["trace"]]
sampling <- read_tsv(snakemake@params[["sampling_details"]]) %>%
  pivot_wider(names_from = param)

burnin <- snakemake@params[["burnin"]]

params_file <- snakemake@output[["xml_paramsf"]]

# XML params
if (length(trace) > 1) {
  # From several trace files, take the median
  config_trace_l <- lapply(trace, function (t) trace_to_config(t, burnin))
  config_trace <- bind_rows(config_trace_l) %>% summarise_all(median) 
  } else {
    config_trace <- trace_to_config(trace, burnin, median_value = snakemake@params[["median"]]) 
  }

config_trace <- config_trace %>%
  bind_cols(sampling)

# Write to json for beast
config_json = toJSON(unbox(fromJSON(toJSON(config_trace))), pretty = T)
write(config_json, params_file)

