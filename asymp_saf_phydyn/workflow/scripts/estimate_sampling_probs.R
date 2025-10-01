# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee Phylodynamics
#        V\ Y /V    Sampling probabilities and rates for phydyn analysis
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

# Load libraries ---------------------------------------------------------------
library(tidyverse)

# Read files
ids <- read_tsv(snakemake@input[["ids"]]) 
active_infections_w <- read_tsv(snakemake@input[["active_infections_w"]]) 
active_infections_d <-  read_tsv(snakemake@input[["active_infections_d"]]) 
  
 
seq_nosymp <- ids %>% filter(deme == "asymp") %>%
  mutate(iso_week = factor(isoweek(date)),
         screening_week = case_when(iso_week == 2 ~ 1,
                                    iso_week == 3 ~ 2,
                                    iso_week == 6 ~ 3,
                                    TRUE ~ NA)) %>%
  count(deme, screening_week) 

# Rho no symptomatics
# 20-29 age category
rho <- left_join(seq_nosymp, active_infections_w, by = "screening_week") %>%
  pivot_wider(names_from = ur, values_from = -c(screening_week, ur, deme, n), names_glue = "{.value}_{ur}") %>%
  mutate(rho = round(n / nosymp_infections_1x, 5), 
         rho_lower = round(n / (nosymp_infections_2.3x * 1.4), 5), # Screening positives / expected active infections # Adding ~38% more for FN https://www.acpjournals.org/doi/10.7326/M20-1495
         rho_upper = round(n / (nosymp_infections_3.1x * 0.5), 5)) %>% # Substracting 50% for FP, very high CT, Walker et al 2021.
  select(screening_week, rho, rho_lower, rho_upper)


# Sampling rate symptomatics
seqs_symp <- ids %>% filter(deme == "symp") %>%
  count(deme)

sampling_rate <- active_infections_d %>%
  filter(stage == "reported", 
         date <= max(ids$date), date >= min(ids$date)) %>%
  summarise(reported_cases = sum(values),
            sp = seqs_symp$n / (reported_cases * 2.7), # mean underreporting
            sp_upper = seqs_symp$n / (reported_cases * 1.0), 
            sp_lower = seqs_symp$n / (reported_cases * 3.1), # max underrerporting
            bur = as.numeric(snakemake@params[["death_rate"]]),
            sr_mean = round(sp * bur, 5),
            sr_upper = round(sp_upper * bur, 5),
            sr_lower = round(sp_lower * bur, 5))

output <- rho %>%
  mutate(rho_lower = min(rho_lower),
         rho_upper = max(rho_upper)) %>%
  pivot_longer(-screening_week, names_to = "param", values_to = "value") %>%
  mutate(param = ifelse(param == "rho", paste0(param, screening_week, "_mean"), param)) %>%
  select(-screening_week) %>% distinct() %>%
  bind_rows(sampling_rate %>% select(sr_mean, sr_upper, sr_lower) %>% 
              pivot_longer(everything(), names_to = "param", values_to = "value"))

write_tsv(output, snakemake@output[["sampling_details"]])

# # debug ------------------------------------------------------------------------
# debugging = TRUE
# setClass(
#   "snakemake_object",
#   contains= "tbl_df",
#   slots = c(input = "character", output = "character", params = "character")
# )
# 
# snakemake <- new("snakemake_object", tibble(),
#                  input = c(ids = "results/data/asymp_symp200_male2029/ids_combined.0.tsv",
#                            active_infections_w = "results/report/active_infections_w.tsv",
#                            active_infections_d = "results/report/active_infections.tsv"),
#                  params = c(death_rate = "36.5"),
#                  output = c(sampling_details ="results_old/report/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029/sampling_details.0.tsv"))




