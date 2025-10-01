# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee Phylodynamics
#        V\ Y /V    Script to get number of active infections in the community 
#    (\   / - \     from BAG data during study period
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


# Functions --------------------------------------------------------------------

get_deconvolved_incidence <- function(reported_cases, delay_list, 
                                      to_future = FALSE, smooth = TRUE) {
  # estimateR deconvolution from reported cases based on delay
  
  if (to_future) {
    reported_cases <- reported_cases %>%
      mutate(date_original = date,
             date = sort(ymd(date_original), decreasing = T)) %>%
      arrange(date)
  }
  
  if (smooth) {
    infections <- get_infections_from_incidence(
      incidence_data = reported_cases$entries,
      delay = delay_list
    ) 
  } else {
    infections <- deconvolve_incidence(
    incidence_data = reported_cases$entries,
    delay = delay_list
    )
  }

  df_new_infections <- infections %>%
    as_tibble() 
  
  if (to_future) {
    df_new_infections <- df_new_infections %>%
      mutate(date = reported_cases$date_original - index_offset) %>%
      select(-index_offset)
  } else {
    df_new_infections <- df_new_infections %>%
      mutate(date = reported_cases$date + index_offset) %>%
      select(-index_offset)
  }
  
  return(df_new_infections)
}

get_delay_gamma_distr <- function(mean_gamma, sd_gamma) {
  # gamma distribution parameters 
  shape_gamma <- mean_gamma^2 / (sd_gamma^2)
  scale_gamma <- (sd_gamma^2) / mean_gamma
  
  delay <- list(name ="gamma", 
                shape = shape_gamma, 
                scale = scale_gamma)
  
  return(delay)
}


# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(estimateR)

# Load BAG datasets ------------------------------------------------------------
cases_bag <- read_tsv(snakemake@input[["cases_byage_d"]])
cases_bag_bysex <- read_tsv(snakemake@input[["tests_bysex_w"]])

# Reported cases 20-29 age group (Issue: army samples also include 18-19 years old)
cases_2029 <- cases_bag %>%
  filter(geoRegion == "CH", ageRange == "20 - 29")  

# Compute proportion of reported cases 20-29 from males
p_cases_bysex <- cases_bag_bysex %>%
  filter(geoRegion == "CH", 
         datum %in% c(202102, 202103, 202106)) %>%
  group_by(sex) %>%
  summarise(entries = sum(entries),
            entries_pos = sum(entries_pos)) %>%
  mutate(p = entries_pos/sum(entries_pos))

cases_2029_male <- cases_2029 %>%
  mutate(entries = entries * p_cases_bysex[p_cases_bysex$sex == "male",]$p)  


# Delay distributions ---------------------------------------------------------
# From Re paper https://elifesciences.org/articles/71345#table2
# Infection to onset of symptoms	M = 5.3,	S = 3.2	(Linton et al., 2020)
# Onset of symptoms to case confirmation	M = 5.5,	S = 3.8	(Bi et al., 2020)
# https://www.nature.com/articles/s41579-022-00822-w#MOESM1
# RT PCR positivity and symptoms (from https://www.nature.com/articles/s41579-022-00822-w)
# 0. Positive any stage: From day 2 after infection to day 19.3 (14 days after symptom onset)
# 1. Presymptomatics: From 2 days after infection (3.3 days before symptom onset) to symptom onset 5.3 days after infection
# 2. Symptomatic: From 5.3 days after infection to 15.3 days after infection (10 days after symptom onset)
# 3. Postsymptomatics (but still positive): From 15.3 days after infection (10 days after symptom onset) to 19.3 days after infections (14 days after symptom onset)

delay_symptoms_start_to_infection <- get_delay_gamma_distr(mean_gamma = 5.3, sd_gamma = 3.2)
delay_symptoms_start_to_positive_start <- get_delay_gamma_distr(mean_gamma = 3.3, sd_gamma = 3.2)
delay_symptoms_start_to_infectious_start <- get_delay_gamma_distr(mean_gamma = 2.3, sd_gamma = 3.2)

delay_report_to_symptoms_start <- get_delay_gamma_distr(mean_gamma = 5.5, sd_gamma = 3.8)

delay_symptoms_end_to_report <- get_delay_gamma_distr(mean_gamma = 4.5, sd_gamma = 3.8)
delay_positive_end_to_report <- get_delay_gamma_distr(mean_gamma = 8.5, sd_gamma = 3.8)
delay_infectious_end_to_report <- get_delay_gamma_distr(mean_gamma = 2.5, sd_gamma = 3.8)


# Deconvolution of reported cases to each 
# of the start and end of stages of infections----------------------------------

# New infections
df_infections_start <- get_deconvolved_incidence(cases_2029_male, 
                                                 list(delay_symptoms_start_to_infection, 
                                                      delay_report_to_symptoms_start), 
                                                 smooth = T) 

# Symptomatic
df_symptoms_start <- get_deconvolved_incidence(cases_2029_male, 
                                               delay_report_to_symptoms_start, 
                                               smooth = T)
df_symptoms_end <- get_deconvolved_incidence(cases_2029_male, 
                                             delay_symptoms_end_to_report, 
                                             to_future = T, smooth = T) 
# Positive to RT-PCR test
df_positive_start <- get_deconvolved_incidence(cases_2029_male, 
                                               list(delay_symptoms_start_to_positive_start, 
                                                    delay_report_to_symptoms_start), 
                                               smooth = T)
df_positive_end <- get_deconvolved_incidence(cases_2029_male, 
                                             delay_positive_end_to_report, 
                                             to_future = T, smooth = T)

# Infectious
df_infectious_start <- get_deconvolved_incidence(cases_2029_male, 
                                                 list(delay_symptoms_start_to_infectious_start, 
                                                      delay_report_to_symptoms_start),
                                                 smooth = T)
df_infectious_end <- get_deconvolved_incidence(cases_2029_male, 
                                               delay_infectious_end_to_report, 
                                               to_future = T, smooth = T)

df_stages_new <- bind_rows(df_infections_start %>% mutate(stage = "infections_start"),
                       df_positive_start %>% mutate(stage = "positive_start"),
                       df_infectious_start %>% mutate(stage = "infectious_start"),
                       df_symptoms_start %>% mutate(stage = "symptoms_start"),
                       df_infectious_end %>% mutate(stage = "infectious_end"),
                       df_symptoms_end %>% mutate(stage = "symptoms_end"),
                       df_positive_end %>% mutate(stage = "positive_end"),
                       cases_2029_male %>% select(values = entries, date) %>% 
                         mutate(stage = "reported_end"))

# Get number of infections active at each stage --------------------------------
# We sum over the new infections in each stage and substract the ones that leave
# that stage.

df_stages <- df_stages_new %>%
  pivot_wider(values_from = "values", names_from = "stage") %>%
  replace(is.na(.), 0) %>%
  mutate(symptoms_new = symptoms_start - symptoms_end,
         symptoms_active = cumsum(symptoms_new),
         presymp_new = positive_start - symptoms_start,
         presymp_active = cumsum(presymp_new),
         positive_new = positive_start - positive_end,
         positive_active = cumsum(positive_new),
         infectious_new = infectious_start - infectious_end,
         infectious_active = cumsum(infectious_new),
         possymp_new = symptoms_end - positive_end,
         possymp_active = cumsum(possymp_new)) %>%
  pivot_longer(infections_start:possymp_active, names_to = c("stage", "time"), 
               names_sep = "_", values_to = "values")

# Save dataframe
write_tsv(df_stages, snakemake@output[["active_infections"]])


# debug ------------------------------------------------------------------------
# setClass(
#   "snakemake_object",
#   contains= "tbl_df",
#   slots = c(input = "character", output = "character", params = "character")
# )
# 
# snakemake <- new("snakemake_object", tibble(),
#                  input = c(screening_weeks = "resources/screening_weeks.tsv",
#                            cases_byage_d = "resources/ext_cases_byage_d.tsv",
#                            tests_bysex_w = "resources/ext_tests_bysex_w.tsv"))


