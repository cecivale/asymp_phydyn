# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee Phylodynamics
#        V\ Y /V    Script for descriptive numbers and figures in first section 
#    (\   / - \     of results
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


library(tidyverse)
library(sf)
library(patchwork)
library(xtable)
library(ggpattern)

set.seed(24)
source("workflow/scripts/plot_opts.R")
source("https://raw.githubusercontent.com/cecivale/talking-to-lapis/refs/heads/main/R/lapis_functions.R")

debugging = FALSE


# This script depends on internal datasets:
#  - Full metadata (from load_data_vineyard.R)
metadata <- read_tsv(snakemake@input[["metadata"]]) 
#  - Screening weeks file 
screening_weeks <- read_tsv(snakemake@input[["screening_weeks"]])
#  - Active infections dataset (calculated in estimateR_infections.R)
active_infections <-  read_tsv(snakemake@input[["active_infections"]])

# and external datasets:
#  - Population
population <- read_tsv(snakemake@input[["population"]]) 
#  - BAG tests: Weekly and daily tests by age, weekly tests by sex
tests_byage_w <- read_tsv(snakemake@input[["tests_byage_w"]])
tests_bysex_w <- read_tsv(snakemake@input[["tests_bysex_w"]])
cases_byage_d <- read_tsv(snakemake@input[["cases_byage_d"]])

# Study period 
study_from <- min(screening_weeks$from) - 10
study_to <- max(screening_weeks$to) + 10
study_period <- interval(ymd(study_from), ymd(study_to))

# ------------------------------------------------------------------------------
# (A) Calculate estimates of number of infections symptomatic and asymptomatic
# in Swiss men 20-29 age group using reported cases and screening positivity
# ------------------------------------------------------------------------------

# 1. Screening positivity ------------------------------------------------------
# We divide the number of positive tests in screening weeks 1, 2, and 3  by the 
# number of total tests reported by the Swiss Armed forces.
# Two approximate values for total tests: from emails and from PDF

#   a. Read files
metadata_screening <- metadata %>%
  filter(army, screening, sim_date) %>%
  mutate(deme = "screening") 

metadata_community <- metadata %>%
  filter(!army, 
         #sex == "Männlich", age_cat10 == "20-29", 
         sim_date, date >= study_from, date <= study_to) %>%
  mutate(deme = "community") 

nrow(metadata_screening)
nrow(metadata_screening %>% filter(qc))

#   b. Calculations
screening_table1_byweek <- metadata_screening %>% 
  count(screening_week, name = "positive_tests") %>%
  left_join(metadata_screening %>% filter(qc) %>%
              count(screening_week, name = "sequences")) %>%
  left_join(screening_weeks, join_by(screening_week)) %>% 
  mutate(positivity = positive_tests / screening_tests,
         positivity_pdf = positive_tests / screening_tests_pdf)

# We assume the positivity calculated for 18-30 years old to be valid for the 20-29 age category
# We cannot calculate positivity 20-29 men because we do not have the data of how many tests were
# done for that demographic group, only the total value.

# Percentage of mens 20-29
sexage_p <- metadata_screening %>% count(sex, age_cat10) %>% mutate(p = n/sum(n))
metadata_screening %>% count(sex) %>% mutate(n/sum(n))
metadata_screening %>% count(age_cat10) %>% mutate(n/sum(n))
metadata_screening %>% filter(sex == "Männlich") %>% count(age_cat10) %>% mutate(n/sum(n))

# In the following calculations, we will use the positivity values based on the 
# test numbers reported in the emails. This is because the PDF may include other 
# tests not related to the screening of new recruits. Moreover, the positivity 
# values from the emails are more consistent — around 2.7% — and we do not expect 
# large fluctuations from week to week.

# Table 1 ----------------------------------------------------------------------
# Two columns: SFA Screening and Community
# Six Rows: RT-PCR tests, Positive RT-PCR tests, RT-PCR positivity, 
#           CT value, Sequences (for analysis), % men 20-29 age group

# Screening
screening_table1 <- screening_table1_byweek %>%
  summarise(category = "SAF screening",
            `RT-PCR tests` = sum(screening_tests),
            `Positive RT-PCR tests` = sum(positive_tests),
            `RT-PCR positivity` = paste0(round(`Positive RT-PCR tests`/`RT-PCR tests` * 100, 2), "%"),
            `% men 20-29 age group` = paste0(sexage_p %>% filter(sex == "Männlich", age_cat10 == "20-29") %>% pull(p) * 100, "%"),
            `Viral Genetic Sequences` = sum(sequences),
  ) %>% bind_cols(metadata_screening %>% 
                    summarise(median_ct = median(ct, na.rm = TRUE),
                              low_ct25 = quantile(ct, 0.25, na.rm = TRUE),
                              high_ct75 = quantile(ct, 0.75, na.rm = TRUE),
                              `CT value (IQR)` = paste0(median_ct, " (", 
                                                        low_ct25, "-", high_ct75, " )"))) %>%
  select(1:4, 10, 5:6)

screening_m2029_table1 <- metadata_screening %>% 
  filter(sex == "Männlich", age_cat10 == "20-29")  %>% 
  count(name = "Positive RT-PCR tests") %>%
  mutate(`% men 20-29 age group` = "100%") %>%
  mutate(category = "SAF screening men 20-29") %>%
  bind_cols(metadata_screening %>% 
              filter(sex == "Männlich", age_cat10 == "20-29", qc) %>%
              count(name = "Viral Genetic Sequences")) %>%
  bind_cols(metadata_screening %>% 
              filter(sex == "Männlich", age_cat10 == "20-29")  %>% 
              summarise(median_ct = median(ct, na.rm = TRUE),
                        low_ct25 = quantile(ct, 0.25, na.rm = TRUE),
                        high_ct75 = quantile(ct, 0.75, na.rm = TRUE),
                        `CT value (IQR)` = paste0(median_ct, " (", 
                                                  low_ct25, "-", high_ct75, " )"))) %>%
  select(1:4, 8)
#left_join(screening_weeks, join_by(screening_week)) %>% 
#mutate(positivity = positive_tests / screening_tests,
#       positivity_pdf = positive_tests / screening_tests_pdf)

# Community
p_community <- cross_join(tests_bysex_w %>% 
                            filter(datum >= 202101, datum <= 202108) %>% 
                            group_by(sex) %>% summarise(entries = sum(entries),
                                                        entries_pos = sum(entries_pos)) %>%
                            ungroup() %>%
                            mutate(sex_p = entries/sum(entries),
                                   sex_pp = entries_pos/sum(entries_pos)) %>% 
                            select(sex, sex_p, sex_pp),
                          tests_byage_w %>% 
                            filter(datum >= 202101, datum <= 202108) %>% 
                            group_by(altersklasse_covid19) %>% 
                            summarise(entries = sum(entries),
                                      entries_pos = sum(entries_pos)) %>%
                            ungroup() %>%
                            mutate(age_p = entries/sum(entries),
                                   age_pp = entries_pos/sum(entries_pos)) %>% 
                            select(altersklasse_covid19, age_p, age_pp)) %>%
  mutate(p = sex_p*age_p, pp = sex_pp*age_pp) %>%
  filter(sex == "male", altersklasse_covid19 == "20 - 29")

community_table1 <- bind_cols(
  tibble(category = "Community"),
  tests_bysex_w %>% 
    filter(datum >= 202101, datum <= 202108) %>%
    summarise(`RT-PCR tests` = sum(entries),
              `Positive RT-PCR tests`  = sum(entries_pos),
              `RT-PCR positivity` = paste0(round(`Positive RT-PCR tests`/`RT-PCR tests` * 100, 2), "%")),
  p_community %>% mutate(`% men 20-29 age group` = paste0(round(pp * 100, 2), "%")) %>% select(`% men 20-29 age group`),
  metadata_community %>% filter(qc) %>% count(name = "Viral Genetic Sequences"),
  metadata_community %>% summarise(median_ct = median(ct, na.rm = TRUE),
                                   low_ct25 = quantile(ct, 0.25, na.rm = TRUE),
                                   high_ct75 = quantile(ct, 0.75, na.rm = TRUE),
                                   `CT value (IQR)` = paste0(median_ct, " (", 
                                                             low_ct25, "-", high_ct75, ")")) %>% select(4))

community_m2029_table1 <- bind_cols(
  tibble(category = "Community men 20-29"),
  tests_byage_w %>% 
    filter(datum >= 202101, datum <= 202108, altersklasse_covid19 == "20 - 29") %>%
    summarise(`RT-PCR tests` = round(sum(entries) * p_community$sex_p),
              `Positive RT-PCR tests`  = round(sum(entries_pos) * p_community$sex_pp),
              `RT-PCR positivity` = paste0(round(`Positive RT-PCR tests`/`RT-PCR tests` * 100, 2), "%")),
  tibble(`% men 20-29 age group` = "100%"),
  metadata_community %>% filter(qc, sex == "Männlich", age_cat10 == "20-29", ) %>% count(name = "Viral Genetic Sequences"),
  metadata_community %>% filter(qc, sex == "Männlich", age_cat10 == "20-29", ) %>% 
    summarise(median_ct = median(ct, na.rm = TRUE),
              low_ct25 = quantile(ct, 0.25, na.rm = TRUE),
              high_ct75 = quantile(ct, 0.75, na.rm = TRUE),
              `CT value (IQR)` = paste0(median_ct, " (", 
                                        low_ct25, "-", high_ct75, ")")) %>% select(4))

community_m_table1 <- bind_cols(
  tibble(category = "Community men"),
  tests_bysex_w %>% 
    filter(datum >= 202101, datum <= 202108, sex == "male") %>%
    summarise(`RT-PCR tests` = sum(entries),
              `Positive RT-PCR tests`  = sum(entries_pos),
              `RT-PCR positivity` = paste0(round(`Positive RT-PCR tests`/`RT-PCR tests` * 100, 2), "%")),
  p_community %>% mutate(`% men 20-29 age group` = paste0(round(age_pp * 100, 2), "%")) %>% select(`% men 20-29 age group`),
  metadata_community %>% filter(qc, sex == "Männlich") %>% count(name = "Viral Genetic Sequences"),
  metadata_community %>% filter(qc, sex == "Männlich") %>% 
    summarise(median_ct = median(ct, na.rm = TRUE),
              low_ct25 = quantile(ct, 0.25, na.rm = TRUE),
              high_ct75 = quantile(ct, 0.75, na.rm = TRUE),
              `CT value (IQR)` = paste0(median_ct, " (", 
                                        low_ct25, "-", high_ct75, ")")) %>% select(4))

community_2029_table1 <- bind_cols(
  tibble(category = "Community 20-29"),
  tests_byage_w %>% 
    filter(datum >= 202101, datum <= 202108, altersklasse_covid19 == "20 - 29") %>%
    summarise(`RT-PCR tests` = sum(entries),
              `Positive RT-PCR tests`  = sum(entries_pos),
              `RT-PCR positivity` = paste0(round(`Positive RT-PCR tests`/`RT-PCR tests` * 100, 2), "%")),
  p_community %>% mutate(`% men 20-29 age group` = paste0(round(sex_pp * 100, 2), "%")) %>% select(`% men 20-29 age group`),
  metadata_community %>% filter(qc, age_cat10 == "20-29", ) %>% count(name = "Viral Genetic Sequences"),
  metadata_community %>% filter(qc, age_cat10 == "20-29", ) %>% 
    summarise(median_ct = median(ct, na.rm = TRUE),
              low_ct25 = quantile(ct, 0.25, na.rm = TRUE),
              high_ct75 = quantile(ct, 0.75, na.rm = TRUE),
              `CT value (IQR)` = paste0(median_ct, " (", 
                                        low_ct25, "-", high_ct75, ")")) %>% select(4))

table1 <- bind_rows(screening_table1, screening_m2029_table1, 
                    #community_2029_table1, community_m_table1,
                    community_table1, community_m2029_table1) %>% 
  mutate_all(as.character) %>% 
  pivot_longer(-category, names_to = "var", values_to = "value") %>%
  pivot_wider(names_from = category)

if (!debugging) { print(xtable(table1, type = "latex"), file = snakemake@output[["table1"]]) }

# 2. Number of no symptomatic active infections in each week -------------------
# Calculated as screening positivity x number of individuals (men, 20-29) without symptoms
#   a. Number of men 20-29 individuals without symptoms = 
#                  Number of Swiss men 20-29 individuals - Number of men 20-29 with symptoms (i.e. that would not have gone to the recruitment)

#     a.1. Number of Swiss men 20-29 individuals
pop_men2029 <- population %>%
  filter(geoRegion == "CH", ageRange == "20 - 29", sex == "male") %>% pull(pop)

#       a.1.1 Percentage pop in the screening (we multiply by 0.75 because 
#             population value is for men 2029 but number of tests includes other age categories and sex)
0.75 * as.numeric(table1[[1,2]]) / pop_men2029

#     a.2. Number of Swiss men 20-29 individuals with symptoms
# From  symptom onset 5.3 days after infection to end of isolation at 15.3 days after infection  (10 days after symptom onset) 
# calculated in estimateR_infections.R script   

active_infections_d <- active_infections %>%
  mutate(iso_week = factor(isoweek(date)),
         screening_week = case_when(iso_week == 2 & year(date) == 2021 ~ 1,
                                    iso_week == 3 & year(date) == 2021 ~ 2,
                                    iso_week == 6 & year(date) == 2021 ~ 3,
                                    TRUE ~ NA)) 

positivity <- as.numeric(gsub("%", "", table1[[3,2]])) / 100

active_infections_w <- active_infections_d %>%
  filter(!is.na(screening_week)) %>%
  group_by(screening_week, stage, time) %>%
  #group_by(stage, time) %>%
  summarise(mean_value_w = mean(values), .groups = "drop") %>%
  # Compute 20-29 year old mean population size without symptoms each screening week
  filter(stage %in% c("symptoms", "presymp", "positive"), time == "active") %>%
  pivot_wider(names_from = c(stage, time), values_from = mean_value_w)

active_reported <- active_infections_w %>%
  mutate(nosymp_pop = pop_men2029 - symptoms_active,
         nosymp_infections = nosymp_pop * positivity,
         asymp_active = nosymp_infections - presymp_active,
         positive_active_nopre = positive_active - presymp_active) 

# 3. Under-reporting -----------------------------------------------------------
# Based on Geneva seroprevalence study https://europepmc.org/article/pmc/pmc8063076)
# We calculate based on our estimates of asymptomatic cases the under-reporting this will imply:

active_ur2.3 <- active_infections_w %>%  
  mutate(across(-screening_week, ~ .x * 2.3),
         nosymp_pop = pop_men2029 - symptoms_active,
         nosymp_infections = nosymp_pop * positivity,
         asymp_active = nosymp_infections - presymp_active,
         positive_active_nopre = positive_active - presymp_active)

active_ur2.7 <- active_infections_w %>%  
  mutate(across(-screening_week, ~ .x * 2.7),
         nosymp_pop = pop_men2029 - symptoms_active,
         nosymp_infections = nosymp_pop * positivity,
         asymp_active = nosymp_infections - presymp_active,
         positive_active_nopre = positive_active - presymp_active)

active_ur3.1 <- active_infections_w %>%  
  mutate(across(-screening_week, ~ .x * 3.1),
         nosymp_pop = pop_men2029 - symptoms_active,
         nosymp_infections = nosymp_pop * positivity,
         asymp_active = nosymp_infections - presymp_active,
         positive_active_nopre = positive_active - presymp_active)

active_w_ur <- bind_rows(
  active_reported %>% mutate(ur = "1x"),
  active_ur2.3 %>% mutate(ur = "2.3x"),
  active_ur2.7 %>% mutate(ur = "2.7x"),
  active_ur3.1 %>% mutate(ur = "3.1x"))

# Proportion of infections being asymptomatic
active_w_ur %>% group_by(ur) %>%
  summarise(median(asymp_active/(asymp_active + positive_active)))

# Number of symp per asymp infection
active_w_ur %>% group_by(ur) %>%
  summarise(median(asymp_active/positive_active))

# Proportion of non symptomatic infections being pre symptomatic
active_w_ur %>% group_by(ur) %>%
  summarise(median(presymp_active/nosymp_infections))

# Prevalence of asymp
active_w_ur %>% group_by(ur) %>%
  summarise(median(asymp_active/pop_men2029))

# Prevalence of symp
active_w_ur %>% group_by(ur) %>%
  summarise(median(positive_active/pop_men2029))

# Interestingly, the screening positivity for the last week is as high as the other weeks, while reported cases are almost half. Why?
# a) Could it be that there is less reporting at the end of the wave? People think that wave is finishing and go less for testing. Check test data 
tests_byage_w   %>%
  filter(altersklasse_covid19 == "20 - 29", datum %in% c(202102, 202103, 202106), geoRegion == "CH")

# It does not look like this is the case, same number of tests more or less were done each week.

# b) Could it be some right censoring at the end of the wave? We see less reported but still many asymp. 
# Is that a signal of effective contact tracing? i.e. symptomatic cases are detected in earlier stages of infection.
# c) Other option is that our screening positivity estimate is too high:
# - They did more tests than said in the email
# - Or non random mixing, so cluster of cases increasing positivity rate

# Save numbers
if (!debugging) { write_tsv(active_w_ur, snakemake@output[["active_infections_w"]]) }


# ------------------------------------------------------------------------------
# (B) Plots for Figure 1
# ------------------------------------------------------------------------------

# 1. Reported cases and estimate in 2029 men group ----------------------------
cases_d <- cases_byage_d %>%
  filter(geoRegion == "CH") %>%
  group_by(date) %>%
  summarise(cases = sum(entries),
            group = "all") %>%
  bind_rows(cases_byage_d %>%
              filter(geoRegion == "CH", ageRange == "20 - 29") %>%
              group_by(date) %>%
              summarise(cases = sum(entries) * p_community$sex_pp,
                        group = "m2029")) %>%
  mutate(study = date %within% study_period)

p1 <- ggplot() +
  geom_bar(data = cases_d %>% filter(date <= study_to + 5, date >= study_from - 5),
           aes(x = date, y = cases, fill = group, alpha = study, size = group), 
           color = "white", stat = "identity") +
  # geom_vline(aes(xintercept = ymd(c(study_from, study_to))),
  #            color = pal_asympsymp[["symp2"]], linetype = 1 ) +
  theme_minimal() +
  scale_y_continuous("Reported SARS-CoV-2 cases") +
  scale_x_date(limits = ymd(c("2020-12-25","2021-03-05"))) +
  scale_fill_manual(name = "", values = c("grey80", pal_asympsymp[["symp"]]),
                    labels = c("Swiss population", "Swiss men 20-29 group (estimated)")) +
  scale_alpha_manual(values = c(0.3, 1), guide = "none") +
  scale_size_manual(values = c(0.15, 0.25), guide = "none") +
  xlab("") +
  theme_asympsymp() +
  theme(legend.position = "top")

p1

# Alpha variant
alpha_p <- read_csv(snakemake@input[["covspectrum_alpha"]])
p2 <- ggplot(alpha_p) +
  geom_line(aes(date, proportion)) +
  geom_ribbon(aes(date, ymin = proportionCILow, ymax = proportionCIHigh), alpha = 0.3) +
  scale_x_date(limits = ymd(c("2020-12-25","2021-03-05")), name = "") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent) +
  theme_asympsymp()
p2

# 3. Screening tests and sequences over time, positivity? ----------------------
# Alpha screening sequences

gisaid_ids <- metadata_screening %>% pull(gisaid_id)

alpha_seqs <- lapis_query_by_id(database = "open",
                                endpoint = "details",
                                samples = gisaid_ids,
                                sample_id = "gisaidEpiIsl",
                                filter = lapis_filter(nextcladeQcOverallScoreFrom = 0,
                                                      nextcladeQcOverallScoreTo = 100,
                                                      variantQuery = "B.1.1.7*")) %>% 
  select(seq_id = genbankAccession, gisaidEpiIsl, genbankAccession, sraAccession, strain,
         date, division, pangoLineage,	nextstrainClade, 
         originatingLab,  database) %>% distinct()


p3 <- metadata_screening %>% mutate(alpha_var = gisaid_id %in% alpha_seqs$gisaidEpiIsl,
                                    qc = ifelse(gisaid_id %in% alpha_seqs$gisaidEpiIsl, T, qc)) %>%
  ggplot() +
  geom_bar(aes(order_date, fill = interaction(qc, alpha_var)), color = "white", size = 0.15) +
  geom_vline(data = metadata_screening %>% filter(qc) %>% group_by(screening_week_bag) %>% summarise(med = median(date)),
             aes(xintercept = med), linetype = 2, linewidth = 0.4) +
  scale_x_date(limits = ymd(c("2020-12-25","2021-03-05")), 
               date_labels = "%b %d", name = "") +
  ylab("Positive RT-PCR\nSAF Screening samples") +
  scale_fill_manual(values = c(pal_asympsymp[["asymp2"]], pal_asympsymp[["asymp"]], pal_asympsymp[["asymp"]]),
                    labels = c("Not sequenced", "Sequenced", "Alpha Sequence"), name = "") +
  theme_asympsymp() +
  theme(legend.position = "bottom")

p3

# 4. Demographics

p4 <- ggplot(bind_rows(metadata_community, metadata_screening) %>% 
               filter(sex != "Unbekannt", qc) %>%
               count(sex, deme, age_cat10) %>% #, qc) %>%
               group_by(deme) %>%
               mutate(p = n/sum(n),
                      p = ifelse(deme == "screening", p, -p))) +
  geom_col_pattern(aes(p, age_cat10, pattern = sex, group = sex, fill = interaction(deme)), #interaction(deme,qc))
                   position = "dodge", width = 1, color = "white",
                   pattern_colour = "white",       # Stripe color
                   pattern_fill = "white",
                   pattern_density = 0.35,         # High density = more lines
                   pattern_spacing = 0.02,        # Tight spacing
                   pattern_angle = 45,             # Vertical stripes
                   pattern_size = 0.01  ) +
  scale_fill_manual(values = c(#pal_asympsymp[["symp2"]], pal_asympsymp[["asymp2"]], 
    pal_asympsymp[["symp"]], pal_asympsymp[["asymp"]]), name = "") +
  scale_pattern_manual(values = c("none", "stripe"), name = "") +
  ylab("Sex and Age category") +
  xlab("Sequences (%)") +
  theme_asympsymp() +
  theme(legend.position = "top")

p4

# 5. Comparison of a/pre/symptomatic infections --------------------------------
p5 <- bind_rows(
  active_reported %>% select(positive_active_nopre, presymp_active, asymp_active) %>%
    pivot_longer(everything(), values_to = "value", names_to = "var") %>%
    mutate(name = "1x"),
  active_ur2.3 %>% select(positive_active_nopre, presymp_active, asymp_active) %>%
    pivot_longer(everything(), values_to = "value", names_to = "var") %>%
    mutate(name = "2.3x"),
  active_ur2.7 %>% select(positive_active_nopre, presymp_active, asymp_active) %>%
    pivot_longer(everything(), values_to = "value", names_to = "var") %>%
    mutate(name = "2.7x"),
  active_ur3.1 %>% select(positive_active_nopre, presymp_active, asymp_active) %>%
    pivot_longer(everything(), values_to = "value", names_to = "var") %>%
    mutate(name = "3.1x")) %>%
  ggplot() +
  geom_bar(aes(value, name, fill = factor(var, levels = c("positive_active_nopre", "presymp_active", "asymp_active" ))),
           stat = "identity", position = "fill",
           color = "white", width = 1) +
  scale_fill_manual(values = c(pal_asympsymp[["symp"]], pal_asympsymp[["symp2"]], pal_asympsymp[["asymp2"]]), name = "") +
  scale_y_discrete(labels = c("None\n(1x)", "Low\n(2.3x)", "Medium\n(2.7x)", "High\n(3.1x)")) +
  scale_x_continuous(labels = scales::percent) +
  ylab("Symptomatic Underreporting Scenario") +
  xlab("Proportion of Infected Population (%)") +
  theme_asympsymp() +
  theme(legend.position = "bottom")

p5


# 6. CT values plot ------------------------------------------------------------
ct_df <- bind_rows(metadata_community %>% mutate(group = "All Community"), 
                   metadata_screening %>% mutate(group = "SAF Screening"),
                   metadata_community %>% filter(age_cat10 == "20-29", sex == "Männlich") %>% mutate(group = "Community\nmen 20-29"), 
                   #metadata_screening %>% filter(age_cat10 == "20-29", sex == "Männlich") %>% mutate(group = "screening_men2029")) 
) %>% filter(!is.na(ct)) %>%
  select(sample_number, qc, ct, deme, age_cat10, sex, group)

test_results <- ggpubr::compare_means(ct ~ group, ct_df,  p.adjust.method = "BH")
my_comparisons <- list( c("All Community", "SAF Screening"), c("All Community", "Community\nmen 20-29"), c("SAF Screening", "Community\nmen 20-29") )

p6 <- ct_df  %>%
  ggplot(aes(group, ct, color = group)) +
  geom_violin(aes(fill = group), alpha = 0.4, position = position_dodge(0.9), color = "transparent") +
  geom_boxplot(aes(color = group, fill = group), width = 0.1,  position = position_dodge(0.9),
               outlier.shape = NA, alpha = 0.6) +
  scale_fill_manual(values = c(pal_asympsymp[["grey"]], pal_asympsymp[["symp"]], pal_asympsymp[["asymp"]])) +
  scale_color_manual(values = c("grey50", pal_asympsymp[["symp"]], pal_asympsymp[["asymp"]])) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, size = 2.5) + 
  theme_asympsymp() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("") +
  xlab("RT-PCR test cycle threshold (CT) value")

p6

# # Subsampled sequences 
# s1 <- read_tsv("results/data/asymp_symp200_men2029/symp/ids_subsampled.0.tsv")
# s2 <- read_tsv("results/data/asymp_symp200_men2029/symp/ids_subsampled.1.tsv")
# s3 <- read_tsv("results/data/asymp_symp200_men2029/symp/ids_subsampled.2.tsv")
# s4 <- read_tsv("results/data/asymp_symp200_men2029/symp/ids_subsampled.3.tsv")
# s5 <- read_tsv("results/data/asymp_symp200_men2029/symp/ids_subsampled.4.tsv")
# seqs_asymp <- read_tsv("results/data/asymp_symp200_men2029/asymp/ids.tsv")
# 
# subsample_metadata <- metadata %>%
#   mutate(s1 = case_when(seq_id %in% s1$seq_id ~ TRUE, TRUE ~ FALSE),
#          s2 = case_when(seq_id %in% s2$seq_id ~ TRUE, TRUE ~ FALSE),
#          s3 = case_when(seq_id %in% s3$seq_id ~ TRUE, TRUE ~ FALSE),
#          s4 = case_when(seq_id %in% s4$seq_id ~ TRUE, TRUE ~ FALSE),
#          s5 = case_when(seq_id %in% s5$seq_id ~ TRUE, TRUE ~ FALSE),
#          s = case_when(seq_id %in% seqs_asymp$seq_id ~ TRUE, TRUE ~ FALSE)) %>%
#   filter(s1 | s2 | s3 | s4 | s5 | s)



# # 4. Screening tests and sequences by canton, map ------------------------------

# Population by canton 
population_canton <- population %>%
  group_by(geoRegion) %>%
  summarise(pop = sum(pop))
population_canton_men2029 <- population %>%
  filter(sex == "male", ageRange == "20 - 29")

# Reported cases by canton
cases_canton <- tests_byage_w %>%
  filter(!str_detect(geoRegion, "CH|CHFL|FL"),
         datum >= sprintf("%d%02d", isoyear(study_from), 
                          isoweek(study_from)),
         datum <= sprintf("%d%02d", isoyear(study_to), 
                          isoweek(study_to))) %>%
  group_by(geoRegion) %>%
  summarise(cases = sum(entries_pos)) %>%
  left_join(population_canton) %>%
  mutate(cases_norm = cases/pop * 100000)

cases_canton_m2029 <- tests_byage_w %>%
  filter(!str_detect(geoRegion, "CH|CHFL|FL"), altersklasse_covid19 == "20 - 29",
         datum >= sprintf("%d%02d", isoyear(study_from), 
                          isoweek(study_from)),
         datum <= sprintf("%d%02d", isoyear(study_to), 
                          isoweek(study_to))) %>%
  group_by(geoRegion) %>%
  summarise(cases2029 = sum(entries_pos),
            cases_men2029 = sum(entries_pos) * p_community$sex_pp) %>%
  left_join(population_canton_men2029) %>%
  mutate(cases_norm_men2029 = cases_men2029/pop * 100000) %>%
  select(kanton = geoRegion, pop, cases_men2029, cases_norm_men2029)

cases_canton_screening <- metadata_screening %>%
  count(name = "cases_screening", kanton) %>%
  left_join(population_canton_men2029, by = c("kanton" = "geoRegion")) %>%
  mutate(cases_norm_screening = cases_screening/pop * 100000) %>%
  select(kanton, cases_screening, cases_norm_screening)

cases_canton_toplot <- cases_canton_m2029 %>%
  arrange(cases_norm_men2029) %>% 
  mutate(cum_cases_men2029 = cumsum(cases_men2029),
         lagcumcases_men2029 = lag(cum_cases_men2029)) %>%
  replace_na(list(lagcumcases_men2029 = 0)) %>%
  left_join(cases_canton_screening %>% arrange(cases_norm_screening) %>% 
              mutate(cum_cases_screening = cumsum(cases_screening),
                     lagcumcases_screening = lag(cum_cases_screening)) %>%
              replace_na(list(lagcumcases_screening = 0)))

# Gradient breaks
breaks_gradient_screening = scales::rescale(c(min(cases_canton_screening$cases_norm_screening), 
                                              quantile(cases_canton_screening$cases_norm_screening, 0.25), 
                                              round(quantile(cases_canton_screening$cases_norm_screening, 0.5)) + 5,
                                              quantile(cases_canton_screening$cases_norm_screening, 0.75),
                                              max(cases_canton_screening$cases_norm_screening)))

breaks_gradient_community = scales::rescale(c(min(cases_canton_m2029$cases_norm_men2029), 
                                              quantile(cases_canton_m2029$cases_norm_men2029, 0.25) - 160, 
                                              round(quantile(cases_canton_m2029$cases_norm_men2029, 0.5)) + 5,
                                              quantile(cases_canton_m2029$cases_norm_men2029, 0.75),
                                              max(cases_canton_m2029$cases_norm_men2029)))

# Map data
ch_shp <- read_sf("resources/map_ch/gadm40_CHE_1.shp") %>%
  mutate(kanton = gsub("CH.", "", HASC_1),
         NAME_1 = case_when(
           NAME_1 == "Basel-Landschaft" ~ "Basel-Land",
           NAME_1 == 	"Genève" ~ "Geneva",
           T ~ NAME_1)) %>%
  mutate(geom = st_union(geometry),
         centroid_lon = st_coordinates(st_centroid(.))[,1],
         centroid_lat = st_coordinates(st_centroid(.))[,2]) %>%
  left_join(cases_canton_toplot) #%>%

# Community map
p7map <- ggplot(data = ch_shp) +
  geom_sf(aes(fill = cases_norm_men2029), color = "transparent") +
  scale_fill_gradientn(colours = colors_gradient_community, values = breaks_gradient_community) + 
  theme_map() +
  theme(legend.position = "top",
        axis.line = element_blank())
p7map

p7a <- ggplot(cases_canton_toplot %>% arrange(cases_norm_men2029)) +
  geom_bar(aes(x = 1, y = cases_men2029, fill = cases_norm_men2029),
           color = "white", stat = "identity", position = "stack", linewidth = 0.5) +
  scale_fill_gradientn(colours = colors_gradient_community, values = breaks_gradient_community) + 
  coord_flip() +
  theme_void() +
  theme(legend.position = "top")
p7a


# bp_s1_data <- subsample_metadata %>% filter(s1) %>%
#   left_join(cases_canton_toplot) %>%
#   filter(!is.na(cases_norm_men2029)) %>% # canton FL?? error for FR?
#   rowwise() %>%
#   mutate(rc = sample(x = lagcumcases_men2029:cum_cases_men2029, size = 1))
# 
# 
# p7b <- ggplot(cases_canton_toplot) +
#   geom_bar(aes(x = 1, y = cases_men2029),
#            color = "transparent", fill = "grey30", stat = "identity", position = "stack") +
#   geom_jitter(data = bp_s1_data, aes(1, rc, color = cases_norm_men2029)) +
#   scale_color_gradientn(colours = colors_gradient_community, values = breaks_gradient_community) + 
#   coord_flip() +
#   theme_void() +
#   theme(legend.position = "top")
# 
# p7b

p8map <- ggplot(data = ch_shp) +
  geom_sf(aes(fill = cases_norm_screening), color = "transparent") +
  scale_fill_gradientn(colours = colors_gradient_screening, values = breaks_gradient_screening,
                       breaks = c(50,100,200,400),
                       labels = c(50,100,200,400),
                       trans = "log1p") +
  theme_map() +
  theme(legend.position = "top",
        axis.line = element_blank())
p8map

p8a <- ggplot(cases_canton_toplot %>% arrange(cases_norm_screening)) +
  geom_bar(aes(x = 1, y = cases_screening, fill = cases_norm_screening),
           color = "white", stat = "identity", position = "stack", linewidth = 0.5) +
  scale_fill_gradientn(colours = colors_gradient_screening, values = breaks_gradient_screening) + 
  coord_flip() +
  theme_void() +
  theme(legend.position = "top")
p8a


# bp_screening_data <- subsample_metadata %>% filter(s) %>%
#   left_join(cases_canton_toplot) %>%
#   filter(!is.na(cases_norm_screening)) %>% # canton FL?? error for FR?
#   rowwise() %>%
#   mutate(rc_screening = sample(x = lagcumcases_screening:cum_cases_screening, size = 1),
#          rc_men_2029 = sample(x = lagcumcases_men2029:cum_cases_men2029, size = 1))
# 
# 
# p8b <- ggplot(cases_canton_toplot) +
#   geom_bar(aes(x = 1, y = cases_screening),
#            color = "transparent", fill = "grey30", stat = "identity", position = "stack") +
#   geom_jitter(data = bp_screening_data, aes(1, rc_screening, color = cases_norm_screening)) +
#   scale_color_gradientn(colours = colors_gradient_screening, values = breaks_gradient_screening) + 
#   coord_flip() +
#   theme_void() +
#   theme(legend.position = "top")
# 
# p8b

p8c <- ggplot(cases_canton_toplot %>% arrange(cases_norm_men2029)) +
  geom_bar(aes(x = 1, y = cases_men2029, fill = cases_norm_screening),
           color = "white", stat = "identity", position = "stack", linewidth = 0.5) +
  scale_fill_gradientn(colours = colors_gradient_screening, values = breaks_gradient_screening) + 
  coord_flip() +
  theme_void() +
  theme(legend.position = "top")
p8c

# p8d <- ggplot(cases_canton_toplot) +
#   geom_bar(aes(x = 1, y = cases_men2029),
#            color = "transparent", fill = "grey30", stat = "identity", position = "stack") +
#   geom_jitter(data = bp_screening_data, aes(1, rc_men_2029, color = cases_norm_screening)) +
#   scale_color_gradientn(colours = colors_gradient_screening, values = breaks_gradient_screening) + 
#   coord_flip() +
#   theme_void() +
#   theme(legend.position = "top")
# 
# p8d

# Arrange plots ----------------------------------------------------------------

A <- p1 / p3 + plot_layout(heights = c(1.5, 1)) 
B1 <- p7map / p8map
#B2 <- ((p7a + theme(legend.position = "none")) / (p7b+ theme(legend.position = "none")) / (p8d+ theme(legend.position = "none"))) 
C <- p4
E <- p5 
D <- p6 + coord_flip()


fig1 <- ((((A / NULL) + plot_layout(heights = c(1.5, 1, 1))  | B1) + plot_layout(widths = c(1.2, 1))) / ((C | (D / E)) + plot_layout(widths = c(1.5, 1)))) + plot_layout(heights = c(1.5,1)) +
  plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")


#Save as PDF
if (!debugging) { 
  ggsave(
    filename = snakemake@output[["fig1"]],
    plot = fig1,         
    width = 174,                       # Full page width
    height = 180,                      # Height aiming for <200mm
    units = "mm",                     
    dpi = 600                          # For vector formats, 300-600 DPI 
  )
}

# # debug ------------------------------------------------------------------------
# debugging = TRUE
# setClass(
#   "snakemake_object",
#   contains= "tbl_df",
#   slots = c(input = "character", output = "character", params = "character")
# )
# 
# snakemake <- new("snakemake_object", tibble(),
#                  input = c(metadata = "results/data/all/metadata.tsv",
#                            screening_weeks = "resources/screening_weeks.tsv",
#                            population = "resources/ext_population.tsv",
#                            tests_byage_w = "resources/ext_tests_byage_w.tsv",
#                            cases_byage_d = "resources/ext_cases_byage_d.tsv",
#                            tests_bysex_w = "resources/ext_tests_bysex_w.tsv",
#                            active_infections = "results/report/active_infections.tsv",
#                           covspectrum_alpha = "resources/AlphaVariantTimeDistributionPlot.csv",
#                            kof_stringencyidx = "resources/ext_kof_stringencyidx.csv"),
#                  output = c(table1 = "results/report/numbers.tex",
#                             fig1 = "results/report/fig1_ggplot.pdf"))
# 
# 

