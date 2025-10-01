# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee Phylodynamics
#        V\ Y /V    Script to get number of sequence to include by canton given cases
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(lubridate)

options("lubridate.week.start" = 1)
set.seed(24)

# 2. Load datasets -------------------------------------------------------------
ids <- read_tsv(snakemake@input[["ids"]])
metadata <- read_tsv(snakemake@input[["metadata"]])

# 3. Get probs for cantonal cases distribution -------------------------------
#setwd("~/Projects/2105-cov-armee/10.published-analysis")
#debug
# ids <- read_tsv("results/data/asymp_symp_male2029/asymp/ids.tsv")
# metadata <- read_tsv("results/data/asymp_symp_male2029/asymp/metadata.tsv")
# cantons <- read_delim("resources/cantons.txt", delim = " ",
#                       col_names = c("geoRegion", "division"))
# screening_weeks <- read_tsv("resources/screening_weeks.tsv")
# cases <- read_tsv("resources/ext_cases_byage_d.tsv")

cantons <- read_delim(snakemake@input[["cantons"]], delim = " ", 
                      col_names = c("geoRegion", "division"))
screening_weeks <- read_tsv(snakemake@input[["screening_weeks"]])
cases <- read_tsv(snakemake@input[["cases"]]) 

prob_canton_2029 <- cases %>%
  filter(!geoRegion %in% c("CH", "CHFL", "FL"), ageRange == "20 - 29") %>%
  left_join(cantons) %>%
  mutate(division = str_replace(division, "d ", "d"),
         division = str_replace(division, "Appenzell A. Rh.", "Appenzell Ausserrhoden")) %>%
  group_by(division) %>%
  summarise(entries = sum(entries), .groups = "drop_last") %>%
  mutate(prob = entries/sum(entries))

max_seqs <- metadata %>% count(division) %>%
  left_join(prob_canton_2029) %>%
  arrange(desc(prob)) %>% slice(1)

groups <-  metadata %>% count(division) %>%
  left_join(prob_canton_2029) %>% 
  mutate(count = ceiling(max(max_seqs$n) * 1 * prob / max(prob))) %>%
  select(division, count)

write_tsv(groups, snakemake@output[["groups"]] )

# prob_canton_2029_week <- cases %>%
#   filter(!geoRegion %in% c("CH", "CHFL", "FL"), ageRange == "20 - 29") %>%
#   mutate(week = floor_date(date, unit = "week")) %>%
#   left_join(cantons) %>%
#   group_by(week, division) %>%
#   summarise(entries = sum(entries), .groups = "drop_last") %>%
#   mutate(prob = entries/sum(entries))
# 
# 
# max_seqs <- df %>% count(week, division) %>%
#   left_join(prob_canton_2029) %>%
#   group_by(week) %>%
#   arrange(desc(prob)) %>% slice(1)
# 
# n_seqs <-  df %>% count(week, division) %>%
#   left_join(prob_canton_2029) %>% 
#   left_join(max_seqs %>% 
#               select(week, max_division = division, max_prob = prob, max_seq = n)) %>%
#   mutate(seqs_toinclude = ceiling(max_seq * 3 * prob / max_prob))
# 
# metadata %>%
#   left_join(groups) %>%
#   group_by(division) %>%
#   sample_n(size = min(n(), count), replace = F)


