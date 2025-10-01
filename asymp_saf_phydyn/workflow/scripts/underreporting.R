# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee Phylodynamics
#        V\ Y /V    Script to explore underreporting data from Stringini et al 
#    (\   / - \     serology studies by age
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

# Underreporting
tests_byage_w <- read_tsv("resources/ext_tests_byage_w.tsv")
tests_byage_d <- read_tsv("resources/ext_cases_byage.tsv")
population <- read_tsv("resources/ext_population.tsv")

tests_byage_d %>%
  filter(date <= "2021-03-05",
         #date >= min(metadata$date))
         date >= "2020-09-03") %>%
  ggplot() +
  geom_point(aes(date, entries, color = ageRange))
geom_histogram(aes(date, entries, fill = ageRange), stat = "identity", position = "dodge")


tests_ge_sw <- tests_byage_d %>%
  filter(date <= "2020-12-08,",
         #date >= min(metadata$date))
         date >= "2020-09-01",
         geoRegion == "GE") %>%
  group_by(ageRange) %>%
  summarise(entries = sum(entries)) 

ggplot(tests_ge_sw, aes(ageRange, entries)) +
           geom_point()

sero_sw <- read_csv("resources/seroprevalence_secondwave_Stringhini2021") %>%
  mutate(from = as.numeric(str_split(age_group, "-", simplify = T)[,1]),
         to = as.numeric(str_split(age_group, "-", simplify = T)[,2])) %>%
  rowwise() %>% 
  mutate(ages = paste(seq(from, to, 1), collapse = ",")) %>%
  tidyr::separate_rows(ages, sep = ',')

sero_fw <- read_csv("resources/seroprevalence_firstwave_Stringhini2020") %>%
  mutate(from = as.numeric(str_split(age_group, "-", simplify = T)[,1]),
         to = as.numeric(str_split(age_group, "-", simplify = T)[,2])) %>%
  rowwise() %>% 
  mutate(ages = paste(seq(from, to, 1), collapse = ",")) %>%
  tidyr::separate_rows(ages, sep = ',')
         
sero_both <- sero_sw %>% left_join(sero_fw, by = "ages") %>%
  mutate(sero_both_mean = seroprevalence_mean - seroprevalence_positive) %>%
  fill(sero_both_mean, .direction = "up") %>%
  mutate(ages = as.numeric(ages),
         ageRange = case_when(
    ages < 9 ~ "0 - 9",
    ages < 19 ~ "10 - 19",
    ages < 29 ~ "20 - 29",
    ages < 39 ~ "30 - 39",
    ages < 49 ~ "40 - 49",
    ages < 59 ~ "50 - 59",
    ages < 69 ~ "60 - 69",
    ages < 79 ~ "70 - 79",
    TRUE ~ "80+"
  )) %>%
  group_by(ageRange) %>%
  summarise(sero_both_mean = mean(sero_both_mean))

pop_ge <- population %>%
  filter(geoRegion == "GE") %>%
  group_by(ageRange) %>%
  summarise(pop = sum(pop))

estim_inf <- sero_both %>% left_join(pop_ge) %>%
  mutate(estim_inf = sero_both_mean * pop / 100) %>%
  left_join(tests_ge_sw) %>%
  mutate(underreporting = estim_inf/entries)

# So it looks like controlling for age does not change the results much, 
# and underreporting from seroprevalence is probably around 2
# How then we explained an underreporting of 5 considering the asymptomatics?

tests_byage_d %>%
  filter(geoRegion == "CH") %>%
  summarise()
  ggplot() +
  geom_bar(aes())
    