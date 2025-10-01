# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee Phylodynamics
#        V\ Y /V    Parse and plot phylodynamic results
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

# 0. Load libraries ------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(treemapify)
library(treedataverse)
library(ggridges)
library(HDInterval)
source("workflow/scripts/plot_opts.R")
source("workflow/scripts/talking_to_beast.R")
source("workflow/scripts/priors.R")

debugging = FALSE

# Parameter names and grouping -------------------------------------------------
names_trace <- c("originBDMMPrime", "f_symp", "p_asymp", 
                 "birthRateSVi0_asymp", "birthRateSVi1_asymp", 
                 "birthRateSVi2_asymp", "birthRateSVi3_asymp",
                 "samplingRateSVi1_symp", "rhoSamplingSVe0asymp",  
                 "rhoSamplingSVe1asymp", "rhoSamplingSVe2asymp",
                 "typeFrequencies.1", "gammaShape", "kappa", 
                 "TreeHeight", "posterior", "likelihood")

names_plots <- c("T[origin]", 
                 "f[phi1]", "p[pop]^a",
                 "lambda[aa*','*e[1]]", "lambda[aa*' '*e[2]]", 
                 "lambda[aa*','*e[3]]",  "lambda[aa*','*e[4]]", 
                 "psi[s*','*2]", "rho[a*','*1]", 
                 "rho[a*','*2]",  "rho[a*','*3]", 
                 "p[origin]^a",   "alpha[Gamma]", 
                 "kappa",  "T[root]",
                 "Posterior", "Tree~likelihood"
)

grouping <- c("Origin~parameters", "Transmission~parameters", "Transmission~parameters",
              "Transmission~parameters", "Transmission~parameters",
              "Transmission~parameters", "Transmission~parameters",
              "Sampling~parameters", "Sampling~parameters",
              "Sampling~parameters", "Sampling~parameters",
              "Origin~parameters", "Phylogenetic~parameters", "Phylogenetic~parameters",
              "Phylogenetic~parameters", "Compound~distributions", "Compound~distributions")

names_df <- tibble(var = names_trace, var_plot = names_plots, grouping = grouping)

# 1. Read trace file -----------------------------------------------------------

burnin <- snakemake@params[["burnin"]]
trace <- read_trace(snakemake@input[["trace"]], 0.1)
MRS <- snakemake@params[["mrs"]]

change_times <- trace %>%
  select(originBDMMPrime, contains("time")) %>%
  mutate(across(everything(), ~ originBDMMPrime - .)) %>%
  slice(1) %>%
  pivot_longer(everything(), names_to = "var", values_to = "age") %>%
  mutate(param = str_split(var, "_|SV|SM|\\.", n = 2, simplify = T)[, 1],
         param = ifelse(param %in% c("originBDMMPrime", "birthRateAmongDemes"), "birthRate", param),
         epoch = str_replace(str_extract(var, c("i[0-9]|e[0-9]")), "e", "i")) %>%
  replace_na(list(epoch = 'i3')) %>%
  select(-var) %>%
  distinct() %>%
  mutate(date = to_date(age, mrs = MRS)) %>%
  arrange(date)

trace_long <- trace %>% select(!contains("time"), -file) %>%
  pivot_longer(-Sample, names_to = "var", values_to = "value") %>%
  mutate(param = str_split(var, "_|SV|SM|\\.", n = 2, simplify = T)[, 1],
         deme = case_when(
           grepl("to_asymp", var) ~ "symp",
           grepl("asymp|\\.1", var) ~ "asymp",
           grepl("to_symp", var) ~ "asymp",
           grepl("symp|\\.2", var) ~ "symp",
           T ~ ""),
         partner = str_split(var, "_", n = 4, simplify = T)[, 4],
         param = ifelse(param == "birthRateAmongDemes", "birthRate", param),
         partner = ifelse(partner == "", deme, partner),
         epoch = str_replace(str_extract(var, c("i[0-9]|e[0-9]")), "e", "i"),
         param2 = case_when(
           param == "birthRate" ~ paste(param, deme, partner, sep = "_"),
           TRUE ~ paste0(param, "_", deme))) %>%
  filter(!param %in% c("migrationRate", "removalProb", "X56", "X55")) %>%
  left_join(change_times) %>%
  replace_na(list(age = 0, date = ymd(MRS))) %>%
  left_join(names_df)

death_rate <- unique(c(trace$deathRateSVasymp, trace$deathRateSVsymp))


prior_samples <- phydyn_asymp_symp_3p %>%
  mutate(param = c("originBDMMPrime", "typeFrequencies", "deathRate", "birthRate", 
                   "f", "p", "samplingRate", "rhoSampling", "gammaShape", "kappa", "Frequencies"),
         upper_limit = as.numeric(str_replace(str_split(density95, " - ", simplify = T)[,2], "\\)", "")))  %>%
  right_join(trace_long %>% select(param, var_plot) %>% distinct() %>% filter(!is.na(var_plot))) %>% 
  select(var_plot, samples, upper_limit, distribution) %>%
  unnest(samples) 


ordered_levels <- trace_long %>%
  filter(!is.na(var_plot)) %>%
  arrange(grouping) %>%
  pull(var_plot) %>%
  unique()
trace_long$var_plot <- factor(trace_long$var_plot, levels = ordered_levels)
prior_samples$var_plot <- factor(prior_samples$var_plot, levels = ordered_levels)

# 2. Supplementary Figures all parameters -----------------------------------------
supfig3 <- ggplot(trace_long %>% filter(!is.na(var_plot))) +
  geom_density(data = prior_samples, aes(samples), fill = "grey90", color = "grey90") + 
  #geom_density(aes(value, fill = grouping, color = grouping, alpha = 0.99)) +
  geom_density_ridges_gradient(aes(value, y = 0, fill = stat(quantile)), quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
  geom_density_ridges_gradient(aes(value, y = 0), fill = "transparent", quantile_lines = TRUE, quantiles = 2, vline_linetype = 1) +
  facet_wrap(~ var_plot, scales = "free", labeller = "label_parsed", ncol = 4) +
  theme_asympsymp() +
  #scale_color_brewer(palette = "Dark2") +
  scale_fill_manual(values = c("transparent", "steelblue", "transparent"), guide = "none") +
  theme(strip.text = element_text(size = 10)) 


#Save as PDF
if (!debugging) { 
  ggsave(
    filename = snakemake@output[["supfig3"]],
    plot = supfig3,
    width = 174,                      
    height = 180,                     
    units = "mm",                      
    dpi = 600                          
  )
}


# 3. Main figure ------------------------------------------------------------------

# A. Summary tree --------------------------------------------------------------
# CCD Tree tip typed tree for swooshy tree
tree <- read.beast(snakemake@input[["ccdtree"]]) %>%
  mutate(type = str_split(label, "\\|", simplify = T)[,2])

# if (!debugging) write.beast(tree, snakemake@output[["tiptyped_tree"]])

options(ignore.negative.edge=TRUE)

t <- ggtree(tree, mrsd = MRS, size = 0.25, color = "grey10") +
  theme_tree2() +
  ggtree::scale_x_ggtree(breaks = seq(2020, 2021.5, 0.25), 
                         labels = decimal2Date(seq(2020, 2021.5, 0.25))) +
  geom_range(range = 'height_0.95_HPD', color = "grey", alpha = 0.5, size = 1) +
  geom_nodepoint(aes(alpha = posterior > 0.95), size = 0.4, color = "grey10") +
  scale_alpha_manual(values = c(0,1))
#geom_tippoint(aes(color = grepl("asymp", label))) +

dftype <- tibble(label = tree@phylo$tip.label, type = str_split(label, "\\|", simplify = T)[,2]) %>%
  column_to_rownames("label")

treeplot_A <- gheatmap(t, dftype, width = 0.02, colnames = FALSE) +
  scale_fill_manual(name = "", labels = c("Aymptomatics", "Symptomatics"), 
                    values = c(pal_asympsymp[["asymp"]], pal_asympsymp[["symp"]])) +
  theme(legend.position = "")  
  


# B. Main parameters -----------------------------------------------------------

# f
plot_B1 <- ggplot(trace_long %>% filter(var == "f_symp") %>% mutate(distr = "posterior") %>% bind_rows(
  prior_samples %>% 
    filter(grepl("Uniform", distribution) | samples < upper_limit) %>%
    filter(var_plot == "f[phi1]") %>% rename(value = samples) %>% mutate(distr = "prior"))
) +
  geom_density_ridges(aes(y = distr , value, fill = distr), color = "white", scale = 7, alpha = 0.8,
                      quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2, lwd = 0.25) +
  geom_density_ridges_gradient(aes(value, y = distr), fill = "transparent", 
                               quantile_lines = TRUE, quantiles = 2, vline_linetype = 1, scale = 7, lwd = 0.25, color = "white")  +
  theme_asympsymp() +
  scale_fill_manual(values = c(pal_asympsymp[["symp"]], pal_asympsymp[["grey"]])) +
  scale_x_continuous(limits = c(0, 15))

plot_B1

# p pop
plot_B2 <- ggplot(trace_long %>% filter(var == "p_asymp") %>% mutate(distr = paste("posterior", deme), type = "asymp") %>% bind_rows(
  trace_long %>% filter(var == "p_asymp") %>% mutate(value = 1 - value) %>% mutate(var ="p_symp", type = "symp") %>% mutate(distr = "posterior")) %>% bind_rows(
  prior_samples %>%
    filter(grepl("Uniform", distribution) | samples < upper_limit) %>%
    filter(var_plot == "p[pop]^a") %>% rename(value = samples) %>% mutate(distr = "prior", type = NA))
) +
  geom_density_ridges(aes(y = distr , value, fill = interaction(distr, type)), color = "white", scale = 7, alpha = 0.8,
                      quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2, lwd = 0.25) +
  geom_density_ridges_gradient(aes(value, y = distr, group = interaction(distr, type)), fill = "transparent", 
                               quantile_lines = TRUE, quantiles = 2, vline_linetype = 1, scale = 7, lwd = 0.25, color = "white")  +
  theme_asympsymp() +
  scale_fill_manual(values = c(pal_asympsymp[["asymp"]], pal_asympsymp[["symp"]]),na.value =  pal_asympsymp[["grey"]]) +
  scale_x_continuous(limits = c(0, 1))

plot_B2

# p infection
prior_probinf <- bind_cols(prior_samples %>% filter(var_plot == "lambda[aa*','*e[1]]") %>% select(birthRateSVi0_asymp = samples),
                           prior_samples %>% filter(var_plot == "p[pop]^a") %>% select(p_asymp = samples),
                           prior_samples %>% filter(var_plot == "f[phi1]") %>% select(f_symp = samples)) %>%
  mutate(p_symp = 1 - p_asymp,
         birthRateSVi0_symp = birthRateSVi0_asymp*f_symp*(1/p_asymp -1),
         inf_asymp0 = birthRateSVi0_asymp / p_asymp,
         inf_symp0 = birthRateSVi0_symp / p_symp,
         inf_r0 = inf_symp0 / inf_asymp0, # should be equal to f
         probinf_asymp = (inf_asymp0 * p_asymp / (inf_asymp0 * p_asymp + inf_symp0 * p_symp))) %>%
         #probinf_symp = (inf_symp0 * p_symp / (inf_symp0 * p_symp + inf_asymp0 * p_asymp))) %>%
  select(probinf_asymp) %>% #, probinf_symp
  pivot_longer(everything(), values_to = "value", names_sep = "_", names_to = c("var", "type"))
                           
trace_probinf <- trace %>%
  mutate(p_symp = 1 - p_asymp,
         inf_asymp0 = birthRateSVi0_asymp / p_asymp,
         inf_symp0 = birthRateSVi0_symp / p_symp,
         inf_r0 = inf_symp0 / inf_asymp0, # should be equal to f
         probinf_asymp = (inf_asymp0 * p_asymp / (inf_asymp0 * p_asymp + inf_symp0 * p_symp)),
         probinf_symp = (inf_symp0 * p_symp / (inf_symp0 * p_symp + inf_asymp0 * p_asymp))) %>%
  select(Sample, f_symp, p_asymp, p_symp, contains("inf_")) 

all(trace_probinf$f_symp - trace_probinf$inf_r0 < 1e-10)
trace_probinf %>% reframe(hpd_probinf_asymp = beastio::getHPD.boa(probinf_asymp, alpha = 0.05),
                          hpd_probinf_symp = beastio::getHPD.boa(probinf_symp, alpha = 0.05),
                          hpd_fsymp = beastio::getHPD.boa(f_symp, alpha = 0.05),
                          hpd_pasymp = beastio::getHPD.boa(p_asymp, alpha = 0.05))

plot_B3 <- ggplot(trace_probinf %>%
         pivot_longer(-Sample,  names_to = c("param", "deme"),
                      names_sep = "_", values_to = "value") %>% 
         filter(param == "probinf") %>% mutate(distr =  paste0("posterior",deme), 
                                               type = deme) %>% bind_rows(prior_probinf %>% mutate(distr = "prior", type = NA))
) +
  geom_density_ridges(aes(y = distr , value, fill = interaction(distr, type), color = interaction(distr, type)), #color = "white", 
                      scale = 7, alpha = 0.6,
                      quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2, lwd = 0.25, color = "white") +
  geom_density_ridges_gradient(aes(value, y = distr, color = interaction(distr, type)), fill = "transparent", 
                               quantile_lines = TRUE, quantiles = 2, vline_linetype = 1, scale = 7, lwd = 0.25, color = "white")  +
  # geom_boxplot(aes(y = distr , value, fill = interaction(distr, type), color = interaction(distr, type)),
  #              width = .3, show.legend = FALSE, outlier.shape = NA, fill = "transparent") + theme_asympsymp() +
  scale_fill_manual(values = c(pal_asympsymp[["asymp"]], pal_asympsymp[["symp"]], pal_asympsymp[["grey"]], pal_asympsymp[["grey"]]),na.value =  pal_asympsymp[["grey"]]) +
  scale_color_manual(values = c(pal_asympsymp[["asymp"]], pal_asympsymp[["symp"]], pal_asympsymp[["grey"]], pal_asympsymp[["grey"]]),na.value =  pal_asympsymp[["grey"]]) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_asympsymp()
  

plot_B3

# C. Reproductive number -------------------------------------------------------
#  Calculate Re for asymp, symp and total population
  
re_df <- trace_long %>%
  select(Sample, param2, epoch, age, date, value) %>%
  filter(str_detect(param2, "birth") | param2 == "p_asymp") %>%
  mutate(epoch = ifelse(param2 == "p_asymp", "i3", epoch)) %>%
  pivot_wider(names_from = param2, values_from = value) %>%
  mutate_at(vars(contains("birth")), list(Re = function(var) var/death_rate)) %>%
  rename_at(vars(contains("Re")), function(name) str_replace(str_replace(name, "birthRate", "Re"), "_Re", "")) %>%
  fill(p_asymp, .direction = "up") %>%
  mutate(p_symp = 1 - p_asymp,
         Re_total = Re_asymp_asymp * p_asymp + Re_asymp_symp * p_asymp + 
                    Re_symp_symp * p_symp + Re_symp_asymp * p_symp,
         Re_symp_total = Re_symp_symp + Re_symp_asymp,
         Re_asymp_total = Re_asymp_symp + Re_asymp_symp) %>%
  pivot_longer(birthRate_asymp_asymp:last_col(), names_to = c("param", "deme", "partner"), 
               names_sep = "_", values_to = "value") 

origin <- trace_long %>% filter(var == "originBDMMPrime") %>%
  reframe(age = as_tibble_col(beastio::getHPD.boa(value))$value,
          epoch = c("low", "median", "high"),
          param = "origin") %>%
  mutate(date = to_date(age, MRS))

re_to_plot <- trace_long %>% filter(var == "originBDMMPrime") %>%
  reframe(age = as_tibble_col(beastio::getHPD.boa(value))$value,
          epoch = paste0("origin_", c("low", "median", "high"))) %>%
  mutate(date = to_date(age, MRS)) %>%
  cross_join(re_df %>% filter(epoch == "i0") %>% select(- c(epoch, age, date))) %>%
  bind_rows(re_df) %>% 
  filter(param == "Re", (deme == "total"| partner == "total")) %>%
  group_by(param, deme, partner, date, epoch) %>%
  reframe(#q = as_tibble_row(quantile(value, probs = seq(0, 1, 0.05))),
          hpd = as_tibble_row(beastio::getHPD.boa(value))) %>%
  mutate(fill_c = paste(deme, partner),
         date = ymd(date))

re_cases <- read_csv(snakemake@input[["cases"]]) %>%
  filter(region == "CHE",
         date <= max(re_to_plot$date),
         date >= min(re_to_plot$date)) 

plot_C <- ggplot(re_to_plot) +
  geom_step(aes(x = date, y = hpd$median, group = interaction(deme, partner), color = fill_c), direction = "vh", linewidth = 0.6) +
  pammtools::geom_stepribbon(aes(x = date, ymin = hpd$lower, ymax = hpd$upper, fill = fill_c, group = interaction(deme, partner)), 
                             direction = "vh", alpha = 0.3) +
  geom_line(data = re_cases %>% filter(data_type == "Confirmed cases", estimate_type == "Cori_step"), 
            aes(date,median_R_mean), color = "#3e0f59") +
  pammtools::geom_stepribbon(data = re_cases %>% filter(data_type == "Confirmed cases", estimate_type == "Cori_step"), 
                             aes(x = ymd(date), ymin = median_R_lowHPD, ymax = median_R_highHPD), 
                             direction = "vh", alpha = 0.2, fill = "#3e0f59") +
  geom_hline(aes(yintercept = 1), linewidth = 0.35, linetype = 2, color = "grey40") +
  scale_fill_manual(name = "", labels = c("Asymptomatics", "Symptomatics", "Total"), 
                    values = c(pal_asympsymp[["asymp"]], pal_asympsymp[["symp"]], pal_asympsymp[["grey"]])) +
  scale_color_manual(name = "", labels = c("Asymptomatics", "Symptomatics", "Total"), 
                     values = c(pal_asympsymp[["asymp"]], pal_asympsymp[["symp"]], pal_asympsymp[["grey"]])) +
  scale_y_continuous(name = "Reproductive number Re", breaks = seq(0, 8 , 2)) +
  scale_x_date(name = "", date_breaks = "month", date_labels = "%b") +
  #coord_cartesian(xlim = ymd(c(to_date(treestats::tree_height(tree@phylo), mrs = MRS), to_date(0, mrs = MRS)))) + 
  coord_cartesian(xlim = ymd(c("2020-02-01", to_date(0, mrs = MRS)))) + 
  theme_asympsymp() +
  theme(panel.grid.minor = element_line(size = 0))

plot_C


# C. Treemap infections --------------------------------------------------------

infections <- trace_probinf %>% select(Sample, probinf_asymp, probinf_symp, p_symp, p_asymp) %>%
  mutate(infections_byasymp_asymp = probinf_asymp * p_asymp,
         infections_byasymp_symp = probinf_asymp * p_symp,
         infections_bysymp_symp = probinf_symp * p_symp,
         infections_bysymp_asymp = probinf_symp * p_asymp) %>%
  select(Sample, infections_byasymp_asymp:infections_bysymp_asymp) %>%
  pivot_longer(-Sample, names_to = c("param", "infector", "infectee"), 
               names_sep = "_", values_to = "value") %>%
  group_by(param, infector, infectee) %>%
  filter(!is.na(value)) %>% # fix in other way NA issue
  reframe(value = as_tibble_col(beastio::getHPD.boa(value))$value,
          hpd = c("low", "median", "high")) 


infections_wide <- infections %>%
  pivot_wider(names_from = hpd, values_from = value) %>%
  # create label with HPD info
  mutate(label = paste0(infectee, " ", infector, "\n", round(median * 100), "% (", 
  round(low * 100), "-", round(high*100), " HPDI 95%)"))

# plot using median as the area
plot_D <- ggplot(infections_wide, 
             aes(area = median, fill = infectee, 
                 label = label, subgroup = infectee)) +
  geom_treemap(colour = "white", size = 2, alpha = 1) +
  geom_treemap_subgroup_border(colour = "white", size = 5) +
  geom_treemap_text(colour = "white", place = "centre",
                    size = 10, grow = F) +
  scale_fill_manual(name = "", labels = c("Asymptomatics", "Symptomatics"), 
                                                            values = c(pal_asympsymp[["asymp"]], pal_asympsymp[["symp"]])) +
  theme_asympsymp()

plot_D


fig3 <- (treeplot_A / plot_C + plot_layout(heights = c(3, 2)) | (plot_B1 / plot_B2 / plot_B3) / plot_D + 
           plot_layout(heights = c(1, 1, 1, 2), ncol = 1)) + plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = "A")
fig3

#Save as PDF
if (!debugging) { 
  ggsave(
    filename = snakemake@output[["output_fig"]],
    plot = fig3,
    width = 174,                      
    height = 150,                     
    units = "mm",                      
    dpi = 600                          
  )
}



# debug ------------------------------------------------------------------------

# setwd("~/Projects/2105-cov-armee/10.published-analysis")
# debugging = TRUE
# 
# setClass(
#   "snakemake_object",
#   contains= "tbl_df",
#   slots = c(input = "character", output = "character", params = "character")
# )
# 
# snakemake <- new("snakemake_object", tibble(),
#                  input = c(
#                    trace = "results/analysis/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029.0.log",
#                    ccdtree = "results/analysis/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029.0.CCD0.tree",
#                    # mcctree = "results/analysis/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029.0.MCC.tree",
#                    cases = "resources/ext_re_cases.csv"),
#                  params = c(burnin = "0.1",
#                             mrs = "2021-02-23"),
#                  output = c(output_fig = "results/report/phydyn_asymp_symp200_male2029.0_ggplot.pdf",
#                             supfig3 = "results/report/suppfig_phydyn_asymp_symp200_male2029.0_ggplot.pdf",
#                             tiptyped_tree = "results/analysis/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029.0.CCD0.tiptyped.tree"))
# 
# 
# snakemake <- new("snakemake_object", tibble(),
#                  input = c(
#                    trace = c("results/analysis/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029.1.log")),
#                  params = c(burnin = "0.1",
#                             mrs = "2021-02-20"),
#                  output = c(output_fig = "results/report/phydyn_asymp_symp200_male2029.1_ggplot.pdf",
#                             supfig3 = "results/report/suppfig_phydyn_asymp_symp200_male2029.1_ggplot.pdf"))
# 
# 
# snakemake <- new("snakemake_object", tibble(),
#                  input = c(
#                    trace = c("results/analysis/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029.2.log")),
#                  params = c(burnin = "0.1",
#                             mrs = "2021-02-23"),
#                  output = c(output_fig = "results/report/phydyn_asymp_symp200_male2029.2_ggplot.pdf",
#                             supfig3 = "results/report/suppfig_phydyn_asymp_symp200_male2029.2_ggplot.pdf"))
# 
# 
# snakemake <- new("snakemake_object", tibble(),
#                  input = c(
#                    trace = c("results/analysis/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029.3.log")),
#                  params = c(burnin = "0.1",
#                             mrs = "2021-02-23"),
#                  output = c(output_fig = "results/report/phydyn_asymp_symp200_male2029.3_ggplot.pdf",
#                             supfig3 = "results/report/suppfig_phydyn_asymp_symp200_male2029.3_ggplot.pdf"))
# 

                      