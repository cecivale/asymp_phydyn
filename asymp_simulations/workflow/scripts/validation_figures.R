# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee Phylodynamics
#        V\ Y /V    Tree stats and Sim inference plot results
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


library(treestats)
library(tidyverse)
library(treedataverse)
library(ape)
library(phangorn)
library(phytools)
library(coda)
source("workflow/scripts/talking_to_beast.R")
source("workflow/scripts/plot_opts.R")

names_trace <- c("originBDMMPrime", "typeFrequencies.2", 
                 "f_type2", "p_type1", 
                 "birthRateSVi0_type1", "birthRateSVi1_type1", 
                 "birthRateSVi2_type1", "birthRateSVi3_type1",
                 "samplingRateSVi1_type2", "rhoSamplingSVe0type1",  
                 "rhoSamplingSVe1type1", "rhoSamplingSVe2type1",
                 "gammaShape", "kappa", 
                 "TreeHeight", "posterior", "likelihood")

names_plots <- c("T[origin]", "p[origin]^s",
                 "f[phi1]", "p[pop]^a",
                 "lambda[aa*','*e[1]]", "lambda[aa*' '*e[2]]", 
                 "lambda[aa*','*e[3]]",  "lambda[aa*','*e[4]]", 
                 "psi[s*','*2]", "rho[a*','*1]", 
                 "rho[a*','*2]",  "rho[a*','*3]", 
                  "alpha[Gamma]", 
                 "kappa",  "T[root]",
                 "Posterior", "Tree~likelihood"
)

grouping <- c("Origin~parameters", "Origin~parameters","Transmission~parameters", "Transmission~parameters",
              "Transmission~parameters", "Transmission~parameters",
              "Transmission~parameters", "Transmission~parameters",
              "Sampling~parameters", "Sampling~parameters",
              "Sampling~parameters", "Sampling~parameters",
               "Phylogenetic~parameters", "Phylogenetic~parameters",
              "Phylogenetic~parameters", "Compound~distributions", "Compound~distributions")

names_df <- tibble(var = names_trace, var_plot = names_plots, grouping = grouping)
names_list <- names_plots
names(names_list) <- names_trace

# 1. Tree Stats ----------------------------------------------------------------

# Define MC and AI stats (based of BATs)
# MC
max_mono_clade <- function(tree, traits, focal_state) {
  max(sapply(unique(tree$edge[,1]), function(node) {
    tips <- extract.clade(tree, node)$tip.label
    if(all(traits[tips] == focal_state)) length(tips) else 0
  }))
}

# AI
association_index <- function(tree, states) {
  if (!all(tree$tip.label %in% names(states)))
    
    states <- as.factor(states[tree$tip.label])
  n_tip <- length(tree$tip.label)
  n_node <- tree$Nnode
  internal_nodes <- (n_tip + 1):(n_tip + n_node)
  
  get_tips <- function(node) phangorn::Descendants(tree, node, "tips")[[1]]
  ai_terms <- sapply(internal_nodes, function(nd) {
    tips <- get_tips(nd)
    m <- length(tips)
    fi <- max(prop.table(table(states[tips])))
    (1 - fi) / (2 * m - 1)
  })
  sum(ai_terms)
}


# debug
# empirical_trees   <- read.beast("../2_asymp_phydyn/results/analysis/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029.0.ds.trees")
# sim_trees <- read.beast(file = "results/siminf_asymp_symp/sim/sim_2types_3p.trees")


# Set seed and read trees
set.seed(24)

empirical_trees   <- read.beast(snakemake@input[["empirical_trees"]])
sim_trees <- read.beast(file = snakemake@input[["sim_trees"]])

# Calculate tree stats for empirical trees
empirical_trees_ds <- empirical_trees[sample(size = 500, seq(1, length(empirical_trees)))]
empiricaltrees_stats_l <- lapply(empirical_trees_ds, function(tree) {
  tree <- tree@phylo
  traits <- str_split(tree$tip.label, "\\|", simplify=T)[,2]
  names(traits) <- tree$tip.label
  traits_factor <- phyDat(as.matrix(traits), type = "USER", levels = unique(traits))
  
  bind_rows(treestats::calc_all_stats(as.phylo(tree), normalize = TRUE)) %>%
    bind_cols(tibble(PS = parsimony(tree, traits_factor, method = "fitch"),
                       MC_symp = max_mono_clade(tree, traits, "symp"),
                       MC_asymp = max_mono_clade(tree, traits, "asymp"),
                       AI = association_index(tree, traits),
                       Total_branch_length <- sum(tree$edge.length)))
  }
)

empiricaltrees_stats <- bind_rows(empiricaltrees_stats_l) %>% 
  mutate(tree = names(empiricaltrees_stats_l), group = "empirical")


# Calculate tree stats for simulated trees
simtrees_stats_l <- lapply(1:length(sim_trees), function(i) {
  tree <- sim_trees[[i]]@phylo
  traits <- as_tibble(sim_trees[[i]]) %>% filter(!is.na(samp)) %>% pull(type)
  names(traits) <- tree$tip.label
  traits_factor <- phyDat(as.matrix(traits), type = "USER", levels = unique(traits))
  
  bind_rows(append(treestats::calc_all_stats(as.phylo(sim_trees[[i]]), normalize = TRUE), c("tree2" = i))) %>%
    bind_cols(tibble(PS = parsimony(tree, traits_factor, method = "fitch"),
                     MC_symp = max_mono_clade(tree, traits, "type2"),
                     MC_asymp = max_mono_clade(tree, traits, "type1"),
                     AI = association_index(tree, traits)))
  }
)

simtrees_stats <- bind_rows(simtrees_stats_l) %>% mutate(tree = names(sim_trees), group = "sim")

# Number of tips
simtrees_tips_l <- lapply(1:length(sim_trees), function(i) {
  sim_trees[[i]] %>% as_tibble() %>% filter(!is.na(samp)) %>% count(samp) %>% 
    mutate(tree = i) %>% pivot_wider(names_from = samp, values_from =  n)
}
)

simtrees_tips <- bind_rows(simtrees_tips_l) %>%
  replace(is.na(.), 0) %>%
  mutate(tips = sample11 + sample12 + sample13 + sample2,
         sample1 = sample11 + sample12 + sample13 ) 

selected_simtrees_stats <- simtrees_stats %>%
  filter(tree2 %in% (simtrees_tips %>% filter(tips >= 50, tips <= 350) %>% pull(tree)))

tree_stats <- bind_rows(empiricaltrees_stats, simtrees_stats) %>%
  pivot_longer(-c(group, tree), names_to = "stat", values_to = "value")

# Select important stats
stats_selection <- c("colless", "sackin", "average_leaf_depth", "mpd",
                      "tree_height", "mean_branch_length", 
                      "var_branch_length", "number_of_lineages",
                     "PS", "MC_symp", "MC_asymp", "AI")
#  Mean pairwise distance (mpd), the mean distance between all pairs of tips in units of branch length 


# pb-value stats
sim_dist <- tree_stats %>%
  filter(group == "sim", stat %in% stats_selection) %>%
  group_by(stat) %>%
  summarise(sim_values = list(value), .groups = "drop")

emp_vals <- tree_stats %>%
  filter(group == "empirical", stat %in% stats_selection) %>%
  select(stat, tree, value)

emp_vals_mean <- emp_vals %>%
  group_by(stat) %>%
  summarise(emp_mean = mean(value), .groups = "drop")

pb_mean <- tree_stats %>%
  filter(group == "sim") %>%
  group_by(stat) %>%
  summarise(sim_values = list(value), N = n(), .groups = "drop") %>%
  inner_join(emp_vals_mean, by = "stat") %>%
  rowwise() %>%
  mutate(
    p_right = (sum(sim_values >= emp_mean) + 0.5) / (N + 1),
    p_left  = (sum(sim_values <= emp_mean) + 0.5) / (N + 1),
    pB_two_sided = pmin(1, 2 * min(p_left, p_right))
  ) %>%
  ungroup()

pb_per_tree <- emp_vals %>%
  group_by(stat) %>%
  inner_join(sim_dist, by = "stat") %>%
  rowwise() %>%
  mutate(
    N = length(sim_values),
    p_right = (sum(sim_values >= value) + 0.5) / (N + 1),
    p_left  = (sum(sim_values <= value) + 0.5) / (N + 1),
    pB_two_sided  = 2 * min(p_left, p_right)
  ) %>%
  ungroup()

pb_summary <- pb_per_tree %>%
  group_by(stat) %>%
  summarise(
    n_empirical = n(),
    frac_signif = sum(pB_two_sided < 0.05)/n_empirical,
    median_pB_two_sided = median(pB_two_sided),
    mean_pB_two_sided = mean(pB_two_sided),
    q05_pB_two_sided    = quantile(pB_two_sided, 0.05),
    q95_pB_two_sided    = quantile(pB_two_sided, 0.95)
  )

# Supplementary figure with all stats
p1 <- ggplot(tree_stats %>% filter(stat %in% stats_selection)) +
  geom_histogram(data = . %>% filter(group == "sim") %>% ungroup %>% group_by(stat) %>%
                   mutate(value = ifelse(value > quantile(value, 0.99), NA, value)), aes(x = value, fill = tree %in% selected_simtrees_stats$tree), 
                 color = "white", alpha = 0.7) +
  geom_point(data = . %>% filter(group == "empirical") %>% group_by(stat) %>% 
               summarise(mean_value = mean(value)),
             aes(x = mean_value, y = -3),
             color = "orange") + 
  geom_errorbar(data = . %>% filter(group == "empirical") %>% group_by(stat) %>% 
               summarise(lower_value = boa::boa.hpd(value, alpha = 0.05)["Lower Bound"],
                         upper_value = boa::boa.hpd(value, alpha = 0.05)["Upper Bound"]), 
               aes(xmin = lower_value, xmax = upper_value, y = -3),
             color = "orange") + 
  geom_text(data = pb_mean, aes(x = emp_mean, y = 40, label = paste("pB = ", round(pB_two_sided, 3)))) +
  facet_wrap(~ stat, scales = "free", ncol = 3) +
  ylab("") +
  xlab("") +
  scale_fill_manual(name = "50 < Number of tips < 350", values = c( "grey80", "steelblue")) +
  theme_asympsymp() +
  theme(legend.position = "top")

p1

ggsave(
    filename = snakemake@output[["supfig_treestats"]],
    plot = p1,         
    width = 174,                       
    height = 180,                      
    units = "mm",                     
    dpi = 600                          
  )


# 2. Simulation inference coverage results -------------------------------------

#debug
# results_siminf_path <- "results/siminf_asymp_symp/inference"
# all_files <- list.files(pattern = "\\.log", path = results_siminf_path, full.names = T)
# trace_files <- all_files[str_detect(all_files, ".model3p_symp_asymp.")]

# Read traces files from each simulation inference
f_trace_l <- lapply(snakemake@input[["trace_files"]], function(file) {
  cat(file)
  full_trace <- read_table(file, comment = "#")
  if (nrow(full_trace) == 0) return()
  chain_length <- nrow(full_trace)
  
  trace <- full_trace[round(0.1*chain_length):chain_length,]
  
  f_trace <- trace %>%
    rename(f_type2 = f_symp, p_type1 = p_asymp) %>%
    pivot_longer(-Sample, names_to = "var", values_to = "value") %>%
    mutate(ansys = file,
           ansys_num = str_split(ansys, "\\.", simplify = T)[,3],
           mcmc_steps = n()) #%>%
  # group_by(var) %>%
  # mutate(ess = effectiveSize(value)) %>%
  # ungroup()
  
  return(f_trace)
})

f_trace <- bind_rows(f_trace_l) %>%
  left_join(names_df)

# Get true values from simulation
json_params <- jsonlite::fromJSON("results/siminf_asymp_symp/sim_2types_3p_params.json")
true_param <- as_tibble(json_params) %>%
  select(originBDMMPrime = sim_time, p_type1, f_type2, 
         birthRateSVi0_type1 = birth_type1_0,  birthRateSVi1_type1 = birth_type1_1, birthRateSVi2_type1 = birth_type1_2, birthRateSVi3_type1 = birth_type1_3,
         samplingRateSVi1_type2 = sampling_type2, rhoSamplingSVe0type1 = rho_type1_0, rhoSamplingSVe1type1 = rho_type1_1, rhoSamplingSVe2type1 = rho_type1_2) %>%
  mutate(typeFrequencies.2 = 1) %>%
  pivot_longer(everything(), names_to = "var", values_to = "true_value")

# Number of tips in each tree to order simulations
tree_nlin <- simtrees_tips %>% mutate(n_tips = sample11 + sample12 + sample13 + sample2,
                                      sample1 = sample11 + sample12 + sample13 ) %>%
  arrange(n_tips) %>% pull(tree)

# Summarise hpd intervals
f_trace_sum1 <- f_trace %>%
  filter(var %in% c("originBDMMPrime", "typeFrequencies.2", "p_type1", "f_type2",
                    "birthRateSVi0_type1", "birthRateSVi1_type1", "birthRateSVi2_type1", "birthRateSVi3_type1",
                    "samplingRateSVi1_type2", "rhoSamplingSVe0type1", "rhoSamplingSVe1type1", "rhoSamplingSVe2type1")) %>%
  left_join(true_param) %>%
  filter(mcmc_steps > 9000) %>%
  group_by(ansys_num, var, var_plot, true_value) %>%
  reframe(hpd95 = as_tibble_row(beastio::getHPD.boa(value, alpha = 0.05, includeMedian = T)),
          hpd90 = as_tibble_row(beastio::getHPD.boa(value, alpha = 0.1, includeMedian = F)),
          hpd50 = as_tibble_row(beastio::getHPD.boa(value, alpha = 0.5, includeMedian = F))) %>%
  unnest(cols = c(hpd95, hpd90, hpd50), names_sep = "_") %>%
  group_by(ansys_num) %>%
  mutate(fill_c = true_value >= hpd95_lower & true_value <= hpd95_upper) 

# Supplementary figure with one simulation per row
suppfig <- f_trace_sum1 %>%
  mutate(var_plot = factor(var_plot, levels = names_plots)) %>%
  ggplot(aes(y = factor(ansys_num, levels = tree_nlin))) +
  geom_segment(aes(x = hpd95_lower, xend = hpd95_upper, color = fill_c), size = 0.5, alpha = 0.4) +
  #geom_segment(aes(x = hpd90_lower, xend = hpd90_upper, color = fill_c), size = 1, alpha = 0.4) +
  #geom_segment(aes(x = hpd50_lower, xend = hpd50_upper, color = fill_c), size = 1, alpha = 0.6) +
  geom_point(aes(x = hpd95_median, fill = fill_c), shape = 22, color = "#456672", size = 1) +
  geom_vline(aes(xintercept = true_value), linetype = 1, linewidth = 1, color = "orange") +
  ylab("simulation") +
  xlab("parameter value") +
  scale_color_manual(values = c("grey", "steelblue")) +
  scale_fill_manual(values = c("grey", "steelblue")) +
  facet_wrap(~var_plot, scales = "free", ncol = 4, labeller = "label_parsed") +
  theme_asympsymp() +
  theme(axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10)
  )

suppfig

ggsave(
  filename = snakemake@output[["supfig_coverage"]],
  plot = suppfig,         
  width = 174,                       
  height = 180,                      
  units = "mm",                     
  dpi = 600                          
)

# Summarise, get Median signed relative error and coverage per parameter
f_trace_sum2 <- f_trace %>%    
  filter(var %in% c("originBDMMPrime", "typeFrequencies.2" ,"p_type1", "f_type2",
                    "birthRateSVi0_type1", "birthRateSVi1_type1", "birthRateSVi2_type1", "birthRateSVi3_type1",
                    "samplingRateSVi1_type2", "rhoSamplingSVe0type1", "rhoSamplingSVe1type1", "rhoSamplingSVe2type1")) %>%
  left_join(true_param) %>%
  left_join(simtrees_tips %>% mutate(ansys_num = as.character(tree))) %>%
  group_by(ansys, ansys_num, var, var_plot, true_value, sample11, sample12, sample13, sample2) %>%
  reframe(hpd95 = as_tibble_row(beastio::getHPD.boa(value, alpha = 0.05, includeMedian = T))) %>%
  unnest(cols = c(hpd95), names_sep = "_") %>%
  group_by(ansys_num, var) %>%
  mutate(fill_c = true_value >= hpd95_lower & true_value <= hpd95_upper,
         n_tips = sample11 + sample12+ sample13 + sample2,
         sample1 = sample11 + sample12+ sample13 ) %>% ungroup() %>% 
  group_by(var, var_plot) %>% 
  mutate(diff_median = (true_value - hpd95_median)/true_value,
         coverage = sum(fill_c)/n())


# Main figure
fig4 <- f_trace_sum2 %>%
  ggplot(aes(y = factor(var, levels = rev(names_trace)), x = diff_median, color = coverage > 0.95)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_vline(aes(xintercept = 0), color = "orange", linetype = 1, size = 1) +
  geom_label(data = f_trace_sum2 %>% select(var, coverage) %>% distinct(), 
             aes(y = factor(var, levels = rev(names_trace)), x = 1, label = round(coverage, 2)),
             size = 3) +
  scale_x_continuous(limits = c(-1.1, 1.1)) +
  scale_y_discrete(labels = function(x) parse(text = names_list[x])) +
  theme_asympsymp() +
  #coord_flip() +
  scale_color_manual(values = c("#54073e", "steelblue")) +
  xlab("Median Relative Error") +
  ylab("") +
  theme(axis.text.y = element_text(size = 10))


fig4

ggsave(
  filename = snakemake@output[["fig_validation"]],
  plot = fig4,         
  width = 100,                       
  height = 100,                      
  units = "mm",                     
  dpi = 600                          
)

