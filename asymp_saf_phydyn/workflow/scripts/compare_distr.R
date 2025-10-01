# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee Phylodynamics
#        V\ Y /V    Compare posterior distribution for data replicates in 
#    (\   / - \     phylodynamic analysis results
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


library(overlapping)
library(tidyverse)
library(patchwork)
library(treemapify)
library(ggnewscale)
source("workflow/scripts/plot_opts.R")
source("workflow/scripts/talking_to_beast.R")

debugging = FALSE

# Function to compare (overlap and JS) each parameter across chains
get_comparison <- function(trace_list, parameter, parameter_min, parameter_max) {
  
  # Extract all chains for this parameter
  ts <- lapply(trace_list, function(chain) chain[[parameter]])
  names(ts) <- paste0("chain ", seq_along(ts) - 1)
  
  # Estimate densities over common range
  ds <- lapply(ts, function(t) density(t, n = 512, from = parameter_min, to = parameter_max))
  w <- ds[[1]]$x[2] - ds[[1]]$x[1]
  ps <- t(sapply(ds, function(d) d$y * w))
  
  # Overlap (pairwise) 
  ov <- overlap(ts, plot = FALSE, pairsOverlap = TRUE, type = "2")
  if (length(trace_list) > 2) { ov_matrix <- matrix(ov$OVPairs)
  } else { ov_matrix <- ov$OV }
  chain_pairs <- combn(names(ts), 2, simplify = FALSE)
  
  ov_df <- map2_dfr(chain_pairs, ov_matrix, function(pair, val) {
    tibble(Var1 = pair[1], Var2 = pair[2], value = val, metric = "overlap")
  })
  
  # # Make symmetric
  ov_df_dup <- ov_df %>%
    rename(Var1_tmp = Var1, Var2_tmp = Var2) %>%
    mutate(Var1 = Var2_tmp, Var2 = Var1_tmp) %>%
    select(-Var1_tmp, -Var2_tmp) %>%
    bind_rows(ov_df)
  
  # Jensen-Shannon Distance
  js <- philentropy::JSD(ps, unit = "log")^0.5
  if (length(trace_list) == 2) {
    js_df_long <- tibble(
      Var1 = names(ts),
      Var2 = rev(names(ts)),
      value = c(js, js),
      metric = c("JS", "JS")
    )
  } else {
    # matrix result
    colnames(js) <- rownames(js) <- names(ts)
    js_df <- as.data.frame(js)
    js_df$Var1 <- rownames(js_df)
    
    
    js_df_long <- js_df %>%
      pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
      filter(Var1 != Var2) %>%
      mutate(metric = "JS")
  }
  
  return(bind_rows(ov_df_dup, js_df_long) %>% mutate(param = parameter))
}
  
# Parameter names and grouping
names_trace <- c("originBDMMPrime", "f_symp", "p_asymp", 
                 "birthRateSVi0_asymp", "birthRateSVi1_asymp", 
                 "birthRateSVi2_asymp", "birthRateSVi3_asymp",
                 "samplingRateSVi1_symp", "rhoSamplingSVe0asymp",  
                 "rhoSamplingSVe1asymp", "rhoSamplingSVe2asymp",
                 "typeFrequencies.1", "gammaShape", "kappa", 
                 "TreeHeight")

names_plots <- c("T[origin]", 
  "f[phi1]", "p[pop]^a",
  "lambda[aa*','*e[1]]", "lambda[aa*' '*e[2]]", 
  "lambda[aa*','*e[3]]",  "lambda[aa*','*e[4]]", 
  "psi[s*','*2]", "rho[a*','*1]", 
  "rho[a*','*2]",  "rho[a*','*3]", 
  "p[origin]^a",   "alpha[Gamma]", 
  "kappa",  "T[root]"
)

grouping <- c("Origin~parameters", "Transmission~parameters", "Transmission~parameters",
              "Transmission~parameters", "Transmission~parameters",
              "Transmission~parameters", "Transmission~parameters",
              "Sampling~parameters", "Sampling~parameters",
              "Sampling~parameters", "Sampling~parameters",
              "Origin~parameters", "Phylogenetic~parameters", "Phylogenetic~parameters",
              "Phylogenetic~parameters")

names_df <- tibble(param = names_trace, param_plot = names_plots, grouping = grouping)

# 1. Read trace files ----------------------------------------------------------
traces_files_inter <- snakemake@input[["inter"]]
traces_files_intra <- snakemake@input[["intra"]]

if (debugging) {
  traces_files_inter <- snakemake@input[1:4]
  traces_files_intra <- snakemake@input[5:7]
}
# traces_files_fprior <- snakemake@input[c(1,8)]

results_l <- lapply(list(traces_files_inter, traces_files_intra), function(traces_files) {

  traces_l <- lapply(traces_files, function(trace_file) {
                     read_trace(trace_file[1], as.numeric(snakemake@params[["burnin"]]))}
                   )

  traces <- bind_rows(traces_l) 

   # 2. Get parameter ranges for all chains combined -----------------
  param_df <- traces %>% 
    select(-file) %>%
    summarise_all(list(mmax  = max, mmin = min), ) %>%
    pivot_longer(everything(), names_to = c("param", "limit"), names_sep = "_m", values_to = "value") %>%
    pivot_wider(values_from = value, names_from = "limit") %>%
    #filter(param == "f_symp")
    filter(param %in% c("originBDMMPrime", "f_symp", "p_asymp", "birthRateSVi0_asymp",
                  "birthRateSVi1_asymp", "birthRateSVi2_asymp", "birthRateSVi3_asymp",
                  "samplingRateSVi1_symp", "rhoSamplingSVe0asymp",  "rhoSamplingSVe1asymp", "rhoSamplingSVe2asymp",
                  "typeFrequencies.1", "gammaShape", "kappa", "TreeHeight"  )) 

  # 3. Compute comparison metrics ----------------------------------------
  comparison <- lapply(1:nrow(param_df), function (i){
    get_comparison(traces_l, param_df$param[i], param_df$min[i], param_df$max[i])
  }) %>% bind_rows()
    
  # 4. Plot heatmap with JS and overlap metrics --------------------------
  p <- ggplot(data = comparison %>% 
                left_join(names_df) %>%
           filter(metric == "overlap" & Var1 > Var2 |
                  metric == "JS" & Var2 > Var1), 
         aes(Var1, Var2)) +
    geom_tile(data = . %>% left_join(names_df) %>% filter(metric == "JS") ,
              aes(fill = cut(value, c(0, 0.05, 0.1, 0.2, 0.4, 1))), 
              color = "white", lwd = 1.5, linetype = 1) +
    scale_fill_brewer(name = "Jensen-Shannon Divergence", type = "seq", 
                      palette = "OrRd", drop = F, na.value = "white") +
    new_scale_fill() +
    geom_tile(data = . %>% filter(metric == "overlap"), 
              aes(fill = cut(value, c(0, 0.5, 0.8, 0.9, 0.95, 1))),
              color = "white", lwd = 1.5, linetype = 1) +
    scale_fill_brewer(name = "Posterior density overlap %", type = "seq", 
                      palette = "YlGn", drop = F, na.value = "white") +
    geom_text(aes(label = round(value, 2),
                  color = ifelse(value < 0.1 | (value >= 0.5 & value <= 0.8), "#6b6b6b", "white")), 
             # color = "white",
              size = 2.5) +
    theme_asympsymp() +
    theme(legend.position = "right",
          strip.text = element_text(size = 10)) +
    facet_wrap(~ grouping + param_plot,  labeller = "label_parsed") +
    scale_color_identity() +
    xlab("") +
    ylab("")
        
  return(list(comparison %>% filter(Var1 > Var2) , p)) 
})

results <- bind_rows(results_l[[1]][[1]] %>% mutate(comparison = "inter-replicate"),
                     results_l[[2]][[1]] %>% mutate(comparison = "intra-replicate")) #,
                     #results_l[[3]][[1]] %>% mutate(comparison = "prior-replicate"))


bp <- ggplot(results %>% left_join(names_df), 
             aes(value, grouping, fill = metric)) +
  geom_violin(aes(color = metric), linewidth = 0.0, alpha = 0.6) + 
  geom_boxplot(width = 0.1, linewidth = 0.4,  outlier.shape = NA, color = "grey30") +
  facet_grid(comparison ~ metric, scales = "free") +
  theme_asympsymp() +
  scale_fill_manual(values = c("orange", "yellowgreen")) +
  scale_color_manual(values = c("orange", "yellowgreen")) +
  #scale_x_continuous(limits = c(0,1)) +
  scale_y_discrete(labels = function(y) parse(text = y)) +
  theme(strip.text = element_text(size = 10)) 

#Save as PDF
if (!debugging) { 
  ggsave(
    filename = snakemake@output[["supfig1"]],
    plot = results_l[[1]][[2]],
    width = 174,                      
    height = 180,                     
    units = "mm",                      
    dpi = 600                          
  )
  
  
  ggsave(
    filename = snakemake@output[["supfig2"]],
    plot = bp,        
    width = 125,                      
    height = 100,                     
    units = "mm",                     
    dpi = 600                         
  )
}



# # debug ------------------------------------------------------------------------
debugging = TRUE
setClass(
  "snakemake_object",
  contains= "tbl_df",
  slots = c(input = "character", output = "character", params = "character")
)

snakemake <- new("snakemake_object", tibble(),
                 input = c("results/analysis/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029.0.log",
                              "results/analysis/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029.1.log",
                              "results/analysis/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029.2.log",
                              "results/analysis/phydyn_asymp_symp_3p_srln/asymp_symp200_male2029.3.log",
                           "results/analysis/phydyn_asymp_symp_3p_srln/chains/asymp_symp200_male2029.0.1.log",
                              "results/analysis/phydyn_asymp_symp_3p_srln/chains/asymp_symp200_male2029.0.2.log",
                              "results/analysis/phydyn_asymp_symp_3p_srln/chains/asymp_symp200_male2029.0.3.log"),
                 params = c(burnin = "0.1"),
                 output = c(supfig1 = "results/report/supfig1_ggplot.pdf",
                            supfig2 = "results/report/supfig2_ggplot.pdf"))



