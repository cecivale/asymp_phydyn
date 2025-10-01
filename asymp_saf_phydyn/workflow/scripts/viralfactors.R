# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee Phylodynamics
#        V\ Y /V    Lineages and viral factors
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


library(tidyverse)
library(xtable)
source("workflow/scripts/plot_opts.R")

# 1. Read and count viral classifications --------------------------------------
metadata <- read_tsv(snakemake@input[["metadata"]]) %>%
  filter(date >= "2021-01-01", date <= "2021-02-24", !army | screening)
# what if we only take sequences from those specific weeks, does it change results? specially B.1.160.16

vlin_2029 <- metadata %>% 
  filter(sex == "Männlich", age_cat10 == "20-29") %>%
  filter(!is.na(pangoLineage)) %>%
  count(pangoLineage, army, screening) %>%
  mutate(group = case_when(army & screening ~ "screening_male2029",
                           !army & !screening ~ "community_male2029")) %>%
  select(-c(screening, army)) %>%
  group_by(group) %>%
  mutate(p = n/sum(n))

vclade_2029 <- metadata %>% 
  filter(sex == "Männlich", age_cat10 == "20-29") %>%
  filter(!is.na(nextstrainClade)) %>%
  count(nextstrainClade, army, screening) %>%
  mutate(group = case_when(army & screening ~ "screening_male2029",
                           !army & !screening ~ "community_male2029")) %>%
  select(-c(screening, army)) %>%
  group_by(group) %>%
  mutate(p = n/sum(n))

vlin_all <- metadata %>% 
  filter(!is.na(pangoLineage)) %>%
  count(pangoLineage, army, screening) %>%
  mutate(group = case_when(army & screening ~ "screening_all",
                           !army & !screening ~ "community_all")) %>%
  select(-c(screening, army)) %>%
  group_by(group) %>%
  mutate(p = n/sum(n))

vclade_all <- metadata %>% 
  filter(!is.na(nextstrainClade)) %>%
  count(nextstrainClade, army, screening) %>%
  mutate(group = case_when(army & screening ~ "screening_all",
                           !army & !screening ~ "community_all")) %>%
  select(-c(screening, army)) %>%
  group_by(group) %>%
  mutate(p = n/sum(n))

vlin <- bind_rows(vlin_2029, vlin_all) %>% filter(group != "screening_male2029")
vclade <- bind_rows(vclade_2029, vclade_all) %>% filter(group != "screening_male2029")
  
# Significantly different? -----------------------------------------------------
vclade_matrix <- vclade %>% select(-p) %>% pivot_wider(names_from = nextstrainClade, values_from = n) %>%
  replace(is.na(.), 0) %>% 
  column_to_rownames("group") %>% as.matrix
vclade_matrix_other <- rowSums(vclade_matrix) - vclade_matrix

vclade_results <- lapply(1:ncol(vclade_matrix), function (i) {
  fisher.test(cbind(vclade_matrix[, i],vclade_matrix_other[, i]), simulate.p.value = TRUE)
})
vclade_pvalues <- tibble(nextstrainClade = colnames(vclade_matrix),
                          pvalue = sapply(vclade_results, function(x) x$p.value),
                          pvalueadj = p.adjust(pvalue, method = "BH"),
                         `p-value (adj)` = paste0(round(pvalue,3), " (", round(pvalueadj,3), ")")) %>%
  select(1,4)


vlin_matrix <- vlin %>% select(-p) %>% pivot_wider(names_from = pangoLineage, values_from = n) %>%
  replace(is.na(.), 0) %>% 
  column_to_rownames("group") %>% as.matrix
vlin_matrix_other <- rowSums(vlin_matrix) - vlin_matrix

vlin_results <- lapply(1:ncol(vlin_matrix), function (i) {
  fisher.test(cbind(vlin_matrix[, i],vlin_matrix_other[, i]), simulate.p.value = TRUE)
})
vlin_pvalues <- tibble(pangoLineage = colnames(vlin_matrix),
                       pvalue = sapply(vlin_results, function(x) x$p.value),
                       pvalueadj = p.adjust(pvalue, method = "BH"),
                       `p-value (adj)` = paste0(round(pvalue,3), " (", round(pvalueadj,3), ")")) %>%
  select(1,4)

# 2. Figure nextclade clades ---------------------------------------------------

ggplot(vclade) +
  geom_bar(aes(group, p, fill = nextstrainClade), stat = "identity")

# 1. Table viral factors -------------------------------------------------------

table_vclade <- vclade %>%
  mutate(`nextstrainClade | pangoLineage` = paste(nextstrainClade, "|-"))%>% 
  pivot_wider(names_from = group, values_from = c(n,p), values_fill = 0) %>%
  arrange(#desc(n_screening_male2029), 
    desc(n_screening_all),
    desc(n_community_male2029),
    desc(n_community_all)) %>%
  mutate_if(is.numeric, coalesce, 0) %>%
  mutate(#`Screening men 20-29` = paste0(round(p_screening_male2029, 4), " (", n_screening_male2029, ")"),
    `SAF Screening` = paste0(round(p_screening_all, 4), " (", n_screening_all, ")"),
    `Community\nmen 20-29` = paste0(round(p_community_male2029, 4), " (", n_community_male2029, ")"),
    `All Community` = paste0(format(round(p_community_all, 4),scientific = F), " (", n_community_all, ")")) %>%
  rowwise %>%
  filter(n_screening_all > 0 | p_community_all > 0.01) %>%
  left_join(vclade_pvalues) %>%
  select_at(c(2, (ncol(.) - 3) : ncol(.)))

lin_ids <- metadata %>% filter(!is.na(pangoLineage)) %>%
  select(nextstrainClade, pangoLineage) %>% distinct() %>% 
  group_by(pangoLineage) %>%
  mutate(`nextstrainClade | pangoLineage` = paste(ifelse(n() > 1,"-",nextstrainClade),
                                                  "|" ,pangoLineage)) %>% ungroup %>%
  select(2,3) %>% distinct()

table_vlin <- vlin  %>% 
  pivot_wider(names_from = group, values_from = c(n,p), values_fill = 0) %>%
  arrange(#desc(n_screening_male2029), 
          desc(n_screening_all),
          desc(n_community_male2029),
          desc(n_community_all)) %>%
  mutate_if(is.numeric, coalesce, 0) %>%
  mutate(#`Screening men 20-29` = paste0(round(p_screening_male2029, 4), " (", n_screening_male2029, ")"),
        `SAF Screening` = paste0(round(p_screening_all, 4), " (", n_screening_all, ")"),
        `Community\nmen 20-29` = paste0(round(p_community_male2029, 4), " (", n_community_male2029, ")"),
         `All Community` = paste0(format(round(p_community_all, 4),scientific = F), " (", n_community_all, ")")) %>%
  #rowwise %>%
  filter(n_screening_all > 1| p_community_all > 0.01) %>%
  left_join(vlin_pvalues) %>%
  left_join(lin_ids) %>%
  select_at(c(ncol(.), (ncol(.) - 4) : (ncol(.) - 1)))

print(xtable(bind_rows(table_vclade,table_vlin)), file = snakemake@output[["table_lin"]], include.rownames = F)

# debug ------------------------------------------------------------------------
# debugging = TRUE

# setClass(
#   "snakemake_object",
#   contains= "tbl_df",
#   slots = c(input = "character", output = "character", params = "character")
# )

# snakemake <- new("snakemake_object", tibble(),
#          input = c(metadata = "results/data/all/metadata.tsv",
#                    screening_weeks = "resources/screening_weeks.tsv"),
#          output = c(table_lin = "results/report/table_lineages.tex"))
