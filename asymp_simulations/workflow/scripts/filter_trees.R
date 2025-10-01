# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee Simulation
#        V\ Y /V    Split trees in one file per tree
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

# Load libraries ---------------------------------------------------------------
library(treeio)
library(stringr)
library(dplyr)
library(readr)


# Read, split, write -----------------------------------------------------------
cat("Reading trees...")
sim_trees <- read.beast(file = snakemake@input[["trees"]])
cat("done!")

if (length(sim_trees) == 1) sim_trees <- list(sim_trees)

trees_l <- lapply(sim_trees, as_tibble)
trees <- bind_rows(trees_l, .id = "column_label")

selected_trees <- trees %>% filter(!is.na(samp)) %>%
  group_by(column_label) %>%
  count(column_label, name = "tips") %>%
  filter(tips >= snakemake@params[["min_tips"]] & tips <= snakemake@params[["max_tips"]]) %>%
         #tips >= 50 & tips <= 350) %>%
  select(column_label)

cat("\n", paste0(nrow(selected_trees), " trees selected from ", length(sim_trees), 
                 " simulated trees with more than ", snakemake@params[["min_tips"]], " and less than ", 
                 snakemake@params[["max_tips"]], " tips.\n" ))

# Extract all time values from a tree with maximum precision
lines <- readLines(snakemake@input[["trees"]])
tree_lines <- grep("^tree STATE_", lines, value = TRUE)

extract_max_time <- function(tree_line) {
  matches <- gregexpr("time=([0-9]+\\.[0-9]+)", tree_line)
  times <- regmatches(tree_line, matches)[[1]]
  numeric_times <- as.numeric(sub("time=", "", times))
  max(numeric_times)
}

max_times <- sapply(tree_lines, extract_max_time)
names(max_times) <- paste0("STATE_", seq_along(max_times) - 1)

# Create tree directory
ifelse(!dir.exists(file.path(snakemake@output[[1]])),
       dir.create(file.path(snakemake@output[[1]])),
       "Directory Exists")

t <- lapply(1:nrow(selected_trees), function(i) {
  tree_state <- selected_trees$column_label[i]
  tree_idx <- gsub("STATE_", "", tree_state)
  mrs <- max_times[tree_state]
  file_name <- paste0(snakemake@params[["out_file"]], tree_idx, ".txt")
  write_lines(sprintf("%.16f", mrs), file = file_name)
  })



# EDA Plot and count, percentage
# df <- trees%>% filter(!is.na(samp))  %>%
#   add_count(column_label, type, name = "typed_tips") %>%
#   add_count(column_label, name = "tips") %>%
#   select(column_label, type, typed_tips, tips) %>% distinct()
# df %>% select(column_label, tips) %>% distinct() %>%
#   ggplot() +
#   geom_histogram(aes(tips, fill = tips >= 50 & tips <= 500))
# 
# df %>% select(column_label, tips) %>% distinct() %>%
#   count(tips >= 50, tips <= 350)
# 
# df %>%  filter(tips >= 50 & tips <= 350, type == "type1") %>%
#   ggplot() +
#   geom_histogram(aes(typed_tips))
# 
# df %>% select(column_label, type, tips, typed_tips) %>% distinct() %>%
#   filter(type == "type1") %>%
#   count(typed_tips >= 20 & tips >= 50 & tips <= 350)



