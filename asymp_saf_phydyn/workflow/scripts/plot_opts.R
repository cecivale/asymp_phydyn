# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project: cov-armee Phylodynamics
#        V\ Y /V    ggplot theme and color palettes
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

library(tidyverse)
library(ggpattern)
library(scales)
library(cowplot)

theme_asympsymp <- function() {
  theme_minimal() +
    theme(#panel.grid.major.y = element_blank(),
          #panel.grid.minor.y = element_blank(),
          #panel.grid.major.x = element_blank(),
          #panel.grid.minor.x = element_blank(),
          legend.position = "none",
          text = element_text(size = 8, color = "black"),
          axis.text = element_text(size = 8, color = "black"),
          axis.title = element_text(size = 8, color = "black", face = "bold"),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_line(size = 0.25),
          panel.grid.major = element_line(size = 0.4),
          panel.grid.minor = element_line(size = 0.25)
          #axis.ticks.length = unit(2, "mm")
          )
}

pal_asympsymp <- c("asymp" = "#264c38", "symp" = "#e28742", 
                       "asymp2" = "#6f9382", "symp2" = "#FBCF96" , 
                       "grey" = "grey70")

colors_gradient_screening <- c("grey60","white", "#6f9382", "#264c38")
colors_gradient_community <- c("grey60","white", "#FBCF96", "#e28742")

