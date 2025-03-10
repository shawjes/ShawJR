# Script to run at start

#### Load packages

#install.packages("devtools")
library(devtools)
# https://github.com/SomaLogic/SomaDataIO
#devtools::install_github("SomaLogic/SomaDataIO")
library(SomaDataIO)
library(ggrepel)
#library(limma)
library(dplyr)
library(tidyr)
library(data.table)
library(broom)
#library(broomExtra)
library(tibble)
library(sjstats)
library(car)
library(lme4)
library(lmerTest)
library(ggplot2)
library(tibble)
library(modelr)
library(tidyverse)
#library(miceadds)
library(ggforce)
require(openxlsx)
library(tidyverse)
#library(caret)
#library(glmnet)
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
#install.packages("reshape")
library(reshape)
library(biomaRt)
library(readxl)    
#library(synapser)

# BiocManager::install("SomaScan.db")
# BiocManager::install("GO.db")
library(dplyr)
library(GO.db)
library(SomaDataIO)
library(SomaScan.db)
library(furrr)
#install.packages("ggfortify")
library(ggfortify)
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations

#### Setting and modifying theme for plots
theme_set(theme_gray(base_size = 12, base_family = "Arial") +
            theme(panel.border = element_rect(colour="black", fill = "transparent"),
                  plot.title = element_text(face="bold", hjust = 0), # lineheight=.8, size=20,
                  axis.text = element_text(color="black", size = 11), #size = 14
                  axis.text.x = element_text(angle = 0, hjust = NULL),
                  strip.background = element_rect(colour="black", fill = "light grey", size = 1), # adjusts facet label borders (if any)
                  panel.background = element_blank(),
                  panel.grid = element_blank()
            ))
RedBlue <- c("#CD3333", "#1874CD")
GrayBlue <- c("grey", "#2b8cbe")

#### Define helper functions
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

`%notin%` <- Negate(`%in%`)

#select <- dplyr::select
#filter <- dplyr::filter

# Setting and modifying theme for plots
theme_set(theme_gray(base_size = 12, base_family = "Arial") +
            theme(panel.border = element_rect(colour="black", fill = "transparent"),
                  plot.title = element_text(face="bold", hjust = 0, size = 12), # lineheight=.8, size=20,
                  axis.text = element_text(color="black", size = 12), 
                  axis.text.x = element_text(angle = 0, hjust = NULL),
                  strip.background = element_rect(colour="black", fill = "light grey", size = 1), # adjusts facet label borders (if any)
                  panel.background = element_blank(),
                  panel.grid = element_blank()
            ));

#####################################################################
##  Define custom color palette for proteomic clock clusters 1-16. ##
#####################################################################
# Lehallier 2019 Cluster  Color
# Cluster 1               Red (Set1[[1]])
# Cluster 2               Medium blue (Set1[[2]])
# Cluster 3               Green (Set1[[3]])
# Cluster 4               Purple (Set1[[4]])
# Cluster 5               Orange (Set1[[5]])
# Cluster 6               Brown (Set1[[7]])
# Cluster 7               Baby pink (Set1[[8]])
# Cluster 8               Dark gray (Set1[[9]])
#1,2,3,4,5,7,8,9
P4C_extended_cluster_colors <- c(
  # Colors for Clock Clusters 1-8 as published in Lehallier et al., 2019:
  brewer.pal(9, "Set1")[[1]], # Red
  brewer.pal(9, "Set1")[[2]], # Medium blue
  brewer.pal(9, "Set1")[[3]], # Green
  brewer.pal(9, "Set1")[[4]], # Purple
  brewer.pal(9, "Set1")[[5]], # Orange
  brewer.pal(9, "Set1")[[7]], # Brown
  brewer.pal(9, "Set1")[[8]], # Baby pink
  brewer.pal(9, "Set1")[[9]], # Dark gray
  # Colors for Clock Clusters 9-18, derived from clustering of remaining Aptamers among D21s: 
  brewer.pal(8, "Set2")[[1]], # seafoam green
  brewer.pal(8, "Set2")[[6]], # gold
  brewer.pal(8, "Set2")[[7]], # camel
  brewer.pal(8, "RdGy")[[4]], # light orange
  brewer.pal(12, "Set3")[[10]], # light purple
  brewer.pal(8, "Set2")[[5]], # lime green,
  brewer.pal(12, "Set3")[[5]], # light blue
  brewer.pal(11, "PiYG")[[2]] # magenta
) 
# IMPORTANT NOTE: We *could* double back and cluster only the proteins associated with Age in D21s, but I want to cluster all proteins not classified by Lehallier et al. 2019, since T21s may have associations with Age that D21s do not.
#length(P4C_extended_cluster_colors)
#[1] 16
#P4C_extended_cluster_colors

