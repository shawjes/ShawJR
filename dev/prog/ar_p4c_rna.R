

#```{r echo=True, include=FALSE}
# For data management
# install.packages('tidyverse')
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# # For visualisation
# install.packages('pheatmap')
# install.packages("DOSE")
# install.packages("enrichplot")
# install.packages("ggupset")

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
#```

#### Run startup function
#```{r echo=True, include=FALSE}
dir.project <- "~/Dropbox/ShawJR/2025/dev"

setwd(paste0(dir.project, "/macro"))
source("ui_init.R")
library(factoextra)
#```

#### Define paths to directories
#```{r}
dir.project <- "~/Dropbox/ShawJR/2025/dev"
dir.donovan2024_source <- paste0(dir.project, "/rawdata/Donovan2024_syn31481952.5")
dir.indata <- paste0(dir.project, "/indata")
dir.aidata <- paste0(dir.project, "/aidata")
dir.ardata <- paste0(dir.project, "/ardata")
dir.ddata <- paste0(dir.project, "/ddata")
dir.docs <- paste0(dir.project, "/docs")
dir.output <- paste0(dir.project, "/output")

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
#```

#### Read Analysis Tracker
#```{r}
# setwd(dir.docs)
# analysis_tracker <- read.xlsx("AnalysisTracker.xlsx") %>%
#   filter(run_model == 1)
# 
# analysis_tracker
#```

#### Read requisite aidata
#```{r}
setwd(dir.aidata)
ai_p4c_rna <- fread("ai_p4c_rna.csv.gz")
ai_p4c_rna_d21_mean_sd <- fread("ai_p4c_rna_d21_mean_sd.csv.gz")
ai_p4c_meta <- fread("ai_p4c_meta.csv.gz")

ai_p4c_rna %>% head()
ai_p4c_rna_d21_mean_sd %>% head()
ai_p4c_meta %>% head()
#```

#```{r}
ai_p4c_rna %>% head()

ai_p4c_rna %>% glimpse()

ai_p4c_rna %>% filter(is.na(Value) | Value == -Inf) %>% head()

ai_p4c_rna$Value %>% summary()
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.000     0.000     0.000     6.695     0.655 25909.685

log10(0.005)
#```

#```{r}
ai_p4c_rna %>% head()
ai_p4c_rna_d21_mean_sd %>% head()
ai_p4c_meta %>% head()
#```


#### ai_p4c_rna01
#```{r}
ai_p4c_rna01 <- ai_p4c_rna %>%
  mutate(In_P4C_RNAseq = 1) %>%
  left_join(ai_p4c_rna_d21_mean_sd,
            by = "Analyte") %>%
  #na.omit() %>%
  group_by(LabID, Analyte) %>%
  mutate(Zlog10_FPKM = (log10(Value + 0.005) - d21_mean_log10FPKM)/d21_sd_log10FPKM) %>%
  ungroup()

ai_p4c_rna01 %>% head(n=100)
ai_p4c_rna01$Value %>% summary()
ai_p4c_rna01$Zlog10_FPKM %>% summary()


#```

#### ar_p4c_rna
#```{r}
ar_p4c_rna <- ai_p4c_rna01 %>%
  left_join(ai_p4c_meta,
            by = "LabID") %>%
  dplyr::select(-c(Value)) %>%
  dplyr::rename(Value = Zlog10_FPKM) %>%
  mutate(Units = "Zlog10_FPKM",
         Omics_type = "Transcriptomics",
         Analyte_type = "Gene expression") %>%
  filter(!is.na(In_P4C_RNAseq)) %>%
  select(LabID, Karyotype, Sex, Age,
         #T21, Female, Karyotype_Sex,
         Omics_type, Analyte_type, Analyte, Value, Units, In_P4C_RNAseq);

setwd(dir.ardata)
fwrite(ar_p4c_rna, "ar_p4c_rna.csv.gz")
#```


#```{r}
ar_p4c_rna %>%
  filter(grepl("AACSP1", Analyte)==TRUE)
#```

