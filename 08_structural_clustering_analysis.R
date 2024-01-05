## -------------------------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(cowplot) 
library(smplot2)
library(ggpubr)
source("code/palettes.R")
source("code/utilityFunctions.R")


## -------------------------------------------------------------------------------------------

# Core gene and unique gene all vs all structural alignments
unigenes_input_files <- list.files("outputs/structural_alignment/unique_genes/",pattern = ".m8", full.names = TRUE) 

core_genes_input_files <- list.files("outputs/structural_alignment/core_genes/",pattern = ".m8", full.names = TRUE) 

# Effector annotation
sequence_anno <- read.csv("fusarium_oxysporum_f._sp._lycopersici_effector_annotation.csv", row.names=1)



## -------------------------------------------------------------------------------------------

# Run 'extractQCdata' function on input files and extract data
unigenes_sa_data <- extractAlignmentData(unigenes_input_files)



## -------------------------------------------------------------------------------------------

# Annotate the QC data with the effector information
unigenes_sa_data_anno <- unigenes_sa_data %>% 
  dplyr::inner_join(sequence_anno[,c("ID", "effector_type_score", "prediction")], join_by(query == ID)) %>%
  dplyr::inner_join(sequence_anno[,c("ID", "prediction")], join_by(target == ID)) %>%
  rename("query_prediction" = "prediction.x", "target_prediction" = "prediction.y", "query_effector_type_score" = "effector_type_score") %>% 
  arrange(desc(alntmscore))


write.csv(unigenes_sa_data_anno, "outputs/unigene_genes_structural_alignment_results.csv")



## -------------------------------------------------------------------------------------------

unigenes_sa_data_anno_filter <- unigenes_sa_data_anno %>%
  filter(alntmscore != 1 & query_prediction != "Non-effector" & alntmscore >= 0.5 & evalue < 0.001) 

write.csv(unigenes_sa_data_anno_filter, "outputs/unigene_genes_structural_alignment_results_filtered.csv")


