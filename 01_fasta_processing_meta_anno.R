## -------------------------------------------------------------------------------------------

library(Biostrings)
library(tidyverse)
library(seqinr)
library(janitor)
library(cowplot)
source("code/palettes.R")



## -------------------------------------------------------------------------------------------

# Filter
unique_effector_sequences <- readAAStringSet("raw_fasta/lycopersici_unique_with_SP.fasta")
core_effector_sequences <- readAAStringSet("raw_fasta/core_Fo_genes_with_SP.fasta")


# Annotation
unique_effector_annotation <- read.csv("annotation/effectorP_unique_genes.csv") %>%
  clean_names() %>%
  mutate("ID" = gsub(" \\|.*", "", .$x_identifier),
         "gene_type" = replicate(nrow(.), "Unique"))

core_effector_annotation <- read.csv("annotation/effectorP_core_genes.csv") %>%
  clean_names() %>%
  mutate("ID" = gsub(" \\|.*", "", .$x_identifier),
         "gene_type" = replicate(nrow(.), "Core"))



## -------------------------------------------------------------------------------------------

pangenes_unique_genes <- read.csv("annotation/getpangenes_unique_genes.csv", header=F)
pangenes_core_genes <- read.csv("annotation/getpangenes_core_genes.csv", header = F)



## -------------------------------------------------------------------------------------------

pangenes_unique_genes$V1 <- gsub("gene-", "", pangenes_unique_genes$V1)
pangenes_core_genes$V1 <- gsub("gene-", "", pangenes_core_genes$V1)

pangenes_core_genes <- pangenes_core_genes %>% 
  separate_longer_delim(c(V1), delim = ",")



## -------------------------------------------------------------------------------------------

unique_effector_sequences_anno <- data.frame(unique_effector_sequences) %>%
  rownames_to_column() %>% mutate(ID = gsub("  \\|.*", "", .$rowname)) %>%
  dplyr::inner_join(unique_effector_annotation)


core_effector_sequences_anno <- data.frame(core_effector_sequences) %>%
  rownames_to_column() %>% mutate(ID = gsub("  \\|.*", "", .$rowname)) %>%
  dplyr::inner_join(core_effector_annotation)




## -------------------------------------------------------------------------------------------

# Unique
unique_effector_sequences <- unique_effector_sequences[order(unique_effector_sequences_anno$ID)]
unique_effector_sequences_anno <- unique_effector_sequences_anno[order(unique_effector_sequences_anno$ID),]


# Core
core_effector_sequences <- core_effector_sequences[order(core_effector_sequences_anno$ID)]
core_effector_sequences_anno <- core_effector_sequences_anno[order(core_effector_sequences_anno$ID),]

#write_delim(data.frame(unique_effector_sequences_anno$ID[duplicated(unique_effector_sequences_anno$unique_effector_sequences)]), "outputs/unique_genes_duplicated_entries.txt")
#write_delim(data.frame(core_effector_sequences_anno$ID[duplicated(core_effector_sequences_anno$core_effector_sequences)]), "outputs/core_genes_duplicated_entries.txt")



## -------------------------------------------------------------------------------------------

#Unique
unique_effector_sequences <- unique_effector_sequences[!duplicated(unique_effector_sequences_anno$unique_effector_sequences)]

unique_effector_sequences_anno <- unique_effector_sequences_anno[!duplicated(unique_effector_sequences_anno$unique_effector_sequences),]


# Core
core_effector_sequences <- core_effector_sequences[!duplicated(core_effector_sequences_anno$core_effector_sequences)]

core_effector_sequences_anno <- core_effector_sequences_anno[!duplicated(core_effector_sequences_anno$core_effector_sequences),]



## -------------------------------------------------------------------------------------------

shared_unique_genes <- intersect(unique_effector_sequences_anno$ID, pangenes_unique_genes$V1)

unique_effector_sequences_anno <- unique_effector_sequences_anno[unique_effector_sequences_anno$ID %in% shared_unique_genes, ]


unique_effector_sequences <- unique_effector_sequences[unique_effector_sequences_anno$rowname]



## -------------------------------------------------------------------------------------------

shared_core_genes <- intersect(core_effector_sequences_anno$ID, pangenes_core_genes$V1)

core_effector_sequences_anno <- core_effector_sequences_anno[core_effector_sequences_anno$ID %in% shared_core_genes, ]


core_effector_sequences <- core_effector_sequences[core_effector_sequences_anno$rowname]



## -------------------------------------------------------------------------------------------

mclm::write_txt(shared_unique_genes, "outputs/unique_genes_IDs.txt", line_glue = "\n") 
mclm::write_txt(shared_core_genes, "outputs/core_genes_IDs.txt", line_glue = "\n") 



## -------------------------------------------------------------------------------------------


names(core_effector_sequences) == core_effector_sequences_anno$rowname
names(unique_effector_sequences) == unique_effector_sequences_anno$rowname



## -------------------------------------------------------------------------------------------

# Create vector of suitable filepaths/file names
core_filenames <- paste0("processed_fasta/core_genes/", core_effector_sequences_anno$ID, ".fasta") 
# Create a list object of the vectors to be passed to pmap
core_fastas <- list(core_effector_sequences, core_effector_sequences_anno$rowname, core_filenames)


# Create vector of suitable filepaths/file names
unique_filenames <- paste0("processed_fasta/unique_genes/", unique_effector_sequences_anno$ID, ".fasta") 
# Create a list object of the vectors to be passed to pmap
unique_fastas <- list(unique_effector_sequences, unique_effector_sequences_anno$rowname, unique_filenames)


# Use pmap to run 'write.fasta' function on the fastas list
pmap(core_fastas, write.fasta)
pmap(unique_fastas, write.fasta)



## -------------------------------------------------------------------------------------------

colnames(unique_effector_sequences_anno)[2] <- "sequences"
colnames(core_effector_sequences_anno)[2] <- "sequences"

effector_anno <- rbind(unique_effector_sequences_anno, core_effector_sequences_anno)
effector_anno$prediction <- gsub(" effector", "", effector_anno$prediction)
effector_anno$prediction <- gsub("Cytoplasmic/apoplastic", "Apoplastic/cytoplasmic", effector_anno$prediction)



effector_anno[effector_anno == "-"] <- NA
effector_anno[,5:7] <- lapply(effector_anno[,5:7], gsub, pattern = "[^0-9.-]", replacement = "")
effector_anno[,5:7] <- lapply(effector_anno[,5:7], as.numeric)

effector_anno$effector_type_score <- rowMeans(effector_anno[,5:7], na.rm = TRUE)


write.csv(effector_anno, "fusarium_oxysporum_f._sp._lycopersici_effector_annotation.csv")



## -------------------------------------------------------------------------------------------

effector_score_plot <- ggplot(effector_anno, aes(x=prediction, y=effector_type_score, fill=prediction)) +
  geom_violin(width=0.8, alpha = 0.3, color = NA) +
  geom_boxplot(width = 0.1, alpha = 0.8, na.rm = TRUE)+
  scale_fill_manual(values = project_palettes$effector) + 
  facet_wrap(~gene_type)+
  theme_cowplot()+
  theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 40,  hjust=1))+
  labs(color="Effector type", fill="Effector type", y="Effector type score") +
  theme(legend.position = "none")

effector_score_plot



## -------------------------------------------------------------------------------------------

summary_table <- effector_anno %>%
  tabyl(gene_type, prediction) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns()


summary_table



## -------------------------------------------------------------------------------------------

ggsave(filename = "figures/01_fasta_processing_meta_anno/effector_score_plot.pdf",
  plot = effector_score_plot,
  width = 7,
  height = 4,
  units = c("in"),
  dpi = 300,
  limitsize = TRUE,
  bg = "white")


write.csv(summary_table, "outputs/effector_type_count_table.csv")



## -------------------------------------------------------------------------------------------


