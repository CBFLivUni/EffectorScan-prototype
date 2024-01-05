## -------------------------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(Biostrings)
library(Peptides)
source("code/palettes.R")
source("code/utilityFunctions.R")

## -------------------------------------------------------------------------------------------

protein_anno <- read.csv("outputs/structure_QC_data.csv", row.names=1)

# Note: need to make sure this file is updated depending on the hits
id_mapping <- read.csv("annotation/uniprot-pdb_id_mapping_unique_genes_hits.csv")



## -------------------------------------------------------------------------------------------


swissprot_hits_data <- processDBQuery(db_hits_file_path = "outputs/structural_alignment/all_unique_genes_vs_swissprot/", protein_anno=protein_anno)

# Note e.value threshold set higher for pdb databases as no hits otherwise
pdb_hits_data <- processDBQuery(db_hits_file_path = "outputs/structural_alignment/all_unique_genes_vs_pdb/", protein_anno=protein_anno, e.value = 0.002)

uniprot50_hits_data <- processDBQuery(db_hits_file_path = "outputs/structural_alignment/all_unique_genes_vs_uniprot50/", protein_anno=protein_anno)

proteome_hits_data <- processDBQuery(db_hits_file_path = "outputs/structural_alignment/all_unique_genes_vs_proteome/", protein_anno=protein_anno)



## -------------------------------------------------------------------------------------------

swissprot_hits <- swissprot_hits_data$hitsDF.filtered
pdb_hits <- pdb_hits_data$hitsDF.filtered
uniprot50_hits <- uniprot50_hits_data$hitsDF.filtered
proteome_hits <- proteome_hits_data$hitsDF.filtered


swissprot_hits$target_ID <- str_match(swissprot_hits$target, "AF-\\s*(.*?)\\s*-F1")[,2]

pdb_hits$target_ID <- gsub("\\..*","", pdb_hits$target)

uniprot50_hits$target_ID <- str_match(uniprot50_hits$target, "AF-\\s*(.*?)\\s*-F1")[,2]

proteome_hits$target_ID <- str_match(proteome_hits$target, "AF-\\s*(.*?)\\s*-F1")[,2]


swissprot_hits <- swissprot_hits %>%
  mutate("database" = replicate(nrow(.), "swissprot"))

pdb_hits <- pdb_hits %>%
  mutate("database" = replicate(nrow(.), "pdb"))

uniprot50_hits <- uniprot50_hits %>%
  mutate("database" = replicate(nrow(.), "uniprot50"))

proteome_hits <- proteome_hits %>%
  mutate("database" = replicate(nrow(.), "proteome"))



## -------------------------------------------------------------------------------------------

all_results <- rbind(swissprot_hits, pdb_hits, uniprot50_hits, proteome_hits)



## -------------------------------------------------------------------------------------------

write.csv(all_results, "outputs/all_unique_genes_structural_hits_all_databases.csv")




## -------------------------------------------------------------------------------------------

all_results_sequence_anno <- all_results %>%
  dplyr::inner_join(protein_anno[,c("ID", "sequences")], join_by(query == ID)) %>% 
  dplyr::inner_join(id_mapping[,c("Entry", "Sequence", "Protein.names", "Gene.Ontology..GO.")], join_by(target_ID == Entry)) %>% 
  dplyr::rename("query_sequence" = "sequences", "target_sequence" = "Sequence", "target_protein_name" = "Protein.names", "target_GO_terms" = "Gene.Ontology..GO.")
  


## -------------------------------------------------------------------------------------------

write.fasta(sequences = as.list(all_results_sequence_anno$target_sequence), names = all_results_sequence_anno$target, file.out = "outputs/unique_genes_hits.fasta")

write.fasta(sequences = as.list(all_results_sequence_anno$query_sequence), names = all_results_sequence_anno$query, file.out = "outputs/unique_genes.fasta")



## -------------------------------------------------------------------------------------------

all_results_sequence_anno$sequence_identity_global <- mapply(globalSequenceAlignment, seq1 = all_results_sequence_anno$query_sequence, seq2 = all_results_sequence_anno$target_sequence)


all_results_sequence_anno$sequence_identity_local <- mapply(globalSequenceAlignment, seq1 = all_results_sequence_anno$query_sequence, seq2 = all_results_sequence_anno$target_sequence, local=T)



## -------------------------------------------------------------------------------------------

all_results_sequence_anno$query_molecular_weight <- all_results_sequence_anno$query_sequence %>%
  map(mw, 
  monoisotopic = FALSE,
  avgScale = "expasy",
  label = "none",
  aaShift = NULL) %>%
  do.call(rbind, .)

all_results_sequence_anno$target_molecular_weight <- all_results_sequence_anno$target_sequence %>%
  map(mw, 
  monoisotopic = FALSE,
  avgScale = "expasy",
  label = "none",
  aaShift = NULL) %>%
  do.call(rbind, .)



## -------------------------------------------------------------------------------------------

tmhelices <- read.csv("annotation/trans_membrane_helices_unique_genes.csv", header=FALSE)

tmhelices <- tmhelices %>% distinct(V1, .keep_all = TRUE)

all_results_sequence_anno <- all_results_sequence_anno %>%
  dplyr::left_join(tmhelices[,c("V1", "V5")], join_by(query == V1)) %>%
  dplyr::rename("query_num_tm_helices" = "V5")

all_results_sequence_anno$query_num_tm_helices <- as.numeric(gsub("PredHel=", "", all_results_sequence_anno$query_num_tm_helices))



## -------------------------------------------------------------------------------------------

tmhelices_hits <- read.csv("annotation/trans_membrane_helices_unique_genes_hits.csv", header=FALSE)

tmhelices_hits <- tmhelices_hits %>% distinct(V1, .keep_all = TRUE)

all_results_sequence_anno <- all_results_sequence_anno %>%
  dplyr::left_join(tmhelices_hits[,c("V1", "V5")], join_by(target == V1)) %>%
  dplyr::rename("target_num_tm_helices" = "V5")

all_results_sequence_anno$target_num_tm_helices <- as.numeric(gsub("PredHel=", "", all_results_sequence_anno$target_num_tm_helices))



## -------------------------------------------------------------------------------------------

signal_peptides <- read.csv("annotation/signalp_unique_genes.csv")

signal_peptides <- signal_peptides %>% distinct(X..ID, .keep_all = TRUE)


all_results_sequence_anno <- all_results_sequence_anno %>%
  dplyr::left_join(signal_peptides[,c("X..ID", "Prediction")], join_by(query == X..ID)) %>%
  dplyr::rename("query_signal_p" = "Prediction")



## -------------------------------------------------------------------------------------------

signal_peptides_hits <- read.csv("annotation/signalp_unique_genes_hits.csv")

signal_peptides_hits <- signal_peptides_hits %>% distinct(X..ID, .keep_all = TRUE)


all_results_sequence_anno <- all_results_sequence_anno %>%
  dplyr::left_join(signal_peptides_hits[,c("X..ID", "Prediction")], join_by(target == X..ID)) %>%
  dplyr::rename("target_signal_p" = "Prediction")



## -------------------------------------------------------------------------------------------

all_results_sequence_anno$query_cysteine_percentage <- (str_count(all_results_sequence_anno$query_sequence, "C"))/(str_length(all_results_sequence_anno$query_sequence)) * 100



## -------------------------------------------------------------------------------------------

all_results_sequence_anno$target_cysteine_percentage <- (str_count(all_results_sequence_anno$target_sequence, "C"))/(str_length(all_results_sequence_anno$target_sequence)) * 100



## -------------------------------------------------------------------------------------------

all_results_sequence_anno <- all_results_sequence_anno %>% arrange(desc(alntmscore))

write.csv(all_results_sequence_anno, "outputs/unique_genes_structural_hits_all_databases_annotated.csv")



## -------------------------------------------------------------------------------------------

all_results_species_filter <- all_results_sequence_anno %>% filter(grepl("Ascomycota|Basidiomycota", taxlineage))


write.csv(all_results_species_filter, "outputs/unique_genes_structural_hits_all_databases_annotated_species_filter.csv")



## -------------------------------------------------------------------------------------------

all_results_species_filter_non_effectors <- all_results_species_filter %>% filter(grepl("Non-effector", query_prediction))


write.csv(all_results_species_filter_non_effectors, "outputs/unique_genes_structural_hits_all_databases_non-effectors.csv")


