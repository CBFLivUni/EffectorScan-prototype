# Function to process log files
extractQCData <-function(input_files){
  
  # Read in and process data
  log_data <- input_files %>% 
    map(read_delim, col_names = c("date", "time", "model", "pLDDT", "pTM"))
  
  names(log_data) <- basename(input_files)
  
  # Remove rows of jobs that have failed and
  # create vector of failed jobs
  failed_jobs <- names(log_data[sapply(log_data, nrow) == 0])
  failed_jobs <- gsub(".txt", "", failed_jobs)
  log_data <- log_data[sapply(log_data, nrow)>0] 
  
  # Bind data
  log_data <- log_data %>% rbindlist(idcol = TRUE) 
  
  # Format PTM and pLDDT columns
  log_data$ID <- gsub(".txt", "", log_data$.id)
  log_data$pLDDT <- as.numeric(gsub("pLDDT=", "", log_data$pLDDT))
  log_data$pTM <- as.numeric(gsub("pTM=", "", log_data$pTM))
  
  # Rename old ID column to be file name instead
  log_data <- log_data %>% rename("file_name" = ".id")
  
  return(list("logData" = log_data, "failedJobs" = failed_jobs))
  
}


# Function to process structural alignment files
extractAlignmentData <-function(input_files, col_names = c("query", "target", "alntmscore", "lddt", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits")){
  
  # Read in and process data
  sa_data <- input_files %>% 
    map(read_delim, col_names = col_names)
  
  names(sa_data) <- basename(input_files)
  
  # Remove (failed) items that have more columns than specified
  sa_data <- sa_data[sapply(sa_data, ncol)==length(col_names)] 
  
  # Bind data
  sa_data <- sa_data %>% rbindlist(idcol = FALSE) 
  
  
  # Remove '.pdb' from query and target column
  sa_data$query <- gsub(".pdb", "", sa_data$query)
  sa_data$target <- gsub(".pdb", "", sa_data$target)
  
  return(sa_data)
  
}



# Function to process database queries for effectors
processDBQuery <-function(db_hits_file_path = "", protein_anno=protein_anno, col_names = c("query", "target", "alntmscore", "lddt", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "taxid", "taxname", "taxlineage"), e.value = 0.001){
  
  # Get list of files from file path
  db_hits_input_files <- list.files(db_hits_file_path, pattern = ".m8", full.names = TRUE) 
  
  # Run the function to extract the alignment data from the file list
  db_hits_data <- extractAlignmentData(db_hits_input_files, col_names = col_names)
  
  # Filter hits to only include fungal species
  db_hits_data <- db_hits_data %>% filter(grepl('Fungi', taxlineage))
  
  # Join to protein annotation to include query effector type, effector score and effector pLDDT
  db_hits_data_anno <- db_hits_data %>% 
    dplyr::inner_join(protein_anno[,c("ID", "effector_type_score", "prediction")], join_by(query == ID)) %>%
    dplyr::rename("query_prediction" = "prediction", "query_effector_type_score" = "effector_type_score") %>% 
    arrange(desc(alntmscore))
  
  # Create filtered version of dataframe that only has
  # TM score >0.5, evalue < 0.001
  db_hits_data_filter <- db_hits_data_anno %>%
    filter(alntmscore >= 0.5 & evalue < e.value) 
  
  return(list(hitsDF = db_hits_data_anno, hitsDF.filtered = db_hits_data_filter))
  
}



# Function to process database queries for effectors
globalSequenceAlignment <-function(seq1 = df$x, seq2 = df$y, local=F){
  
  data(BLOSUM62)
  
  if(local){
    
    global_align <- pairwiseAlignment(as.character(seq1),
                                      as.character(seq2),
                                      substitutionMatrix = BLOSUM62, 
                                      gapOpening = -2,
                                      gapExtension = -8, 
                                      scoreOnly = FALSE,
                                      type="local")
    
  } else {
    
    global_align <- pairwiseAlignment(as.character(seq1),
                                      as.character(seq2),
                                      substitutionMatrix = BLOSUM62, 
                                      gapOpening = -2,
                                      gapExtension = -8, 
                                      scoreOnly = FALSE)
    
    
  }
  
  
  sequence_identity <- pid(global_align)/100
  return(sequence_identity)
  
}
