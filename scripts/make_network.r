# Author: Jake Harvey jakekharvey@gmail.com
# Make a binary (or probabilistic) host-parasite network

library(DESeq2)
library(tidyverse) # Must load after DESeq2 because of conflict with DESeq2 dependancy
library(vegan)

get_name <- function(x){
  n <- filter(samples, sample_code == substring(x, 2))
  if(nrow(n) == 1){
    return(n$odonate_spp)
  } else {
    return(NA)
  }
}
get_name <- Vectorize(get_name)

get_num_mites <- function(x){
  n <- filter(samples, sample_code == substring(x, 2))
  if(nrow(n) == 1){
    return(n$mite_num)
  } else {
    return(NA)
  }
}
get_num_mites <- Vectorize(get_num_mites)

remove_singletons <- function(mites){
  
  m <- mutate(mites, across(is.numeric, ceiling))
  
  keep <- m %>%
    select(2:last_col()) %>%
    pivot_longer(1:last_col()) %>%
    group_by(name) %>%
    summarize(count = sum(value)) %>%
    filter(count > 1)
  
  mites <- select(mites, odonate_spp | keep$name)
  return(mites)
}

samples <- read_csv("datasets_primary/mite_samples.csv")

samples <- samples %>% 
  rename(c(sample_code = "Sample Code", odonate_spp = "Odonate Species",
           mite_num = "Est Num. of Mites")) %>%
  mutate(odonate_spp = str_replace_all(odonate_spp, " ", "_"))

make_network <- function(seq_file) {
  
  path <- paste("datasets_derived/sequencing/", seq_file, sep = "")
  mites <- read_csv(path)
  
  controls <- mites %>%
    select(contains(c("neg" , "phylum" , "class" , "order" , "family" , "genus" , "species")))
  
  write_csv(controls, "datasets_derived/controls.csv")
  
  # Get rid of junk cols, only need the asv counts
  mites <- mites %>%
    select(!contains(c("sequence", "phylum" , "class" , "order" , "family" , "genus" , "species"))) %>%
    select(2:last_col())
  
  mites <- varianceStabilizingTransformation(as.matrix(mites) + 1, fitType = "local")
  
  mites <- mites %>% 
    as.data.frame() %>%
    select(matches("(X[[:digit:]]{4}_[[:digit:]]{1,2})")) %>% #Exclude the control samples
    t() %>%
    as.matrix()
  
  mites <- mites %>%
    as.data.frame() %>%
    rownames_to_column("odonate_spp") %>%
    mutate(across(is.numeric, Vectorize(function(x) if(x > 0){1}else{0}))) %>%# filter pos detec. from vst
    filter(rowSums(across(where(is.numeric))) <= get_num_mites(odonate_spp)) # Remove samples with more sequences than mites 
  
  mites <- mites %>%
    mutate(odonate_spp = get_name(odonate_spp)) %>% # Getting the spp for each sample code
    group_by(odonate_spp) %>%
    summarise(across(everything(), sum)) %>% #Sum s/t number indicates number of times an association is detected
    mutate_if(is.numeric, funs(./sum(.))) %>%
    select_if(~ !any(is.nan(.))) %>% # Previous op causes div 0 if there are no detections.
    remove_singletons()
  
  out_path <- paste("datasets_derived/", str_extract(seq_file, "[:graph:]+(?=\\.)"),
                    "_prob_network.csv", sep = "")
  write_csv(mites, out_path)
  
  mites <- mutate(mites, across(is.numeric, ceiling))
  out_path <- paste("datasets_derived/", str_extract(seq_file, "[:graph:]+(?=\\.)"),
                    "_bin_network.csv", sep = "")
  write_csv(mites, out_path)
  
}

make_network("mite_sequences_annotated.csv")
make_network("mite_sequences_otu90_annotated.csv")
make_network("mite_sequences_otu97_annotated.csv")
