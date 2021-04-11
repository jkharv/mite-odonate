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
mites <- read_csv("datasets_derived/sequencing/mite_sequences_annotated.csv")

samples <- samples %>% 
  rename(c(sample_code = "Sample Code", odonate_spp = "Odonate Species")) %>%
  mutate(odonate_spp = str_replace_all(odonate_spp, " ", "_"))

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
  mutate(across(is.numeric, Vectorize(function(x) if(x > 0){1}else{0}))) %>% # filter pos detec. from vst
  mutate(odonate_spp = get_name(odonate_spp)) %>% # Getting the spp for each sample code
  group_by(odonate_spp) %>%
  summarise(across(everything(), sum)) %>% #Sum s/t number indicates number of times an association is detected
  mutate_if(is.numeric, funs(./sum(.))) %>%
  select_if(~ !any(is.nan(.))) %>% # Previous op causes div 0 if there are no detections.
  remove_singletons()

write_csv(mites, "datasets_derived/prob_network.csv")

mites <- mutate(mites, across(is.numeric, ceiling))
write_csv(mites, "datasets_derived/bin_network.csv")
