# Author: Jake Harvey jakekharvey@gmail.com
# Make a binary (or probabilistic) host-parasite network

library(tidyverse)

get_name <- function(x){
  n <- filter(samples, sample_code == substring(x, 2))
  if(nrow(n) == 1){
    return(n$odonate_spp)
  } else {
    return(NA)
  }
}
get_name <- Vectorize(get_name)

samples <- read_csv("datasets_primary/mite_samples.csv")
mites <- read_csv("datasets_derived/sequencing/mite_sequences_annotated.csv")

samples <- samples %>% 
  rename(sample_code   = "Sample Code", 
         odonate_spp            = "Odonate Species", 
         sampling_site          = "Sampling Site", 
         mites_no               = "Est Num. of Mites", 
         odonates_no            = "Num. of Odonates", 
         removal_date           = "Mite Removal Date", 
         removed_by             = "Mites Removed By", 
         sample_no              = "Sample No.", 
         extraction_date        = "DNA Extraction date", 
         kit_used               = "Kit Used", 
         preservation_treatment = "Preservation treatment", 
         comments               = "Comments") %>%
  mutate(odonate_spp = str_replace_all(odonate_spp, " ", "_")) %>%
  mutate(sampling_site = str_replace_all(sampling_site, " ", "_"))

controls <- mites %>%
  select(contains("neg") | "phylum" | "class" | "order" | "family" | "genus" | "species")

write_csv(controls, "datasets_derived/controls.csv")

mites <- mites %>%
  select(matches("(X[[:digit:]]{4}_[[:digit:]]{1,2})")) %>% #Exclude the control samples
  t() %>%
  as.data.frame() %>%
  rownames_to_column("odonate_spp") %>%
  mutate(odonate_spp = get_name(odonate_spp)) %>% # Getting the spp for each sample code
  mutate(across(is.numeric, Vectorize(function(x) if(x>1000){1}else{0}))) %>% #Threshold, reads -> presence
  group_by(odonate_spp) %>%
  summarise(across(everything(), sum)) %>% #Sum s/t number indicates number of times an association is detected
  mutate_if(is.numeric, funs(./sum(.))) 

write_csv(mites, "datasets_derived/prob_network.csv")

mites <- mutate(mites, across(is.numeric, ceiling))
write_csv(mites, "datasets_derived/bin_network.csv")