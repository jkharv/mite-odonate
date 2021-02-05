# Author: Jacob Harvey - jakekharvey@gmail.com

library(tidyverse)

a2015 <- read_csv("datasets_primary/2015_data.csv", na = c("", "NA", "N/A"))
a2019 <- read_csv("datasets_primary/2019_data.csv", na = c("", "NA", "N/A"))
a2020 <- read_csv("datasets_primary/2020_data.csv", na = c("", "NA", "N/A"))

a2015 <- a2015 %>%
  mutate(species = paste(Genus, Species, sep = "_")) %>%
  rename(Mite_Abundance = "Mite Load") %>%
  filter(species != "Lestes sp" & species != "Aeshna sp") %>%
  drop_na(c(Species, Mite_Abundance)) %>%
  select(c(species, Mite_Abundance))

# Julie used old taxonomy without Gomphus split up. I'm using Gomphus in the analysis.
# Gomphurus, Hylogomphus, Phanogomphus, and Stenogomphurus -> Gomphus
a2020 <-  a2020 %>%
  mutate_if(is.character, str_replace_all, 
            pattern = "Gomphurus|Hylogomphus|Phanogomphus|Stenogomphurus", 
            replacement = "Gomphus") %>%
  mutate(species = paste(Genus, Species, sep = "_")) %>%
  drop_na(c(Species, Mite_Abundance)) %>%
  select(c(species, Mite_Abundance))

a2019 <-  a2019 %>%
  mutate_if(is.character, str_replace_all, 
            pattern = "Gomphurus|Hylogomphus|Phanogomphus|Stenogomphurus", 
            replacement = "Gomphus") %>%
  mutate(species = paste(Genus, Species, sep = "_")) %>%
  drop_na(c(Species, Mite_Abundance)) %>%
  select(c(species, Mite_Abundance))
  
# Julie's data used a standard sampling protocol.
abundances <- a2015 %>%
  group_by(species) %>%
  tally() %>%
  rename(abundance = n)

all_data <- union_all(union_all(a2015, a2019), a2020)

all_data <- all_data %>%
  mutate(abundance = 1) %>%
  group_by(species) %>%
  summarise_all(sum) %>%
  mutate(prevalence = Mite_Abundance/abundance)

write_csv(abundances, "datasets_derived/odonate_abundances.csv")
write_csv(abundances, "datasets_derived/mite_prevalences.csv")