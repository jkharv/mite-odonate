# Author: Jake Harvey - jakekharvey@gmail.com
#
# Take all the disparate CSV files from different projects and combine them
# all into one. Data cleaning too.

library(tidyverse)
library(lubridate)

a2015 <- read_csv("datasets_primary/2015_data.csv", 
                  na = c("", "NA", "N/A", "n/a", "NA.2", "NA.1"))
a2019 <- read_csv("datasets_primary/2019_data.csv", 
                  na = c("", "NA", "N/A", "n/a"))
a2020 <- read_csv("datasets_primary/2020_data.csv", 
                  na = c("", "NA", "N/A", "n/a"))
masses <- read_csv("datasets_primary/2019-2020-mass.csv", 
                   na = c("", "NA", "N/A", "n/a"))
tania_ids <- read_csv("datasets_primary/tania_id_mess.csv", 
                      na = c("", "NA", "N/A", "n/a"))

# 2019 uses a different date format and the year is left implicit.
a2019$Date <- a2019$Date %>%
  paste("_2019", sep = "") %>%
  mdy()

# 2019 also says Lake rather than site
a2019 <- rename(a2019, Site = Lake)

# 2019 also says Lake rather than site
a2020 <- rename(a2020, Date = "Date (YMD)")

# Combining both years of the immune study.
both <- intersect(names(a2019), names(a2020))
a2019 <- select(a2019, all_of(both))
a2020 <- select(a2020, all_of(both))
immune_study <- union(a2019, a2020)

# Unfortunately, some misguided person decided to delete all the cols they
# didn't immediately need in the mass dataset so we can't unambiguously
# decide which specimens in each dataset represent the same specimen.

# Produce a hopefully unique ID with all the info we have left in the mass data.
immune_study <- unite(immune_study, temp_id, 
                      c(Site, Specimen_Number, Species, Sex), sep="_",
                      remove = FALSE)

# Produce the same ID in the mass dataset
masses <- masses %>%
  #Some had .2 .3 ... appended to make them unique. Must remove.
  separate(Sequence_num, into = c("seq_num", "suffix"), "\\.") %>%
  unite(temp_id, c(seq_num, Species, Sex), sep="_", remove = FALSE)

# We have only 32 which are ambiguous with this ID.
counts_masses <- count(masses, temp_id)
ambiguous_ids <- filter(counts_masses, n > 1)
masses_unambiguous <-  masses %>%
  filter(!(temp_id %in% ambiguous_ids$temp_id)) %>%
  select(c("temp_id", "Weight(g)"))

# Join all the mass measurements which can unambiguously be associated with the
# specimen they belong to.
immune_study <- immune_study %>%
  left_join(masses_unambiguous, by = "temp_id", keep = FALSE)
# I cooulldd do those 64 specimens manually, but I won't

immune_study$Identifier <- "Maggie Blondeau"

# Now we can join it to the 2015 data.

# Checking which names need to be changed
intersect(names(a2015), names(immune_study))

# Renaming 2015 columns which need renaming.
a2015 <- a2015 %>%
  rename(Site = SiteName) %>%
  rename(Mite_Abundance = `Mite Load`) %>%
  rename(Specimen_Number = SpecimenNumber) %>%
  rename(Identifier = Author) %>%
  rename(Mass = "Mass (g)") %>%
  rename(Note = "Notes")


# Both datasets need to have Date type for the date col.
a2015$Date <- dmy(a2015$Date)

# Renaming columns in immune_study which need renaming.
immune_study <- immune_study %>%
  rename(Mass = "Weight(g)")



odonates <- full_join(a2015, immune_study)

# Now we can do some cleaning.
odonates <- odonates %>% 
  drop_na(Species) %>%
  filter(Species != "sp") %>%
  mutate(Species = str_replace(Species, "obstrusum", "obtrusum")) %>%
  mutate(Species = str_replace(Species, "obtusum", "obtrusum")) %>%
  # We'll use the old taxonomy here cause it's simpler to change and the 
  # phylogeny is already using it.
  mutate(Genus = str_replace(Genus, "Phanogomphus", "Gomphus")) %>%
  mutate(Genus = str_replace(Genus, "Agrigomphus", "Gomphus")) %>%
  # Someone put Gomphus rather than Gomphidae in a bunch of places.
  mutate(Family = str_replace(Family, "Gomphus", "Gomphidae")) %>%
  mutate(Family = str_replace(Family, "coenagrionidae", "Coenagrionidae")) %>%
  mutate(Family = str_replace(Family, "Coenagionidae", "Coenagrionidae")) %>%
  mutate(Family = str_replace(Family, "Corduligastridae", "Cordulegastridae")) %>%
  unite(Genus_Species, c(Genus, Species), sep = "_", remove = FALSE)

cols <- setdiff(names(odonates), 
                c("temp_id", "X48", "Mites Removed?", "G_species", "Life stage"))

odonates <- odonates %>%
  select(all_of(cols)) %>%
  rename_all(tolower) %>%
  rename_with(function(x) str_replace_all(x, " ", "_")) %>%
  rename(water_temperature = temperature) %>%
  relocate(c(date, site, sitecode, station, latitude, longitude, 
             specimen_number, water_temperature, ph, suborder, family, genus, 
             species, genus_species, sex, mass, mite_abundance, identifier))

write_csv(odonates, "datasets_derived/odonates_merged.csv")
