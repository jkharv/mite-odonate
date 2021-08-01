# Author: Jake Harvey - jakekharvey@gmail.com
#
# Take all the disparate CSV files from different projects and combine them
# all into one. Data cleaning too.

library(tidyverse)
library(lubridate)

a2015 <- read_csv("datasets_primary/2015_data.csv", 
                  na = c("", "NA", "N/A", "n/a", "NA.2", "NA.1"))

immune_resp <- read_csv("datasets_primary/Blondeau.BIOL490.Data.csv")
immune_resp$Identifier <- "Maggie Blondeau"

# Checking which names need to be changed
intersect(names(a2015), names(immune_resp))

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
immune_resp$Date <- parse_date_time(immune_resp$Year, "Y")

# Renaming columns in immune_study which need renaming.
immune_resp <- immune_resp %>%
  rename(Mass = "Weight")

odonates <- full_join(a2015, immune_resp)

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

odonates <- odonates %>%
  select(Date, Site, Latitude, Longitude, SubOrder, Family, Genus, Species,
         Mass, Mite_Abundance, Melanization, Control) %>%
  rename_with(tolower)

write_csv(odonates, "datasets_derived/odonates_merged.csv")

# We'll also output some site level summary info to make maps with.
sites <- odonates %>%
  mutate(year = year(date)) %>%
  select(c(latitude, longitude, site, mite_abundance, year)) %>%
  mutate(mite_abundance = if_else(is.na(mite_abundance), 0, mite_abundance)) %>%
  drop_na(site, latitude, longitude)

site_parasitism <- sites %>%
  select(site, mite_abundance) %>%
  group_by(site) %>%
  summarise(mite_total = sum(mite_abundance))

sites <- distinct(sites, site, .keep_all = TRUE)

sites <- left_join(sites, site_parasitism, by = "site") %>%
  select(!mite_abundance)

write_csv(sites, "datasets_derived/sites.csv")
