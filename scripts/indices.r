#Author: Jake Harvey (jakekharvey@gmail.com)

library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(picante)
library(ape)

# This script is used on both the observed and predicted data, so I feed in the
# data from a command line argument. This is all orchestrated in the makefile.
cmd_args <- commandArgs(trailingOnly = TRUE)
file_path <- cmd_args

# For debugging.
file_path <- "datasets_derived/predicted_network.csv"

network <- read_csv(file_path)
network <- rename(network, species = odonate_spp)

phylo_path <- "datasets_primary/phylogeny/FinalOdonatetree.tre"
phylo <- read.nexus(phylo_path)

# Filtering out anything in our network that we don't have in our phylogeny.
# We caught like four species that weren't in the Arrowsmith et al. phylogeny.
keep <- intersect(phylo$tip.label, network$species)
network <- filter(network, species %in% keep)

# Reformat our network into a Picante sample matrix format. It's a bit strange
# to convert to a Phylocom sample format only to use sample2matrix(), but for
# some reason Picante didn't like my sample matrix that I created directly.
x <- network %>%
  pivot_longer(2:last_col()) %>%
  rename(mite = name) %>%
  arrange(mite) %>%
  relocate(mite, value, species) %>%
  sample2matrix()

# Calculate the mean pairwise phylogenetic distance and the standardized effect
# size of the MPD against a null model. The null model we're using randomly
# selects the same number of hosts from our sample pool 999 times and does MPD.
ses_mpd <- ses.mpd(x, cophenetic.phylo(phylo), null.model = "sample.pool",
                   abundance.weighted = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(mite = rowname, ses_mpd = mpd.obs.z, mpd = mpd.obs) %>%
  select(mite, ses_mpd, mpd)

# Output to the correct file depending on whether this is the observed or
# predicted data.
if(grepl("predicted", file_path)){
  write_csv(ses_mpd, "datasets_derived/mite_pd_predicted.csv")
} else {
  write_csv(ses_mpd, "datasets_derived/mite_pd.csv")
}
