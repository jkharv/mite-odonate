# Author: Jacob Harvey - jakekharvey@gmail.com
# Take the network in edge-list form add join some information about the mites
# degree of specialization to it (ses and mpd)

library(dplyr)
library(readr)
library(tidyr)

# This script is used on both the observed and predicted data, so I feed in the
# data from a command line argument. This is all orchestrated in the makefile.
cmd_args <- commandArgs(trailingOnly = TRUE)
network_path <- cmd_args[1]
mites_path <- cmd_args[2]

# For debugging.
network_path <- "datasets_derived/predicted_network.csv"
mites_path <- "datasets_derived/mite_pd_predicted.csv"

network <- read_csv(network_path)
mites <- read_csv(mites_path)

# Turn the adjaceny matrix into an edge list, and join the mite data to it.
mite_interactions <- network %>% 
  pivot_longer(2:last_col()) %>%
  rename(mite = name, interaction = value) %>%
  left_join(mites, by = "mite")

# Output to the correct file, whether this is predicted or observed data.
if(grepl("predicted", mites_path)){
  write_csv(mite_interactions, "datasets_derived/mite_summaries_predicted.csv")
} else {
  write_csv(mite_interactions, "datasets_derived/mite_summaries.csv")
}
