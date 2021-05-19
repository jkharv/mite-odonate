
library(tidyverse)
library(GGally)
library(network)
library(bipartite)

net <- read_csv("datasets_derived/mite_sequences_otu90_annotated_bin_network.csv")

net <- net %>%
  column_to_rownames("odonate_spp") %>%
  as.matrix() 

net <- net[,order(colSums(net))]
net <- net[order(rowSums(net)),]

#  as.network(matrix.type = "bipartite")

plotweb(net, method = "normal")
