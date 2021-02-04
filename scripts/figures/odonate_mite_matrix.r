# Author: Jake Harvey - jakekharvey@gmail.com

library(tidyverse)
library(ggtree)
library(ape)

bin_network <- read_csv("datasets_derived/bin_network.csv")

# Filtering out species that we don't have in our phylogeny.
bin_network <- bin_network %>%
  filter(odonate_spp != "Enallagma_annexum") %>%
  filter(odonate_spp != "Lestes_forcipatus") %>%
  filter(odonate_spp != "Lestes_inaequalis") %>%
  filter(odonate_spp != "Lestes_vigilax")

phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")
trimmed_phylo <- keep.tip(phylo, bin_network$odonate_spp)

tree <- ggtree(trimmed_phylo) +
  geom_tiplab() + 
  xlim_tree(450)
print(tree)

x <- column_to_rownames(as.data.frame(bin_network), "odonate_spp")
x <- as.matrix(x)
x <- x[,order(colSums(x))]

gheatmap(tree, x, offset = 175)
