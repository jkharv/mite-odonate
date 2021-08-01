# Author: Jake Harvey -- jakekharvey@gmail.com

library(tidyverse)
library(ggtree)
library(ape)
library(phytools)
library(phylogram)

net <- read_csv("datasets_derived/mite_sequences_annotated_bin_network.csv")

odonate_phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")
keep <- intersect(odonate_phylo$tip.label, net$odonate_spp)
odonate_phylo <- keep.tip(odonate_phylo, keep)

net <- net %>%
  # A couple of the species we have are not in the phylogeny we have.
  filter(odonate_spp %in% keep)
write_csv(net, "datasets_derived/mite_net_for_plotting.csv")

net <- net %>%  
  column_to_rownames("odonate_spp") %>%
  as.matrix()

mite_phylo <- read.tree("datasets_derived/mite_tree.tre")
keep <- intersect(mite_phylo$tip.label, colnames(net))
mite_phylo <- keep.tip(mite_phylo, keep)

# Sort the rows/cols to be in the same order as the phylogeny
net <- net[,mite_phylo$tip.label]
net <- net[odonate_phylo$tip.label,]

svg("figures/phylo_heatmap.svg", width = 15, height = 10)
heatmap(net, Rowv = as.dendrogram(odonate_phylo), 
        Colv = as.dendrogram(mite_phylo),
        scale = "none")
dev.off()
