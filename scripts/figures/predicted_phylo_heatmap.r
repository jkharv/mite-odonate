# Author: Jake Harvey -- jakekharvey@gmail.com

library(tidyverse)
library(ggtree)
library(ape)
library(phytools)
library(phylogram)

net <- read_csv("datasets_derived/predicted_network.csv")

odonate_phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")
keep <- intersect(odonate_phylo$tip.label, net$species)
odonate_phylo <- keep.tip(odonate_phylo, keep)

net <- net %>%  
  column_to_rownames("species") %>%
  as.matrix()

mite_phylo <- read.tree("datasets_derived/mite_tree.tre")
keep <- intersect(mite_phylo$tip.label, colnames(net))
mite_phylo <- keep.tip(mite_phylo, keep)

# Sort the rows/cols to be in the same order as the phylogeny
net <- net[,mite_phylo$tip.label]
net <- net[odonate_phylo$tip.label,]

svg("figures/predicted_phylo_heatmap.svg", width = 15, height = 10)
heatmap(net, Rowv = as.dendrogram(odonate_phylo), 
        Colv = as.dendrogram(mite_phylo))
dev.off()


