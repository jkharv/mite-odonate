
library(tidyverse)
library(GGally)
library(ggtree)
library(network)
library(bipartite)
library(ape)
library(phytools)

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

odonate_tree <- ggtree(odonate_phylo, ladderize = TRUE) +
  geom_tiplab() + 
  xlim_tree(500)
ggsave("figures/odonate_phylogeny.svg", odonate_tree, width = 25, height = 25, units ="cm")

mite_tree <- ggtree(mite_phylo, ladderize = TRUE) +
  geom_tiplab()
ggsave("figures/mite_phylogeny.svg", mite_tree, width = 25, height = 25, units = "cm")

# Get the labels from the phylo in the order which they are graphed. Don't completely
# know how this works, just copied from an answer from Guangchuang Yu
d <- fortify(odonate_phylo)
d = subset(d, isTip)
odonate_labs <- with(d, label[order(y, decreasing=T)])

d <- fortify(mite_phylo)
d = subset(d, isTip)
mite_labs <- with(d, label[order(y, decreasing=T)])

# Sort the network so it graphs lined up with the phylogenies.
net <- net[,mite_labs]
net <- net[odonate_labs,]

plotweb(net, method = "normal", text.rot = 90, col.interaction = "black")
