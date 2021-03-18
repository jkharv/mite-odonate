# Author: Jake Harvey jakekharvey@gmail.com

library(tidyverse)
library(ggtree)
library(ape)

mites <- read_csv("datasets_derived/mite_summaries.csv")
phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")

# I have data for some odonates which aren't in the phylogeny.
# These are left out of all the analysis.
keep <- intersect(phylo$tip.label, mites$species)
mites <- filter(mites, species %in% keep)
phylo <- keep.tip(phylo, keep)

tree <- ggtree(phylo) +
  geom_tiplab() + 
  xlim_tree(450)

scale_plot <- facet_plot(tree, panel = "Phylogenetic Scale", data = mites,
                         geom = geom_point, aes(x = mite_scale)) +
                         xlim_expand(c(25, 100), panel = "Phylogenetic Scale") + 
                         theme_tree2()
print(scale_plot)
