# Author: Jake Harvey -- jakekharvey@gmail.com

library(tidyverse)
library(ape)
library(ggtree)

mites <- read_csv("datasets_derived/mite_phylo_scale_asv.csv")
phylo <- read.tree("datasets_derived/mite_tree.tre")

# Trim the mite phylogeny to only the sequences which made it 
# through data cleaning steps.
keep <- intersect(phylo$tip.label, mites$mite)
mites <- filter(mites, mite %in% keep)
phylo <- keep.tip(phylo, keep)

tree <- ggtree(phylo, ladderize = TRUE) +
  geom_tiplab() #+
  #xlim_tree(0.5)
plot(tree)

scale_plot <- facet_plot(tree, panel = "Phylogenetic Scale", data = mites,
                         geom = geom_jitter, aes(x = simpson), width = 0) +
  xlim_expand(c(25, 100), panel = "Phylogenetic Scale")
ggsave("figures/mite_scale_phylo_plot.svg", scale_plot,
       height = 30, width = 30, units = "cm")
