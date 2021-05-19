# Author: Jake Harvey jakekharvey@gmail.com

library(tidyverse)
library(ggtree)
library(ape)
  
scale_plot <- function(mites_path){
  mites <- read_csv(mites_path)
  phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")
  
  # I have data for some odonates which aren't in the phylogeny.
  # These are left out of all the analysis.
  keep <- intersect(phylo$tip.label, mites$species)
  mites <- filter(mites, species %in% keep)
  phylo <- keep.tip(phylo, keep)
  
  tree <- ggtree(phylo, ladderize = TRUE) +
    geom_tiplab() + 
    xlim_tree(500) +
    theme_tree2()
  
  scale_plot <- facet_plot(tree, panel = "Phylogenetic Scale", data = mites,
                           geom = geom_jitter, aes(x = simpson), width = 0) +
                           xlim_expand(c(25, 100), panel = "Phylogenetic Scale")
  
  if(grepl("asv", mites_path)){
    ggsave("figures/pscale_phylo_plot_asv.svg", scale_plot, width = 8, height = 8)
    return()
  }
  if(grepl("otu90", mites_path)){
    ggsave("figures/pscale_phylo_plot_otu90.svg", scale_plot, width = 8, height = 8)
    return()
  }
  if(grepl("otu97", mites_path)){
    ggsave("figures/pscale_phylo_plot_otu97.svg", scale_plot, width = 8, height = 8)
    return()
  }
}

scale_plot("datasets_derived/mite_summaries_asv.csv")
scale_plot("datasets_derived/mite_summaries_otu90.csv")
scale_plot("datasets_derived/mite_summaries_otu97.csv")