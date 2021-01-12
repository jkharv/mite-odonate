# Author: Jake Harvey jakekharvey@gmail.com
#
# TODO. Filtering out odonates not in my phylogeny leveas some mites with zeroes
# all the way through, creating NaNs as scales.
#
# TODO. GAM can't fit on the using the full phylogeny because there isn't enough
# species smapled to estimate all the leaves and internal nodes. Can only estimate
# as many nodes as I have data points. Must order tips based on distance from sampled
# tips and discard as many as needed to allow the GAM to fit.
#
# TODO. Separate the analysis from the algorithms, put them in their own file.
#
# TODO. More permanent solution to MRFtools compatibility issue with Gratia.
#

library(MRFtools)
library(ape)
library(tidyverse)
library(mgcv)

bin_network <- read_csv("datasets_derived/bin_network.csv")

# Filtering out species that we don't have in our phylogeny.
bin_network <- bin_network %>%
  filter(odonate_spp != "Enallagma_annexum") %>%
  filter(odonate_spp != "Lestes_forcipatus") %>%
  filter(odonate_spp != "Lestes_inaequalis") %>%
  filter(odonate_spp != "Lestes_vigilax")

phylo_path <- "datasets_primary/phylogeny/FinalOdonatetree.tre"
phylo <- read.nexus(phylo_path)

# Partial application to fix all the parameters but mite.
# Also vectorize cause mutate needs this.
partial_scale <- Vectorize(function(m){combinational_phylo_scale(bin_network, phylo, m)})

mites <- tibble(mite = colnames(bin_network)[-1]) %>%
  mutate(combo_phylo_scale = partial_scale(mite))

write_csv(mites, "datasets_derived/mite_phylo_scale.csv")

#---------------- Combinations Method ----------------------

combinational_phylo_scale <- function(network, phylo, mite){
  
  parasitized_only <- filter(network , !!mite > 0)
  
  combinations <- cross2(parasitized_only$odonate_spp, parasitized_only$odonate_spp)
  node_list <- data.frame(node = numeric(), coef = numeric())
  
  for(pair in combinations){
    
    pair <- unlist(pair)
    mrca <- getMRCA(phylo, pair)
    coef <- parasitized_only[[which(parasitized_only$odonate_spp == pair[1]), mite]] +
            parasitized_only[[which(parasitized_only$odonate_spp == pair[2]), mite]]
    node_list <- add_to_df(node_list, mrca, coef)
  }
  
  #node_depths <- node.depth.edgelength(phylo)
  node_depths <- node.depth.edgelength(phylo)
  node_depths <- max(node_depths) - node_depths 
  
  z <- node_list %>%
    mutate(weight = coef / sum(node_list$coef)) %>%
    mutate(depth = node_depths[node]) %>%
    mutate(weighted_depth = weight * depth)
  
  phylo_scale <- sum(z$weighted_depth)
  
  return(phylo_scale)
}

add_to_df <- function(df, node, coef){
  
  if(node %in% df$node){
    x <- df[which(df$node == node), "coef"]
    df[which(df$node == node), "coef"] <- x + coef
  } else {
    df <- bind_rows(df, data.frame("node" = node, "coef" = coef))
  }
  return(df)
}

#----------------- MRF GAM METHOD --------------------------

# Phylo trimmed to only the species I've got in the dataset.
trimmed_phylo <- keep.tip(phylo, bin_network$odonate_spp)
plot(trimmed_phylo)

pmatrix <- mrf_penalty(trimmed_phylo)

#get_edf <- function(){
# interaction_presence ~ s(odonate_sp, bs= phylo_pen) )
bin_network$odonate_spp <- factor(bin_network$odonate_spp, levels = rownames(trimmed_phylo))
levels(bin_network$odonate_spp)
dim(pmatrix)
rownames(pmatrix)
x <- gam(V38 ~ s(odonate_spp, bs = 'mrf', xt = list(penalty = pmatrix), k = 3), 
         data = bin_network, 
         method = "REML",
         drop.unused.levels = FALSE,
         family = binomial)
summary(x)$edf
#}
