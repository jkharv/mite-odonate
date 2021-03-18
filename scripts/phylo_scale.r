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

#---------------- Combinations Method ----------------------

combinational_phylo_scale <- function(network, phylo, mite){
  
  parasitized_only <- filter(network , !!sym(mite) > 0)
  
  combinations <- cross2(parasitized_only$odonate_spp, 
                         parasitized_only$odonate_spp)
  node_list <- data.frame(node = numeric(), coef = numeric())
  # For some reason idk why I get different results from getMRCA than fram mrca
  mrcas <- mrca(phylo)
  
  for(pair in combinations){
    
    pair <- unlist(pair)
    mrca <- mrcas[pair[1], pair[2]]
    coef <- parasitized_only[[which(parasitized_only$odonate_spp == pair[1]), 
                              mite]] +
            parasitized_only[[which(parasitized_only$odonate_spp == pair[2]), 
                              mite]]
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
