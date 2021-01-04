# Author: Jake Harvey jakekharvey@gmail.com
# Use the mrf methods to get phylogenetic scale of specialization for each mite. 

library(devtools)
devtools::install_github("jkharv/MRFtools")
devtools::install("~/Code/MRFtools")

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


v1scale <- combinational_phylo_scale(bin_network, phylo, "V36")

#---------------- Combinations Method ----------------------

combinational_phylo_scale <- function(network, phylo, mite){

  parasitized_only <- filter(network , !!sym(mite) > 0)
  
  combinations <- cross2(parasitized_only$odonate_spp, parasitized_only$odonate_spp)
  node_list <- data.frame(node = numeric(), coef = numeric())
  
  for(pair in combinations){
    
    pair <- unlist(pair)
    mrca <- getMRCA(phylo, pair)
    coef <- parasitized_only[[which(parasitized_only$odonate_spp == pair[1]), mite]] +
            parasitized_only[[which(parasitized_only$odonate_spp == pair[2]), mite]]
    node_list <- add_to_df(node_list, mrca, coef)
  }
  
  node_depths <- node.depth.edgelength(phylo)
  
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

#get_edf <- function(){ # I'm thinking update.formula is the way to go.
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
