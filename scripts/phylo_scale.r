# Author: Jake Harvey jakekharvey@gmail.com
# Use the mrf methods to get phylogenetic scale of specialization for each mite. 

library(devtools)
devtools::install_github("jkharv/MRFtools")
devtools::install("~/Code/MRFtools")

library(MRFtools)
library(ape)
library(tidyverse)
library(mgcv)

bin_network <- read_csv("datasets_derived/bin_network_copy.csv")

phylo_path <- "datasets_primary/phylogeny/FinalOdonatetree.tre"
phylo <- read.nexus(phylo_path)

# Phylo trimmed to only the species I've got in the dataset.
trimmed_phylo <- keep.tip(phylo, bin_network$odonate_spp)
plot(trimmed_phylo)

pmatrix <- mrf_penalty(trimmed_phylo)
# interaction_presence ~ s(odandate_sp, bs= phylo_pen) )
bin_network$odonate_spp <- factor(bin_network$odonate_spp)
x <- gam(V1 ~ s(odonate_spp, bs = 'mrf', xt = list(penalty = pmatrix), k = 3), 
         data = bin_network, 
         method = "REML",
         drop.unused.levels = FALSE,
         family = binomial)
summary(x)
rownames(pmatrix)
levels(bin_network$odonate_spp)
str(pmatrix)
