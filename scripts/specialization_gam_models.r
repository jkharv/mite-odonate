# Author: Jacob Harvey - jakekharvey@gmail.com  

library(MRFtools)
library(ape)
library(tidyverse)
library(mgcv)
library(car)

odonates <- read_csv("datasets_derived/odonate_summaries.csv")

odonates$species <- factor(odonates$species)
odonates <- drop_na(odonates, mite_scale_mean)
odonates <- mutate(odonates, sr = specialist / generalist)

phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")
# Phylogeny is missing many species. I have mite data for four of the missing.
keep <- intersect(phylo$tip.label, odonates$species)
odonates <- filter(odonates, species %in% keep)
phylo <- keep.tip(phylo, keep)

pmatrix <- mrf_penalty(phylo)

spec_ratio <- cbind(odonates$specialist, odonates$generalist)

phylo_model <- gam(spec_ratio ~ s(species, bs = 'mrf', xt = list(penalty = pmatrix), k = 10), 
         data = odonates, method = "REML", family = binomial)
summary(phylo_model)

r_model <- gam(phylo_model$residuals ~ log(mass), data = odonates, method = "REML", family = "gaussian")
summary(r_model)

mass_model <- gam(spec_ratio ~ log(mass), 
         data = odonates, method = "REML", family = binomial)
summary(mass_model)
