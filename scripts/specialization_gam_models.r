# Author: Jacob Harvey - jakekharvey@gmail.com  

library(MRFtools)
library(ape)
library(tidyverse)
library(mgcv)
library(car)

mites <- read_csv("datasets_derived/mite_summaries.csv")

mites$species <- factor(mites$species)

phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")
# Phylogeny is missing many species. I have mite data for four of the missing.
keep <- intersect(phylo$tip.label, mites$species)
mites <- filter(mites, species %in% keep)
phylo <- keep.tip(phylo, keep)
plot(phylo)

pmatrix <- mrf_penalty(phylo)
qqPlot(mites$mite_scale)
hist(mites$mite_scale)

m <- gam(sqrt(rank(mite_scale)) ~ s(species, bs = 'mrf', xt = list(penalty = pmatrix)) + s(prevalence) + log(mass), 
         data = mites, method = "REML", family = gaussian)
summary(m)


x<-ggplot(aes(x = mite_scale, y = prevalence), data = mites) +
  geom_point() + 
  scale_y_log10() + 
  geom_smooth()
print(x)
