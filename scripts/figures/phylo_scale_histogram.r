# Author: Jake Harvey - jakekharvey@gmail.com

library(tidyverse)

mites <- read_csv("datasets_derived/mite_phylo_scale.csv")

#Filter out the zeros, there are a bunch of them.
mites <- filter(mites, combo_phylo_scale > 0.01) 

hist <- ggplot(mites, aes(combo_phylo_scale)) +
  geom_histogram()
print(hist)
