# Author: Jake Harvey jakekharvey@gmail.com

library(tidyverse)

indices <- read_csv("datasets_derived/mite_phylo_scale.csv")

p <- ggplot(data = indices, aes(x = num_host, y = combo_phylo_scale)) +
    geom_point()
print(p)
