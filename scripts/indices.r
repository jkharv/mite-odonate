#Author: Jake Harvey (jakekharvey@gmail.com)

library(picante)

source("scripts/phylo_scale.r")

bin_network <- read_csv("datasets_derived/bin_network.csv")

# Filtering out species that we don't have in our phylogeny.
bin_network <- bin_network %>%
  filter(odonate_spp != "Enallagma_annexum") %>%
  filter(odonate_spp != "Lestes_forcipatus") %>%
  filter(odonate_spp != "Lestes_inaequalis") %>%
  filter(odonate_spp != "Lestes_vigilax")

phylo_path <- "datasets_primary/phylogeny/FinalOdonatetree.tre"
phylo <- read.nexus(phylo_path)

#~~~~~~~~~~~~~~~~~~~Phylo scale~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Partial application to fix all the parameters but mite.
# Also vectorize cause mutate needs this.
partial_scale <- Vectorize(function(m){combinational_phylo_scale(bin_network, phylo, m)})

mites <- tibble(mite = colnames(bin_network)[-1]) %>%
  mutate(combo_phylo_scale = partial_scale(mite))

write_csv(mites, "datasets_derived/mite_phylo_scale.csv")

#~~~~~~~~~~~~~~~~~~Specialization indices~~~~~~~~~~~~~~~~~


x<-tibble(sample = rep("A", 31), mite = bin_network$V1, odonate = bin_network$odonate_spp)

x <- bin_network %>%
  pivot_longer(2:last_col()) %>%
  rename(mite = name) %>%
  arrange(mite) %>%
  relocate(mite, value, odonate_spp)

# Ok this is literally *exact* code from picante sample2matrix which would not run
# and gave some inscrutable error. Somehow it runs perfectly not in a function. Fuck R.
colnames(x) <- c("plot","abund","id")
y <- tapply(x$abund, list(x$plot, x$id), sum)
y[is.na(y)] <- 0
x<-as.data.frame(y)

div <- raoD(x, phylo)

# getMRCA(phylo, c("Leucorrhinia_hudsonica", "Leucorrhinia_hudsonica"))
# node_depths <- node.depth.edgelength(phylo)
# node_depths <- max(node_depths) - node_depths 
# node_depths[370]
# c<-mrca(phylo)
# pair <-  c("Leucorrhinia_hudsonica", "Leucorrhinia_hudsonica")
# c[pair[1], pair[2]]