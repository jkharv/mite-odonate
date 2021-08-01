# Author: Jake Harvey - jakekharvey@gmail.com

library(tidyverse)
library(phangorn)
library(phylotools)

mite_sequences <-read_csv("datasets_derived/sequencing/mite_sequences_annotated.csv") 

seqs <- mite_sequences %>%
  select(c(mite, sequence)) %>%
  rename(seq.text = sequence) %>%
  rename(seq.name = mite)

dat2fasta(seqs, "datasets_derived/sequencing/mite_sequences.fasta")

in_file  <- "datasets_derived/sequencing/mite_sequences.fasta"
out_file <- "datasets_derived/sequencing/mite_sequences_aligned.fasta"
# Aligning sequences using MUSCLE
muscle_cmd <- paste("muscle -in ", in_file, " -out ", out_file, sep = "")
system(muscle_cmd, wait = TRUE)

seqs_aligned <- read.phyDat(out_file, format = "FASTA")

# Create a neighbour joining tree as a base to start off from.
mites_dm  <- dist.ml(seqs_aligned)
mites_nj  <- NJ(mites_dm)

tree_init <- pml(mites_nj, data  = seqs_aligned)

ml_tree <- optim.pml(tree_init, model="GTR", optInv=TRUE, optGamma=TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))

ml_tree_bs <- bootstrap.pml(ml_tree, bs = 500, optNni = TRUE, 
                            control = pml.control(trace = 0))
#plotBS(midpoint(ml_tree$tree), ml_tree_bs, p = 50, type="p")

ml_tree <- midpoint(ml_tree$tree)

write.tree(ml_tree, file = "datasets_derived/mite_tree.tre")
