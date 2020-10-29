# Author: Jacob Harvey jakekharvey@gmail.com

library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16", force = TRUE)
library(dada2)

seqs <- read.csv("datasets_primary/sequencing/seqtab.nochim_MITE5.csv")
mite_ref <- "datasets_derived/sequencing/reference_db.fasta"

seqs_tax <- assignTaxonomy(seqs$ASV.sequence, mite_ref)

