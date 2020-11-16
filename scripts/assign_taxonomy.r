# Author: Jacob Harvey jakekharvey@gmail.com

if(!require(dada2)){
  library("devtools")
  devtools::install_github("benjjneb/dada2", ref="v1.16")
}

library(dada2)
library(tidyverse)

seqs <- read.csv("datasets_primary/sequencing/seqtab.nochim_MITE5.csv")
seqs <- rename(seqs, sequence = ASV.sequence)
mite_ref <- "datasets_derived/sequencing/reference_db.fasta"

seqs_tax <- assignTaxonomy(seqs$sequence, mite_ref, tryRC = TRUE)
seqs_tax <- as.data.frame(seqs_tax)
seqs_tax <- rownames_to_column(seqs_tax)
seqs_tax <- rename(seqs_tax, sequence = rowname, phylum = Kingdom, class = Phylum, 
                   order = Class, family = Order, genus = Family, species = Genus)

s <- left_join(seqs, seqs_tax, by = "sequence")
# only interested in mite sequences.
mites <- filter(s, order == "Trombidiformes")

outfile <- "datasets_derived/sequencing/annotated_mite_sequences.csv"
write_csv(mites, outfile)






