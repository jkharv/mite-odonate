# Author: Jacob Harvey jakekharvey@gmail.com

library(dada2)
library(tidyverse)

file_name <- "mite_sequences.csv"
seqs <- read.csv(paste("datasets_derived/sequencing/", file_name, sep = ""))
seqs <- rename(seqs, sequence = ASV.sequence)
mite_ref <- "datasets_derived/sequencing/reference_db.fasta"

message("Assigning taxonomy to sequences.")
seqs_tax <- assignTaxonomy(seqs$sequence, mite_ref, tryRC = TRUE, multithread = 8)
seqs_tax <- as.data.frame(seqs_tax)
seqs_tax <- rownames_to_column(seqs_tax)
seqs_tax <- rename(seqs_tax, sequence = rowname, phylum = Kingdom, class = Phylum, 
                   order = Class, family = Order, genus = Family, species = Genus)

message("Filtering out non-mite sequences")
s <- left_join(seqs, seqs_tax, by = "sequence")
# only interested in mite sequences.
mites <- filter(s, order == "Trombidiformes")


outfile <- paste("datasets_derived/sequencing/", 
                 str_extract(file_name, "[:graph:]+(?=\\.)"), 
                 "_annotated.csv", 
                 sep="")
write_csv(mites, outfile)
message("Done")






