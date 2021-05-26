# Author: Jacob Harvey jakekharvey@gmail.com

library(dada2)
library(tidyverse)
library(phylotools)

assign_taxonomy <- function(file_name){
    
    seqs <- read_csv(paste("datasets_derived/sequencing/", file_name, sep = ""))
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
    # Put the common mites at the top, less likely to get cleaned out and cause gaps
    # in the mite numbering.
    mites <- mites[order(rowSums(mites[3:(ncol(mites)-6)]), decreasing = TRUE),]
    mites <- rename(mites, mite = X1)
    mites <- mutate(mites, c = "M")
    mites <- mutate(mites, n = 1:nrow(mites))
    mites <- unite(mites, mite, c(c, n), sep = "")
    mites <- relocate(mites, mite)
    
    outfile <- paste("datasets_derived/sequencing/", 
                     str_extract(file_name, "[:graph:]+(?=\\.)"), "_annotated.csv", 
                     sep="")
    write_csv(mites, outfile)

    
    outfile <- paste("datasets_derived/sequencing/", 
                     str_extract(file_name, "[:graph:]+(?=\\.)"), 
                     "_mites.fasta", 
                     sep="")
    fasta <- select(mites, c(mite, sequence))
    fasta <- rename(fasta, seq.name = mite, seq.text = sequence)
    dat2fasta(fasta, outfile)
    
    message("Done")
}

assign_taxonomy("mite_sequences.csv")
assign_taxonomy("mite_sequences_otu97.csv")
assign_taxonomy("mite_sequences_otu90.csv")
