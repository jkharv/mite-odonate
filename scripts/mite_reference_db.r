library(bold)
library(ape)
library(taxize)


odonate_ref <- read.FASTA("odonata_reference_seqs.fas")
mite_ref    <- read.FASTA("arachnida_reference_seqs.fas")
refseq <- c(odonate_ref, mite_ref)

sppnames <- lapply(names(refseq), FUN = function(x) unlist(strsplit(as.character(x), "\\|"))[2])

taxonomies <- classification(get_boldid(as.character(sppnames), division = "Animalia"), db = "bold")









get_taxonomy <- function(seq){ # 
  
  sppname <- unlist(strsplit(names(seq), "\\|"))[2]
  spec <- classification(get_boldid(sppname))
  taxonomy <- cbind(spec$phylum_name, spec$class_name, spec$order_name, spec$family_name, spec$genus_name, spec$species_name)
  taxonomy <- apply(taxonomy, 1, paste, collapse=";")
}
