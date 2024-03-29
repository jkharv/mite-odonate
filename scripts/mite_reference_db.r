library(taxize)
library(memoise)
library(tidyverse)
library(phylotools)

main <- function() {

odonate_ref <- read.fasta("datasets_primary/sequencing/odonata_reference_seqs.fas")
mite_ref    <- read.fasta("datasets_primary/sequencing/arachnida_reference_seqs.fas")
refseq <- bind_rows(odonate_ref, mite_ref)
rm(odonate_ref); rm(mite_ref) # These are kinda, big, I don't really want to keep them around...

# DADA2 doesn't support any non ACTG characters in the sequence. Remove - used for alignment
refseq$seq.text <- lapply(refseq$seq.text, function(x) str_replace_all(x, "-", ""))

outfile <- "datasets_derived/sequencing/reference_db.fasta"
dat2fasta(refseq, outfile)
message("FASTA files loaded")

sppnames <- lapply(refseq$seq.name, FUN = function(x) unlist(strsplit(as.character(x), "\\|"))[2])
limit = length(sppnames) 
#limit = 5 # One-stop for doing this on a limited set for debugging purposes


taxonomies <- data.frame(phylum = character(0), class = character(0), order = character(0), 
                         family = character(0), genus = character(0), species = character(0))
message("Making taxonomy requests")
i <- 1
while (i <= limit){
  t <- taxonomy_request(sppnames, i)
  t <- c(phylum = t$phylum, class = t$class, order = t$order, family = t$family, genus = t$genus, species = t$species)
  taxonomies <- bind_rows(taxonomies, t)
  i <- i + 1
  if(i %% 1000 == 0){
    m <- paste(as.character(i/length(sppnames[1:limit]) * 100), "% of taxonomies retrieved")
    message(m)
  }
}
message("Taxonomy information Retrieved")

taxonomies <- unite(taxonomies, taxonomy, sep = ";", na.rm = TRUE)
taxonomies <- modify(taxonomies, function(x) paste(x, ";", sep = "")) # DADA2 needs each taxonomy line to end with ;
# Replace spaces with _ in species names, probably useful later on... (can't do this earlier, doesn't match with BOLD)
taxonomies <- modify(taxonomies, function(x) str_replace_all(x, " ", "_"))

reftable <- data.frame(refseq$seq.name, taxonomies)
rename.fasta(infile = outfile, ref_table = reftable, outfile = outfile)
message("Done")

} # main

cs <- function(x){ # separate out just the part that makes a network request. 
  classification(get_boldid(x, division = "Animalia"), db = "bold")
}
cache <- cache_filesystem("taxon_cache.rcache")
classify <- memoise(cs, cache = cache) # Get a memoized version to save time on duplicate requests (there are many).

taxonomy_request <- function(sppnames, i){ # tries to get taxonomy info from bold in batches
  
  sppname <- sppnames[i]
  t <- "Not NULL"
  attempt <- 0
  repeat{
    attempt <- attempt + 1
    t <- classify(as.character(sppname))
    if(attempt == 5 | !is.na(t)){
      break
    }
  }
  if(attempt > 4){
    return(data.frame(phylum = NA, class = NA, order = NA, family = NA, genus = NA, species = NA))
  }
  taxonomy <- cbind(t)
  return(taxonomy)
}

main()