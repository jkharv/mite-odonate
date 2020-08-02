# Author: Jake Harvey - jakekharvey@gmail.com

library(taxize)

args <- commandArgs(trailingOnly = TRUE)
input_file  <- args[2] 
output_file <- paste("datasets_derived/annotated_sequences/",args[1], ".csv", sep = "")

blast <- read.csv(input_file, header = FALSE)

# The CSV must be sorted. Should be already if it's straight from BLAST.
seqs <- blast[!duplicated(blast$V1),][c("V1", "V2", "V3")]
colnames(seqs) <- c("Sequence", "Accession_number", "Percent_indentity")

# API call to NCBI, get taxon ids for the accession numbers.
# And do it silently, I don't have a key.
suppressMessages(results <- invisible(genbank2uid(seqs$Accession_number, batch_size = 1)))


# Get the taxonomy from each accession number.
suppressMessages(res <- classification(results, db = "ncbi", return_id = FALSE))
res <- cbind(res)

taxonomy <- subset(res, TRUE, select = c("superkingdom", "kingdom","phylum","class", "order", "family","genus", "species"))

write.csv(cbind(seqs, taxonomy), output_file, row.names = FALSE)