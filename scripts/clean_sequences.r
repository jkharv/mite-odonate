# Author: Tonia de Bellis / JULY 6 2020
#
# Cleans the raw sequences and produces and ASV and an OTU table.
#
library(tidyverse)
library(reshape2)
library(dada2)

library(phyloseq)
library(picante)

library(DECIPHER)
library(Biostrings)
library(speedyseq)
ncores = 3

# Reference
# https://benjjneb.github.io/dada2/tutorial.html
# https://alexiscarter.github.io/metab/Dada_script.html#assignation-taxonomique-taxonomy-assignment
# https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#construct_phylogenetic_tree

#load("dada2_ARG_AMF_clean_MAY25.20.RData")
 
# set file paths
path <- "datasets_primary/sequencing/raw_fastq"   # set to the directory containing the fastq files 

list.files(path) # You should see the names of the fastq files.
#"sample1-miteCOI_S343_L001_R1_001.fastq.gz" 
#plate4-H08-AMF_S321_L001_R1_001.fastq.gz"   

# The sort function ensures forward/reverse reads are in the same order.
fnFs <- sort(list.files(path, pattern="L001_R2_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="L001_R1_001.fastq.gz", full.names = TRUE))
head (fnFs)

basename(fnFs)
# "sample1-miteCOI_S343_L001_R2_001.fastq.gz"
#need to become "sample1-miteCOI"

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# use 1 here because want name only in part 1
# sample1-miteCOI_(part 1) S343_(pt 2) L001_(pt 3) R2(pt 4) _001.fastq.gz (part 5)
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) %>%
  str_subset("(20[:digit:]{2}-[:digit:]{1,2})|(Negative[:graph:]*)|(PCR[:graph:]*)") %>%
  str_extract("(20[:digit:]{2}-[:digit:]{1,2}|(Negative-[:digit:]{1})|(PCR-neg))") %>%
  str_replace_all("-", "_")
sample.names # Correct format for other scripts to recognize

if(interactive()){
  # Plot quality profile of forward reads
  plotQualityProfile(fnFs[1:9]) # nice
  #nice to at least ~ 250 QS30
  # Plot quality profile of reverse reads
  plotQualityProfile(fnRs[1:9]) # very nice ~280
}

# must set filtered file folder path
filt_path.MITES <- file.path(path, "dada2-filtered_MITES") # Place filtered files in filtered/subdirectory
filtFs.MITES <- file.path(filt_path.MITES, paste0(sample.names, "_F_filt.fastq"))
filtRs.MITES <- file.path(filt_path.MITES, paste0(sample.names, "_R_filt.fastq"))

length(filtFs.MITES)#94
length(filtRs.MITES)#94

# need to trim around 20 from start of both reads (gets rid of primer)
# F = 26   / R = 26

out.MITES <- filterAndTrim(fnFs, filtFs.MITES, fnRs, filtRs.MITES, 
                             trimLeft = c(26,26), truncLen=c(250,280),
                             maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                             compress=TRUE, multithread= ncores, verbose=TRUE)  
out.MITES

#Filtering visualization
filvis <- cbind((out.MITES[,2]/out.MITES[,1])*100) # Percentage filtered sequence / non-filtered sequence
filvis_disc <- cbind(out.MITES, filvis) # combines out and pourc
filvis_disc

(mean(out.MITES[,2])/mean(out.MITES[,1]))*100 #94.7 Mean percentage 

# Learn error rates
# multithread argument = number of cores to use
errF <- learnErrors(filtFs.MITES, multithread = ncores)
errR <- learnErrors(filtRs.MITES, multithread = ncores)

if(interactive()){
  # sanity check - visualize error rates
  plotErrors(errF, nominalQ=TRUE)
  plotErrors(errR, nominalQ=TRUE)
}

# Dereplicate the filtered fastq files
# Combines all identical sequencing reads into into unique sequences with a corresponding abundance. 
derepFs <- derepFastq(filtFs.MITES, verbose=TRUE)
derepRs <- derepFastq(filtRs.MITES, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err = errF,  pool = "pseudo", multithread = ncores)
dadaRs <- dada(derepRs, err = errR,  pool = "pseudo", multithread = ncores)

# e.g. inspect results
dadaFs[[1]] #77 ASVs
dadaRs[[1]] #84 ASVs


# Merge the dereplicated forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                      minOverlap = 10, 
                      maxMismatch = 0, returnRejects=F, verbose=TRUE)
# inspect results
head(mergers[[1]])
max(mergers[[1]]$nmatch) #  Largest overlap 224
min(mergers[[1]]$nmatch) #  Smallest overlap 164

# make ASVs table (ASVs are stored in the seqtab object.)
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # 5 568
table(nchar(getSequences(seqtab))) 
seqtab[,1]

rownames(seqtab)

if(interactive()){
  hist(nchar(getSequences(seqtab)),xlab="Size", ylab="Frequency", main = "ASVs length", xlim=c(240,375), ylim=c(0,800)) 
}
# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=40, verbose=TRUE)
# Identified  bimeras out of  input sequences.
dim(seqtab.nochim) # 5 177 
sum(seqtab.nochim)/sum(seqtab) # 0.88

# Length of the non-chimeric sequences
hist(nchar(getSequences(seqtab.nochim)),xlab="Size", ylab="Frequency", main = "Non-chimeric ASVs length", xlim=c(240,375), ylim=c(0,300)) 

#transform the ASVs occurrences in presence / absence which will allow to quantify the number of ASVs per sample.
seqtab.nochim.bin <- ifelse(seqtab.nochim>0,1,0) 

#following table shows how many sequences were eliminated at each step.
getN <- function(x) sum(getUniques(x))
track <- data.frame(Input=as.numeric(out.MITES[,1]), # input
                    Filtered=as.numeric(out.MITES[,2]), # filtered
                    "Filt//In"=as.numeric(round(((out.MITES[,2]/out.MITES[,1])*100),2)),# % (Filtered / Input)
                    Merge = as.numeric(sapply(mergers, getN)), # Merged 
                    "Mer//In"=as.numeric(round(((sapply(mergers, getN)/out.MITES[,1])*100),2)),# % (Merged / Input)
                    Nonchim = as.numeric(rowSums(seqtab.nochim)),# Non-chimeric                       
                    "Nonchim//In"=as.numeric(round(((rowSums(seqtab.nochim)/out.MITES[,1])*100),2)),# % (Non-chimeric / Input)
                    ASV = as.numeric(rowSums(seqtab.nochim.bin))) # Number of ASVs per sample 
rownames(track) <- sample.names # Row names
head(track) 
print(track)

# graph it
gtrack<- track[,c(1,2,4,6)]
gtrack$ID <- rownames(gtrack)

lgtrack <- melt(gtrack, id.vars="ID")
bar_track <- ggplot(lgtrack ,aes(x=ID, y=as.numeric(value), fill=variable)) +
  geom_bar(stat="identity", position = "identity") + 
  theme_classic() + # Theme
  theme(axis.ticks.length=unit(0.3,"cm")) + # Ticks size
  theme(axis.text.x = element_text(angle=45) , legend.title = element_blank())+ # Changes the x labels orientation & delete legend title
  scale_x_discrete(name ="Sample ID", limits=rownames(track))+ # Changes x-axis title & sorts the x label names
  scale_y_continuous(name="Abundance", breaks=seq(from = 0, to = 1000, by = 100))+ #Changes y-axis title & sets the y breaks.
  ggtitle("Track")# Main title
bar_track

#try Steve track code again
getN <- function(x) sum(getUniques(x))
track_S2 <- cbind(out.MITES, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track_S2) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track_S2) <- sample.names
head(track_S2)

# assign taxonomoy --> must blast files?
#save.image("dada2_MITES_JUL6.20.RData")

seqtab_asv_out <- as.data.frame(t(seqtab.nochim))
seqtab_asv_out <- rownames_to_column(seqtab_asv_out, var = "ASV.sequence")
  
write.csv(seqtab_asv_out, file = "datasets_derived/sequencing/mite_sequences.csv")

### Maybe try to Clump ASVs to OTUs -------------------------
# https://github.com/benjjneb/dada2/issues/947

# OTUs at 97%
asv_sequences <- colnames(seqtab.nochim)
sample_names <- rownames(seqtab.nochim)
dna <- Biostrings::DNAStringSet(asv_sequences)

nproc <- ncores # Increase to use multiple processors
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
clusters <- DECIPHER::IdClusters(
                                 d, 
                                 method = "complete",
                                 cutoff = 0.03, # use `cutoff = 0.03` for a 97% OTU 
                                 processors = nproc
                                )

## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
clusters <- clusters %>%
  add_column(sequence = asv_sequences)

seqs <- seqtab.nochim %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sequence") %>%
  full_join(clusters, by = "sequence")

otu_tab <- seqs %>%
  group_by(cluster) %>%
  summarize(across(.cols = contains("_"), sum)) %>%
  left_join(select(seqs, c(sequence, cluster)), by = "cluster") %>%
  distinct(cluster, .keep_all = TRUE) %>%
  relocate(c(cluster, sequence)) %>%
  select(2:last_col()) %>%
  rename("sequence" = "ASV.sequence")

write.csv(otu_tab, file = "datasets_derived/sequencing/mite_sequences_otu97.csv")

# OTUs at 90%
asv_sequences <- colnames(seqtab.nochim)
sample_names <- rownames(seqtab.nochim)
dna <- Biostrings::DNAStringSet(asv_sequences)

nproc <- ncores # Increase to use multiple processors
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
clusters <- DECIPHER::IdClusters(
  d, 
  method = "complete",
  cutoff = 0.1, # use `cutoff = 0.03` for a 97% OTU 
  processors = nproc
)

## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
clusters <- clusters %>%
  add_column(sequence = asv_sequences)

seqs <- seqtab.nochim %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sequence") %>%
  full_join(clusters, by = "sequence")

otu_tab <- seqs %>%
  group_by(cluster) %>%
  summarize(across(.cols = contains("_"), sum)) %>%
  left_join(select(seqs, c(sequence, cluster)), by = "cluster") %>%
  distinct(cluster, .keep_all = TRUE) %>%
  relocate(c(cluster, sequence)) %>%
  select(2:last_col()) %>%
  rename("sequence" = "ASV.sequence")

write.csv(otu_tab, file = "datasets_derived/sequencing/mite_sequences_otu90.csv")
