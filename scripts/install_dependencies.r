library(remotes)

remotes::install_github("benjjneb/dada2", ref="v1.16")
install.packages("BiocManager")
install.packages("picante")
BiocManager::install("DECIPHER")

remotes::install_github("jkharv/MRFtools")

BiocManager::install("phyloseq")
BiocManager::install("DESeq2")
BiocManager::install("muscle")
install.packages("phangorn")
remotes::install_github("mikemc/speedyseq")
remotes::install_github("YuLab-SMU/ggtree")
