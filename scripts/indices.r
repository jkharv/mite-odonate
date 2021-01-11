#Author: Jake Harvey (jakekharvey@gmail.com)
#
#
#

library(picante)

x<-tibble(sample = rep("A", 31), mite = bin_network$V1, odonate = bin_network$odonate_spp)
x<sample2matrix(x)

ph_pd(x, phylo)
typeof(x$odonate)
