library(ggplot2)

maggieRaw <- read.xlsx("Documents/Lessard Lab/Mites-Odonates/Datasets/Maggie's Data.xlsx", sheet = 2, na.strings = "N/A")
immuneData <- subset(maggieRaw, Immune.Response > 0)
immuneData <-immuneData[c(8:12,47)]

# Spss average immune responses
sppImm <- aggregate(immuneData[6], by = list(immuneData$G_species), FUN = mean)
ggplot(sppImm)
