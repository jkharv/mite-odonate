# Maggie's data
spplk <- c(1,5,0,1,1,2,4,1,1,0,2,1,2,4,0,3,2,1,2,4,2,0,3,5,5)
mean(spplk)
qqnorm(spplk)
qqline(spplk)
t.test(spplk)


# Julie's data
spplk2 <- read.csv("julies_only_parasitized.csv")
mean(spplk2$spss)
qqnorm(spplk2$spss)
qqline(spplk2$spss)
t.test(spplk2$spss) # 2.82 * 40 = 113 samples