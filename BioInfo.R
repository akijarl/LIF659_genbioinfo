setwd("C:/Users/aki/Desktop/LIF659_genbioinfo/")
list.files()

library(vcfR)

x <- read.vcfR("PopGen_project/sample.vcf.gz")

x@meta
x@fix
x@gt

x <- read.vcfR("PopGen_project/Algea_sample.vcf")

queryMETA(x)

x@meta
x@fix
x@gt

dp <- extract.info(x, element = "DP", as.numeric=TRUE)

summary(dp)
min(dp,na.rm = T)
max(dp,na.rm = T)

hist(dp, main="Depth distribution", xlab="DP")
abline(v=mean(dp,na.rm = T)+(2*sd(dp,na.rm = T)),col="red")
abline(v=mean(dp,na.rm = T)-(2*sd(dp,na.rm = T)),col="red")

dp <- extract.gt(x, element = "DP", as.numeric=TRUE)

summary(dp)
str(dp)

queryMETA(x, element = "DP")

hist(dp, main="Depth distribution", xlab="DP")
abline(v=mean(dp,na.rm = T)+(2*sd(dp,na.rm = T)),col="red")
abline(v=mean(dp,na.rm = T)-(2*sd(dp,na.rm = T)),col="red")

library(ggplot2)

ggplot(dp)+
  geom_bar(aes(dp))

colMeans(dp,na.rm = T)

barplot(colMeans(dp),ylab = "Depth")

geno <- extract.gt(x) # Character matrix Containing the genotypes
position <- getPOS(x) # Positions in bp
chromosome <- getCHROM(x) # Chromosome information
pos_loc<-paste(chromosome,position,sep="_")

require(hierfstat)
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno), dimnames = list(pos_loc,colnames(geno)))
G[geno %in% c("0/0")] <- 11
G[geno %in% c("0/1", "1/0")] <- 12
G[geno %in% c("1/1")] <- 22
G[geno %in% NA] <- 9

pop<-c(rep("AKR",7),rep("GAR",8),rep("GRO",8),rep("HEL",8),rep("POR",8))

G2<-data.frame(t(G))
colnames(G2)<-pos_loc
G2$Pop<-pop
G2<-G2[,c(ncol(G2),1:(ncol(G2)-1))]

AlStat<-basic.stats(G2)
str(AlStat)

perLoc<-AlStat$perloc        
