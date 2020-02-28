##YOU HAVE TO INSTALL THOSE LIBRARIES, IF YOU HAVE THEM JUST IGNORE THIS STEP
install.packages("BiocManager")
BiocManager::install("SNPRelate")
BiocManager::install("gdsfmt")

##SET THE PWD
setwd("~/Escritorio/ProyectoUni5/bin")

##LOAD THE PACKAGES RELATED WITH SNPRELATE
library("ggplot2")
library("BiocManager")
library("gdsfmt")
library("SNPRelate")
library("ape")

###DOWNLOAD THE FILE 
##CREAT DATA GDS FROM PLILNK
snpgdsBED2GDS("../data/lobosplink.bed",
              "../data/lobosplink.fam",
              "../data/lobosplink.bim", 
              out.gdsfn="../data/lobos.gds", 
              option = snpgdsOption(Z=38)) # 138 cromosomes

## CHECK ABSTRACT
snpgdsSummary("../data/lobos.gds")

##WORK WITH THE .GDDS DOCUMENTS
genofile <- snpgdsOpen("../data/lobos.gds")

##CHEK snp.ids
head(read.gdsn(index.gdsn(genofile,"snp.id")))

###CHECK THE NAMES OF SAMPLES
sample.id <- read.gdsn(index.gdsn(genofile,"sample.id"))
sample.id

##CREAT A PCA
#PCA
pca <- snpgdsPCA(genofile, num.thread = 2)

##CALCULATE THE PERCENT OF VARIATION OF THE FIRST COMPONENTS
pc.percent <- pca$varprop*100
head(round(pc.percent,2))
