#dependencies
library(magrittr)
library(survival)
library(survminer)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(org.Hs.eg.db)
library(KEGG.db)
library(fgsea)
library(enrichplot)
library(clusterProfiler)
library(magrittr)
library(dummies) 
library(formula.tools)
library(foreach)
library(gdata)
library(gtools)
library(gmodels)
library(MASS)
library(gplots)
library(formula.tools)
library(Hmisc) 
library(Matrix)
library(dplyr)
library(GSVA)
library(doParallel)

#install and load FPSOmics
#devtools::install_github('Yelab2020/FPSOmics')
library(FPSOmics)

#Calculate ferroptosis score
print('FPS')
data(m1_input_mRNA,package='FPSOmics') #load example mRNA matrix
Dataset <- read.table("...txt",  #load your own mRNA matrix
                      header=TRUE, sep="\t", na.strings="NA", dec=".",
                      strip.white=TRUE)
Dataset <- data.frame(Dataset, row.names = 1)
FPS_score=FPSOmics::FPS(Dataset)
head(FPS_score)
