################
### 1. Setup
################

#BiocManager::install("edgeR")

library(edgeR)

counts <- read.table("....txt")
head(counts)

#Create DGEList object
d0 <- DGEList(counts)

################
### 2. Preprocessing
################

d0 <- calcNormFactors(d0)
d0

#Filter low-expressed genes
keep <- rowSums(cpm(d0)>=5) >=1
d1 <- d0[keep, , keep.lib.sizes=FALSE]

counts.per.milion.limma<-cpm(d1$counts)
write.table(counts.per.milion.limma, file = "....txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#Derive experiment information from the sample names
snames <- colnames(counts) # Sample names
snames

variant <- substr(snames, 1, nchar(snames) - 3) 
line <- substr(snames, 4, nchar(snames) - 0)
variant
line

group <- interaction(variant)

plotMDS(d1, col = as.numeric(group))

################
### 3. Voom transformation and calculation of variance weights
################

#Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)
mm <- model.matrix(~0 + group)
#The above specifies a model where each coefficient corresponds to a group mean

#voom
y <- voom(d1, mm, plot = T)

################
### 4. Fitting linear models in limma
################

#lmFit fits a linear model using weighted least squares for each gene:
fit <- lmFit(y, mm)
head(coef(fit))

#Specify which groups to compare:
contr <- makeContrasts(group___ - group____, levels = colnames(coef(fit)))
contr

#Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "logFC", n = Inf, confint=FALSE, adjust.method="BH")
head(top.table, 30) 
top.table$Gene <- rownames(top.table)
