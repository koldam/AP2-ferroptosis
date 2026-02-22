library(pathview)

data <- read.delim("....txt", row.names = 1)

pv.out <- pathview(gene.data = data, pathway.id = "____", 
                   species = "hsa", gene.idtype = "symbol", 
                   out.suffix = "out", same.layer=T, limit = list(gene=1.5),
                   bins = list(gene=20), low = "blue", mid = "gray", high = "red",
                   res = 600, key.pos = "topright")