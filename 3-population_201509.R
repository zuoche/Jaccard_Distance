setwd("H:/plink-1.07-dos/")

library(permute)
library(lattice)
library(vegan

ab <- c("CCP-1 cit","Vim60-75 cit","Vim2-17 cit", "Fib36-52", "Fib573", "Fib 591", "CEP-1", "ptm12", "ptm13", "ptm15", "ptm35", "ptm36", "ptm37", "Fib alpha621-635 cit", "Fib alpha36-50 cit", "Fib beta60-74 cit")
ab15 <- c("CCP-1 cit","Vim60-75 cit","Vim2-17 cit", "Fib36-52", "Fib573", "Fib 591", "CEP-1", "ptm12", "ptm13", "ptm35", "ptm36", "ptm37", "Fib alpha621-635 cit", "Fib alpha36-50 cit", "Fib beta60-74 cit")


pheno <- read.table("Merged_pheno1-35_caseonly.txt")
eira_wtccc <- read.table("Eira_wtccc.fam")
narac <- read.table("US.fam")
discovery <- subset(pheno, pheno[,2] %in% eira_wtccc[,2])
replication <- subset(pheno, pheno[,2] %in% narac[,2])
discovery16 <- discovery[,3:18]
replication16 <- replication[,3:18]
row.names(discovery16) <- ab
row.names(replication16) <- ab
#jaccard distance for eira+wtccc
tmp <- replace(discovery16, discovery16==1, 0)
discovery16 <- replace(tmp, tmp==2, 1)
discovery16 <- t(discovery16)
discovery_dist <- vegdist(discovery16, method="jaccard")
#jaccard distance for narac
tmp <- replace(replication16, replication16==1, 0)
replication16 <- replace(tmp, tmp==2, 1)
replication16 <- t(replication16)
replication_dist <- vegdist(replication16, method="jaccard")
replication15 <- replication16[-10,]
replication_dist_new <- vegdist(replication15, method="jaccard")

#discovery_dist & replication_dist are 2 lists of 120 values
#(15+14+13+12+11+10+1+9+8+7+6+5+4+3+2+1)
#write to .xlsx and filling the transposed cells


library(pheatmap)
pheatmap(discovery16, scale="none", clustering_distance_col=discovery_dist, cluster_rows=F, cluster_cols=T, fontsize=15, fontsize_rows=3)
pheatmap(replication15, scale="none", clustering_distance_col=replication_dist, cluster_rows=F, cluster_cols=T, fontsize=15, fontsize_rows=3)
pheatmap(replication15, scale="none", clustering_distance_col=replication_dist_new, cluster_rows=F, cluster_cols=T, fontsize=15, fontsize_rows=3)
