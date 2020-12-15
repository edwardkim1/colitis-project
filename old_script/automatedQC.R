############################################################
## Quality control and clustering of Martin 2019 scRNA-seq
## edward_kim@college.harvard.edu - December 2019
#############################################################
​
##############################
# 0 - Load librairies
##############################
library(dplyr)
library(ggplot2)
library(future)
library(DropletUtils) # for empty droplets
library(Seurat)
library(DoubletFinder)


##############################
# 1 - Definitions & Settings
##############################
'%ni%' = Negate('%in%')
all.genes <- rownames(s2_postQC) # do after importing counts matrix
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
options(future.globals.maxSize = 11000 * 1024^2)

​############################## 
# 2 - Source file 
##############################
source("clustering_functions.R")
source("htesting_functions.R")

##############################
# 3 - Remove empty droplets
##############################
setwd("../data/martin/")
input.dirs <- dir(pattern="_")
lapply(input.dirs,remove_emptyDroplets) # remove empty droplets

## Diagnostics ##
#table(Limited=e.out$Limited, Significant=is.cell)
#plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
#    xlab="Total UMI count", ylab="-Log Probability")

##############################
# 4 - Apply QC on each sample
##############################
setwd("../data/martin/qc_files")
dirname <- dir("../",pattern="GSM")
mclapply(dirname,qc)

##############################
# 5 - Data integration
##############################
dirname <- dir("../",pattern="GSM")
cells <- mclapply(paste(dirname,"postQC.rds", sep="_"),readRDS)
#### Merging cells excluding patient 5 and 16
s2.wo5 <- merge(cells[[3]], c(cells[[4]], cells[[5]], cells[[6]], cells[[7]], cells[[8]], cells[[9]], cells[[10]], cells[[11]], cells[[12]], cells[[13]], cells[[14]], cells[[15]], cells[[16]], cells[[17]], cells[[18]], cells[[19]], cells[[20]]))
s2 <- s2.wo5[['RNA']]@counts
s2.qc1 <- filter_mito(s2, save=T, filename="integrated_mito_stats.rds")
s2.qc2 <- CreateSeuratObject(s2.qc1,min.cells=10)
s2.qc2$orig.ident <- s2.wo5$orig.ident
s2_postQC <- subset(s2.qc2, subset = nFeature_RNA > 200)
## Clustering steps
s2_postQC <- cluster1(s2_postQC, n.var.features = 2000, parallel=T, n.cores=40)
s2_postQC <- cluster2(s2_postQC,npc=20, k.param=100, resolution=0.4, min.dist=0.5, parallel=T, n.cores=40)
save_figures(s2_postQC)
## Differential gene analysis
plan("multiprocess", workers = 40)
s3_2_regress.markers <- FindAllMarkers(s3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
top10 <- s3_2_regress.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


##############################
# 6 - Subset T cells
##############################
s3 <- subset(s2_postQC,idents= c(1,2,3,5,6,8,10,12,15,17,18))
cd.counts <- Matrix::colSums(s3[['RNA']]@counts[c('CD3E','CD3G','CD3D','CD4','CD8A','CD8B','CD28','ICOS','CD2','CD7'),])
s3_2 <- s3[,cd.counts>0]

s3_2 <- cluster1(s3_2, n.var.features = 2000, parallel=T, n.cores=40)
s3_2 <- cluster2(s3_2, npc=15, k.param=100, resolution=0.5, min.dist=0.5, parallel=T, n.cores=40)
save_figures(s3_2)

##############################
# 7 - Remove ILC's
##############################
s3_3 <- subset(s3_2,idents=c(0,1,2,3,4,5,8,9))
s3_3 <- cluster1(s3_3, n.var.features = 3000, parallel=T, n.cores=40)
s3_3 <- cluster2(s3_3, npc=12, k.param=100, resolution=0.6, min.dist=0.3, parallel=T, n.cores=40)
s3_3 <- label.SingleR(s3_3)
save_figures(s3_3)

df.prop.s3_3 <- wilcoxon.df(s3_3)

pdf("s3v3_barplot_sem_error_9.pdf", width = 10, height = 10)
df.prop.summary <- data_summary(df.prop.s3_3, varname="proportion", groupnames=c("cluster", "status"))
p <- ggplot(df.prop.summary, aes(x=cluster, y=proportion, fill=status)) + 
	geom_bar(stat="identity", position=position_dodge()) +
	geom_errorbar(aes(ymin=proportion-sem, ymax=proportion+sem), width=.2, position=position_dodge(.9)) +
	labs(y = "Fraction of total cells") + scale_fill_brewer(palette="Paired") + theme_minimal()
p
dev.off()

s3_3.pvals <- wilcoxon.p.values(0:9,df.prop.s3_3)

##############################
# 8 - Subset CD8+ from Ts
##############################
sCD8 <- subset(s3_3,subset= CD8A > 1 & CD4 < 1)
sCD8 <- cluster1(sCD8, 2000, parallel=T, n.cores=40)
#sCD8 <- cluster2(sCD8, npc=15, k.param=100, resolution=0.4,min.dist=0.3, parallel=T, n.cores=40) #v1
sCD8 <- cluster2(sCD8, npc=20, k.param=100, resolution=0.5,min.dist=0.3, parallel=T, n.cores=40) #v2
sCD8 <- label.SingleR(sCD8)
save_figures(sCD8)
df.prop.sCD8 <- wilcoxon.df(sCD8)

pdf("sCD8_barplot_sem_error.pdf", width = 10, height = 10)
df.prop.summary <- data_summary(df.prop.sCD8, varname="proportion", groupnames=c("cluster", "status"))
p <- ggplot(df.prop.summary, aes(x=cluster, y=proportion, fill=status)) + 
	geom_bar(stat="identity", position=position_dodge()) +
	geom_errorbar(aes(ymin=proportion-sem, ymax=proportion+sem), width=.2, position=position_dodge(.9)) +
	labs(y = "Fraction of total cells") + scale_fill_brewer(palette="Paired") + theme_minimal()
p
dev.off()

sCD8.pvals <- wilcoxon.p.values(c(0,1,2,3),df.prop.sCD8)

##############################
# 9 - Subset CD4+ from Ts
##############################
sCD4 <- subset(s3_3,subset= CD4 > 1 & CD8A < 1)
sCD4 <- cluster1(sCD4, 2000, parallel=T, n.cores=40)
sCD4 <- cluster2(sCD4, npc=20, k.param=100, resolution=0.5,min.dist=0.3, parallel=T, n.cores=40)
sCD4 <- label.SingleR(sCD4)
save_figures(sCD4)
df.prop.sCD4 <- wilcoxon.df(sCD4)

pdf("sCD4_barplot_sem_error.pdf", width = 10, height = 10)
df.prop.summary <- data_summary(df.prop.sCD4, varname="proportion", groupnames=c("cluster", "status"))
p <- ggplot(df.prop.summary, aes(x=cluster, y=proportion, fill=status)) + 
	geom_bar(stat="identity", position=position_dodge()) +
	geom_errorbar(aes(ymin=proportion-sem, ymax=proportion+sem), width=.2, position=position_dodge(.9)) +
	labs(y = "Fraction of total cells") + scale_fill_brewer(palette="Paired") + theme_minimal()
p
dev.off()

sCD4.pvals <- wilcoxon.p.values(c(0,1,2),df.prop.sCD4)




