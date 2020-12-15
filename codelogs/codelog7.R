######################################################################
## 
## Analysis of all T cells with PD1
## 
## edward_kim@college.harvard.edu - August 2020
######################################################################

##############################
# 0 - Load librairies
##############################
library(dplyr)
library(ggplot2)
library(Seurat)
library(scran)
library(pheatmap)

##############################
# 1 - Definitions & Settings
##############################
'%ni%' = Negate('%in%')
all.genes <- rownames(s2_postQC) # do after importing counts matrix
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
options(future.globals.maxSize = 11000 * 1024^2)
setwd("~/1.projects/PD1")

############################## 
# 2 - Source file 
##############################
source("scripts/Stat111functions.R")
source("scripts/new_clustering_functions.R")

#################################
## Basic Inference: onlyTcells ##
#################################
# Metadata Key: 
# sampletype = patient ID
# colitis = Colitis | Control | No-Colitis
# colitis2 = Colitis(aPD1) | Colitis (combo) | Control | 
# 			 No-Colitis (aPD1) | No-colitis (combo)
# condition = Batch1 | Batch2 | Batch3 | Batch4 | BatchPD1
# orig.ident = CD3_Tcell
# others: Phase, S.Score, G2M.Score
# will create: colitis2
# 
# PD-1 info:
# colitis+: "p101519_CD3","p020620_CD45"
# no colitis: "DFCI-1294_CD3","DFCI-1545_CD3","DFCI-1618_CD3"
s <- readRDS("saved_objects/onlyTcells.rds")

# Viz 1: UMAP - group by colitis
p.status <- DimPlot(s, reduction = "umap", label=F, group.by="colitis")
ggsave("figures/onlyTcells_group_colitis.pdf")
# Viz 2: UMAP - split by colitis
p.status.split <- DimPlot(s, reduction = "umap", label=F, split.by="colitis", ncol=3)
ggsave("figures/onlyTcells_split_colitis.pdf", width= 12, height= 4, units= "in")
# Viz 3: UMAP - group by colitis2
p.colitis2 <- DimPlot(s, reduction = "umap", label=F, group.by="colitis2")
ggsave("figures/onlyTcells_group_colitis2.pdf")
# Viz 4: cluster proportions
p.prop <-



DimPlot(seurat.object, reduction = "umap", label=F, split.by="status")
ggsave("figures/onlyTcells_split_status.pdf")

#####################################
## Basic Inference: onlyTcells1550 ##
#####################################
# Metadata Key: 
# sampletype = patient ID
# colitis = Colitis | Control | No-Colitis
# colitis2 = Colitis(aPD1) | Colitis (combo) | Control | 
# 			 No-Colitis (aPD1) | No-colitis (combo)
# condition = Batch1 | Batch2 | Batch3 | Batch4 | BatchPD1
# orig.ident = CD3_Tcell
# others: Phase, S.Score, G2M.Score
# will create: colitis2
# 
# PD-1 info:
# colitis+: "p101519_CD3","p020620_CD45"
# no colitis: "DFCI-1294_CD3","DFCI-1545_CD3","DFCI-1618_CD3"
s <- readRDS("saved_objects/onlyTcells1550_labeled.rds")

## Viz 1: UMAP - group by colitis and colitis2
p <- DimPlot(s, reduction = "umap", label=F, group.by="colitis")
q <- DimPlot(s, reduction = "umap", label=F, group.by="colitis2")
p+q
ggsave("figures/onlyTcells1550_group_colitis.pdf",width= 14, height= 7, units= "in")
## Viz 2: UMAP - split by colitis2
p.status.split <- DimPlot(s, reduction = "umap", label=F, split.by="colitis2", ncol=3)
ggsave("figures/onlyTcells1550_split_colitis2.pdf", width= 12, height= 8, units= "in")
## Viz 3: UMAP - split by colitis
p.status.split <- DimPlot(s, reduction = "umap", label=F, split.by="colitis", ncol=3)
ggsave("figures/onlyTcells1550_split_colitis.pdf", width= 12, height= 5, units= "in")

## Viz 4a: cluster proportions
table(s$seurat_clusters,s$colitis2) %>% unclass() %>% prop.table()
data.frame(cluster =s$seurat_clusters, condition= s$colitis) %>%
	ggplot() + geom_bar(
		mapping = aes(x= cluster, fill= condition),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/onlyTcells1550_cluster_prop.pdf")

## Viz 4b: cluster proportions colitis 2
data.frame(cluster =s$seurat_clusters, condition= s$colitis2) %>%
	ggplot() + geom_bar(
		mapping = aes(x= cluster, fill= condition),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/onlyTcells1550_cluster_prop_colitis2.pdf")

## Viz 4c: cluster proportions by condition colitis 2
reordered.conditions <- factor(s$colitis2, levels = c("Colitis (aPD1)" ,"No-Colitis (aPD1)" , "Colitis (combo)" ,  "No-Colitis (combo)", "Control"))
data.frame(cluster =s$seurat_clusters, condition= reordered.conditions) %>%
	ggplot() + geom_bar(
		mapping = aes(x= condition, fill= cluster),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/onlyTcells1550_cluster_prop_colitis2_by_condition.pdf")

## Viz 4d: cluster proportions by condition colitis 2 (combotherapy)
clusters <- subset(s, colitis2== "Colitis (combo)" | colitis2== "No-Colitis (combo)" )$seurat_clusters
conditions <- subset(s, colitis2== "Colitis (combo)" | colitis2== "No-Colitis (combo)" )$colitis2
data.frame(cluster = clusters, condition= conditions) %>%
	ggplot() + geom_bar(
		mapping = aes(x= cluster, fill= condition),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/onlyTcells1550_cluster_prop_combotherapy.pdf")

## Viz 4e: cluster proportions by condition colitis 2 (monotherapy)
clusters <- subset(s, colitis2== "Colitis (aPD1)"| colitis2== "No-Colitis (aPD1)")$seurat_clusters
conditions <- subset(s, colitis2== "Colitis (aPD1)"| colitis2== "No-Colitis (aPD1)")$colitis2
data.frame(cluster = clusters, condition= conditions) %>%
	ggplot() + geom_bar(
		mapping = aes(x= cluster, fill= condition),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/onlyTcells1550_cluster_prop_monotherapy.pdf")

## Viz 4f: cluster proportions  (per patient analysis for monotherapy)
clusters <- subset(s, colitis2== "Colitis (aPD1)"| colitis2== "No-Colitis (aPD1)")$seurat_clusters
patient <- subset(s, colitis2== "Colitis (aPD1)"| colitis2== "No-Colitis (aPD1)")$sampletype
data.frame(cluster = clusters, patient= patient) %>%
	ggplot() + geom_bar(
		mapping = aes(x= patient, fill= cluster),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/onlyTcells1550_cluster_prop_monotherapy_patient.pdf")

## Viz 4g: DA
library(edgeR)
x <- subset(s, colitis2== "Colitis (aPD1)"| colitis2== "No-Colitis (aPD1)")
abundances <- table(x$seurat_clusters, x$sampletype) %>% unclass()
extra.info <- x@meta.data[match(colnames(abundances),x$sampletype),]
y.ab <- DGEList(abundances,samples=extra.info)

keep <- filterByExpr(y.ab, group=y.ab$samples$colitis2) #just for procedural reasons; no clusters are actually removed
y.ab <- y.ab[keep,]

design <- model.matrix(~factor(colitis2),y.ab$samples)
design

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)

pdf("figures/onlyTcells1550_BVC.pdf")
plotBCV(y.ab, cex=1)
dev.off()

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
summary(fit.ab$df.prior)

pdf("figures/onlyTcells1550_QLDisp.pdf")
plotQLDisp(fit.ab, cex=1)
dev.off()

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)


## Viz 5: FeaturePlots of various genes
FeaturePlot(s, features= c('CD3E','CD3G','CD3D','CD28'))
ggsave("figures/onlyTcells1550_Features_Tcell1.pdf", width= 12, height= 12, units= "in")
FeaturePlot(s, features= c('CD4','CD8A','CD8B','ICOS'))
ggsave("figures/onlyTcells1550_Features_Tcell2.pdf", width= 12, height= 12, units= "in")
FeaturePlot(s, features= c('HAVCR2','LAG3','CTLA4','PDCD1'))
ggsave("figures/onlyTcells1550_Features_exhaustion.pdf", width= 12, height= 12, units= "in")
FeaturePlot(s, features= c('KLRB1','NKG7','TOX','TOX2'))
ggsave("figures/onlyTcells1550_Features_Misc1.pdf", width= 12, height= 12, units= "in")
FeaturePlot(s, features= c('IL17A','IL10','IFNG','TGFB1'))
ggsave("figures/onlyTcells1550_Features_cytokines.pdf", width= 12, height= 12, units= "in")




## Pearson's Chi-squared test
# Calculating the control's exepected distribution
control <- subset(s, colitis=="Control")$seurat_clusters %>% table()
no.colitis <- subset(s, colitis=="No-Colitis")$seurat_clusters %>% table()
colitis <- subset(s, colitis=="Colitis")$seurat_clusters %>% table()
chi_squared_test <- function(theoretical.counts, observed.counts) {
	expected.prop <- theoretical.counts/sum(theoretical.counts)
	N.observed <- sum(observed.counts)
	expected.counts <- N.observed*expected.prop
	chi.sq.stat <- sum((observed.counts - expected.counts)^2 / expected.counts)
	deg.of.freedom <- nrow(observed.counts) - 1
	p.value<- pchisq(chi.sq.stat,df=deg.of.freedom, lower.tail=F)
	print(chi.sq.stat)
	print(deg.of.freedom)
	print(p.value)
}
chi_squared_test(control,no.colitis)
chi_squared_test(control,colitis)




