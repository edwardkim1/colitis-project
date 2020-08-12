######################################################################
## Integration, QC, and clustering of 
## all checkpoint colitis T cells
## edward_kim@college.harvard.edu - June 2020
######################################################################
​
##############################
# 0 - Load librairies
##############################
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(Seurat)
library(sctransform)
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
#source("/scripts/new_clustering_functions.R")
load("data/ssuo_CD3/seurat.object.RData")

####################################
# 3 - Preparing data for integration
####################################
# Metadata Key: 
# sampletype = patient ID
# colitis = Colitis | Control | No-Colitis
# condition = Batch1 | Batch2 | Batch3 | Batch4 | BatchPD1
# orig.ident = CD3_Tcell
# others: Phase, S.Score, G2M.Score
# will create: colitis2
# 
# PD-1 info:
# colitis+: "p101519_CD3","p020620_CD45"
# no colitis: "DFCI-1294_CD3","DFCI-1545_CD3","DFCI-1618_CD3"

## Modifying T cell data
Tcell$colitis2 <- ifelse(Tcell$colitis == "Colitis", "Colitis (combo)",Tcell$colitis)
Tcell$colitis2 <- ifelse(Tcell$colitis2 == "No-Colitis", "No-Colitis (combo)",Tcell$colitis2)

# p101519
s <- readRDS("saved_objects/p101519_CD3_3.RDS")
s$sampletype <- "p101519"
s$colitis <- "Colitis"
s$condition <- "BatchPD1"
s$orig.ident <- "CD3_Tcell"
s$colitis2 <- "Colitis (aPD1)"
# p020620
s1 <- readRDS("saved_objects/p020620_CD3subset.RDS")
s1$sampletype <- "p020620"
s1$colitis <- "Colitis"
s1$condition <- "BatchPD1"
s1$orig.ident <- "CD3_Tcell"
s1$colitis2 <- "Colitis (aPD1)"
# DFCI-1294_CD3
s2 <- readRDS("saved_objects/DFCI-1294_CD3_3.RDS")
s2$sampletype <- "DFCI1294"
s2$colitis <- "No-Colitis"
s2$condition <- "BatchPD1"
s2$orig.ident <- "CD3_Tcell"
s2$colitis2 <- "No-Colitis (aPD1)"
# DFCI-1545_CD3
s3 <- readRDS("saved_objects/DFCI-1545_CD3_3.RDS")
s3$sampletype <- "DFCI1545"
s3$colitis <- "No-Colitis"
s3$condition <- "BatchPD1"
s3$orig.ident <- "CD3_Tcell"
s3$colitis2 <- "No-Colitis (aPD1)"
# DFCI-1618_CD3
s4 <- readRDS("saved_objects/DFCI-1618_CD3_3.RDS")
s4$sampletype <- "DFCI1618"
s4$colitis <- "No-Colitis"
s4$condition <- "BatchPD1"
s4$orig.ident <- "CD3_Tcell"
s4$colitis2 <- "No-Colitis (aPD1)"

####################################
# 4 - Merging data
####################################
s <- merge(Tcell, c(s,s1,s2,s3,s4))
# removing unapplicable meta.data columns
m.names <- colnames(s@meta.data)
columns.to.remove <- m.names[grepl("^pANN*|^DF.c*|seurat_clusters|^RNA*|cells.to.remove",m.names)]
for(i in columns.to.remove) {
  s[[i]] <- NULL
}
#33538x81352 to 18188x81352
x <- s[["RNA"]]@counts
nnz_row <- tabulate(x@i + 1)
keep <- rownames(x)[nnz_row>10]
s <- subset(s, features = keep)
saveRDS(s, "saved_objects/allTcells.RDS")

####################################
# 5 - Clustering
####################################
s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")
stats <- filter_stats(s, save=T, filename="saved_objects/allTcells_filterstats.RDS")
s$cells.to.remove <- stats$cells.to.remove
saveRDS(s,"saved_objects/allTcells_1.RDS")
s <- s[,!s$cells.to.remove]
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
nPCs <- 20
s <- ScaleData(s, features = all.genes, vars.to.regress="percent.mt") %>% RunPCA(features=VariableFeatures(s)) %>% RunUMAP(dims = 1:nPCs) %>% FindNeighbors(dims = 1:nPCs) %>% FindClusters(resolution = 0.5)

saveRDS(s,"saved_objects/allTcells_2.RDS")

##############################
# 6 - Diagnostic Plots
##############################
## Figures from QC
stats <- readRDS("saved_objects/allTcells_filterstats.RDS")

# Mito% Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre mito%", length(stats$mt.pre)),rep("post mito%", length(stats$mt.post))), percent.mito = c(stats$mt.pre,stats$mt.post))) + geom_violin(aes(x= x, y=percent.mito))
ggsave("figures/allTcells_MITOviolin.pdf")

# Log(No. of Features) Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre lnf", length(stats$lnf.pre)),rep("post lnf", length(stats$lnf.post))), log.num.feats = c(stats$lnf.pre,stats$lnf.post))) + geom_violin(aes(x= x, y=log.num.feats))
ggsave("figures/allTcells_LNFviolin.pdf")

# Log.no.features vs. mito% Figure
p <- readRDS("saved_objects/allTcells_1.RDS") %>% FeatureScatter(feature1="nFeature_RNA", feature2="percent.mt", group.by="cells.to.remove")
ggsave("figures/allTcells_LNFxMITO.pdf")

# Mito% Model Fitting Figure
p <- plot_kde(stats$mt.pre, kernel = "gaussian", bw=0.5, lab.x = "mitochondrial gene percentage (per cell)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$mt.median,stats$mt.mad)}, 
	support=seq(0,100,0.1))
p + geom_vline(xintercept= stats$mt.lim, color = "red") +
	geom_vline(xintercept= stats$mt.median-3*stats$mt.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$mt.median+3*stats$mt.mad, linetype = "dotted")
ggsave("figures/allTcells_MITO-KDE-Normal.pdf")

# Log.no.features Model Fitting Figure
p <- plot_kde(stats$lnf.pre, kernel = "gaussian", bw=0.01, lab.x = "Log(No. of Features)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$lnf.median,stats$lnf.mad)})
p + geom_vline(xintercept= stats$lnf.lim, color = "red") +
	geom_vline(xintercept= stats$lnf.median-3*stats$lnf.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$lnf.median+3*stats$lnf.mad, linetype = "dotted")
ggsave("figures/allTcells_LNF-KDE-Normal.pdf")


## Figures from Object 2
s <- readRDS("saved_objects/allTcells_2.RDS")

# Variable Features Save & Figure
VariableFeatures(s)%>%saveRDS(file="saved_objects/allTcells_hvf.RDS")
p <- VariableFeaturePlot.Tcells(s) %>% LabelPoints(points=head(VariableFeatures(s),10), repel= TRUE)
ggsave("figures/allTcells_VariableFeatures.pdf")

# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 30)
ggsave("figures/allTcells_elbowplot.pdf")

# Evaling cluster stability: Bootstrap
nPCs<- 20
myClusterFUN <- function(seurat.object) {
    g <- seurat.object %>% FindNeighbors(dims = 1:nPCs, verbose=F) %>% FindClusters(resolution = 0.5, verbose=F)
    g$seurat_clusters
}
originals <- s$seurat_clusters
set.seed(111)
coassign <- bootstrapCluster(s, FUN=myClusterFUN, clusters=originals, iterations=50)
pdf("figures/allTcells_ClusterStability_50iter.pdf")
pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE, color=rev(viridis::magma(100)))
dev.off()

# Clustering Figure
p <- DimPlot(s)
ggsave("figures/allTcells_UMAP.pdf", width= 10, height= 10, units= "in")


#Data Frame with information about number of cells removed
x <- data.frame(
	"Orig.count" = length(stats$mt.pre),
	"Mito.U.lnf.removed" =  stats$total.remove,
	"Mito.upper.thres" = stats$mt.lim,
	"Feats.lower.thres" = exp(stats$lnf.lim),
	"Doublet_hi" = table(s$DF_hi.lo)[[1]],
	"Doublet_lo" = table(s$DF_hi.lo)[[2]],
	"Singlet" = table(s$DF_hi.lo)[[3]],
	"Final_count" = length(stats$mt.pre)-stats$total.remove-table(s$DF_hi.lo)[[1]]
	)


####################################
# *** Clustering (SCTransform) ***
####################################
s <- readRDS("saved_objects/allTcells_1.RDS")
library(future)
plan(multicore)
s <- SCTransform(s, vars.to.regress = "percent.mt",verbose=F)
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TR.V")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
VariableFeatures(s) <- head(VariableFeatures(s),1500)
s <- RunPCA(s, features=VariableFeatures(s))
nPCs <- 30
s <- RunUMAP(s, dims = 1:nPCs) %>% FindNeighbors(dims = 1:nPCs) %>% FindClusters(resolution= 0.4)
saveRDS(s, "saved_objects/allTcells_SCT_hvf1500.RDS")

##############################
# *** Diagnostic Plots ***
##############################
## Figures from Object 2
s <- readRDS("saved_objects/allTcells_SCT_hvf1500.RDS")

# Variable Features Save & Figure
VariableFeatures(s)%>%saveRDS(file="saved_objects/allTcells_SCT_hvf.RDS")
p <- VariableFeaturePlot.Tcells.SCT(s, assay="SCT") %>% LabelPoints(points=head(VariableFeatures(s),10), repel= TRUE)
ggsave("figures/allTcells_SCT_VariableFeatures_2.pdf")

# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 40)
ggsave("figures/allTcells_SCT_elbowplot_2.pdf")

# Evaling cluster stability: Bootstrap
nPCs<- 30
myClusterFUN <- function(seurat.object) {
    g <- seurat.object %>% FindNeighbors(dims = 1:nPCs, verbose=F) %>% FindClusters(resolution = 0.5, verbose=F)
    g$seurat_clusters
}
originals <- s$seurat_clusters
set.seed(111)
coassign <- bootstrapCluster(s, FUN=myClusterFUN, clusters=originals, iterations=100)
dim(coassign)
pdf("figures/allTcells_SCT_ClusterStability_100iter.pdf")
pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE, color=rev(viridis::magma(100)))
dev.off()

# Clustering Figure
p <- DimPlot(s)
ggsave("figures/allTcells_SCT_UMAP_2.pdf", width= 10, height= 10, units= "in")


#Data Frame with information about number of cells removed
x <- data.frame(
	"Orig.count" = length(stats$mt.pre),
	"Mito.U.lnf.removed" =  stats$total.remove,
	"Mito.upper.thres" = stats$mt.lim,
	"Feats.lower.thres" = exp(stats$lnf.lim),
	"Doublet_hi" = table(s$DF_hi.lo)[[1]],
	"Doublet_lo" = table(s$DF_hi.lo)[[2]],
	"Singlet" = table(s$DF_hi.lo)[[3]],
	"Final_count" = length(stats$mt.pre)-stats$total.remove-table(s$DF_hi.lo)[[1]]
	)

####################################
# 5 - Removing bad clusters
####################################
## Feature plot
s <- readRDS("saved_objects/p020620_CD45_3.RDS")
FeaturePlot(s, features= c('CD3E','CD3G','CD3D','CD28'))
ggsave("figures/p020620_CD45_TcellFeatures.pdf", width= 12, height= 12, units= "in")
FeaturePlot(s, features= c('CD4','CD8A','CD8B','ICOS'))
ggsave("figures/p020620_CD45_TcellFeatures2.pdf", width= 12, height= 12, units= "in")



s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")
stats <- filter_stats(s, save=T, filename="saved_objects/allTcells_filterstats.RDS")
s$cells.to.remove <- stats$cells.to.remove
saveRDS(s,"saved_objects/allTcells_1.RDS")
s <- s[,!s$cells.to.remove]
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
nPCs <- 20
s <- ScaleData(s, features = all.genes, vars.to.regress="percent.mt") %>% RunPCA(features=VariableFeatures(s)) %>% RunUMAP(dims = 1:nPCs) %>% FindNeighbors(dims = 1:nPCs) %>% FindClusters(resolution = 0.5)

saveRDS(s,"saved_objects/allTcells_2.RDS")
