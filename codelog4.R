############################################################
## Quality control and clustering of PD1 samples
## edward_kim@college.harvard.edu - June 2020
#############################################################
​

##############################
# 0 - Load librairies
##############################
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(Seurat)
library(sctransform)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("scran")
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

​############################## 
# 2 - Source file 
##############################
source("scripts/Stat111functions.R")
#source("/scripts/new_clustering_functions.R")
load("data/ssuo_CD3/seurat.object.RData")

##############################
# 3 - Apply QC on each sample
##############################
DFCI1294 <- Read10X(data.dir = "data/DFCI-1294-CD3/outs/filtered_feature_bc_matrix")
s <- CreateSeuratObject(counts = DFCI1294, project = "DFCI-1294_CD3", min.cells = 0, min.features = 0) %>% PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")
stats <- filter_stats(s, save=T, filename="saved_objects/DFCI-1294_CD3_filterstats.RDS")
s$cells.to.remove <- stats$cells.to.remove
saveRDS(s,"saved_objects/DFCI-1294_CD3_1.RDS")
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
s <- label_doublets(s,nPCs=nPCs, save.pK.figure=T, filename= "figures/DFCI-1294_CD3_pK_ggplot.pdf")

saveRDS(s,"saved_objects/DFCI-1294_CD3_2.RDS")
s <- subset(s, subset = DF_hi.lo != "Doublet_hi")
saveRDS(s,"saved_objects/DFCI-1294_CD3_3.RDS")


##############################
# 4 - Diagnostic Plots
##############################
## Figures from QC
stats <- readRDS("saved_objects/DFCI-1294_CD3_filterstats.RDS")

# Mito% Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre mito%", length(stats$mt.pre)),rep("post mito%", length(stats$mt.post))), percent.mito = c(stats$mt.pre,stats$mt.post))) + geom_violin(aes(x= x, y=percent.mito))
ggsave("figures/DFCI-1294_CD3_MITOviolin.pdf")

# Log(No. of Features) Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre lnf", length(stats$lnf.pre)),rep("post lnf", length(stats$lnf.post))), log.num.feats = c(stats$lnf.pre,stats$lnf.post))) + geom_violin(aes(x= x, y=log.num.feats))
ggsave("figures/DFCI-1294_CD3_LNFviolin.pdf")

# Log.no.features vs. mito% Figure
p <- readRDS("saved_objects/DFCI-1294_CD3_1.RDS") %>% FeatureScatter(feature1="nFeature_RNA", feature2="percent.mt", group.by="cells.to.remove")
ggsave("figures/DFCI-1294_CD3_LNFxMITO.pdf")

# Mito% Model Fitting Figure
p <- plot_kde(stats$mt.pre, kernel = "gaussian", bw=0.5, lab.x = "mitochondrial gene percentage (per cell)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$mt.median,stats$mt.mad)}, 
	support=seq(0,100,0.1))
p + geom_vline(xintercept= stats$mt.lim, color = "red") +
	geom_vline(xintercept= stats$mt.median-3*stats$mt.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$mt.median+3*stats$mt.mad, linetype = "dotted")
ggsave("figures/DFCI-1294_CD3_MITO-KDE-Normal.pdf")

# Log.no.features Model Fitting Figure
p <- plot_kde(stats$lnf.pre, kernel = "gaussian", bw=0.01, lab.x = "Log(No. of Features)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$lnf.median,stats$lnf.mad)})
p + geom_vline(xintercept= stats$lnf.lim, color = "red") +
	geom_vline(xintercept= stats$lnf.median-3*stats$lnf.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$lnf.median+3*stats$lnf.mad, linetype = "dotted")
ggsave("figures/DFCI-1294_CD3_LNF-KDE-Normal.pdf")


## Figures from Object 2
s <- readRDS("saved_objects/DFCI-1294_CD3_2.RDS")

# Variable Features Save & Figure
VariableFeatures(s)%>%saveRDS(file="saved_objects/DFCI-1294_CD3_hvf.RDS")
p <- VariableFeaturePlot.Tcells(s) %>% LabelPoints(points=head(VariableFeatures(s),10), repel= TRUE)
ggsave("figures/DFCI-1294_CD3_VariableFeatures.pdf")

# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 30)
ggsave("figures/DFCI-1294_CD3_elbowplot.pdf")

# Evaling cluster stability: Bootstrap
nPCs<- 20
myClusterFUN <- function(seurat.object) {
    g <- seurat.object %>% FindNeighbors(dims = 1:nPCs, verbose=F) %>% FindClusters(resolution = 0.5, verbose=F)
    g$seurat_clusters
}

originals <- s$seurat_clusters
set.seed(111)
coassign <- bootstrapCluster(s, FUN=myClusterFUN, clusters=originals, iterations=1000)
dim(coassign)

pdf("figures/DFCI-1294_CD3_ClusterStability_1000iter.pdf")
pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE, color=rev(viridis::magma(100)))
dev.off()

# Clustering and Doublet Detection Figure
p <- DimPlot(s)
q <- DimPlot(s, reduction = "umap", group.by="DF_hi.lo")+ scale_colour_manual(values=c("red","yellow","gray"))
p + q
ggsave("figures/DFCI-1294_CD3_UMAP_doublets.pdf", width= 12, height= 6, units= "in")

#Vector with information about number of cells removed
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
x


#########################
### DFCI 1545 #######
#########################
DFCI1545 <- Read10X(data.dir = "data/DFCI-1545-CD3/outs/filtered_feature_bc_matrix")
s <- CreateSeuratObject(counts = DFCI1545, project = "DFCI-1545_CD3", min.cells = 0, min.features = 0) %>% PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")
stats <- filter_stats(s, save=T, filename="saved_objects/DFCI-1545_CD3_filterstats.RDS")
s$cells.to.remove <- stats$cells.to.remove
saveRDS(s,"saved_objects/DFCI-1545_CD3_1.RDS")
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
s <- label_doublets(s,nPCs=nPCs, save.pK.figure=T, filename= "figures/DFCI-1545_CD3_pK_ggplot.pdf")

saveRDS(s,"saved_objects/DFCI-1545_CD3_2.RDS")
s <- subset(s, subset = DF_hi.lo != "Doublet_hi")
saveRDS(s,"saved_objects/DFCI-1545_CD3_3.RDS")


##############################
# 4 - Diagnostic Plots
##############################
## Figures from QC
stats <- readRDS("saved_objects/DFCI-1545_CD3_filterstats.RDS")

# Mito% Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre mito%", length(stats$mt.pre)),rep("post mito%", length(stats$mt.post))), percent.mito = c(stats$mt.pre,stats$mt.post))) + geom_violin(aes(x= x, y=percent.mito))
ggsave("figures/DFCI-1545_CD3_MITOviolin.pdf")

# Log(No. of Features) Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre lnf", length(stats$lnf.pre)),rep("post lnf", length(stats$lnf.post))), log.num.feats = c(stats$lnf.pre,stats$lnf.post))) + geom_violin(aes(x= x, y=log.num.feats))
ggsave("figures/DFCI-1545_CD3_LNFviolin.pdf")

# Log.no.features vs. mito% Figure
p <- readRDS("saved_objects/DFCI-1545_CD3_1.RDS") %>% FeatureScatter(feature1="nFeature_RNA", feature2="percent.mt", group.by="cells.to.remove")
ggsave("figures/DFCI-1545_CD3_LNFxMITO.pdf")

# Mito% Model Fitting Figure
p <- plot_kde(stats$mt.pre, kernel = "gaussian", bw=0.5, lab.x = "mitochondrial gene percentage (per cell)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$mt.median,stats$mt.mad)}, 
	support=seq(0,100,0.1))
p + geom_vline(xintercept= stats$mt.lim, color = "red") +
	geom_vline(xintercept= stats$mt.median-3*stats$mt.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$mt.median+3*stats$mt.mad, linetype = "dotted")
ggsave("figures/DFCI-1545_CD3_MITO-KDE-Normal.pdf")

# Log.no.features Model Fitting Figure
p <- plot_kde(stats$lnf.pre, kernel = "gaussian", bw=0.01, lab.x = "Log(No. of Features)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$lnf.median,stats$lnf.mad)})
p + geom_vline(xintercept= stats$lnf.lim, color = "red") +
	geom_vline(xintercept= stats$lnf.median-3*stats$lnf.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$lnf.median+3*stats$lnf.mad, linetype = "dotted")
ggsave("figures/DFCI-1545_CD3_LNF-KDE-Normal.pdf")


## Figures from Object 2
s <- readRDS("saved_objects/DFCI-1545_CD3_2.RDS")

# Variable Features Save & Figure
VariableFeatures(s)%>%saveRDS(file="saved_objects/DFCI-1545_CD3_hvf.RDS")
p <- VariableFeaturePlot.Tcells(s) %>% LabelPoints(points=head(VariableFeatures(s),10), repel= TRUE)
ggsave("figures/DFCI-1545_CD3_VariableFeatures.pdf")

# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 30)
ggsave("figures/DFCI-1545_CD3_elbowplot.pdf")

# Clustering and Doublet Detection Figure
p <- DimPlot(s)
q <- DimPlot(s, reduction = "umap", group.by="DF_hi.lo")+ scale_colour_manual(values=c("red","yellow","gray"))
p + q
ggsave("figures/DFCI-1545_CD3_UMAP_doublets.pdf", width= 12, height= 6, units= "in")

#Vector with information about number of cells removed
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
x

#########################
### DFCI 1618 #######
#########################
DFCI1618 <- Read10X(data.dir = "data/DFCI-1618-CD3/outs/filtered_feature_bc_matrix")
s <- CreateSeuratObject(counts = DFCI1618, project = "DFCI-1618_CD3", min.cells = 0, min.features = 0) %>% PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")
stats <- filter_stats(s, save=T, filename="saved_objects/DFCI-1618_CD3_filterstats.RDS")
s$cells.to.remove <- stats$cells.to.remove
saveRDS(s,"saved_objects/DFCI-1618_CD3_1.RDS")
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
s <- label_doublets(s,nPCs=nPCs, save.pK.figure=T, filename= "figures/DFCI-1618_CD3_pK_ggplot.pdf")

saveRDS(s,"saved_objects/DFCI-1618_CD3_2.RDS")
s <- subset(s, subset = DF_hi.lo != "Doublet_hi")
saveRDS(s,"saved_objects/DFCI-1618_CD3_3.RDS")


##############################
# 4 - Diagnostic Plots
##############################
## Figures from QC
stats <- readRDS("saved_objects/DFCI-1618_CD3_filterstats.RDS")

# Mito% Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre mito%", length(stats$mt.pre)),rep("post mito%", length(stats$mt.post))), percent.mito = c(stats$mt.pre,stats$mt.post))) + geom_violin(aes(x= x, y=percent.mito))
ggsave("figures/DFCI-1618_CD3_MITOviolin.pdf")

# Log(No. of Features) Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre lnf", length(stats$lnf.pre)),rep("post lnf", length(stats$lnf.post))), log.num.feats = c(stats$lnf.pre,stats$lnf.post))) + geom_violin(aes(x= x, y=log.num.feats))
ggsave("figures/DFCI-1618_CD3_LNFviolin.pdf")

# Log.no.features vs. mito% Figure
p <- readRDS("saved_objects/DFCI-1618_CD3_1.RDS") %>% FeatureScatter(feature1="nFeature_RNA", feature2="percent.mt", group.by="cells.to.remove")
ggsave("figures/DFCI-1618_CD3_LNFxMITO.pdf")

# Mito% Model Fitting Figure
p <- plot_kde(stats$mt.pre, kernel = "gaussian", bw=0.5, lab.x = "mitochondrial gene percentage (per cell)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$mt.median,stats$mt.mad)}, 
	support=seq(0,100,0.1))
p + geom_vline(xintercept= stats$mt.lim, color = "red") +
	geom_vline(xintercept= stats$mt.median-3*stats$mt.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$mt.median+3*stats$mt.mad, linetype = "dotted")
ggsave("figures/DFCI-1618_CD3_MITO-KDE-Normal.pdf")

# Log.no.features Model Fitting Figure
p <- plot_kde(stats$lnf.pre, kernel = "gaussian", bw=0.01, lab.x = "Log(No. of Features)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$lnf.median,stats$lnf.mad)})
p + geom_vline(xintercept= stats$lnf.lim, color = "red") +
	geom_vline(xintercept= stats$lnf.median-3*stats$lnf.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$lnf.median+3*stats$lnf.mad, linetype = "dotted")
ggsave("figures/DFCI-1618_CD3_LNF-KDE-Normal.pdf")


## Figures from Object 2
s <- readRDS("saved_objects/DFCI-1618_CD3_2.RDS")

# Variable Features Save & Figure
VariableFeatures(s)%>%saveRDS(file="saved_objects/DFCI-1618_CD3_hvf.RDS")
p <- VariableFeaturePlot.Tcells(s) %>% LabelPoints(points=head(VariableFeatures(s),10), repel= TRUE)
ggsave("figures/DFCI-1618_CD3_VariableFeatures.pdf")

# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 30)
ggsave("figures/DFCI-1618_CD3_elbowplot.pdf")

# Evaling cluster stability: Bootstrap
nPCs<- 20
myClusterFUN <- function(seurat.object) {
    g <- seurat.object %>% FindNeighbors(dims = 1:nPCs, verbose=F) %>% FindClusters(resolution = 0.5, verbose=F)
    g$seurat_clusters
}
originals <- s$seurat_clusters
set.seed(111)
coassign <- bootstrapCluster(s, FUN=myClusterFUN, clusters=originals, iterations=100)
dim(coassign)
pdf("figures/DFCI-1618_CD3_ClusterStability_100iter.pdf")
pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE, color=rev(viridis::magma(100)))
dev.off()

# Clustering and Doublet Detection Figure
p <- DimPlot(s)
q <- DimPlot(s, reduction = "umap", group.by="DF_hi.lo")+ scale_colour_manual(values=c("red","yellow","gray"))
p + q
ggsave("figures/DFCI-1618_CD3_UMAP_doublets.pdf", width= 12, height= 6, units= "in")

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

#########################
###    p101519-CD3  #####
#########################
p101519 <- Read10X(data.dir = "data/p101519-CD3/outs/filtered_feature_bc_matrix")
s <- CreateSeuratObject(counts = p101519, project = "p101519_CD3", min.cells = 0, min.features = 0) %>% PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")
stats <- filter_stats(s, save=T, filename="saved_objects/p101519_CD3_filterstats.RDS")
s$cells.to.remove <- stats$cells.to.remove
saveRDS(s,"saved_objects/p101519_CD3_1.RDS")
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
s <- label_doublets(s,nPCs=nPCs, save.pK.figure=T, filename= "figures/p101519_CD3_pK_ggplot.pdf")

saveRDS(s,"saved_objects/p101519_CD3_2.RDS")
s <- subset(s, subset = DF_hi.lo != "Doublet_hi")
saveRDS(s,"saved_objects/p101519_CD3_3.RDS")


##############################
# 4 - Diagnostic Plots
##############################
## Figures from QC
stats <- readRDS("saved_objects/p101519_CD3_filterstats.RDS")

# Mito% Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre mito%", length(stats$mt.pre)),rep("post mito%", length(stats$mt.post))), percent.mito = c(stats$mt.pre,stats$mt.post))) + geom_violin(aes(x= x, y=percent.mito))
ggsave("figures/p101519_CD3_MITOviolin.pdf")

# Log(No. of Features) Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre lnf", length(stats$lnf.pre)),rep("post lnf", length(stats$lnf.post))), log.num.feats = c(stats$lnf.pre,stats$lnf.post))) + geom_violin(aes(x= x, y=log.num.feats))
ggsave("figures/p101519_CD3_LNFviolin.pdf")

# Log.no.features vs. mito% Figure
p <- readRDS("saved_objects/p101519_CD3_1.RDS") %>% FeatureScatter(feature1="nFeature_RNA", feature2="percent.mt", group.by="cells.to.remove")
ggsave("figures/p101519_CD3_LNFxMITO.pdf")

# Mito% Model Fitting Figure
p <- plot_kde(stats$mt.pre, kernel = "gaussian", bw=0.5, lab.x = "mitochondrial gene percentage (per cell)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$mt.median,stats$mt.mad)}, 
	support=seq(0,100,0.1))
p + geom_vline(xintercept= stats$mt.lim, color = "red") +
	geom_vline(xintercept= stats$mt.median-3*stats$mt.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$mt.median+3*stats$mt.mad, linetype = "dotted")
ggsave("figures/p101519_CD3_MITO-KDE-Normal.pdf")

# Log.no.features Model Fitting Figure
p <- plot_kde(stats$lnf.pre, kernel = "gaussian", bw=0.01, lab.x = "Log(No. of Features)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$lnf.median,stats$lnf.mad)})
p + geom_vline(xintercept= stats$lnf.lim, color = "red") +
	geom_vline(xintercept= stats$lnf.median-3*stats$lnf.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$lnf.median+3*stats$lnf.mad, linetype = "dotted")
ggsave("figures/p101519_CD3_LNF-KDE-Normal.pdf")


## Figures from Object 2
s <- readRDS("saved_objects/p101519_CD3_2.RDS")

# Variable Features Save & Figure
VariableFeatures(s)%>%saveRDS(file="saved_objects/p101519_CD3_hvf.RDS")
p <- VariableFeaturePlot.Tcells(s) %>% LabelPoints(points=head(VariableFeatures(s),10), repel= TRUE)
ggsave("figures/p101519_CD3_VariableFeatures.pdf")

# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 30)
ggsave("figures/p101519_CD3_elbowplot.pdf")

# Evaling cluster stability: Bootstrap
nPCs<- 20
myClusterFUN <- function(seurat.object) {
    g <- seurat.object %>% FindNeighbors(dims = 1:nPCs, verbose=F) %>% FindClusters(resolution = 0.5, verbose=F)
    g$seurat_clusters
}
originals <- s$seurat_clusters
set.seed(111)
coassign <- bootstrapCluster(s, FUN=myClusterFUN, clusters=originals, iterations=100)
dim(coassign)
pdf("figures/p101519_CD3_ClusterStability_100iter.pdf")
pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE, color=rev(viridis::magma(100)))
dev.off()

# Clustering and Doublet Detection Figure
p <- DimPlot(s)
q <- DimPlot(s, reduction = "umap", group.by="DF_hi.lo")+ scale_colour_manual(values=c("red","yellow","gray"))
p + q
ggsave("figures/p101519_CD3_UMAP_doublets.pdf", width= 12, height= 6, units= "in")

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

#########################
###    p020620-CD45 (npc 20, hvf filtered)   #####
#########################
p020620 <- Read10X(data.dir = "data/p020620-CD45/outs/filtered_feature_bc_matrix")
s <- CreateSeuratObject(counts = p020620, project = "p020620_CD45", min.cells = 0, min.features = 0) %>% PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")
stats <- filter_stats(s, save=T, filename="saved_objects/p020620_CD45_filterstats.RDS")
s$cells.to.remove <- stats$cells.to.remove
saveRDS(s,"saved_objects/p020620_CD45_1.RDS")
s <- s[,!s$cells.to.remove]
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
nPCs <- 15
s <- ScaleData(s, features = all.genes, vars.to.regress="percent.mt") %>% RunPCA(features=VariableFeatures(s)) %>% RunUMAP(dims = 1:nPCs) %>% FindNeighbors(dims = 1:nPCs) %>% FindClusters(resolution = 0.5)
s <- label_doublets(s,nPCs=nPCs, save.pK.figure=T, filename= "figures/p020620_CD45_pK_ggplot.pdf")

saveRDS(s,"saved_objects/p020620_CD45_2.RDS")
s <- subset(s, subset = DF_hi.lo != "Doublet_hi")
saveRDS(s,"saved_objects/p020620_CD45_3.RDS")


##############################
# 4 - Diagnostic Plots
##############################
## Figures from QC
stats <- readRDS("saved_objects/p020620_CD45_filterstats.RDS")

# Mito% Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre mito%", length(stats$mt.pre)),rep("post mito%", length(stats$mt.post))), percent.mito = c(stats$mt.pre,stats$mt.post))) + geom_violin(aes(x= x, y=percent.mito))
ggsave("figures/p020620_CD45_MITOviolin.pdf")

# Log(No. of Features) Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre lnf", length(stats$lnf.pre)),rep("post lnf", length(stats$lnf.post))), log.num.feats = c(stats$lnf.pre,stats$lnf.post))) + geom_violin(aes(x= x, y=log.num.feats))
ggsave("figures/p020620_CD45_LNFviolin.pdf")

# Log.no.features vs. mito% Figure
p <- readRDS("saved_objects/p020620_CD45_1.RDS") %>% FeatureScatter(feature1="nFeature_RNA", feature2="percent.mt", group.by="cells.to.remove")
ggsave("figures/p020620_CD45_LNFxMITO.pdf")

# Mito% Model Fitting Figure
p <- plot_kde(stats$mt.pre, kernel = "gaussian", bw=0.5, lab.x = "mitochondrial gene percentage (per cell)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$mt.median,stats$mt.mad)}, 
	support=seq(0,100,0.1))
p + geom_vline(xintercept= stats$mt.lim, color = "red") +
	geom_vline(xintercept= stats$mt.median-3*stats$mt.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$mt.median+3*stats$mt.mad, linetype = "dotted")
ggsave("figures/p020620_CD45_MITO-KDE-Normal.pdf")

# Log.no.features Model Fitting Figure
p <- plot_kde(stats$lnf.pre, kernel = "gaussian", bw=0.01, lab.x = "Log(No. of Features)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$lnf.median,stats$lnf.mad)})
p + geom_vline(xintercept= stats$lnf.lim, color = "red") +
	geom_vline(xintercept= stats$lnf.median-3*stats$lnf.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$lnf.median+3*stats$lnf.mad, linetype = "dotted")
ggsave("figures/p020620_CD45_LNF-KDE-Normal.pdf")


## Figures from Object 2
s <- readRDS("saved_objects/p020620_CD45_2.RDS")

# Variable Features Save & Figure
VariableFeatures(s)%>%saveRDS(file="saved_objects/p020620_CD45_hvf.RDS")
p <- VariableFeaturePlot.Tcells(s) %>% LabelPoints(points=head(VariableFeatures(s),10), repel= TRUE)
ggsave("figures/p020620_CD45_VariableFeatures.pdf")

# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 30)
ggsave("figures/p020620_CD45_elbowplot.pdf")

# Evaling cluster stability: Bootstrap
nPCs<- 20
myClusterFUN <- function(seurat.object) {
    g <- seurat.object %>% FindNeighbors(dims = 1:nPCs, verbose=F) %>% FindClusters(resolution = 0.5, verbose=F)
    g$seurat_clusters
}
originals <- s$seurat_clusters
set.seed(111)
coassign <- bootstrapCluster(s, FUN=myClusterFUN, clusters=originals, iterations=100)
dim(coassign)
pdf("figures/p020620_CD45_ClusterStability_100iter.pdf")
pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE, color=rev(viridis::magma(100)))
dev.off()

# Clustering and Doublet Detection Figure
p <- DimPlot(s)
q <- DimPlot(s, reduction = "umap", group.by="DF_hi.lo")+ scale_colour_manual(values=c("red","yellow","gray"))
p + q
ggsave("figures/p020620_CD45_UMAP_doublets.pdf", width= 12, height= 6, units= "in")

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

#########################
###    p020620-CD45 (npc 15)   #####
#########################
s <- readRDS("saved_objects/p020620_CD45_1.RDS")
s <- s[,!s$cells.to.remove]
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
#hvf <- VariableFeatures(s)
#'%ni%' = Negate('%in%')
#all.genes <- rownames(s)
#trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
#igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
#VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
nPCs <- 15
s <- ScaleData(s, features = all.genes, vars.to.regress="percent.mt") %>% RunPCA(features=VariableFeatures(s)) %>% RunUMAP(dims = 1:nPCs) %>% FindNeighbors(dims = 1:nPCs) %>% FindClusters(resolution = 0.5)
s <- label_doublets(s,nPCs=nPCs, save.pK.figure=T, filename= "figures/p020620_CD45v2_pK_ggplot.pdf")

saveRDS(s,"saved_objects/p020620_CD45v2_2.RDS")
s <- subset(s, subset = DF_hi.lo != "Doublet_hi")
saveRDS(s,"saved_objects/p020620_CD45v2_3.RDS")


##############################
# 4 - Diagnostic Plots
##############################
## Figures from Object 2
s <- readRDS("saved_objects/p020620_CD45v2_2.RDS")

# Variable Features Save & Figure
VariableFeatures(s)%>%saveRDS(file="saved_objects/p020620_CD45v2_hvf.RDS")
p <- VariableFeaturePlot(s) %>% LabelPoints(points=head(VariableFeatures(s),10), repel= TRUE)
ggsave("figures/p020620_CD45v2_VariableFeatures.pdf")

# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 30)
ggsave("figures/p020620_CD45v2_elbowplot.pdf")

# Evaling cluster stability: Bootstrap
nPCs<- 20
myClusterFUN <- function(seurat.object) {
    g <- seurat.object %>% FindNeighbors(dims = 1:nPCs, verbose=F) %>% FindClusters(resolution = 0.5, verbose=F)
    g$seurat_clusters
}
originals <- s$seurat_clusters
set.seed(111)
coassign <- bootstrapCluster(s, FUN=myClusterFUN, clusters=originals, iterations=100)
dim(coassign)
pdf("figures/p020620_CD45_ClusterStability_100iter.pdf")
pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE, color=rev(viridis::magma(100)))
dev.off()

# Clustering and Doublet Detection Figure
p <- DimPlot(s)
q <- DimPlot(s, reduction = "umap", group.by="DF_hi.lo")+ scale_colour_manual(values=c("red","yellow","gray"))
p + q
ggsave("figures/p020620_CD45v2_UMAP_doublets.pdf", width= 12, height= 6, units= "in")

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

#############################
####   Subset T cells    ####
#############################
#'CD3E','CD3G','CD3D','CD28','CD4','CD8A','CD8B','ICOS','CD2','CD7'
s <- readRDS("saved_objects/p020620_CD45_3.RDS")
FeaturePlot(s, features= c('CD3E','CD3G','CD3D','CD28'))
ggsave("figures/p020620_CD45_TcellFeatures.pdf", width= 12, height= 12, units= "in")
FeaturePlot(s, features= c('CD4','CD8A','CD8B','ICOS'))
ggsave("figures/p020620_CD45_TcellFeatures2.pdf", width= 12, height= 12, units= "in")
#found clusters 2,3,5,6 are strongly T cells. Not sure what 11 is.
s <- subset(s, idents= c(2,3,5,6))
saveRDS(s,"saved_objects/p020620_CD3subset.RDS")

#############################
#### Data Frame to Latex ####
#############################
get.info <- function(name) {
	stats <- readRDS(paste("saved_objects/",name,"_filterstats.RDS", sep=""))
	s <- readRDS(paste("saved_objects/",name,"_2.RDS", sep=""))
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
	rownames(x) <- name
	x
}
stats <- readRDS("saved_objects/p020620_CD45_filterstats.RDS")
s <- readRDS("saved_objects/p020620_CD45v2_2.RDS")
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
rownames(x) <- "p020620_CD45v2"
y <- rbind(get.info("DFCI-1294_CD3"),
	get.info("DFCI-1545_CD3"),
	get.info("DFCI-1618_CD3"),
	get.info("p101519_CD3"),
	get.info("p020620_CD45"),
	x)

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
xtable(y)

