##############################
# 0 - Load librairies
##############################
library(dplyr)
library(ggplot2)
library(Seurat)
library(future)
library(future.apply)
plan("multicore", workers = 10)
#library(scran)
#library(pheatmap)

##############################
# 1 - Definitions & Settings
##############################
'%ni%' = Negate('%in%')
plan("multicore")
options(future.globals.maxSize = 11000 * 1024^2)

############################## 
# 2 - Source file 
##############################
source("scripts/Stat111functions.R")
source("scripts/new_clustering_functions.R")

######################################
# 3 - (intra-Martin) choose CD3
######################################
# s2.rds is 20650 by 67790; corresponds to s2v2 umap3 in figures/CD_martin
#DimPlot(s, label=T)
#ggsave("figures/martin_umap_labeling.pdf")
#DimPlot(s, label=T)
#ggsave("figures/martin_umap.pdf")
s <- readRDS("saved_objects/CD_martin_final/s2.rds")

s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")
stats <- filter_stats(s, save=T, filename="saved_objects/martin_collected_filterstats.RDS")
#custom removal based on LNFxMITO plot
s$cells.to.remove <-  s$nFeature_RNA > exp(median(stats$lnf.pre)+3*mad(stats$lnf.pre))
# 869 cells removed; nFeature_RNA > 2177
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
s <- ScaleData(s, features = all.genes, vars.to.regress=c("percent.mt","patient.id")) %>% RunPCA(features=VariableFeatures(s))
s <-  cluster_umap(s, npc=15, k.param=120, resolution=0.5, min.dist=0.3)
DimPlot(s, label=T)
ggsave("figures/martin_npc15_k120_umap.pdf")
# Do the rest
s <- ScaleData(s, features = all.genes, vars.to.regress=c("percent.mt")) %>% RunPCA(features=VariableFeatures(s))
s <-  cluster_umap(s, npc=15, k.param=100, resolution=0.5, min.dist=0.3)
DimPlot(s, label=T)
ggsave("figures/martin_nopatientregress_npc15_k100_umap.pdf")


######################################
# 3 - labeling (percent.mt and patient.id)
######################################
library(SingleR)
library(BiocParallel)
martin <- read.csv(file="data/martin_categories.csv", header=TRUE,sep=",") %>% process_martin()
pred.Martin <- SingleR(method='cluster',  
		test=s[['RNA']]@data, 
		clusters= s$seurat_clusters,
		ref=martin, 
		labels=colnames(martin),
		BPPARAM=MulticoreParam(10))
s$SingleR.Martin <- s$seurat_clusters
levels(s$SingleR.Martin) <- pred.Martin$labels

ref <- MonacoImmuneData()
pred.Monaco <- SingleR(method="cluster",
		test=s[['RNA']]@data,
		clusters= s$seurat_clusters, 
		ref=ref, 
		labels=ref$label.fine, 
		BPPARAM=MulticoreParam(10))
s$SingleR.Monaco <- s$seurat_clusters
levels(s$SingleR.Monaco) <- pred.Monaco$labels

require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco")
ggarrange(p1,p3, ncol=2, nrow=1)
ggsave("figures/martin_nopatientregress_npc15_k120_labeling.pdf", width= 16, height= 6, units= "in")

######################################
# 3 - labeling (percent.mt only)
######################################
library(SingleR)
library(BiocParallel)
martin <- read.csv(file="data/martin_categories.csv", header=TRUE,sep=",") %>% process_martin()
pred.Martin <- SingleR(method='cluster',  
		test=s[['RNA']]@data, 
		clusters= s$seurat_clusters,
		ref=martin, 
		labels=colnames(martin),
		BPPARAM=MulticoreParam(10))
s$SingleR.Martin <- s$seurat_clusters
levels(s$SingleR.Martin) <- pred.Martin$labels

ref <- MonacoImmuneData()
pred.Monaco <- SingleR(method="cluster",
		test=s[['RNA']]@data,
		clusters= s$seurat_clusters, 
		ref=ref, 
		labels=ref$label.fine, 
		BPPARAM=MulticoreParam(10))
s$SingleR.Monaco <- s$seurat_clusters
levels(s$SingleR.Monaco) <- pred.Monaco$labels

require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco")
ggarrange(p1,p3, ncol=2, nrow=1)
ggsave("figures/martin_npc15_k100_labeling.pdf", width= 16, height= 6, units= "in")

######################################
# 3 - (intra-Martin) choose CD3 from only percent.mt regress
######################################
FeaturePlot(s, features= c('CD3E','CD3G','CD3D','CD28'))
ggsave("figures/martin15100_TcellFeatures.pdf", width= 12, height= 12, units= "in")
FeaturePlot(s, features= c('CD4','CD8A','CD8B','ICOS'))
ggsave("figures/martin15100_TcellFeatures2.pdf", width= 12, height= 12, units= "in")
# Feature plot shows that clusters 0,2,3,5,6 are T cells

s <- subset(s,idents= c(0,2,3,5,6))
saveRDS(s, "saved_objects/martinTcells.rds") #38604 cells
######################################
# 3 - Recluster the CD3 cells from only percent.mt regress
######################################
s <- readRDS("saved_objects/martinTcells.rds")
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
s <- ScaleData(s, features = all.genes, vars.to.regress=c("percent.mt")) %>% RunPCA(features=VariableFeatures(s))
s <-  cluster_umap(s, npc=15, k.param=50, resolution=0.5, min.dist=0.3)
s1 <-  cluster_umap(s, npc=15, k.param=100, resolution=0.5, min.dist=0.3)
s2 <-  cluster_umap(s, npc=15, k.param=200, resolution=0.5, min.dist=0.3)

p <- ElbowPlot(s)
ggsave("figures/martin_CD3_elbowplot.pdf")

p <- DimPlot(s)
ggsave("figures/martin1550_CD3_umap.pdf")
p <- DimPlot(s1)
ggsave("figures/martin15100_CD3_umap.pdf")
p <- DimPlot(s1, group.by="patient.id")
ggsave("figures/martin15100_CD3_groupby_patient.pdf")
p <- DimPlot(s2)
ggsave("figures/martin15200_CD3_umap.pdf")
p <- DimPlot(s2, group.by="patient.id")
ggsave("figures/martin15200_CD3_groupby_patient.pdf")

inflamed <-c("_1$","_3$","_6$","_7$","_10$","_12$","_14$","_16$","_18$")
toMatch <- paste(inflamed,collapse="|")
s1$status <- ifelse(grepl(colnames(s), pattern=toMatch),"inflamed","non-inflamed")
p <- DimPlot(s1, group.by="status")
ggsave("figures/martin15100_CD3_groupby_status.pdf")

######################################
# 3 - Labeling
######################################
load('data/ssuo_CD3/seurat.object.RData')
Tcell <- process_ssuo_Tcells(Tcell)
martin <- read.csv(file="data/martin_categories.csv", header=TRUE,sep=",") %>% process_martin()

require(SingleR)
require(BiocParallel)
pred.CPI <- SingleR(method='cluster',
	test=as.matrix(s[['RNA']]@data), 
	clusters= s$seurat_clusters,
	ref=as.matrix(Tcell[['RNA']]@data), 
	labels=Tcell$names,
	de.method="wilcox",
	BPPARAM=MulticoreParam(10))
s$SingleR.Luoma <- s$seurat_clusters
levels(s$SingleR.Luoma) <- pred.CPI$labels

pred.Martin <- SingleR(method='cluster',  
	test=s[['RNA']]@data, 
	clusters= s$seurat_clusters,
	ref=martin, 
	labels=colnames(martin),
	BPPARAM=MulticoreParam(10))
s$SingleR.Martin <- s$seurat_clusters
levels(s$SingleR.Martin) <- pred.Martin$labels

ref <- MonacoImmuneData()
pred.Monaco <- SingleR(method="cluster",
	test=s[['RNA']]@data,
	clusters= s$seurat_clusters, 
	ref=ref, 
	labels=ref$label.fine, 
	BPPARAM=MulticoreParam(10))
s$SingleR.Monaco <- s$seurat_clusters
levels(s$SingleR.Monaco) <- pred.Monaco$labels

#saveRDS(list(pred.CPI,pred.Martin,pred.Monaco),"saved_objects/label_info_for_cdcpi1550.rds")
require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin") + labs(title="Martin reference")
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma") + labs(title="Luoma-Suo reference")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco") + labs(title="Monaco reference")
ggarrange(p1,p2,p3, ncol=3, nrow=1)
ggsave("figures/martin15200_labeling.pdf", width= 24, height= 6, units= "in")

saveRDS(s,"saved_objects/martinTcells15200.rds")

######################################
# 3 - Label for before reclustering
######################################
load('data/ssuo_CD3/seurat.object.RData')
Tcell <- process_ssuo_Tcells(Tcell)
martin <- read.csv(file="data/martin_categories.csv", header=TRUE,sep=",") %>% process_martin()

require(SingleR)
require(BiocParallel)
pred.CPI <- SingleR(method='cluster',
	test=as.matrix(s[['RNA']]@data), 
	clusters= s$seurat_clusters,
	ref=as.matrix(Tcell[['RNA']]@data), 
	labels=Tcell$names,
	de.method="wilcox",
	BPPARAM=MulticoreParam(10))
s$SingleR.Luoma <- s$seurat_clusters
levels(s$SingleR.Luoma) <- pred.CPI$labels

pred.Martin <- SingleR(method='cluster',  
	test=s[['RNA']]@data, 
	clusters= s$seurat_clusters,
	ref=martin, 
	labels=colnames(martin),
	BPPARAM=MulticoreParam(10))
s$SingleR.Martin <- s$seurat_clusters
levels(s$SingleR.Martin) <- pred.Martin$labels

ref <- MonacoImmuneData()
pred.Monaco <- SingleR(method="cluster",
	test=s[['RNA']]@data,
	clusters= s$seurat_clusters, 
	ref=ref, 
	labels=ref$label.fine, 
	BPPARAM=MulticoreParam(10))
s$SingleR.Monaco <- s$seurat_clusters
levels(s$SingleR.Monaco) <- pred.Monaco$labels

#saveRDS(list(pred.CPI,pred.Martin,pred.Monaco),"saved_objects/label_info_for_cdcpi1550.rds")
require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin") + labs(title="Martin reference")
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma") + labs(title="Luoma-Suo reference")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco") + labs(title="Monaco reference")
ggarrange(p1,p2,p3, ncol=3, nrow=1)
ggsave("figures/martin_CD3_before_reclustering_labeling.pdf", width= 24, height= 6, units= "in")


######################################
# 3 - Recluster the CD3 cells from only percent.mt regress without any regress out.
######################################
s <- readRDS("saved_objects/martinTcells.rds")
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
s <- ScaleData(s, features = all.genes) %>% RunPCA(features=VariableFeatures(s))
s <-  cluster_umap(s, npc=15, k.param=100, resolution=0.5, min.dist=0.3)

p <- ElbowPlot(s)
ggsave("figures/martin_CD3_noregress_elbowplot.pdf")

p <- DimPlot(s)
ggsave("figures/martin15100_CD3_noregress_umap.pdf")
p <- DimPlot(s, group.by="patient.id")
ggsave("figures/martin15100_CD3_noregress_groupby_patient.pdf")

