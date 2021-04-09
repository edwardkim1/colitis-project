############################################################
## Clustering UC Chang
## edward_kim@college.harvard.edu - December 2020
#############################################################

##############################
# 0 - Load librairies
##############################
library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(future)
library(future.apply)
library(hdf5r)
library(parallel)
#library(fgsea)
#library(harmony)
#library(scran)
#library(pheatmap)

##############################
# 1 - Definitions & Settings
##############################
'%ni%' = Negate('%in%')
plan("multicore", workers = 10)
options(future.globals.maxSize = 20000 * 1024^2)

############################## 
# 2 - Source file 
##############################
source("scripts/Stat111functions.R")
source("scripts/new_clustering_functions.R")

####################################
# UC Chang
####################################
#df <- read.table(file = "data/UC_chang/GSM3576396_C9_R_cell-gene_UMI_table.tsv", sep = "\t", header = TRUE)
#df2 <- read.table(file = "data/UC_chang/GSM3576399_U4_R_cell-gene_UMI_table.tsv", sep = "\t", header = TRUE)

pbmc.files <- dir("data/UC_chang/",pattern="pBMC_cell-gene")
rectum.files <- dir("data/UC_chang/",pattern="R_cell-gene")

lapply(paste("data/UC_chang/",pbmc.files,sep=""), 
	function(x) {read.table(x, sep= "\t", header=T, row.names=1)}) -> df.pbmc
lapply(paste("data/UC_chang/",rectum.files,sep=""), 
	function(x) {read.table(x, sep= "\t", header=T, row.names=1)}) -> df.rec
lapply(df.pbmc,nrow) %>% unlist()
lapply(df.rec, nrow) %>% unlist()

lapply(df.pbmc[3:9],nrow) %>% unlist() +
lapply(df.rec[4:10], nrow) %>% unlist()


# create Seurat objects
proj.names <- paste("U",1:7, sep="")
create.pbmc <- function(x) {
	df.pbmc[[x]] %>% as.matrix() %>% t() %>% CreateSeuratObject(project = proj.names[x-2], min.cells = 0, min.features = 0)
}

s.list <- lapply(3:9, create.pbmc)

####################################
# Naive integration with cells from all patients
####################################
s <- merge(s.list[[1]],c(s.list[[2]],s.list[[3]],s.list[[4]],s.list[[5]],s.list[[6]],s.list[[7]]))

########################################################
# 3 - Add patient annotations and inflammed annotations
########################################################
library(stringr)
# rename sampletype to sample.id and merge with orig.ident material
s$sample.id <- s$orig.ident
s$orig.ident.backup <- s$orig.ident
#### mutate orig.ident and edit sample.id
s$orig.ident <- "UC_chang"
#### add patient.id
s$patient.id <- s$sample.id
# add inflamed id for Crohns
s$colitis <- "UC_inflamed"
# remove genes with less than 10
x <- s[["RNA"]]@counts
nnz_row <- tabulate(x@i + 1)
keep <- rownames(x)[nnz_row>10]
s <- subset(s, features = keep)
#10105
#############################
# 4 - perform clustering
#############################
s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")

s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
s <- ScaleData(s, features = all.genes, vars.to.regress=c("percent.mt")) %>% RunPCA(features=VariableFeatures(s))
s <-  cluster_umap(s, npc=20, k.param=120, resolution=0.5, min.dist=0.3)
DimPlot(s, label=T)
# make directory: UC_chang_013021
ggsave("figures/UC_chang_013021/uc20120_umap.pdf")
DimPlot(s, group.by="sample.id" , label=F)
ggsave("figures/UC_chang_013021/uc20120_persample.pdf")

######################################################
# 3 - CD3E and CD4/CD8 and other QC metrics per cluster
######################################################
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/uc20120_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/UC_chang_013021/uc20120_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/uc20120_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/uc20120_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the IGKC and JCHAIN
FeaturePlot(s, features= c('IGKC','JCHAIN'), reduction="umap")
ggsave("figures/UC_chang_013021/uc20120_IGKC-JCHAIN_featureplot.pdf", width= 12, height= 7, units= "in")
# violin plot of IGKC per cluster
p <- VlnPlot(s,features = "IGKC", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave(paste("figures/UC_chang_013021/uc20120_VlnPlot_IGKC_by_cluster.pdf", sep=""), width= 7, height= 7, units= "in")

####################################
# 5 - Subset T cells
####################################
# violin plot of CD3E per cluster
s1 <- s
s <- subset(s1, idents=c("0","1","2","3","5","7"))
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 
## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
s <- ScaleData(s, features = all.genes, vars.to.regress=c("percent.mt")) %>% RunPCA(features=VariableFeatures(s))
s <-  cluster_umap(s, npc=20, k.param=100, resolution=0.5, min.dist=0.3)
DimPlot(s, label=T)
ggsave("figures/UC_chang_013021/ucT20100_umap.pdf")
DimPlot(s, group.by="sample.id" , label=F)
ggsave("figures/UC_chang_013021/ucT20100_persample.pdf")

####################################
# 5 - Redo clustering with anchors
####################################
s.list <- SplitObject(s, split.by="sample.id")
# variable features
s.features <- SelectIntegrationFeatures(object.list = s.list, nfeatures = 3000)
'%ni%' = Negate('%in%')
trvgenes <- s.features[grepl(x=s.features, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- s.features[grepl(x=s.features, pattern = "^IG.V")]
s.features <- s.features[s.features %ni% c(trvgenes,igvgenes)]
# Finding the integration anchors and integrating
s.anchors <- FindIntegrationAnchors(object.list = s.list, 
	normalization.method= "LogNormalize",
	anchor.features= s.features,
	dims = 1:20)
s <- IntegrateData(anchorset = s.anchors)

DefaultAssay(s) <- "integrated"
# Do the rest; no need to regress out percent.mt because those genes are not in the list
s <- ScaleData(s) %>% RunPCA()
# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 20)
ggsave("figures/UC_chang_013021/ucTa20100_anchor_elbowplot.pdf")

s <- FindNeighbors(s, dims = 1:15, k.param = 100) 
# names(s@graphs)
s <- FindClusters(s,resolution = 0.5, graph.name ="integrated_snn") 
s <- RunUMAP(s,min.dist = 0.3, graph="integrated_snn")

#visualize
p <- DimPlot(s, label=T)
ggsave("figures/UC_chang_013021/ucTa20100_PCA.pdf")
p <- DimPlot(s, reduction="umap", label=T)
ggsave("figures/UC_chang_013021/ucTa20100_UMAP.pdf")
p <- DimPlot(s, reduction="umap", group.by="sample.id" , label=F)
ggsave("figures/UC_chang_013021/ucTa20100_persample.pdf")

# Differential Expression
DE_heatmap(s,"figures/UC_chang_013021/ucTa20100_DE_avgexp.pdf")

######################################################
# 3 - CD3E and CD4/CD8 and other QC metrics per cluster
######################################################
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/ucTa20100_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/UC_chang_013021/ucTa20100_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/ucTa20100_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/ucTa20100_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the IGKC and JCHAIN
FeaturePlot(s, features= c('IGKC','JCHAIN'), reduction="umap")
ggsave("figures/UC_chang_013021/ucTa20100_IGKC-JCHAIN_featureplot.pdf", width= 12, height= 7, units= "in")

# violin plot of IGKC per cluster
p <- VlnPlot(s,features = "IGKC", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave(paste("figures/UC_chang_013021/ucTa20100_VlnPlot_IGKC_by_cluster.pdf", sep=""), width= 7, height= 7, units= "in")

# save RDS
saveRDS(s, "saved_objects/UC_chang_013021/ucTa20100.rds")

######################################################
# 3 - Integration with CDCPI reference
######################################################
s1 <- readRDS("saved_objects/UC_chang_013021/ucTa20100.rds")
s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
DefaultAssay(s1) <- "RNA"
s.list <- SplitObject(s1, split.by="sample.id")
s.list <- c(s2,s.list)
# variable features
s.features <- SelectIntegrationFeatures(object.list = s.list, nfeatures = 2000)
'%ni%' = Negate('%in%')
trvgenes <- s.features[grepl(x=s.features, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- s.features[grepl(x=s.features, pattern = "^IG.V")]
s.features <- s.features[s.features %ni% c(trvgenes,igvgenes)]
# Finding the integration anchors and integrating
s.anchors <- FindIntegrationAnchors(object.list = s.list, 
	normalization.method= "LogNormalize",
	reference = 1,
	anchor.features= s.features,
	dims = 1:50)
s <- IntegrateData(anchorset = s.anchors, dims=1:50)

DefaultAssay(s) <- "integrated"
# Do the rest; no need to regress out percent.mt because those genes are not in the list
s <- ScaleData(s) %>% RunPCA()
# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 20)
ggsave("figures/UC_chang_013021/ucTa20100_anchor_elbowplot.pdf")

s <- FindNeighbors(s, dims = 1:15, k.param = 100) 
# names(s@graphs)
s <- FindClusters(s,resolution = 0.5, graph.name ="integrated_snn") 
s <- RunUMAP(s,min.dist = 0.3, graph="integrated_snn")

#visualize
p <- DimPlot(s, label=T)
ggsave("figures/UC_chang_013021/uccpi20100_PCA.pdf")
p <- DimPlot(s, reduction="umap", label=T)
ggsave("figures/UC_chang_013021/uccpi20100_UMAP.pdf")
p <- DimPlot(s, reduction="umap", group.by="sample.id" , label=F)
ggsave("figures/UC_chang_013021/uccpi20100_persample.pdf")

# Differential Expression
DE_heatmap(s,"figures/UC_chang_013021/uccpi20100_DE_avgexp.pdf")

# save integration
saveRDS(s, "saved_objects/UC_chang_013021/uccpi_pbmc.rds")


# did not execute below:

# distribution by sample
data.frame(cluster = s$seurat_clusters, patient= s$sample.id) %>%
        ggplot() + geom_bar(
                mapping = aes(x= patient, fill= cluster),
                position = "fill"
        ) + labs(y="proportion")+ theme_classic() + coord_flip()


ggsave("figures/UC_chang_013021/martin_naive_npc20_k120_distpersample.pdf")


####################################
# 5 - Cluster Annotations
####################################
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

require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin") + labs(title="Martin reference")
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma") + labs(title="Luoma-Suo reference")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco") + labs(title="Monaco reference")
ggarrange(p1,p2,p3, ncol=3, nrow=1)
ggsave("figures/UC_chang_013021/martin_naive_npc20_k120_singleR.pdf", width= 24, height= 6, units= "in")


# save RDS
saveRDS(s, paste("saved_objects/CD_martin_qc_", date, "/martin_20120.rds", sep=""))



