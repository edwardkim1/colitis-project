############################################################
## Mapping lab IBD on the CPI dataset
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
# cluster IBD (negative and positive)
###################################
# mkdir CD_luoma_011621
# CD45-
date <- "011621"
dirname <- "p1089neg"
input.directory <- "data/CD_luoma/p1089neg-GEX-Pool7/outs/filtered_feature_bc_matrix"
qc_CD(dirname, date, "luoma", input.directory)
save_figures_CD(dirname, date, "luoma", input.directory)
# CD45+
dirname <- "p1089pos"
input.directory <- "data/CD_luoma/p1089pos-GEX-Pool6/outs/filtered_feature_bc_matrix"
qc_CD(dirname, date, "luoma", input.directory, custom.lnf.lim = 6.7)
save_figures_CD(dirname, date, "luoma", input.directory)

# load the pos and neg IBD samples
ibdpos <- readRDS("saved_objects/CD_luoma_qc_011621/p1089pos_3.RDS")
p <- DimPlot(ibdpos, reduction="umap", label=T) #+ ggtitle("p1089pos")
ggsave("figures/CD_luoma_011621/p1089pos_umap25pcs.pdf")
# Differential Expression
DE_heatmap(ibdpos,"figures/CD_luoma_011621/p1089pos_DE_avgexp.pdf")

####################################
# Subset T cells
###################################
ibdpos <- readRDS("saved_objects/CD_luoma_qc_011621/p1089pos_3.RDS")
s <- subset(ibdpos, idents=c("0","2","4","6","9"))
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
ggsave("figures/CD_luoma_011621/p1089t20100_umap.pdf")
# Differential Expression
DE_heatmap(s,"figures/CD_luoma_011621/p1089t20100_DE_avgexp.pdf")
# save RDS
saveRDS(s,"saved_objects/CD_luoma_qc_011621/p1089t20100.RDS")

####################################
# Merge with existing cdcpi cluster
###################################
# new conda session with Seurat v4:
remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
remotes::install_github("jlmelville/uwot")
remotes::install_github("mojaveazure/seurat-disk")


library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(patchwork)

# load cdcpi
s.a <- readRDS("saved_objects/cdcpi3_anchor.rds")
# load the pos and neg IBD samples
ibdpos <- readRDS("saved_objects/CD_luoma_qc_011621/p1089t20100.RDS")

# add metadata to the IBD samples
ibdpos$orig.ident <- "CD_luoma"
ibdpos$sample.id <- "p1089pos"
ibdpos$patient.id <- "p1089"
ibdpos$colitis <- "CD_inflamed"
# edit one orig.ident detail of cdcpi
s.a$orig.ident <- ifelse(s.a$orig.ident == "CD3_Tcell","CPI_colitis",s.a$orig.ident)

anchors <- FindTransferAnchors(
  reference = s.a,
  query = ibdpos,
  reference.assay = "integrated",
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

s.a <- RunUMAP(s.a, return.model=TRUE,dims=1:25)
p <- DimPlot(s.a, reduction="umap", label=T, raster=F)
ggsave("figures/CD_luoma_011621_refmap/cdcpi3_umap25pcs.pdf")

s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
s2 <- RunUMAP(s2, return.model=TRUE,dims=1:25)
p <- DimPlot(s2, reduction="umap", label=T, raster=T)
ggsave("figures/CD_luoma_011621_refmap/cpi_umap25pcs.pdf")

ibdpos.postMQ <- MapQuery(
  anchorset = anchors,
  query = ibdpos,
  reference = s.a,
  refdata = s.a$seurat_clusters,
  reference.reduction = "pca",
  reduction.model = "umap"
)

anchors <- FindTransferAnchors(
  reference = s2,
  query = ibdpos,
  reference.assay = "RNA",
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

ibdpos.postMQ2 <- MapQuery(
  anchorset = anchors,
  query = ibdpos,
  reference = s2,
  refdata = s2$seurat_clusters,
  reference.reduction = "pca",
  reduction.model = "umap"
)

#visualizing the mapping results
p1 = DimPlot(ibdpos.postMQ, reduction = "ref.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(ibdpos.postMQ2, reduction = "ref.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2
ggsave("figures/CD_luoma_011621_refmap/ibdpos_cpirefUMAP.pdf", width= 12, height= 7, units= "in")

# plotting CPI cells just to check labels
load('data/ssuo_CD3/seurat.object.RData')
Tcell <- process_ssuo_Tcells(Tcell)
p <- DimPlot(Tcell, reduction="umap", label=T, raster=T)
ggsave("figures/CD_luoma_011621_refmap/cpi_ssuo_orig_umap25pcs.pdf")

# mapping onto just CPI cells
Tcell <- RunUMAP(Tcell, return.model=TRUE,dims=1:25)
p <- DimPlot(Tcell, reduction="umap", label=T, raster=T)
ggsave("figures/CD_luoma_011621_refmap/cpi_ssuo_umap25pcs.pdf")

anchors <- FindTransferAnchors(
  reference = Tcell,
  query = ibdpos,
  reference.assay = "RNA",
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

ibdpos.postMQ3 <- MapQuery(
  anchorset = anchors,
  query = ibdpos,
  reference = Tcell,
  refdata = Tcell$seurat_clusters,
  reference.reduction = "pca",
  reduction.model = "umap"
)

p1 = DimPlot(ibdpos.postMQ3, reduction = "ref.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
ggsave("figures/CD_luoma_011621_refmap/ibdpos_cpirefUMAP.pdf", width= 7, height= 7, units= "in")

# Plot barplot comparing cluster proportions
plot.propbar(group1 = as.character(ibdpos.postMQ3$predicted.id),
		group2 = as.character(Tcell$seurat_clusters[Tcell$colitis == "Colitis"]),
		factor.levels = 0:12,
		group.names = c("MGH CD","CPI colitis"),
		filename = "figures/CD_luoma_011621_refmap/ibdpos_cpi_barplot.pdf")

# violin plot of predicted id score per cluster
ibdpos.postMQ3$predicted.id <- factor(ibdpos.postMQ3$predicted.id, levels=0:12)
p <- VlnPlot(ibdpos.postMQ3,features = "predicted.id.score", pt.size= 0, group.by="predicted.id") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_luoma_011621_refmap/ibdpos_cpissuo_map_scoredist.pdf", width= 7, height= 7, units= "in")

#merge reference and query [[got error]] ************
Tcell$id <- 'reference'
ibdpos.postMQ3$id <- 'query'
refquery <- merge(Tcell, ibdpos.postMQ3)
#refquery[["pca"]] <- merge(Tcell[["pca"]], ibdpos.postMQ3[["ref.pca"]]) # got error
#refquery <- RunUMAP(refquery, reduction = 'pca', dims = 1:50)
#DimPlot(refquery, group.by = 'id')
#ggsave("figures/CD_luoma_011621_refmap/ibdpos_cpissuo_merged_umap.pdf")


# DE per cluster
DE_volcano <- function(seurat.object, ident.1, group.by,subset.ident, FCcutoff,drawConnectors=T, title, filename) {
	require(EnhancedVolcano)
	s.markers <- FindMarkers(seurat.object, ident.1 =ident.1 ,group.by=group.by, subset.ident= subset.ident, min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
	head(s.markers, n=20) %>% rownames() -> labels
	up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
	down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
	labels <- unique(c(labels, up.labels,down.labels))
	EnhancedVolcano(s.markers,
	    lab = rownames(s.markers),
	    x = 'avg_log2FC',
	    y = 'p_val_adj',
	  #  xlim = c(-20, 20),
	    title = title,
	    pCutoff = 0.05,
	    FCcutoff = FCcutoff,
	    pointSize = 3.0,
	    labSize = 3.0,
	   selectLab = labels,
	    drawConnectors = drawConnectors,
	    widthConnectors = 0.2,
	    colConnectors = 'grey30')
	ggsave(filename, width= 10, height= 8, units= "in")
	return(s.markers)
}

Idents(refquery) <- c(as.character(Tcell$seurat_clusters),as.character(ibdpos.postMQ3$predicted.id))
library(EnhancedVolcano)
markers.out <- lapply(0:12, function(x) DE_volcano(refquery, ident.1 ="query" ,group.by= "id", subset.ident= x, FCcutoff = 0.8,
	title = paste('Cluster',x, 'MGH CD vs. CPI colitis'),
	filename= paste("figures/CD_luoma_011621_refmap/DE_cluster_",x,".pdf" ,sep="")))










