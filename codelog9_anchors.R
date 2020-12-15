##############################
# 0 - Load librairies
##############################
library(dplyr)
library(ggplot2)
library(Seurat)
library(future)
library(future.apply)
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

############################## 
# 3 - Integration
##############################
s1 <- readRDS("saved_objects/martinTcells.rds")
s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
s.list <- list(CPI=s2,CD=s1)
# variable features
s.features <- SelectIntegrationFeatures(object.list = s.list, nfeatures = 2000)
'%ni%' = Negate('%in%')
all.genes <- unique(c(rownames(s.list$CPI), rownames(s.list$CD)))
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
s.features <- s.features[s.features %ni% c(trvgenes,igvgenes)]
# Finding the integration anchors and integrating
reference.number <- 1
s.anchors <- FindIntegrationAnchors(object.list = s.list, 
	reference = reference.number, 
	normalization.method= "LogNormalize",
	anchor.features= s.features,
	dims = 1:30)
s <- IntegrateData(anchorset = s.anchors)

#### Observation
# length(unique(c(rownames(s.list$CPI), rownames(s.list$CD))))
# > 24385
# is different from
# 26385 x 119063 integrated matrix (26385 # of features)

####################################
# 4 - Adding/Editing meta data
####################################
# add inflamed id for Crohns
inflammed <-c("_1$","_3$","_6$","_7$","_10$","_12$","_14$","_16$","_18$")
toMatch <- paste(inflammed,collapse="|")
status <- ifelse(grepl(colnames(s), pattern=toMatch),"CD_inflamed","CD_uninflamed")
s$colitis <- ifelse(is.na(s$colitis), status ,s$colitis)
s$colitis <- ifelse(s$colitis =="Colitis", "CPI_colitis", s$colitis)
s$colitis <- ifelse(s$colitis =="No-Colitis", "CPI_no-colitis", s$colitis)
s$colitis2 <- ifelse(is.na(s$colitis2), status ,s$colitis2)
# rename patient.id to CD.patient.id
s$CD.patient.id <- s$patient.id
s$patient.id <- NULL
# rename sampletype to sample.id and merge with orig.ident material
s$sample.id <- s$sampletype
s$sampletype <- NULL
s$sample.id <- ifelse(is.na(s$sample.id), s$orig.ident ,s$sample.id)
table(s$sample.id)
# removing unapplicable meta.data columns
m.names <- colnames(s@meta.data)
columns.to.remove <- m.names[grepl("^pANN*|^DF.c*|seurat_clusters|^RNA*|cells.to.remove|Single*|percent.mito",m.names)]
for(i in columns.to.remove) {
  s[[i]] <- NULL
}
#26385x119063 to 23202x119063
x <- s[["RNA"]]@counts
nnz_row <- tabulate(x@i + 1)
keep <- rownames(x)[nnz_row>10]
s <- subset(s, features = keep)

####################################
# 5 - Clustering
####################################
DefaultAssay(s) <- "integrated"
# Do the rest; no need to regress out percent.mt because those genes are not in the list
s <- ScaleData(s) %>% RunPCA()
# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 30)
ggsave("figures/cdcpi_anchor_elbowplot.pdf")

s <- FindNeighbors(s, dims = 1:15, k.param = 50) 
# names(s@graphs)
s <- FindClusters(s,resolution = 0.5, graph.name ="integrated_snn") 
s <- RunUMAP(s,min.dist = 0.3, graph="integrated_snn")

#visualize
p <- DimPlot(s, label=T)
ggsave("figures/cdcpi_anchor_PCA.pdf")

p <- DimPlot(s, reduction="umap", label=T)
ggsave("figures/cdcpi_anchor_UMAP.pdf")

## Viz 1: UMAP - group by colitis and colitis2
p <- DimPlot(s, reduction = "umap", label=F, group.by="colitis")
q <- DimPlot(s, reduction = "umap", label=F, group.by="colitis2")
p+q
ggsave("figures/cdcpi_anchor_group_colitis.pdf",width= 16, height= 6, units= "in")

## Viz 4a: cluster proportions
table(s$seurat_clusters,s$colitis2) %>% unclass() %>% prop.table()
data.frame(cluster =s$seurat_clusters, condition= s$colitis) %>%
	ggplot() + geom_bar(
		mapping = aes(x= cluster, fill= condition),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/cdcpi_anchor_cluster_prop.pdf")

## Viz 4c: cluster proportions by condition colitis 2
reordered.conditions <- factor(s$colitis2, levels = c("Colitis (aPD1)" ,"No-Colitis (aPD1)" , "Colitis (combo)" ,  "No-Colitis (combo)", "Control", "CD_inflamed","CD_uninflamed"))
data.frame(cluster =s$seurat_clusters, condition= reordered.conditions) %>%
	ggplot() + geom_bar(
		mapping = aes(x= condition, fill= cluster),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/cdcpi_anchor_cluster_prop_colitis2_by_condition.pdf")

####################################
# 5 - Cluster Annotations
####################################
load('data/ssuo_CD3/seurat.object.RData')
Tcell <- process_ssuo_Tcells(Tcell)
martin <- read.csv(file="data/martin_categories.csv", header=TRUE,sep=",") %>% process_martin()

require(SingleR)
require(BiocParallel)
pred.CPI <- SingleR(method='cluster',
	test=as.matrix(s[['integrated']]@data), 
	clusters= s$seurat_clusters,
	ref=as.matrix(Tcell[['RNA']]@data), 
	labels=Tcell$names,
	de.method="wilcox",
	BPPARAM=MulticoreParam(10))
s$SingleR.Luoma <- s$seurat_clusters
levels(s$SingleR.Luoma) <- pred.CPI$labels

pred.Martin <- SingleR(method='cluster',  
	test=s[['integrated']]@data, 
	clusters= s$seurat_clusters,
	ref=martin, 
	labels=colnames(martin),
	BPPARAM=MulticoreParam(10))
s$SingleR.Martin <- s$seurat_clusters
levels(s$SingleR.Martin) <- pred.Martin$labels

ref <- MonacoImmuneData()
pred.Monaco <- SingleR(method="cluster",
	test=s[['integrated']]@data,
	clusters= s$seurat_clusters, 
	ref=ref, 
	labels=ref$label.fine, 
	BPPARAM=MulticoreParam(10))
s$SingleR.Monaco <- s$seurat_clusters
levels(s$SingleR.Monaco) <- pred.Monaco$labels

saveRDS(list(pred.CPI,pred.Martin,pred.Monaco),"saved_objects/label_info_for_cdcpi1550.rds")


require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin") + labs(title="Martin reference")
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma") + labs(title="Luoma-Suo reference")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco") + labs(title="Monaco reference")
ggarrange(p1,p2,p3, ncol=3, nrow=1)
ggsave("figures/cdcpi_anchor_labeling.pdf", width= 24, height= 6, units= "in")

saveRDS(s, "saved_objects/cdcpi_anchor.rds")

