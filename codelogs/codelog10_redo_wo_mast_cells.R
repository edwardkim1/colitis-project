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
# 5 - Subsetting T cells
##############################
s <- readRDS("saved_objects/smillie_unfiltered.RDS")
Tclusters <- rownames(table(s$Cluster))[c(1:4,7:9,12,20,23)] #accidentally included mast cells, so removed
Idents(s) <- "Cluster"
s <- subset(s, idents = Tclusters)
sum(s$cells.to.remove) # 3540 cells to remove

# Plot showing the distribution of cells to remove by cluster
p <- VlnPlot(s,features = "cells.to.remove", pt.size= 0, group.by="Cluster") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/smillie_VlnPlot_cells2remove_by_cluster.pdf", width= 7, height= 7, units= "in")

# Plot showing the distribution of %mito by cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="Cluster") + theme(axis.text.x = element_text(angle = 90), axis.title.x=element_blank(), legend.position = "none")
ggsave("figures/smillie_VlnPlot_percentmt_by_cluster.pdf", width= 7, height= 7, units= "in")

##############################
# 6 - Clustering T cells
##############################
DefaultAssay(s) <- "RNA"
s <- s[,!s$cells.to.remove]
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
s <- ScaleData(s, features = all.genes, vars.to.regress=c("percent.mt")) %>% RunPCA(features=VariableFeatures(s))
s <-  cluster_umap(s, npc=15, k.param=100, resolution=0.5, min.dist=0.3)
# vizualize; comparison
p <- DimPlot(s,label=T, group.by="seurat_clusters")
q <- DimPlot(s,label=T, group.by="Cluster")
p+q
ggsave("figures/smillie_Tcell_UMAP.pdf", width= 17, height= 7, units= "in")

saveRDS(s,"saved_objects/smillie_CD3.RDS")

##############################
# 7 - Integration
##############################
#s <- readRDS("saved_objects/smillie_CD3.RDS")
s1 <- readRDS("saved_objects/martinTcells.rds")
s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
s.list <- list(CPI=s2,CD=s1,UC=s)
# variable features
s.features <- SelectIntegrationFeatures(object.list = s.list, nfeatures = 3000)
'%ni%' = Negate('%in%')
all.genes <- unique(c(rownames(s.list$CPI), rownames(s.list$CD),rownames(s.list$UC)))
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

####################################
# 8 - Adding/Editing meta data (must be done in the specified order)
####################################
library(stringr)
#### removing unapplicable meta.data columns
m.names <- colnames(s@meta.data)
columns.to.remove <- m.names[grepl("^pANN*|^DF.c*|seurat_clusters|^RNA*|cells.to.remove|Single*|percent.mito",m.names)]
for(i in columns.to.remove) {
  s[[i]] <- NULL
}
#### rename sampletype to sample.id and merge with orig.ident and Sample
s$sample.id <- s$sampletype
s$sample.id <- ifelse(s$orig.ident == "UC_smillie", s$Sample ,s$sample.id)
s$sample.id <- ifelse(grepl(s$orig.ident,pattern="^GSM*"), s$orig.ident, s$sample.id)
s$sampletype <- NULL
s$Sample <- NULL
#### mutate orig.ident
s$orig.ident <- ifelse(grepl(s$orig.ident,pattern="^GSM*"), "CD_martin", s$orig.ident)
#### merge patient.id and Subject
s$patient.id <- str_replace(s$patient.id,"Patient ","CD_P")
s$patient.id <- ifelse(s$orig.ident == "UC_smillie",str_glue("UC_{s$Subject}"), s$patient.id)
s$patient.id <- ifelse(s$orig.ident == "CD3_Tcell",str_glue("CPI_{s$sample.id}"),s$patient.id)
s$Subject <- NULL
#### update colitis and colitis2 with CD and UC info
s$colitis <- str_glue("CPI_{str_to_lower(s$colitis)}")
inflamed <-c("_1$","_3$","_6$","_7$","_10$","_12$","_14$","_16$","_18$")
toMatch <- paste(inflamed,collapse="|")
status <- ifelse(grepl(colnames(s), pattern=toMatch),"CD_inflamed","CD_non-inflamed")
s$colitis <- ifelse(s$orig.ident == "CD_martin", status ,s$colitis)
status2 <- str_replace(str_to_lower(s$Health),"healthy", "control")
s$colitis <- ifelse(s$orig.ident == "UC_smillie", str_glue("UC_{status2}"),s$colitis)
s$colitis2 <- ifelse(s$orig.ident != "CD3_Tcell", s$colitis, s$colitis2)
s$colitis <- ifelse(s$colitis == "UC_control" | s$colitis == "CPI_control", "Control", s$colitis)
s$Health <- NULL
#### Remove more meta.data columns that are unnecessary
s$S.Score <- NULL
s$G2M.Score <- NULL
s$nGene <- NULL
####################################
# 9 - Clustering
####################################
DefaultAssay(s) <- "integrated"
# Do the rest; no need to regress out percent.mt because those genes are not in the list
s <- ScaleData(s) %>% RunPCA()
# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 40)
ggsave("figures/cduccpi_anchor_elbowplot.pdf")

s <- FindNeighbors(s, dims = 1:15, k.param = 50) 
# names(s@graphs)
s <- FindClusters(s,resolution = 0.5, graph.name ="integrated_snn") 
s <- RunUMAP(s,min.dist = 0.3, graph="integrated_snn")

#visualize
p <- DimPlot(s, label=T)
ggsave("figures/cduccpi_anchor_PCA.pdf")

p <- DimPlot(s, reduction="umap", label=T)
ggsave("figures/cduccpi_anchor_UMAP.pdf")

## Viz 1: UMAP - group by colitis and colitis2
p <- DimPlot(s, reduction = "umap", label=F, group.by="colitis")
q <- DimPlot(s, reduction = "umap", label=F, group.by="colitis2")
p+q
ggsave("figures/cduccpi_anchor_group_colitis.pdf",width= 16, height= 6, units= "in")

## Viz 4a: cluster proportions
table(s$seurat_clusters,s$colitis2) %>% unclass() %>% prop.table()
data.frame(cluster =s$seurat_clusters, condition= s$colitis) %>%
	ggplot() + geom_bar(
		mapping = aes(x= cluster, fill= condition),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/cduccpi_anchor_cluster_prop.pdf")

## Viz 4c: cluster proportions by condition colitis 2
reordered.conditions <- factor(s$colitis2, levels = c("Colitis (aPD1)" ,"No-Colitis (aPD1)" , "Colitis (combo)" ,  "No-Colitis (combo)", "CD_inflamed","CD_non-inflamed", "UC_inflamed", "UC_non-inflamed", "CPI_control","UC_control"))
data.frame(cluster =s$seurat_clusters, condition= reordered.conditions) %>%
	ggplot() + geom_bar(
		mapping = aes(x= condition, fill= cluster),
		position = "fill"
	) + labs(y="proportion")+ theme_classic() + theme(axis.text.x = element_text(angle = 50,hjust=1))
ggsave("figures/cduccpi_anchor_cluster_prop_colitis2_by_condition.pdf")

####################################
# 10 - Cluster Annotations
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

#s3 <- readRDS("saved_objects/smillie_filtered.RDS")
pred.Smillie <- SingleR(method="cluster",
	test=s[['integrated']]@data,
	clusters=s$seurat_clusters, 
	ref=s3[['RNA']]@data, 
	labels=s3$Cluster,
	de.method="wilcox", 
	BPPARAM=MulticoreParam(10))
s$SingleR.Smillie <- s$seurat_clusters
levels(s$SingleR.Smillie) <- pred.Smillie$labels

saveRDS(list(pred.CPI,pred.Martin,pred.Monaco,pred.Smillie),"saved_objects/label_info_for_cduccpi1550.rds")


require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin") + labs(title="Martin reference")
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma") + labs(title="Luoma-Suo reference")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco") + labs(title="Monaco reference")
p4 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Smillie") + labs(title="Smillie reference")
ggarrange(p1,p2,p3,p4, ncol=2, nrow=2)
ggsave("figures/cduccpi_anchor_labeling.pdf", width= 16, height= 12, units= "in")

saveRDS(s, "saved_objects/cduccpi_anchor.rds")

##################################################
# 10 - Clustering check with reverse annotating
##################################################
## Get annotations from Luoma-Suo cells
require(stringr)
require(future.apply)
# Finding the correct mapping: turns out it is the "-1_1, -2_1, etc." endings
# that match the Suo Tcells
colnames(s2) %>% str_sub(17) %>% table()
colnames(Tcell) %>% str_sub(17) %>% table()
colnames(s2)[as.logical(cells)] %>% str_sub(17) %>% table()
Tcell$names %>% class()
# executing the addtion
x <- Tcell$names
names(x) <- str_glue("{names(x)}_1")
s2$old.Tcell.annotations <- x
# graph showing clusters: onlyTcells (CPI)
p <- DimPlot(s2,reduction="umap",group.by="old.Tcell.clusters")
ggsave("figures/apply_old_labels_to_onlyTcells1550.pdf")
# Repeating the search: turns out it is the "-1_1, -2_1, etc." endings
# that match the Suo Tcells
colnames(s)[which(s$orig.ident=="CD3_Tcell")] %>% str_sub(17) %>% table()
#executing the addition
x <- Tcell$names
names(x) <- str_glue("{names(x)}_1")
s$old.CPI.annotations <- x
# graph showing clusters: cduccpi (CPI+CD+UC)
p <- DimPlot(s,reduction="umap",group.by="old.CPI.annotations")
ggsave("figures/apply_old_CPI_labels_to_cduccpi1550.pdf", width= 8, height= 7, units= "in" )
p <- DimPlot(s,reduction="umap",group.by="old.CPI.annotations", cells=!is.na(s$old.CPI.annotations))
ggsave("figures/apply_old_CPI_labels_to_cduccpi1550_noNA.pdf", width= 8, height= 7, units= "in" )

# Get annotations from Smillie cells
s$old.UC.annotations <- s3$Cluster
# graph showing clusters: cduccpi (CPI+CD+UC)
p <- DimPlot(s,reduction="umap",group.by="old.UC.annotations")
ggsave("figures/apply_old_UC_labels_to_cduccpi1550.pdf", width= 8, height= 7, units= "in" )
p <- DimPlot(s,reduction="umap",group.by="old.UC.annotations", cells=!is.na(s$old.UC.annotations),order=sample(colnames(s)[which(!is.na(s$old.UC.annotations))]))
ggsave("figures/apply_old_UC_labels_to_cduccpi1550_noNA.pdf", width= 8, height= 7, units= "in" )






