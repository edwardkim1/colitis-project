##############################################
## cluster labeling: onlyTcells nPC 15 k 50 ##
##############################################
library(dplyr)
library(SingleR)
library(BiocParallel)
source("scripts/new_clustering_functions.R")

s <- readRDS("saved_objects/onlyTcells_npc15_k50.rds")

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

saveRDS(list(pred.CPI,pred.Martin,pred.Monaco),"saved_objects/label_info_for_onlyTcells1550.rds")

saveRDS(s,"saved_objects/onlyTcells1550_labeled.rds")

require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin")
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco")
ggarrange(p1,p2,p3, ncol=3, nrow=1)
ggsave("figures/onlyTcells_npc15_k50_labeling.pdf", width= 24, height= 6, units= "in")
# Diagnostic plot



