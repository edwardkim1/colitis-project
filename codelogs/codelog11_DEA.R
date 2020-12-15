############################################################
## Differential expression by cluster
## edward_kim@college.harvard.edu - Aug 2020
#############################################################

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

##################################################
# 2 - Differential expression of previous clusters 
##################################################
s1 <- readRDS("saved_objects/martinTcells.rds")
s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
s3 <- readRDS("saved_objects/smillie_CD3.RDS")


DE_heatmap <- function(seurat.object, filename, n.cores=10) {
	require(Seurat)
	require(future)
	plan(multicore, workers= n.cores)
	s_regress.markers <- FindAllMarkers(seurat.object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
	top5 <- s_regress.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
	# change spaces in cluster names to underscores
	orig.levels <- levels(seurat.object)
	Idents(seurat.object) <- gsub(pattern = " ", replacement = "_", x = Idents(seurat.object))
	orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
	levels(seurat.object) <- orig.levels
	# store cluster averages in a Seurat object
	cluster.averages <- AverageExpression(seurat.object, return.seurat = TRUE)
	# visualize
	p <- DoHeatmap(cluster.averages,features = top5$gene, size = 3, draw.lines = FALSE, raster=F)
	ggsave(filename)
}

DE_heatmap(s2,"figures/avgexp_onlyTcells1550.pdf")
DE_heatmap(s3,"figures/avgexp_smillie1550.pdf")
DE_heatmap(s1,"figures/avgexp_martin1550.pdf")

unlist(TopFeatures(s[["pca"]], balanced = TRUE))


