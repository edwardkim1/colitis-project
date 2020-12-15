######################################################################
## 
## Analysis of all T cells with PD1
## 
## edward_kim@college.harvard.edu - June 2020
######################################################################
​
# load merged PD1 and colitis T cells (filtered again)
s <- readRDS("saved_objects/allTcells_2.RDS") 
options(future.globals.maxSize = 11000 * 1024^2)

################################
## Differential expr analysis ##
################################
plan(multicore)
s_regress.markers <- FindAllMarkers(s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
top5 <- s_regress.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
p <- DoHeatmap(s, features = top5$gene) + NoLegend()
ggsave("figures/allTcells_heatmap.pdf", width = 30, height = 20, units="in")

######################
## Cluster labeling ##
######################
# Check if the Luoma reference clusters are labeled properly:
# load('data/ssuo_CD3/seurat.object.RData')
# p <- DimPlot(Tcell, reduction = "umap", label=T, group.by="seurat_clusters")
# Tcell$names <- Tcell$seurat_clusters
# levels(Tcell$names) <- c("Type 1 cytokines Trm", "CD8 TRM", "IEL: CD8/gd", "Cytotoxic effector", "LP Trm", "Type 3 cytokines Trm", "Treg", "Naive/CM", "TFH","Th1 effector", "Cycling", "5", "MAIT")
# q <- DimPlot(Tcell, reduction = "umap", label=T, group.by="names")
# p + q
# ggsave("figures/ssuo_CD3_check_labeling.pdf", width= 12, height= 6, units= "in")
library(SingleR)
library(BiocParallel)
load('data/ssuo_CD3/seurat.object.RData')
martin <- read.csv(file="data/martin_categories.csv", header=TRUE,sep=",") %>% process_martin()
Tcell <- process_ssuo_Tcells(Tcell)
s <- label_singleR(s,save.diagnostics=T, save.name="saved_objects/allTcells_SingleR_diagnostics.RDS")
p <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin")
q <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma")
p+q
ggsave("figures/allTcells_labeling.pdf", width= 12, height= 6, units= "in")

########################
## Remove non T cells ##
## and relabel        ##
########################
slotNames(s[["RNA"]])
umap <- as.data.frame(s@reductions$umap@cell.embeddings)
ggplot(umap) + geom_histogram(aes(x=UMAP_1), binwidth=0.5)
ggsave("figures/allTcells_UMAP_1_hist.pdf")
nonT <- rownames(subset(umap, UMAP_1 < -7 | UMAP_1 > 9.5)) # get non T cells based on UMAP coords
s <- subset(s, cells=nonT, invert=T)
p <- DimPlot(s, reduction = "umap", label=T)
ggsave("figures/allTcells_UMAP_post-nonT-removal.pdf")
saveRDS(s,"saved_objects/allTcells_post-nonT-removal.RDS")

## Recluster ##
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
p <- DimPlot(s, reduction = "umap", label=T)
ggsave("figures/onlyTcells_UMAP_recluster.pdf")
# Elbow plot
p <- ElbowPlot(s, ndims= 30)
ggsave("figures/onlyTcells_elbowplot.pdf")
# cluster v2
nPCs <- 15
s1 <- RunUMAP(s,dims = 1:nPCs) %>% FindNeighbors(dims = 1:nPCs) %>% FindClusters(resolution = 0.5)
p <- DimPlot(s1, reduction = "umap", label=T)
ggsave("figures/onlyTcells_UMAP_recluster2.pdf")
# cluster v3
nPCs <- 10
s1 <- RunUMAP(s,dims = 1:nPCs) %>% FindNeighbors(dims = 1:nPCs) %>% FindClusters(resolution = 0.5)
p <- DimPlot(s1, reduction = "umap", label=T)
ggsave("figures/onlyTcells_UMAP_recluster3.pdf")

###################################
## Experimenting with clustering ##
###################################
# Vary the k.param
require(ggpubr)
p1 <- cluster_umap(s, npc=15, k.param=20, resolution=0.5, min.dist=0.3) %>% DimPlot(reduction = "umap", label=T) + labs(title="k.param 20")
p2 <- cluster_umap(s, npc=15, k.param=50, resolution=0.5, min.dist=0.3) %>% DimPlot(reduction = "umap", label=T) + labs(title="k.param 50")
p3 <- cluster_umap(s, npc=15, k.param=100, resolution=0.5, min.dist=0.3) %>% DimPlot(reduction = "umap", label=T) + labs(title="k.param 100")
p4 <- cluster_umap(s, npc=15, k.param=500, resolution=0.5, min.dist=0.3) %>% DimPlot(reduction = "umap", label=T) + labs(title="k.param 500")
ggarrange(p1,p2,p3, ncol=2, nrow=2)
ggsave("figures/onlyTcells_vary_kparam.pdf", width= 12, height= 12, units= "in")
# Vary the npcs
p1 <- cluster_umap(s, npc=10, k.param=100, resolution=0.5, min.dist=0.3) %>% DimPlot(reduction = "umap", label=T) + labs(title="10 PCs")
p2 <- cluster_umap(s, npc=15, k.param=100, resolution=0.5, min.dist=0.3) %>% DimPlot(reduction = "umap", label=T) + labs(title="15 PCs")
p3 <- cluster_umap(s, npc=30, k.param=100, resolution=0.5, min.dist=0.3) %>% DimPlot(reduction = "umap", label=T) + labs(title="30 PCs")
p4 <- cluster_umap(s, npc=40, k.param=100, resolution=0.5, min.dist=0.3) %>% DimPlot(reduction = "umap", label=T) + labs(title="50 PCs")
ggarrange(p1,p2,p3,p4, ncol=2, nrow=2)
ggsave("figures/onlyTcells_vary_npcs.pdf", width= 12, height= 12, units= "in")
## Vary k.param again
p1 <- cluster_umap(s, npc=15, k.param=120, resolution=0.5, min.dist=0.3) %>% DimPlot(reduction = "umap", label=T) + labs(title="k.param 120")
p2 <- cluster_umap(s, npc=15, k.param=200, resolution=0.5, min.dist=0.3) %>% DimPlot(reduction = "umap", label=T) + labs(title="k.param 200")
ggarrange(p1,p2, ncol=2, nrow=1)
ggsave("figures/onlyTcells_vary_kparam2.pdf", width= 12, height= 6, units= "in")


##################################
## Final clustering: onlyTcells ##
##################################
s <- cluster_umap(s, npc=10, k.param=200, resolution=0.5, min.dist=0.3)
saveRDS(s,"saved_objects/onlyTcells.rds")

s <- cluster_umap(s, npc=15, k.param=50, resolution=0.5, min.dist=0.3)
saveRDS(s,"saved_objects/onlyTcells_npc15_k50.rds")

####################
## DE: onlyTcells ##
####################
s_regress.markers <- FindAllMarkers(s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
top5 <- s_regress.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
p <- DoHeatmap(s, features = top5$gene) + NoLegend()
ggsave("figures/onlyTcells_heatmap.pdf", width = 30, height = 20, units="in")
ggsave("figures/onlyTcells_heatmapv2.pdf")

################################
## DE: onlyTcells npc 15 k 50 ##
################################
s_regress.markers <- FindAllMarkers(s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
top5 <- s_regress.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
p <- DoHeatmap(s, features = top5$gene, raster=FALSE) + NoLegend()
ggsave("figures/onlyTcells1550_heatmap.pdf")
#ggsave("figures/onlyTcells1550_heatmap.pdf", width = 30, height = 20, units="in")

##################################
## cluster labeling: onlyTcells ##
##################################
library(SingleR)
library(BiocParallel)

# ref <- MonacoImmuneData()
# pred <- SingleR(test=s[['RNA']]@data, ref=ref, labels=ref$label.fine, BPPARAM=MulticoreParam()) #takes too long unfortunately

# #### Several Visualizations for single cell SingleR:
# ## Viz1: Heatmap of single cell assignments
# plotScoreHeatmap(pred)
# ggsave("figures/onlyTcells_labeling_Monaco_viz1.pdf")
# ## Viz2: Heatmap of cluster-level assignments based on number of cells
# library(pheatmap)
# tab <- table(pred$pruned.labels, s$seurat_clusters)
# pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))
# ggsave("figures/onlyTcells_labeling_Monaco_viz2.pdf")
# # conclusion: takes too long!!!

##### cpi reference
##### martin reference
load('data/ssuo_CD3/seurat.object.RData')
Tcell <- process_ssuo_Tcells(Tcell)
martin <- read.csv(file="data/martin_categories.csv", header=TRUE,sep=",") %>% process_martin()
s <- label_singleR(s,save.diagnostics=T, save.name="saved_objects/onlyTcells_SingleR_diagnostics.rds")
require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin")
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco")
ggarrange(p1,p2,p3, ncol=3, nrow=1)
ggsave("figures/onlyTcells_labeling.pdf", width= 12, height= 3, units= "in")
# Diagnostic plot






