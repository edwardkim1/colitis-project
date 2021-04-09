############################################################
## Continuation of codelog18.2 - subsetting B cells
## edward_kim@college.harvard.edu - Feb 2021
#############################################################

##############################
# 0 - Load librairies
##############################
library(dplyr)
library(stringr)
library(Seurat)
library(future)
library(future.apply)
library(hdf5r)
library(parallel)
library(ggplot2)
library(gridExtra)
library(grid)
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
options(ggrepel.max.overlaps = Inf)

############################## 
# 2 - Source file 
##############################
source("scripts/Stat111functions.R")
source("scripts/new_clustering_functions.R")

####################################
# 5 - Subset B cells
####################################
s <- readRDS("saved_objects/UC_chang_021321/uccpirec-allCD45.rds")
# B cells vizualization feature plots
DefaultAssay(s) <- "RNA"
FeaturePlot(s, c("IGHA1","IGHG1"))
ggsave("figures/UC_chang_021321/uccpirec20100_IGHA1-G1.pdf", width= 16, height= 8, units= "in")

DefaultAssay(s) <- "RNA"
FeaturePlot(s, c("CD19","CD79A"))
ggsave("figures/UC_chang_021321/uccpirec20100_CD19_CD79A.pdf", width= 16, height= 8, units= "in")
# subsetting
s1 <- s
s <- subset(s1, idents=c("1","2","4","8","14"))

# begin integration
s.list <- SplitObject(s, split.by="patient.id")
##### add this part for re-evaluating the variable features
s.list <- lapply(X = s.list, FUN = function(x) {
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# variable features
s.features <- SelectIntegrationFeatures(object.list = s.list, nfeatures = 2000)
'%ni%' = Negate('%in%')
trvgenes <- s.features[grepl(x=s.features, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- s.features[grepl(x=s.features, pattern = "^IG.V")]
s.features <- s.features[s.features %ni% c(trvgenes,igvgenes)]
# Finding the integration anchors and integrating
s.anchors <- FindIntegrationAnchors(object.list = s.list, 
        normalization.method= "LogNormalize",
        #reference = 1:16,
        anchor.features= s.features,
        dims = 1:30)
s <- IntegrateData(anchorset = s.anchors, dims=1:30)

DefaultAssay(s) <- "integrated"
# Do the rest; no need to regress out percent.mt because those genes are not in the list
s <- ScaleData(s) %>% RunPCA()
# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 20)
ggsave("figures/UC_chang_021321/uccpiB2_elbowplot.pdf")

s <- FindNeighbors(s, dims = 1:15, k.param = 100)  
# names(s@graphs)
s <- FindClusters(s,resolution = 0.5, graph.name ="integrated_snn") 
s <- RunUMAP(s,min.dist = 0.5, graph="integrated_snn")

#visualize
p <- DimPlot(s, reduction = "pca",label=T)
ggsave("figures/UC_chang_021321/uccpiB2_PCA.pdf")
p <- DimPlot(s, reduction="umap", label=T)
ggsave("figures/UC_chang_021321/uccpiB2_UMAP.pdf")
p <- DimPlot(s, reduction="umap", group.by="patient.id" , label=F)
ggsave("figures/UC_chang_021321/uccpiB2_persample.pdf", width=14, height=10, units= "in")

# Differential Expression
DefaultAssay(s) <- "RNA"
out <- DE_heatmap(s)
ggsave("figures/UC_chang_021321/uccpiB2_DE_avgexp.pdf",out$plot, width=7, height=9, units= "in")
DefaultAssay(s) <- "integrated"
out <- DE_heatmap(s)
ggsave("figures/UC_chang_021321/uccpiB2_DE_avgexp_integrated.pdf",out$plot, width=7, height=9, units= "in")

####################################
# 5 - Feature Plots and Vln Plots
####################################
# B cells vizualization feature plots
DefaultAssay(s) <- "RNA"
FeaturePlot(s, c("IGHA1","IGHA2", "IGHD", "IGHG1", "IGHG2", "CD79A"),ncol=3)
ggsave("figures/UC_chang_021321/uccpiB2_IGfeatures.pdf", width= 16, height= 8, units= "in")
# Feature plots from the UC Chang study
DefaultAssay(s) <- "RNA"
FeaturePlot(s, c("PRDM1","XBP1", "CD79A","IGHM", "CD19", "CD3E"),ncol=3)
ggsave("figures/UC_chang_021321/uccpiB2_Bcellfeatures.pdf", width= 16, height= 8, units= "in")
# Vln plots of IGHM expression
p <- VlnPlot(s,features = "IGHM", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_021321/uccpiB2_VlnPlot_IGHM_by_cluster.pdf", width= 7, height= 7, units= "in")
# Vln plots of IGHM expression
p <- VlnPlot(s,features = "CD27", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_021321/uccpiB2_VlnPlot_CD27_by_cluster.pdf", width= 7, height= 7, units= "in")
# Vln plots of PDCD1LG2 expression
p <- VlnPlot(s,features = "PDCD1LG2", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_021321/uccpiB2_VlnPlot_PDCD1LG2_by_cluster.pdf", width= 7, height= 7, units= "in")
# Vln plots of CD45 expression
p <- VlnPlot(s,features = "PTPRC", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_021321/uccpiB2_VlnPlot_CD45RB_PTPRC_by_cluster.pdf", width= 7, height= 7, units= "in")
# Vln plots of CD20 (MS4A1) expression
p <- VlnPlot(s,features = "MS4A1", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_021321/uccpiB2_VlnPlot_CD20_PTPRC_by_cluster.pdf", width= 7, height= 7, units= "in")

####################################
# 5 - CPI Labels
####################################
temp <- s[,s$orig.ident=="CPI_colitis"]
p <- DimPlot(temp, reduction="umap", group.by="anno2" , label=F)
ggsave("figures/UC_chang_021321/uccpiB2_labels_CPIanno2.pdf", width=9, height=7, units= "in")
# highlight certain cells
Idents(temp) <- temp$anno2
p <- DimPlot(temp, reduction="umap", group.by="anno2",cells.highlight=WhichCells(temp, idents = "16: Cycling T"), label=F)
ggsave("figures/UC_chang_021321/uccpiB2_labels_highlight_cyclingT.pdf", width=9, height=7, units= "in")
p <- DimPlot(temp, reduction="umap", group.by="anno2",cells.highlight=WhichCells(temp, idents = "4: Mast"), label=F)
ggsave("figures/UC_chang_021321/uccpiB2_labels_highlight_Mast.pdf", width=9, height=7, units= "in")

####################################
# 5 - Single R: Monaco
####################################
library(SingleR)
library(BiocParallel)
ref <- MonacoImmuneData()
pred.Monaco <- SingleR(method="cluster",
        test=s[['integrated']]@data,
        clusters= s$seurat_clusters, 
        ref=ref, 
        labels=ref$label.fine, 
        BPPARAM=MulticoreParam(10))
s$SingleR.Monaco1 <- s$seurat_clusters
levels(s$SingleR.Monaco1) <- pred.Monaco$labels

ref <- MonacoImmuneData()
pred.Monaco <- SingleR(method="cluster",
        test=s[['RNA']]@data,
        clusters= s$seurat_clusters, 
        ref=ref, 
        labels=ref$label.fine, 
        BPPARAM=MulticoreParam(10))
s$SingleR.Monaco2 <- s$seurat_clusters
levels(s$SingleR.Monaco2) <- pred.Monaco$labels


require(ggpubr)
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco1") + labs(title="Monaco (Integration-corrected)")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco2") + labs(title="Monaco (RNA)")
ggarrange(p2,p3, ncol=2, nrow=1)
ggsave("figures/UC_chang_021321/uccpiB2_anchor_labeling.pdf", width= 18, height= 6, units= "in")


####################################
# 5 - Cluster 7 IgG vs IgA Coexpression
####################################
# Correlation between IgHA expression and IgHB total
DefaultAssay(s) <- "RNA"
FeatureScatter(s,feature1="IGHA1", feature2="IGHG1")
ggsave(paste("figures/UC_chang_021321/uccpiB2_IGHA1xIGHG1_total.pdf", sep=""), width= 7, height= 7, units= "in")
# Correlation between IgHA expression and IgHB in cluster 7
DefaultAssay(s) <- "RNA"
FeatureScatter(s[,s$seurat_clusters=="7"],feature1="IGHA1", feature2="IGHG1")
ggsave(paste("figures/UC_chang_021321/uccpiB2_IGHA1xIGHG1_cluster7.pdf", sep=""), width= 7, height= 7, units= "in")

####################################
# Adding Annotations
####################################
bcell.anno <- c("(0,2,4,8) IgA Plasma B", 
        "(1) MHCII+ Memory B",
        "(0,2,4,8) IgA Plasma B",
        "(3) Naive",
        "(0,2,4,8) IgA Plasma B",
        "(5) Early Plasma B",
        "(6) IL32+ T cell / Doublet",
        "(7) IgG Plasma B",
        "(0,2,4,8) IgA Plasma B",
        "(9) Cycling B",
        "(10) MHCII+ Memory")
s$bcell.anno <- s$seurat_clusters
levels(s$bcell.anno) <- bcell.anno

# save RDS
saveRDS(s, "saved_objects/UC_chang_021321/uccpirec-Bcell_subset.rds")

####################################
# 5 - DA figure
####################################
p <- plot.propbar.errors(annotations = s$bcell.anno,
        sample.id = s$sample.id,
        conditions = s$colitis3,
        clusters= names(table(s$bcell.anno)),
        group.names= c("Control","UC-Control","No-Colitis","Colitis","UC"))
ggsave("figures/UC_chang_021321/uccpiB2_barplot.pdf", width= 12, height= 7, units= "in")

#### Wilcoxon Test (Unpaired) ####
wilcoxon.p.values <- function(clusters, df.prop) {
        p.values <- unlist(lapply(clusters,function(x) {wilcox.test(proportion ~ status, data = subset(df.prop, cluster==x), paired = FALSE)$p.value}))
        names(p.values) <- unique(df.prop$cluster)
        return(p.values)
}
# get p-values using wilcoxon test
df.s <- wilcoxon.df(s, s$SingleR.Luoma2)
s.pvals <- wilcoxon.p.values(names(table(s$SingleR.Luoma2)),subset(df.s, status=="Colitis" | status=="UC"))
