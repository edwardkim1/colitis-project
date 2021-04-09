##############################
# 0 - Load librairies
##############################
library(dplyr)
library(stringr)
library(Seurat)
library(future)
#library(future.apply)
#library(hdf5r)
#library(parallel)
library(ggplot2)
#library(gridExtra)
#library(grid)
#library(fgsea)
#library(harmony)
#library(scran)
#library(pheatmap)

##############################
# 1 - Definitions & Settings
##############################
'%ni%' = Negate('%in%')
plan(multicore, workers = 10)
options(future.globals.maxSize = 20000 * 1024^2)
options(ggrepel.max.overlaps = Inf)
key.colitis.colors = c("#a3d2ca", "#ffcda3", "#5eaaa8","#056676","#db6400")

############################## 
# 2 - Source file 
##############################
source("scripts/Stat111functions.R")
source("scripts/new_clustering_functions.R")


############################## 
# 2 - Integration (failed)
# ##############################
# s <- readRDS("saved_objects/UC_chang_021321/uccpirec-allCD45.rds")
# ibdpos <- readRDS("saved_objects/CD_luoma_qc_011621/p1089pos_3.RDS")

# # meta data
# ibdpos$sample.id <- ibdpos$orig.ident
# ibdpos$orig.ident <- "MGH_CD"
# ibdpos$patient.id <- "MGH CD"
# ibdpos$colitis <- "CD"
# ibdpos$colitis2 <- "CD"
# ibdpos$colitis3 <- "CD"

# DefaultAssay(s) <- "integrated"
# s.list <- list(s,ibdpos)
# # variable features
# s.features <- SelectIntegrationFeatures(object.list = s.list, nfeatures = 3000)
# '%ni%' = Negate('%in%')
# trvgenes <- s.features[grepl(x=s.features, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
# igvgenes <- s.features[grepl(x=s.features, pattern = "^IG.V")]
# s.features <- s.features[s.features %ni% c(trvgenes,igvgenes)]
# # Finding the integration anchors and integrating
# s.anchors <- FindIntegrationAnchors(object.list = s.list, 
#         normalization.method= "LogNormalize",
#         reference = 1,
#         anchor.features= s.features,
#         dims = 1:50)
# s <- IntegrateData(anchorset = s.anchors, dims=1:50, new.assay.name="integrated2")

# DefaultAssay(s) <- "integrated2"
# # Do the rest; no need to regress out percent.mt because those genes are not in the list
# s <- ScaleData(s) %>% RunPCA()
# # PCA Elbow Plot Figure
# p <- ElbowPlot(s, ndims= 20)
# ggsave("figures/MGHCD_031621/cduccpi_elbowplot.pdf")

# s <- FindNeighbors(s, dims = 1:15, k.param = 100)  
# # names(s@graphs)
# s <- FindClusters(s,resolution = 0.5, graph.name ="integrated_snn") 
# #s <- RunUMAP(s,min.dist = 0.4, graph="integrated_snn")
# s <- RunUMAP(s, dims=1:25)

# #visualize
# p <- DimPlot(s, reduction = "pca",label=T)
# ggsave("figures/MGHCD_031621/cduccpi_PCA.pdf")
# p <- DimPlot(s, reduction="umap", label=T)
# ggsave("figures/MGHCD_031621/cduccpi_UMAP.pdf")
# p <- DimPlot(s, reduction="umap", group.by="patient.id" , label=F)
# ggsave("figures/MGHCD_031621/cduccpi_persample.pdf", width=14, height=10, units= "in")

# # Differential Expression
# DefaultAssay(s) <- "RNA"
# out <- DE_heatmap(s)
# ggsave("figures/MGHCD_031621/cduccpi_DE_avgexp.pdf",out$plot, width=4, height=7, units= "in")
# DefaultAssay(s) <- "integrated2"

# # save object
# saveRDS(s, "saved_objects/MGHCD_031621/cduccpi.rds")

############################## 
# 2 - Reference Mapping (all CD45)
##############################
# mapping onto just CPI cells
s <- readRDS("saved_objects/UC_chang_021321/uccpirec-allCD45.rds")
ibdpos <- readRDS("saved_objects/CD_luoma_qc_011621/p1089pos_3.RDS")

DefaultAssay(s) <- "integrated"
s <- RunUMAP(s, return.model=TRUE,dims=1:25)
p <- DimPlot(s, reduction="umap", label=T, raster=T)
ggsave("figures/MGHCD_031621/uccpi_umap25pcs.pdf" ,width= 6, height= 6)

temp.names <- names(table(s$SingleR.Luoma2))
s$SingleR.Luoma2 <- factor(s$SingleR.Luoma2, levels=temp.names[c(2,3,7,1,6,5,4)])
p <- DimPlot(s, reduction="umap", group.by="SingleR.Luoma2", label=T, raster=T) + NoLegend()
ggsave("figures/MGHCD_031621/uccpi_umap25pcs.pdf" ,width= 6, height= 6)
DimPlot(s, reduction= "umap", label=FALSE, group.by="colitis3", ncol=1, cols=key.colitis.colors)
ggsave("figures/MGHCD_031621/UCCPI_UMAP_groupby_dataset.pdf", width= 7, height= 6, units= "in")

anchors <- FindTransferAnchors(
  reference = s,
  query = ibdpos,
  reference.assay = "integrated",
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

ibdpos.postMQ3 <- MapQuery(
  anchorset = anchors,
  query = ibdpos,
  reference = s,
  refdata = s$SingleR.Luoma2,
  reference.reduction = "pca",
  reduction.model = "umap"
)

ibdpos.postMQ3$predicted.id <- ifelse(ibdpos.postMQ3$predicted.id=="lgA.plasma.B","Plasma.B", ibdpos.postMQ3$predicted.id)
temp.names <- names(table(ibdpos.postMQ3$predicted.id))
ibdpos.postMQ3$predicted.id <- factor(ibdpos.postMQ3$predicted.id, levels=temp.names[c(2,3,7,1,6,5,4)])

p1 = DimPlot(ibdpos.postMQ3, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
ggsave("figures/MGHCD_031621/ibdpos_cpirefUMAP_predicted.pdf", width= 6, height= 6, units= "in")

# violin plot of predicted id score per cluster
p <- VlnPlot(ibdpos.postMQ3,features = "predicted.id.score", pt.size= 0, group.by="predicted.id") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/MGHCD_031621/ibdpos_cpissuo_map_scoredist.pdf", width= 4, height= 4, units= "in")


# DA figures
combined.annotations <- c(as.character(ibdpos.postMQ3$predicted.id), as.character(s$SingleR.Luoma2))
combined.sample.id <- c(rep("MGH CD", ncol(ibdpos.postMQ3)), s$sample.id)
combined.condition <- c(rep("MGH CD", ncol(ibdpos.postMQ3)), as.character(s$colitis3))
key.colitis.colors.CD <- c(key.colitis.colors, "#C7D3D4FF")

p <- plot.propbar.errors(annotations = combined.annotations,
        sample.id = combined.sample.id,
        conditions = combined.condition,
        clusters= names(table(s$SingleR.Luoma2))[c(2,3,7,1,6,5,4)],
        group.names= c("CPI-Control", "UC-Control" , "CPI-No-Colitis", "CPI-Colitis", "UC", "MGH CD"))
p + scale_fill_manual(values=key.colitis.colors.CD)
ggsave("figures/MGHCD_031621/uccpi_CD45_barplot.pdf", width= 10, height= 5, units= "in")

saveRDS(ibdpos.postMQ3, "saved_objects/MGHCD_031621/MGHCD_refmap.rds")





############################## 
# 2 - Reference Mapping (T cells)
##############################
ibdposT <- readRDS("saved_objects/CD_luoma_qc_011621/p1089t20100.RDS")
sT <- readRDS("saved_objects/UC_chang_031321/uccpiTa_withUCcontrols.rds")
DefaultAssay(sT) <- "integrated"
# plotting CPI cells just to check labels
sT <- RunUMAP(sT, return.model=TRUE,dims=1:25)
p <- DimPlot(sT, reduction="umap", group.by="seurat_clusters", label=T, raster=T)
ggsave("figures/MGHCD_T_031621/uccpiT_umap25pcs.pdf")
DimPlot(sT, reduction= "umap", label=FALSE, group.by="colitis3", ncol=1, cols=key.colitis.colors)
ggsave("figures/MGHCD_T_031621/UCCPI_UMAP_groupby_dataset.pdf", width= 7, height= 6, units= "in")

# mapping
anchors <- FindTransferAnchors(
  reference = sT,
  query = ibdposT,
  reference.assay = "integrated",
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

postMQ4 <- MapQuery(
  anchorset = anchors,
  query = ibdposT,
  reference = sT,
  refdata = sT$seurat_clusters,
  reference.reduction = "pca",
  reduction.model = "umap"
)
saveRDS(postMQ4, "saved_objects/MGHCD_T_031621/MGHCD_T_refmap.rds")
postMQ4 <- readRDS("saved_objects/MGHCD_T_031621/MGHCD_T_refmap.rds")

postMQ4$predicted.id <- factor(postMQ4$predicted.id, levels=0:11)

p1 = DimPlot(postMQ4, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
ggsave("figures/MGHCD_T_031621/ibdpos_cpirefUMAP_predicted.pdf", width= 6, height= 6, units= "in")

p1 = DimPlot(postMQ4, reduction = "ref.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
ggsave("figures/MGHCD_T_031621/ibdpos_cpirefUMAP_predicted_orig.pdf", width= 6, height= 6, units= "in")


# violin plot of predicted id score per cluster
p <- VlnPlot(postMQ4,features = "predicted.id.score", pt.size= 0, group.by="predicted.id") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/MGHCD_T_031621/ibdpos_cpissuo_map_scoredist.pdf", width= 4, height= 5, units= "in")


# DA figures
postMQ4$predicted.id <- factor(postMQ4$predicted.id)
levels(postMQ4$predicted.id) <- c("CD4 Trm", "CD8 Trm","LP Trm","IEL:CD8/gd","Cytotoxic effector",
	"Naive/CM","TFH","Treg","Th1 effector","Cycling","10","MAIT")

combined.annotations <- c(as.character(postMQ4$predicted.id), as.character(sT$names))
combined.sample.id <- c(rep("MGH CD", ncol(postMQ4)), sT$sample.id)
combined.condition <- c(rep("MGH CD", ncol(postMQ4)), as.character(sT$colitis3))
key.colitis.colors.CD <- c(key.colitis.colors, "#C7D3D4FF")

p <- plot.propbar.errors(annotations = combined.annotations,
        sample.id = combined.sample.id,
        conditions = combined.condition,
        clusters= levels(postMQ4$predicted.id),
        group.names= c("CPI-Control", "UC-Control" , "CPI-No-Colitis", "CPI-Colitis", "UC", "MGH CD"))
p + scale_fill_manual(values=key.colitis.colors.CD) +
	scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
ggsave("figures/MGHCD_T_031621/uccpi_CD45_barplot.pdf", width= 10, height= 5, units= "in")








