############################################################
## Redo UCCPI with UC controls
## edward_kim@college.harvard.edu - Feb 2021
#############################################################

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
plan("multicore", workers = 10)
options(future.globals.maxSize = 20000 * 1024^2)
options(ggrepel.max.overlaps = Inf)
key.colitis.colors = c("#a3d2ca", "#ffcda3", "#5eaaa8","#056676","#db6400")

############################## 
# 2 - Source file 
##############################
source("scripts/Stat111functions.R")
source("scripts/new_clustering_functions.R")

############################## 
# 3 - Load datasets
##############################
s <- readRDS("saved_objects/UC_chang_031321/ucrec20120.rds")
# Featureplot to identify T cells
FeaturePlot(s, c("CD3E","CD8A","CD4"),ncol=3)
ggsave("figures/UC_chang_031321/ucrec_featureplot_Tgenes.pdf", width= 16, height= 5, units= "in")
# subset T cells
s.t <- subset(s, idents=c(0,3,5,6,7))

load('data/ssuo_CD3/seurat.object.RData')
Tcell <- process_ssuo_Tcells(Tcell)
Tcell$sample.id <- Tcell$sampletype
Tcell$patient.id <- Tcell$sampletype
###########################################
# 3 - adding metadata for the CPI Tcells
###########################################
library(stringr)
# add patient.id and sample.id for CPI
s$orig.ident <- ifelse(s$orig.ident=="CD3_Tcell", "CPI_colitis", s$orig.ident)
s$sample.id <- ifelse(s$orig.ident=="CPI_colitis", s$sampletype, s$sample.id)
s$patient.id <- ifelse(s$orig.ident=="CPI_colitis", s$sampletype, s$patient.id)
s$sampletype <- NULL

### condense sample.id for CPI
# for CPI colitis
Tcell$sample.id[Tcell$colitis=="Colitis"] %>% 
table() %>% 
rownames() -> CPI_C.samples # "DFCI551" "p1032"   "p1059"   "p369"    "p728"    "p806"    "p896" "p980"
replace <- paste("CT",1:length(CPI_C.samples),sep="")
for(i in 1:length(CPI_C.samples)) {
        Tcell$sample.id <- ifelse(Tcell$patient.id==CPI_C.samples[i],replace[i], Tcell$sample.id)
}
# for CPI no colitis
Tcell$sample.id[Tcell$colitis=="No-Colitis"] %>% 
table() %>% 
rownames() -> CPI_NC.samples #"p57"    "p590"   "p667"   "p685"   "p931-1" "p949-1"
replace <- paste("NC",1:length(CPI_NC.samples),sep="")
for(i in 1:length(CPI_NC.samples)) {
        Tcell$sample.id <- ifelse(Tcell$patient.id==CPI_NC.samples[i],replace[i], Tcell$sample.id)
}
# for Controls
Tcell$sample.id[Tcell$colitis=="Control"] %>% 
table() %>% 
rownames() -> Control.samples #"p031919" "p062519" "p1001"   "p1069"   "p1073"   "p1074"   "p1075" "p1076"
replace <- paste("C",1:length(Control.samples),sep="")
for(i in 1:length(Control.samples)) {
        Tcell$sample.id <- ifelse(Tcell$patient.id==Control.samples[i],replace[i], Tcell$sample.id)
}

###########################################
# 3 - begin integration
###########################################
# begin integration
s.list <- SplitObject(s.t, split.by="sample.id")
##### add this part for re-evaluating the variable features
s.list <- lapply(X = s.list, FUN = function(x) {
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

s.list <- c(Tcell,s.list)
# variable features
s.features <- SelectIntegrationFeatures(object.list = s.list, nfeatures = 3000)
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
ggsave("figures/UC_chang_031321/uccpirec20100_elbowplot.pdf")

s <- FindNeighbors(s, dims = 1:15, k.param = 100)  
# names(s@graphs)
s <- FindClusters(s,resolution = 0.5, graph.name ="integrated_snn") 
#s <- RunUMAP(s,min.dist = 0.4, graph="integrated_snn")
s <- RunUMAP(s, dims=1:25)

#visualize
p <- DimPlot(s, reduction = "pca",label=T)
ggsave("figures/UC_chang_031321/uccpirec20100_PCA.pdf")
p <- DimPlot(s, reduction="umap", label=T)
ggsave("figures/UC_chang_031321/uccpirec20100_UMAP.pdf")
p <- DimPlot(s, reduction="umap", group.by="patient.id" , label=F)
ggsave("figures/UC_chang_031321/uccpirec20100_persample.pdf", width=14, height=10, units= "in")

# Differential Expression
DefaultAssay(s) <- "RNA"
out <- DE_heatmap(s)
ggsave("figures/UC_chang_031321/uccpirec20100_DE_avgexp.pdf",out$plot, width=4, height=7, units= "in")
DefaultAssay(s) <- "integrated"

# save object
saveRDS(s, "saved_objects/UC_chang_031321/uccpiTa_withUCcontrols.rds")

library(Seurat)
setwd("1.projects/colitis/")
s <- readRDS("saved_objects/UC_chang_031321/uccpiTa_withUCcontrols.rds")



######################################################
# 3 - QC metrics per cluster
######################################################
s$orig.ident <- ifelse(s$orig.ident == "CD3_Tcell","CPI_colitis",s$orig.ident)
DefaultAssay(s) <- "RNA"
s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")

######################################################
# 3 - CD3E and CD4/CD8 and other QC metrics per cluster
######################################################
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_031321/uccpirec20100_VlnPlot_CD3E_by_cluster.pdf", width= 4, height= 4, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/UC_chang_031321/uccpirec20100_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_031321/uccpirec20100_VlnPlot_percentmito_by_cluster.pdf", width= 4, height= 4, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_031321/uccpirec20100_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 4, height= 4, units= "in")
# feature plots the IGKC and JCHAIN
FeaturePlot(s, features= c('IGKC','JCHAIN'), reduction="umap")
ggsave("figures/UC_chang_031321/uccpirec20100_IGKC-JCHAIN_featureplot.pdf", width= 7, height= 4, units= "in")
# violin plot of IGKC per cluster
p <- VlnPlot(s,features = "IGKC", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave(paste("figures/UC_chang_031321/uccpirec20100_VlnPlot_IGKC_by_cluster.pdf", sep=""), width= 4, height= 4, units= "in")
# distribution by sample
data.frame(cluster = s$seurat_clusters, patient= s$sample.id) %>%
        ggplot() + geom_bar(
                mapping = aes(x= patient, fill= cluster),
                position = "fill"
        ) + labs(y="proportion")+ theme_classic() + coord_flip()
ggsave("figures/UC_chang_031321/uccpirec20100_cluster_prop_per_sample.pdf",width= 4, height= 6, units= "in")
# visualize only previous cluster annotations from CPI dataset
s2 <- s[,s$orig.ident=="CPI_colitis"]
s2$names <- as.character(Tcell$names)
Idents(s2) <- s2$names
p <- DimPlot(s2, reduction="umap", label=T)
ggsave("figures/UC_chang_031321/uccpirec20100_UMAP_CPI_reference.pdf",width= 7, height= 6, units= "in")


######################
# DA Figure
######################
# s$colitis3 <- s$colitis
# s$colitis3 <- ifelse(s$colitis=="Control" & s$orig.ident=="UC_chang", "UC-Control",s$colitis3)
# s$colitis3 <- ifelse(s$colitis3== "Colitis", "CPI-Colitis", s$colitis3)
# s$colitis3 <- ifelse(s$colitis3== "Control", "CPI-Control", s$colitis3)
# s$colitis3 <- ifelse(s$colitis3== "No-Colitis", "CPI-No-Colitis", s$colitis3)
# s$colitis3 <- factor(s$colitis3, levels= c("CPI-Control", "UC-Control" , "CPI-No-Colitis", "CPI-Colitis", "UC"))
# s$orig.ident <- ifelse(s$orig.ident== "UC_chang", "UC_Boland", s$orig.ident)
# saveRDS(s, "saved_objects/UC_chang_031321/uccpiTa_withUCcontrols.rds")

s$names <- factor(s$seurat_clusters)
levels(s$names) <- c("CD4 Trm", "CD8 Trm","LP Trm","IEL:CD8/gd","Cytotoxic effector",
	"Naive/CM","TFH","Treg","Th1 effector","Cycling","10","MAIT")

p <- plot.propbar.errors(annotations = s$names,
        sample.id = s$patient.id,
        conditions = s$colitis3,
        clusters= levels(s$names),
        group.names= c("CPI-Control", "UC-Control" , "CPI-No-Colitis", "CPI-Colitis", "UC"))
p + scale_fill_manual(values=key.colitis.colors) + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1)) + 
	scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
ggsave("figures/UC_chang_031321/uccpi_Tcells_barplot.pdf", width= 10, height= 5, units= "in")


# DA analysis with qlglm
s2 <- subset(s, colitis3 == "CPI-Colitis" | colitis3 == "UC")
s2$sample.id <- droplevels(s2$sample.id)
s2$UCoverCPI <- s2$colitis3 == "UC"
design <- function(y) model.matrix(~factor(UCoverCPI),y)
out <- DA_analysis(s2, design, title = "UC_chang_031321/uccpiT_DA", annotations=s2$names)
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)


# DA analysis with qlglm - controls
s2 <- subset(s, colitis3 == "CPI-Control" | colitis3 == "UC-Control")
s2$sample.id <- droplevels(s2$sample.id)
s2$UCoverCPI <- s2$colitis3 == "UC-Control"
design <- function(y) model.matrix(~factor(UCoverCPI),y)
out <- DA_analysis(s2, design, title = "UC_chang_031321/uccpiT_controls_DA", annotations=s2$names)
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)







