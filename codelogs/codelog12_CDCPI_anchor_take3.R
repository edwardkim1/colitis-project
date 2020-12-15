############################################################
## Anchor redo with the CD (15100) + CPI -> cdcpi3
## edward_kim@college.harvard.edu - October 2020
#############################################################

##############################
# 0 - Load librairies
##############################
library(dplyr)
library(ggplot2)
library(Seurat)
library(future)
library(future.apply)
library(harmony)
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
s1 <- readRDS("saved_objects/CD_martin_qc_100720/martin_naive_CD3_15100.rds")
s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
s.list <- list(CPI=s2,CD=s1)
# variable features
s.features <- SelectIntegrationFeatures(object.list = s.list, nfeatures = 2000)
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

#### Observation
# length(unique(c(rownames(s.list$CPI), rownames(s.list$CD))))
# > 24385
# is different from
# 26385 x 119063 integrated matrix (26385 # of features)

####################################
# 4 - Adding/Editing meta data
####################################
library(stringr)
# rename sampletype to sample.id and merge with orig.ident material
s$sample.id <- s$sampletype
s$sampletype <- NULL
s$sample.id <- ifelse(is.na(s$sample.id), s$orig.ident ,s$sample.id)
#### mutate orig.ident and edit sample.id
s$orig.ident <- ifelse(grepl(s$orig.ident,pattern="^GSM*"), "CD_martin", s$orig.ident)
s$sample.id <- ifelse(s$orig.ident=="CD_martin", str_glue("CD_{str_sub(s$sample.id,-3)}"), s$sample.id)
s$sample.id <- ifelse(s$orig.ident=="CD3_Tcell", str_glue("CPI_{s$sample.id}"), s$sample.id)
#### add and merge patient.id
s <- add_patient_id_all(s)
s$patient.id <- str_replace(s$patient.id,"Patient ","CD_P")
s$patient.id <- ifelse(s$orig.ident == "CD3_Tcell",str_glue("CPI_{s$sample.id}"),s$patient.id)
# add inflamed id for Crohns
inflamed <-c("CD__69$","CD_122$","CD_128$","CD_138$","CD_158$","CD_181$","CD_187$","CD_193$","CD_190$","CD_196$","CD_209$")
toMatch <- paste(inflammed,collapse="|")
status <- ifelse(grepl(s$sample.id, pattern=toMatch),"CD_inflamed","CD_uninflamed")
s$colitis <- ifelse(s$orig.ident == "CD_martin", status ,s$colitis)
s$colitis <- ifelse(s$colitis =="Colitis", "CPI_colitis", s$colitis)
s$colitis <- ifelse(s$colitis =="No-Colitis", "CPI_no-colitis", s$colitis)
s$colitis2 <- ifelse(s$orig.ident == "CD_martin", status ,s$colitis2)
## Convert colitis2 and sample.id to factors
s$colitis2 <- factor(s$colitis2, levels = c("CD_inflamed","CD_uninflamed", "Colitis (aPD1)" ,"No-Colitis (aPD1)" , "Colitis (combo)" ,  "No-Colitis (combo)", "Control"))
conditions <- levels(s$colitis2)
order <- vector(mode="character", length=0)
for(i in 1:length(conditions)) {
	ids <- as.character(s$sample.id)
	order <- c(order, rownames(table(ids[s$colitis2 == conditions[i]])))
}
s$sample.id <- factor(s$sample.id, levels = order)
# removing unapplicable meta.data columns
m.names <- colnames(s@meta.data)
columns.to.remove <- m.names[grepl("^pANN*|^DF.c*|seurat_clusters|^RNA*|cells.to.remove|Single*|percent.mito",m.names)]
for(i in columns.to.remove) {
  s[[i]] <- NULL
}

####################################
# 4 - Clustering
####################################
DefaultAssay(s) <- "integrated"
# Do the rest; no need to regress out percent.mt because those genes are not in the list
s <- ScaleData(s) %>% RunPCA()
# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 30)
ggsave("figures/cdcpi3/cdcpi3_anchor_elbowplot.pdf")

s <- FindNeighbors(s, dims = 1:15, k.param = 100) 
# names(s@graphs)
s <- FindClusters(s,resolution = 0.5, graph.name ="integrated_snn") 
s <- RunUMAP(s,min.dist = 0.3, graph="integrated_snn")

#visualize
p <- DimPlot(s, label=T)
ggsave("figures/cdcpi3/cdcpi3_anchor_PCA.pdf")

p <- DimPlot(s, reduction="umap", label=T)
ggsave("figures/cdcpi3/cdcpi3_anchor_UMAP.pdf")


#################################################################
# 10 - Clustering check with comparison to within-batch clusters
#################################################################
library(pheatmap)
s<- readRDS("saved_objects/cdcpi3_anchor.rds")
s1 <- readRDS("saved_objects/CD_martin_qc_100720/martin_naive_CD3_15100.rds")
s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
# For the first batch (adding +10 for a smoother color transition
# from zero to non-zero counts for any given matrix entry).
tab <- table(paste("after", s$seurat_clusters[s$orig.ident=="CD_martin"]),
    paste("before", s1$seurat_clusters))
heatCD <- pheatmap(matrix.sort.no.diag(log10(tab+10)), cluster_row=FALSE, cluster_col=FALSE,
    main="CD comparison", silent=TRUE)

# For the second batch.
tab <- table(paste("after", s$seurat_clusters[s$orig.ident=="CD3_Tcell"]),
    paste("before", s2$seurat_clusters))
heatCPI <- pheatmap(matrix.sort.no.diag(log10(tab+10)), cluster_row=FALSE, cluster_col=FALSE,
    main="CPI comparison", silent=TRUE)

pdf("figures/cdcpi3/cdcpi3_cluster_comparison.pdf")
gridExtra::grid.arrange(heatCD[[4]], heatCPI[[4]])
dev.off()

####################################
# 5 - Basic Vizualizations
####################################
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/cdcpi3/cdcpi3_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of CD4 per cluster
p <- VlnPlot(s,features = "CD4", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/cdcpi3/cdcpi3_VlnPlot_CD4_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of CD8 per cluster
p <- VlnPlot(s,features = "CD8A", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/cdcpi3/cdcpi3_VlnPlot_CD8_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/cdcpi3/cdcpi3_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# distribution by sample
data.frame(cluster = s$seurat_clusters, patient= s$sample.id) %>%
	ggplot() + geom_bar(
		mapping = aes(x= patient, fill= cluster),
		position = "fill"
	) + labs(y="proportion")+ theme_classic() + coord_flip()
ggsave("figures/cdcpi3/cdcpi3_cluster_prop_per_sample.pdf",width= 7, height= 7, units= "in")
# differential gene analysis
DE_heatmap(s,"figures/cdcpi3/cdcpi3_DE_avgexp.pdf")

####################################
# DA analysis: CD inflamed vs CD non-inflamed
####################################
library(edgeR)
x <- subset(s, orig.ident == "CD_martin")
x$sample.id <- droplevels(x$sample.id)
x$CD_inf <- x$colitis == "CD_inflamed"
design <- function(y) model.matrix(~factor(CD_inf),y)
out <- DA_analysis(x, design, title = "cdcpi2_CD_inf_vs_uninf")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

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

#saveRDS(list(pred.CPI,pred.Martin,pred.Monaco),"saved_objects/label_info_for_cdcpi21550.rds")


require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin") + labs(title="Martin reference")
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma") + labs(title="Luoma-Suo reference")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco") + labs(title="Monaco reference")
ggarrange(p1,p2,p3, ncol=3, nrow=1)
ggsave("figures/cdcpi3/cdcpi3_anchor_labeling.pdf", width= 24, height= 6, units= "in")

saveRDS(s, "saved_objects/cdcpi3_anchor.rds")


