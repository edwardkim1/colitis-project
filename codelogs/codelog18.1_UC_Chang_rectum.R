############################################################
## Clustering UC Chang (Rectum samples)
## edward_kim@college.harvard.edu - Jan 2021
#############################################################

####################################
# UC Chang
####################################
rectum.files <- dir("data/UC_chang/",pattern="R_cell-gene")
df.rec <- lapply(paste("data/UC_chang/",rectum.files,sep=""), 
	function(x) {read.table(x, sep= "\t", header=T, row.names=1)}) 

# create Seurat objects
proj.names <- paste("U",1:7, sep="")
create.rec <- function(x) {
	df.rec[[x]] %>% as.matrix() %>% t() %>% CreateSeuratObject(project = proj.names[x-3], min.cells = 0, min.features = 0)
}

s.list <- lapply(4:10, create.rec)

####################################
# Naive integration with cells from all patients
####################################
s <- merge(s.list[[1]],c(s.list[[2]],s.list[[3]],s.list[[4]],s.list[[5]],s.list[[6]],s.list[[7]]))

########################################################
# 3 - Add patient annotations and inflammed annotations
########################################################
library(stringr)
# rename sampletype to sample.id and merge with orig.ident material
s$sample.id <- s$orig.ident
s$orig.ident.backup <- s$orig.ident
#### mutate orig.ident and edit sample.id
s$orig.ident <- "UC_chang"
#### add patient.id
s$patient.id <- s$sample.id
# add inflamed id for Crohns
s$colitis <- "UC_inflamed"
# remove genes with less than 10
x <- s[["RNA"]]@counts
nnz_row <- tabulate(x@i + 1)
keep <- rownames(x)[nnz_row>10]
s <- subset(s, features = keep)
#10104
#############################
# 4 - perform clustering
#############################
s <- PercentageFeatureSet(s, pattern = "^MT.", col.name = "percent.mt")

s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
s <- ScaleData(s, features = all.genes, vars.to.regress=c("percent.mt")) %>% RunPCA(features=VariableFeatures(s))
s <-  cluster_umap(s, npc=20, k.param=120, resolution=0.5, min.dist=0.3)
DimPlot(s, label=T)
# make directory: UC_chang_013021
ggsave("figures/UC_chang_013021/ucrec20120_umap.pdf")
DimPlot(s, group.by="sample.id" , label=F)
ggsave("figures/UC_chang_013021/ucrec20120_persample.pdf")

######################################################
# 3 - CD3E and CD4/CD8 and other QC metrics per cluster
######################################################
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/ucrec20120_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/UC_chang_013021/ucrec20120_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/ucrec20120_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/ucrec20120_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the IGKC and JCHAIN
FeaturePlot(s, features= c('IGKC','JCHAIN'), reduction="umap")
ggsave("figures/UC_chang_013021/ucrec20120_IGKC-JCHAIN_featureplot.pdf", width= 12, height= 7, units= "in")

# violin plot of IGKC per cluster
p <- VlnPlot(s,features = "IGKC", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave(paste("figures/UC_chang_013021/ucrec20120_VlnPlot_IGKC_by_cluster.pdf", sep=""), width= 7, height= 7, units= "in")

####################################
# 5 - Subset T cells
####################################
# violin plot of CD3E per cluster
s1 <- s
s <- subset(s1, idents=c("1","2","4","5","7"))
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
ggsave("figures/UC_chang_013021/ucrecT20100_umap.pdf")
DimPlot(s, group.by="sample.id" , label=F)
ggsave("figures/UC_chang_013021/ucrecT20100_persample.pdf")

# Differential Expression
DE_heatmap(s,"figures/UC_chang_013021/ucrecT20100_DE_avgexp.pdf")

######################################################
# 3 - CD3E and CD4/CD8 and other QC metrics per cluster
######################################################
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/ucrecT20100_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/UC_chang_013021/ucrecT20100_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/ucrecT20100_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/ucrecT20100_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the IGKC and JCHAIN
FeaturePlot(s, features= c('IGKC','JCHAIN'), reduction="umap")
ggsave("figures/UC_chang_013021/ucrecT20100_IGKC-JCHAIN_featureplot.pdf", width= 12, height= 7, units= "in")

# violin plot of IGKC per cluster
p <- VlnPlot(s,features = "IGKC", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave(paste("figures/UC_chang_013021/ucrecT20100_VlnPlot_IGKC_by_cluster.pdf", sep=""), width= 7, height= 7, units= "in")

# save RDS
saveRDS(s, "saved_objects/UC_chang_013021/ucrecT20100.rds")

######################################################
# 3 - Integration with CDCPI reference
######################################################
s1 <- readRDS("saved_objects/UC_chang_013021/ucrecT20100.rds")
s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
DefaultAssay(s1) <- "RNA"
s.list <- SplitObject(s1, split.by="sample.id")
s.list <- c(s2,s.list)
# variable features
s.features <- SelectIntegrationFeatures(object.list = s.list, nfeatures = 2000)
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
ggsave("figures/UC_chang_013021/uccpirec20100_elbowplot.pdf")

s <- FindNeighbors(s, dims = 1:15, k.param = 100) 
# names(s@graphs)
s <- FindClusters(s,resolution = 0.5, graph.name ="integrated_snn") 
s <- RunUMAP(s,min.dist = 0.3, graph="integrated_snn")

#visualize
p <- DimPlot(s, label=T)
ggsave("figures/UC_chang_013021/uccpirec20100_PCA.pdf")
p <- DimPlot(s, reduction="umap", label=T)
ggsave("figures/UC_chang_013021/uccpirec20100_UMAP.pdf")
p <- DimPlot(s, reduction="umap", group.by="sample.id" , label=F)
ggsave("figures/UC_chang_013021/uccpirec20100_persample.pdf")

# Differential Expression
DE_heatmap(s,"figures/UC_chang_013021/uccpirec20100_DE_avgexp.pdf")

######################################################
# 3 - CD3E and CD4/CD8 and other QC metrics per cluster
######################################################
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/uccpirec20100_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/UC_chang_013021/uccpirec20100_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/uccpirec20100_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/uccpirec20100_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the IGKC and JCHAIN
FeaturePlot(s, features= c('IGKC','JCHAIN'), reduction="umap")
ggsave("figures/UC_chang_013021/uccpirec20100_IGKC-JCHAIN_featureplot.pdf", width= 12, height= 7, units= "in")
# violin plot of IGKC per cluster
p <- VlnPlot(s,features = "IGKC", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave(paste("figures/UC_chang_013021/uccpirec20100_VlnPlot_IGKC_by_cluster.pdf", sep=""), width= 7, height= 7, units= "in")
# distribution by sample
data.frame(cluster = s$seurat_clusters, patient= s$sample.id) %>%
        ggplot() + geom_bar(
                mapping = aes(x= patient, fill= cluster),
                position = "fill"
        ) + labs(y="proportion")+ theme_classic() + coord_flip()
ggsave("figures/UC_chang_013021/uccpirec20100_cluster_prop_per_sample.pdf",width= 7, height= 7, units= "in")


###########################################
# 3 - adding metadata for the CPI portion
###########################################
library(stringr)
# add patient.id and sample.id for CPI
s$orig.ident <- ifelse(s$orig.ident=="CD3_Tcell", "CPI_colitis", s$orig.ident)
s$sample.id <- ifelse(s$orig.ident=="CPI_colitis", s$sampletype, s$sample.id)
s$patient.id <- ifelse(s$orig.ident=="CPI_colitis", s$sampletype, s$patient.id)
s$sampletype <- NULL

### condense sample.id for CPI
# for CPI colitis
s$sample.id[s$colitis=="Colitis"] %>% 
table() %>% 
rownames() -> CPI_C.samples # DFCI551 p020620 p101519   p1032   p1059    p369    p728    p806    p896    p980
replace <- paste("CPI_C",1:length(CPI_C.samples),sep="")
for(i in 1:length(CPI_C.samples)) {
        s$sample.id <- ifelse(s$patient.id==CPI_C.samples[i],replace[i], s$sample.id)
}
# for CPI no colitis
s$sample.id[s$colitis=="No-Colitis"] %>% 
table() %>% 
rownames() -> CPI_NC.samples # DFCI1294 DFCI1545 DFCI1618      p57     p590     p667     p685   p931-1 p949-1
replace <- paste("CPI_NC",1:length(CPI_NC.samples),sep="")
for(i in 1:length(CPI_NC.samples)) {
        s$sample.id <- ifelse(s$patient.id==CPI_NC.samples[i],replace[i], s$sample.id)
}
# for Controls
s$sample.id[s$colitis=="Control"] %>% 
table() %>% 
rownames() -> Control.samples # DFCI1294 DFCI1545 DFCI1618      p57     p590     p667     p685   p931-1 p949-1
replace <- paste("C",1:length(Control.samples),sep="")
for(i in 1:length(Control.samples)) {
        s$sample.id <- ifelse(s$patient.id==Control.samples[i],replace[i], s$sample.id)
}
# add levels
s$sample.id <- factor(s$sample.id, levels = c("C1","C2","C3","C4","C5","C6","C7","C8",
"CPI_NC1", "CPI_NC2", "CPI_NC3", "CPI_NC4", "CPI_NC5", "CPI_NC6", "CPI_NC7", "CPI_NC8" , "CPI_NC9",
"CPI_C1", "CPI_C2",  "CPI_C3",  "CPI_C4" , "CPI_C5" , "CPI_C6",  "CPI_C7",  "CPI_C8",  "CPI_C9" ,"CPI_C10",
"U1","U2","U3","U4","U5","U6","U7"))

## rename UC_inflamed to just UC
s$colitis <- ifelse(s$colitis == "UC_inflamed", "UC", s$colitis)


###########################################
# 3 - Redo select figures
###########################################
# distribution by sample
data.frame(cluster = s$seurat_clusters, patient= s$sample.id) %>%
        ggplot() + geom_bar(
                mapping = aes(x= patient, fill= cluster),
                position = "fill"
        ) + labs(y="proportion",x="sample") theme_classic() + coord_flip()
ggsave("figures/UC_chang_013021/uccpirec20100_cluster_prop_per_sample2.pdf",width= 7, height= 7, units= "in")


####################################
# 5 - Cluster Annotations
####################################
load('data/ssuo_CD3/seurat.object.RData')
Tcell <- process_ssuo_Tcells(Tcell)

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

ref <- MonacoImmuneData()
pred.Monaco <- SingleR(method="cluster",
        test=s[['integrated']]@data,
        clusters= s$seurat_clusters, 
        ref=ref, 
        labels=ref$label.fine, 
        BPPARAM=MulticoreParam(10))
s$SingleR.Monaco <- s$seurat_clusters
levels(s$SingleR.Monaco) <- pred.Monaco$labels

require(ggpubr)
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma") + labs(title="Luoma-Suo reference")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco") + labs(title="Monaco reference")
ggarrange(p2,p3, ncol=2, nrow=1)
ggsave("figures/UC_chang_013021/uccpirec20100_anchor_labeling.pdf", width= 18, height= 6, units= "in")

###########################################
# 3 - Remove Cluster 12 and redo PCA/UMAP
###########################################
s1 <- s
s <- subset(s1, idents=c("12"),invert=T)
# Do the rest; no need to regress out percent.mt because those genes are not in the list
s <- ScaleData(s) %>% RunPCA()
# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 20)
ggsave("figures/UC_chang_013021/uccpirec2_elbowplot.pdf")

s <- FindNeighbors(s, dims = 1:15, k.param = 100) 
# names(s@graphs)
s <- FindClusters(s,resolution = 0.5, graph.name ="integrated_snn") 
s <- RunUMAP(s,min.dist = 0.3, graph="integrated_snn")

#visualize
p <- DimPlot(s, reduction="pca", label=T)
ggsave("figures/UC_chang_013021/uccpirec2_PCA.pdf")
p <- DimPlot(s, reduction="umap", label=T)
ggsave("figures/UC_chang_013021/uccpirec2_UMAP.pdf")
p <- DimPlot(s, reduction="umap", group.by="sample.id" , label=F)
ggsave("figures/UC_chang_013021/uccpirec2_persample.pdf")

# Differential Expression
DE_heatmap(s,"figures/UC_chang_013021/uccpirec2_DE_avgexp.pdf")

######################################################
# 3 - CD3E and CD4/CD8 and other QC metrics per cluster
######################################################
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/uccpirec2_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/UC_chang_013021/uccpirec2_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/uccpirec2_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_013021/uccpirec2_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the IGKC and JCHAIN
FeaturePlot(s, features= c('IGKC','JCHAIN'), reduction="umap")
ggsave("figures/UC_chang_013021/uccpirec2_IGKC-JCHAIN_featureplot.pdf", width= 12, height= 7, units= "in")
# violin plot of IGKC per cluster
p <- VlnPlot(s,features = "IGKC", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave(paste("figures/UC_chang_013021/uccpirec2_VlnPlot_IGKC_by_cluster.pdf", sep=""), width= 7, height= 7, units= "in")
# distribution by sample
data.frame(cluster = s$seurat_clusters, patient= s$sample.id) %>%
        ggplot() + geom_bar(
                mapping = aes(x= patient, fill= cluster),
                position = "fill"
        ) + labs(y="proportion")+ theme_classic() + coord_flip()
ggsave("figures/UC_chang_013021/uccpirec2_cluster_prop_per_sample.pdf",width= 7, height= 7, units= "in")


####################################
# 5 - Cluster Annotations
####################################
load('data/ssuo_CD3/seurat.object.RData')
Tcell <- process_ssuo_Tcells(Tcell)

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

ref <- MonacoImmuneData()
pred.Monaco <- SingleR(method="cluster",
        test=s[['integrated']]@data,
        clusters= s$seurat_clusters, 
        ref=ref, 
        labels=ref$label.fine, 
        BPPARAM=MulticoreParam(10))
s$SingleR.Monaco <- s$seurat_clusters
levels(s$SingleR.Monaco) <- pred.Monaco$labels

require(ggpubr)
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma") + labs(title="Luoma-Suo reference")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco") + labs(title="Monaco reference")
ggarrange(p2,p3, ncol=2, nrow=1)
ggsave("figures/UC_chang_013021/uccpirec2_anchor_labeling.pdf", width= 18, height= 6, units= "in")


#save RDS
saveRDS(s, "saved_objects/UC_chang_013021/uccpirec2.rds")


####################################
# 5 - DE per cluster (CPI/UC)
####################################
get_labels <- function(s.markers, FCcutoff = 0.8) {
        head(s.markers, n=20) %>% rownames() -> labels
        up.labels <- rownames(s.markers[s.markers$avg_log2FC > FCcutoff,])
        down.labels <- rownames(s.markers[s.markers$avg_log2FC < -FCcutoff,])
        labels <- unique(c(labels, up.labels,down.labels))
        return(labels)
}
lapply(markers.out,get_labels)


DE_volcano <- function(seurat.object, ident.1, group.by,subset.ident, FCcutoff = 0.8,drawConnectors=T, title, xlim,legendPosition = "top") {
        require(EnhancedVolcano)
        s.markers <- FindMarkers(seurat.object, ident.1 =ident.1 ,group.by=group.by, subset.ident= subset.ident, min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
        head(s.markers, n=20) %>% rownames() -> labels
        up.labels <- rownames(s.markers[s.markers$avg_log2FC > FCcutoff,])
        down.labels <- rownames(s.markers[s.markers$avg_log2FC < -FCcutoff,])
        labels <- unique(c(up.labels,down.labels))
        p <- EnhancedVolcano(s.markers,
            lab = rownames(s.markers),
            x = 'avg_log2FC',
            y = 'p_val_adj',
            xlim = xlim,
            title = title,
            pCutoff = 0.05,
            FCcutoff = FCcutoff,
            pointSize = 3.0,
            labSize = 4.0,
            selectLab = labels,
            drawConnectors = drawConnectors,
            widthConnectors = 0.2,
            colConnectors = 'grey30',
            border="full",
            borderWidth = 1.0,
            borderColour = 'black',
            colAlpha = 1,
            gridlines.major = FALSE,
            gridlines.minor= FALSE,
            legendPosition = legendPosition)
        return(list(markers = s.markers, plot=p))
}


# compare the UC/CPI datasets
out <- lapply(0:11, function(x) DE_volcano(s, ident.1 ="UC_chang" ,group.by= "orig.ident", subset.ident= x, FCcutoff = 0.8, xlim = c(-2, 2),
        title = paste('Cluster',x-1, 'UC vs. CPI dataset'),
        legendPosition= "none"))

library(gridExtra)
library(grid)
options(ggrepel.max.overlaps = Inf)
cluster.names <- c("Type 1 cytokines Trm",
        "Cytotoxic effector",
        "CD8 Trm",
        "LP Trm",
        "IEL: CD8/gd",
        "TFH",
        "Treg",
        "Naive/CM",
        "Th1 effector",
        "Cycling",
        "10",
        "MAIT")
plist <- lapply(1:12, function(x) {out[[x]]$plot + theme(axis.title.x=element_blank(), axis.title.y=element_blank())} + labs(subtitle = paste("Cluster ",x-1,": ",cluster.names[x],sep=""),title=NULL) )
p <- grid.arrange(grobs= plist, ncol=4,
        bottom= textGrob(label= expression(paste("Log"[2]*" fold change")),
                gp = gpar(col = "black", fontsize = 20)), 
        left=textGrob(label= expression(paste("-Log"[10]*italic(P))), rot=90, 
                gp = gpar(col = "black", fontsize = 20)),
        top = textGrob(label= expression(bold("UC/CPI Per Cluster Dataset Differential Expression")), just = "center",
                gp = gpar(col = "black", fontsize = 24)))

ggsave(paste("figures/UC_chang_013021/uccpirec2_DE_dataset_comparison.pdf" ,sep=""),p, width= 8*4, height= 8*3, units= "in")



# compare UC/CPI colitis conditions
s2 <- s[,s$colitis=="Colitis" | s$colitis=="UC"]
out2 <- lapply(0:11, function(x) DE_volcano(s2, ident.1 ="UC_chang" ,group.by= "orig.ident", subset.ident= x, FCcutoff = 0.8, xlim = c(-2, 2),
        title = paste('Cluster',x-1, 'UC vs. CPI dataset'),
        legendPosition= "none"))
plist <- lapply(1:12, function(x) {out2[[x]]$plot + theme(axis.title.x=element_blank(), axis.title.y=element_blank())} + labs(subtitle = paste("Cluster ",x-1,": ",cluster.names[x],sep=""),title=NULL) )
p <- grid.arrange(grobs= plist, ncol=4,
        bottom= textGrob(label= expression(paste("Log"[2]*" fold change")),
                gp = gpar(col = "black", fontsize = 20)), 
        left=textGrob(label= expression(paste("-Log"[10]*italic(P))), rot=90, 
                gp = gpar(col = "black", fontsize = 20)),
        top = textGrob(label= expression(bold("UC/CPI-colitis Per Cluster Differential Expression")), just = "center",
                gp = gpar(col = "black", fontsize = 24)))

ggsave(paste("figures/UC_chang_013021/uccpirec2_DE_colitis_comparison.pdf" ,sep=""),p, width= 8*4, height= 8*3, units= "in")



####################################
# 5 - DA figure
####################################
p <- plot.propbar.errors(seurat.object= s,
        annotations = s$seurat_clusters,
        clusters= 0:11,
        group.names= c("Control","No-Colitis","Colitis","UC"))
ggsave("figures/UC_chang_013021/uccpirec2_barplot.pdf", width= 12, height= 7, units= "in")

#### Wilcoxon Test (Unpaired) ####
wilcoxon.p.values <- function(clusters, df.prop) {
        p.values <- unlist(lapply(clusters,function(x) {wilcox.test(proportion ~ status, data = subset(df.prop, cluster==x), paired = FALSE)$p.value}))
        names(p.values) <- unique(df.prop$cluster)
        return(p.values)
}
# get p-values using wilcoxon test
s.pvals <- wilcoxon.p.values(0:11,subset(df.s, status=="Colitis" | status=="UC"))

##############################################
# DA analysis: CPI colitis vs CPI no-colitis
##############################################
s2 <- subset(s, colitis == "Colitis" | colitis == "UC")
s2$sample.id <- droplevels(s2$sample.id)
s2$UCoverCPI <- s2$colitis == "UC"
design <- function(y) model.matrix(~factor(UCoverCPI),y)
out <- DA_analysis(s2, design, title = "UC_chang_013021/uccpirec2_DA")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)






