############################################################
## Clustering UC Chang (Rectum samples) all cells
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
# UC Chang
####################################
rectum.files <- dir("data/UC_chang/",pattern="R_cell-gene")
df.rec <- lapply(paste("data/UC_chang/",rectum.files,sep=""), 
	function(x) {read.table(x, sep= "\t", header=T, row.names=1)}) 

# create Seurat objects
# Controls 9-16 are UC
proj.names <- c(paste("C",9:11,sep="") ,paste("U",1:7, sep=""),paste("C",12:16,sep="") )

create.rec <- function(x) {
	temp <- df.rec[[x]] %>% as.matrix() %>% t() 
        rownames(temp) <- gsub("\\.", "-", rownames(temp))
        CreateSeuratObject(temp, project = proj.names[x], min.cells = 0, min.features = 0)
}

s.list <- lapply(1:15, create.rec)

####################################
# Naive integration with cells from all patients
####################################
s <- merge(s.list[[1]],c(s.list[[2]],s.list[[3]],s.list[[4]],s.list[[5]],s.list[[6]],s.list[[7]],s.list[[8]],s.list[[9]],s.list[[10]],s.list[[11]],s.list[[12]],s.list[[13]],s.list[[14]],s.list[[15]]))

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
s$colitis <- ifelse(grepl("C",s$sample.id), "Control", "UC")
# remove genes with less than 10
x <- s[["RNA"]]@counts
nnz_row <- tabulate(x@i + 1)
keep <- rownames(x)[nnz_row>10]
s <- subset(s, features = keep)
#10105
#############################
# 4 - perform clustering
#############################
s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")

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
ggsave("figures/UC_chang_021321/ucrec20120_umap.pdf")
DimPlot(s, group.by="sample.id" , label=F)
ggsave("figures/UC_chang_021321/ucrec20120_persample.pdf")

# differential expression
DE_heatmap(s,"figures/UC_chang_021321/ucrec20120_DE_avgexp.pdf")

# save RDS
saveRDS(s, "saved_objects/UC_chang_021321/ucrec20120.rds")


#######################################################
# 3 - Load CD45 cells and use as integration reference
#######################################################
load('data/ssuo_CD45/seurat.object.new.RData')
#DimPlot(Tcell, label=T, reduction = "umap")
#ggsave("figures/UC_chang_021321/ssuoCD45_umap.pdf")
# Prepare the metadata before integration
Tcell$sample.id <- Tcell$sampletype
Tcell$patient.id <- Tcell$patientid

# begin integration
s.list <- SplitObject(s, split.by="sample.id")
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
ggsave("figures/UC_chang_021321/uccpirec20100_elbowplot.pdf")

s <- FindNeighbors(s, dims = 1:15, k.param = 100)  
# names(s@graphs)
s <- FindClusters(s,resolution = 0.5, graph.name ="integrated_snn") 
s <- RunUMAP(s,min.dist = 0.4, graph="integrated_snn")

#visualize
p <- DimPlot(s, reduction = "pca",label=T)
ggsave("figures/UC_chang_021321/uccpirec20100_PCA.pdf")
p <- DimPlot(s, reduction="umap", label=T)
ggsave("figures/UC_chang_021321/uccpirec20100_UMAP.pdf")
p <- DimPlot(s, reduction="umap", group.by="patient.id" , label=F)
ggsave("figures/UC_chang_021321/uccpirec20100_persample.pdf", width=14, height=10, units= "in")

# Differential Expression
DefaultAssay(s) <- "RNA"
out <- DE_heatmap(s)
ggsave("figures/UC_chang_021321/uccpirec20100_DE_avgexp.pdf",out$plot, width=7, height=9, units= "in")
DefaultAssay(s) <- "integrated"

######################################################
# 3 - QC metrics per cluster
######################################################
s$orig.ident <- ifelse(s$orig.ident == "10X_Tcell","CPI_colitis",s$orig.ident)
DefaultAssay(s) <- "RNA"
s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")

# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_021321/uccpirec20100_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_021321/uccpirec20100_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")

####################################
# 5 - Cluster Annotations (CD45+)
####################################
load('data/ssuo_CD45/seurat.object.new.RData')

require(SingleR)
require(BiocParallel)
pred.CPI <- SingleR(method='cluster',
        test=as.matrix(s[['integrated']]@data), 
        clusters= s$seurat_clusters,
        ref=as.matrix(Tcell[['RNA']]@data), 
        labels=Tcell$anno2,
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
ggsave("figures/UC_chang_021321/uccpirec20100_anchor_labeling.pdf", width= 18, height= 6, units= "in")

# second labeling with the RNA as default
pred.CPI <- SingleR(method='cluster',
        test=as.matrix(s[['RNA']]@data), 
        clusters= s$seurat_clusters,
        ref=as.matrix(Tcell[['RNA']]@data), 
        labels=Tcell$anno2,
        de.method="wilcox",
        BPPARAM=MulticoreParam(10))
s$SingleR.Luoma <- s$seurat_clusters
levels(s$SingleR.Luoma) <- pred.CPI$labels

ref <- MonacoImmuneData()
pred.Monaco <- SingleR(method="cluster",
        test=s[['RNA']]@data,
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
ggsave("figures/UC_chang_021321/uccpirec20100_anchor_labeling_RNA.pdf", width= 18, height= 6, units= "in")

# third labeling figure
pred.CPI <- SingleR(method='cluster',
        test=as.matrix(s[['integrated']]@data), 
        clusters= s$seurat_clusters,
        ref=as.matrix(Tcell[['RNA']]@data), 
        labels=Tcell$banno1,
        de.method="wilcox",
        BPPARAM=MulticoreParam(10))
s$SingleR.Luoma2 <- s$seurat_clusters
levels(s$SingleR.Luoma2) <- pred.CPI$labels

p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma2") + labs(title="Luoma-Suo reference")
ggsave("figures/UC_chang_021321/uccpirec20100_anchor_labeling2.pdf", width= 8, height= 6, units= "in")

# fourth labeling figure with the UC annotations
#uc.meta <- read.csv("data/UC_chang/GSE125527_cell_metadata.csv")
#cell.names <- colnames(s[,s$orig.ident=="UC_chang"])
#cell.names <- word(cell.names,1,sep = "\\_")


#s$patient.id2 <- s$patient.id
#replace <- uc.meta$patient_assignment[which(cell.names %in% uc.meta$cell_id)]
#replace <- c(rep(NA, sum(s$orig.ident=="CPI_colitis")), replace)

#s$patient.id2 <- if(s$orig.ident=="UC_chang", )


# save RDS
saveRDS(s, "saved_objects/UC_chang_021321/uccpirec-allCD45.rds")

####################################
# 5 - DA figure
####################################
p <- plot.propbar.errors(annotations = s$SingleR.Luoma2,
        sample.id = s$sample.id,
        conditions = s$colitis,
        clusters= names(table(s$SingleR.Luoma2))[c(4,2,1,3,7,5,6)],
        group.names= c("Control","No-Colitis","Colitis","UC"))
ggsave("figures/UC_chang_021321/uccpirec20100_barplot.pdf", width= 12, height= 7, units= "in")

#### Wilcoxon Test (Unpaired) ####
wilcoxon.p.values <- function(clusters, df.prop) {
        p.values <- unlist(lapply(clusters,function(x) {wilcox.test(proportion ~ status, data = subset(df.prop, cluster==x), paired = FALSE)$p.value}))
        names(p.values) <- unique(df.prop$cluster)
        return(p.values)
}
# get p-values using wilcoxon test
df.s <- wilcoxon.df(s, s$SingleR.Luoma2)
s.pvals <- wilcoxon.p.values(names(table(s$SingleR.Luoma2)),subset(df.s, status=="Colitis" | status=="UC"))

##############################################
# DE analysis: UC dataset vs CPI dataset
##############################################
out <- DE_volcano(s, ident.1 ="UC_chang" ,assay = "RNA",group.by= "orig.ident", subset.ident= NULL, FCcutoff = 0.8, xlim = c(-2, 2),
        title = paste('UC dataset / CPI dataset'),
        legendPosition= "none")

ggsave("figures/UC_chang_021321/uccpirec20100_DE_dataset_comparison.pdf",out$p + labs(subtitle = NULL), width= 8, height= 8, units= "in")

##############################################
# DA analysis
##############################################
s2 <- subset(s, colitis == "Colitis" | colitis == "UC")
s2$sample.id <- droplevels(s2$sample.id)
s2$UCoverCPI <- s2$colitis == "UC"
design <- function(y) model.matrix(~factor(UCoverCPI),y)
out <- DA_analysis(s2, design, title = "UC_chang_021321/uccpirec2_DA", annotations = s2$SingleR.Luoma2)
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

##############################################
# Reference Mapping
##############################################
s <- readRDS("saved_objects/UC_chang_021321/ucrec20120.rds")

Tcell <- RunUMAP(Tcell, return.model=TRUE,dims=1:25)
p <- DimPlot(Tcell, reduction="umap", label=F)
ggsave("figures/UC_chang_021321/ssuoCD45_new_umap.pdf", width= 9, height= 7, units= "in")

anchors <- FindTransferAnchors(
  reference = Tcell,
  query = s,
  reference.assay = "RNA",
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

s.postMQ <- MapQuery(
  anchorset = anchors,
  query = s,
  reference = Tcell,
  refdata = Tcell$banno1,
  reference.reduction = "pca",
  reduction.model = "umap"
)

#visualizing the mapping results
p1 = DimPlot(s.postMQ, reduction = "ref.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
ggsave("figures/UC_chang_021321/UC_cpirefUMAP.pdf", width= 7, height= 7, units= "in")

# Plot barplot comparing cluster proportions
s.postMQ$annotations <- s.postMQ$predicted.id
Tcell$annotations <- Tcell$banno1
s2 <- merge(Tcell,s.postMQ)

p <- plot.propbar.errors(seurat.object= s2,
        annotations = s2$annotations,
        clusters= names(table(s2$annotations))[c(5,1,2,3,4,6,7)],
        group.names= c("Control","No-Colitis","Colitis","UC"))
ggsave("figures/UC_chang_021321/uccpi_ref_barplot.pdf", width= 12, height= 7, units= "in")

# violin plot of predicted id score per cluster
s.postMQ$predicted.id <- factor(s.postMQ$predicted.id, levels=names(table(s.postMQ$predicted.id))[c(5,1,2,3,4,6,7)])
p <- VlnPlot(s.postMQ,features = "predicted.id.score", pt.size= 0, group.by="predicted.id") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/UC_chang_021321/uccpi_ref_predictedID_vlnplot.pdf", width= 7, height= 7, units= "in")

#################################################################
# 10 - Clustering check with comparison to within-batch clusters
#################################################################
library(pheatmap)
library(grid)
# For the first batch (adding +1 for a smoother color transition
# from zero to non-zero counts for any given matrix entry).
int.labels <- s1[,s1$orig.ident== "UC_chang"]$SingleR.Luoma2
rf.labels <- s.postMQ$predicted.id

tab <- table(int.labels,rf.labels)
heatmap <- pheatmap(matrix.sort.no.diag(tab), cluster_row=FALSE, cluster_col=FALSE,
    main="Annotation Comparison", silent=TRUE)

pdf("figures/UC_chang_021321/uccpi_ref_v_integration_comparison.pdf")
heatmap
dev.off()


#log10
heatmap <- pheatmap(matrix.sort.no.diag(log10(tab+1)), cluster_row=FALSE, cluster_col=FALSE,
    main="Annotation Comparison", silent=TRUE)

pdf("figures/UC_chang_021321/uccpi_ref_v_integration_comparison_log10.pdf")
heatmap
dev.off()


####################################################
# 10 - Comparing Controls (Integrated dataset)
####################################################
s$colitis3 <- ifelse(s$colitis == "Control" & s$orig.ident=="UC_chang", "UC-Control", s$colitis)

p <- plot.propbar.errors(annotations = s$SingleR.Luoma2,
        sample.id = s$sample.id,
        conditions = s$colitis3,
        clusters= names(table(s$SingleR.Luoma2))[c(4,2,1,3,7,5,6)],
        group.names= c("Control","UC-Control","No-Colitis","Colitis","UC"))
ggsave("figures/UC_chang_021321/uccpirec2_barplot.pdf", width= 12, height= 7, units= "in")

##############################################
# fGSEA dataset comparison
##############################################
# install.packages("msigdbr")
# BiocManager::install("org.Hs.eg.db")
library(fgsea)
library(data.table)
library(ggplot2)

s <- readRDS("saved_objects/UC_chang_021321/uccpirec-allCD45.rds")

# pulling the hallmarks geneset from the msigdbr
library(msigdbr)
hallmarks <- msigdbr(species = "Homo sapiens", category = "H")

# getting the genes
out <- DE_volcano(s, ident.1 ="UC_chang" ,assay = "RNA",group.by= "orig.ident", subset.ident= NULL, FCcutoff = 0.8, xlim = c(-3, 3),
        title = paste('UC dataset / CPI dataset'),
        legendPosition= "none")

# compute fgseaRes
fRes <- get_fgseaRes(msigdbr.gs= hallmarks, seurat4.markers= out$markers)
# no hallmarks matched

# plot the fgsea
fRes <- list(fRes)
anno.index <- 1
topPathwaysUp <- fRes[[anno.index]]$result %>% 
                   subset(ES > 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10, pval) %>%
                   pull(pathway)
topPathwaysDown <- fRes[[anno.index]]$result %>% 
                   subset(ES < 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10,pval) %>%
                   pull(pathway)
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf("figures/UC_chang_021321/uccpirec2_fgsea_dataset_comparison.pdf",width=12,height=5)
plotGseaTable(fRes[[anno.index]]$pathways[topPathways], fRes[[anno.index]]$ranks, fRes[[anno.index]]$result, gseaParam=0.5)
dev.off()


##############################################
# fGSEA per broad cell type annotation
##############################################
library(fgsea)
library(data.table)
library(ggplot2)

s <- readRDS("saved_objects/UC_chang_021321/uccpirec-allCD45.rds")

# pulling the hallmarks geneset from the msigdbr
library(msigdbr)
hallmarks <- msigdbr(species = "Homo sapiens", category = "H")

# getting the genes
# subset UC colitis and CPI-colitis
s2 <- s[,s$colitis=="UC" | s$colitis=="Colitis"]
broad.clusters <- s$SingleR.Luoma2 %>% table() %>% names()

DefaultAssay(s2) <- "RNA"
Idents(s2) <- s2$SingleR.Luoma2
out <- lapply(broad.clusters, function(x) DE_volcano(s2, ident.1 ="UC_chang" , assay = "RNA", group.by= "orig.ident", subset.ident= x, FCcutoff = 0.8, xlim = c(-3, 3),
        title = paste('Cluster ',x, ': UC/CPI-colitis', sep=""),
        legendPosition= "none"))

# plot the DE
plist <- lapply(1:7, function(x) {out[[x]]$plot + theme(axis.title.x=element_blank(), axis.title.y=element_blank())} + labs(subtitle = paste("Cell Type: ",broad.clusters[x],sep=""),title=NULL) )
p <- grid.arrange(grobs= plist, ncol=3,
        bottom= textGrob(label= expression(paste("Log"[2]*" fold change")),
                gp = gpar(col = "black", fontsize = 20)), 
        left=textGrob(label= expression(paste("-Log"[10]*italic(P))), rot=90, 
                gp = gpar(col = "black", fontsize = 20)),
        top = textGrob(label= expression(bold("UC/CPI-Colitis Differential Expression")), just = "center",
                gp = gpar(col = "black", fontsize = 24)))

ggsave(paste("figures/UC_chang_021321/uccpirec2_DE_banno2_comparison.pdf" ,sep=""),p, width= 8*4, height= 8*3, units= "in")

# compute fgseaRes
fRes <- lapply(1:7, function(x) get_fgseaRes(msigdbr.gs= hallmarks, seurat4.markers= out[[x]]$markers))

# plot the fRes

#collapsedPathways <- collapsePathways(fRes[order(pval)][padj < 0.01], 
#                                      examplePathways, exampleRanks)
#mainPathways <- fRes[pathway %in% collapsedPathways$mainPathways][
#                         order(-NES), pathway]

# IgA Plasma B cells
topPathways <- fRes[[2]]$pathway
pdf("figures/UC_chang_021321/uccpirec2_IgA_fgsea.pdf",width=10,height=3)
plotGseaTable(fRes[[2]]$pathways[topPathways], fRes[[2]]$ranks, fRes[[2]]$result, gseaParam=0.5)
dev.off()

# CD8 T cells
topPathways <- fRes[[2]]$pathway
pdf("figures/UC_chang_021321/uccpirec2_CD8_fgsea.pdf",width=10,height=3)
plotGseaTable(fRes[[3]]$pathways[topPathways], fRes[[3]]$ranks, fRes[[3]]$result, gseaParam=0.5)
dev.off()

# Myeloid cells
topPathwaysUp <- fRes[[5]]$result %>% 
                   subset(ES > 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10, pval) %>%
                   pull(pathway)
topPathwaysDown <- fRes[[5]]$result %>% 
                   subset(ES < 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10,pval) %>%
                   pull(pathway)
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf("figures/UC_chang_021321/uccpirec2_myeloid_fgsea.pdf",width=10,height=5)
plotGseaTable(fRes[[5]]$pathways[topPathways], fRes[[5]]$ranks, fRes[[5]]$result, gseaParam=0.5)
dev.off()

# IgA Plasma B cells Take 2 (with fgsea without the 0.8 constraint)
anno.index <- 2
topPathwaysUp <- fRes[[anno.index]]$result %>% 
                   subset(ES > 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10, pval) %>%
                   pull(pathway)
topPathwaysDown <- fRes[[anno.index]]$result %>% 
                   subset(ES < 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10,pval) %>%
                   pull(pathway)
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf("figures/UC_chang_021321/uccpirec2_IgA_fgsea_2.pdf",width=10,height=5)
plotGseaTable(fRes[[anno.index]]$pathways[topPathways], fRes[[anno.index]]$ranks, fRes[[anno.index]]$result, gseaParam=0.5)
dev.off()

# CD8 T cells Take 2 (with fgsea without the 0.8 constraint)
anno.index <- 3
topPathwaysUp <- fRes[[anno.index]]$result %>% 
                   subset(ES > 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10, pval) %>%
                   pull(pathway)
topPathwaysDown <- fRes[[anno.index]]$result %>% 
                   subset(ES < 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10,pval) %>%
                   pull(pathway)
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf("figures/UC_chang_021321/uccpirec2_CD8_fgsea_2.pdf",width=10,height=5)
plotGseaTable(fRes[[anno.index]]$pathways[topPathways], fRes[[anno.index]]$ranks, fRes[[anno.index]]$result, gseaParam=0.5)
dev.off()

# CD4 T cells (with fgsea without the 0.8 constraint)
anno.index <- 1
topPathwaysUp <- fRes[[anno.index]]$result %>% 
                   subset(ES > 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10, pval) %>%
                   pull(pathway)
topPathwaysDown <- fRes[[anno.index]]$result %>% 
                   subset(ES < 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10,pval) %>%
                   pull(pathway)
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf("figures/UC_chang_021321/uccpirec2_CD4_fgsea_2.pdf",width=10,height=5)
plotGseaTable(fRes[[anno.index]]$pathways[topPathways], fRes[[anno.index]]$ranks, fRes[[anno.index]]$result, gseaParam=0.5)
dev.off()

# B cells (with fgsea without the 0.8 constraint)
anno.index <- 4
topPathwaysUp <- fRes[[anno.index]]$result %>% 
                   subset(ES > 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10, pval) %>%
                   pull(pathway)
topPathwaysDown <- fRes[[anno.index]]$result %>% 
                   subset(ES < 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10,pval) %>%
                   pull(pathway)
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf("figures/UC_chang_021321/uccpirec2_B_fgsea_2.pdf",width=10,height=5)
plotGseaTable(fRes[[anno.index]]$pathways[topPathways], fRes[[anno.index]]$ranks, fRes[[anno.index]]$result, gseaParam=0.5)
dev.off()

# Mast cells (with fgsea without the 0.8 constraint)
anno.index <- 6
topPathwaysUp <- fRes[[anno.index]]$result %>% 
                   subset(ES > 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10, pval) %>%
                   pull(pathway)
topPathwaysDown <- fRes[[anno.index]]$result %>% 
                   subset(ES < 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10,pval) %>%
                   pull(pathway)
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf("figures/UC_chang_021321/uccpirec2_Mast_fgsea_2.pdf",width=10,height=5)
plotGseaTable(fRes[[anno.index]]$pathways[topPathways], fRes[[anno.index]]$ranks, fRes[[anno.index]]$result, gseaParam=0.5)
dev.off()

# ILCs (with fgsea without the 0.8 constraint)
anno.index <- 7
topPathwaysUp <- fRes[[anno.index]]$result %>% 
                   subset(ES > 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10, pval) %>%
                   pull(pathway)
topPathwaysDown <- fRes[[anno.index]]$result %>% 
                   subset(ES < 0 & padj < 0.05) %>%
                   arrange(pval) %>%
                   top_n(10,pval) %>%
                   pull(pathway)
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf("figures/UC_chang_021321/uccpirec2_ILCs_fgsea_2.pdf",width=10,height=5)
plotGseaTable(fRes[[anno.index]]$pathways[topPathways], fRes[[anno.index]]$ranks, fRes[[anno.index]]$result, gseaParam=0.5)
dev.off()



##############################################
# fGSEA per broad cell type annotation (figure)
##############################################
library(fgsea)
library(data.table)
library(ggplot2)

s <- readRDS("saved_objects/UC_chang_021321/uccpirec-allCD45.rds")

# pulling the hallmarks geneset from the msigdbr
library(msigdbr)
hallmarks <- msigdbr(species = "Homo sapiens", category = "H")

# getting the genes
# subset UC colitis and CPI-colitis
s2 <- s[,s$colitis=="UC" | s$colitis=="Colitis"]
broad.clusters <- s$SingleR.Luoma2 %>% table() %>% names()

DefaultAssay(s2) <- "RNA"
Idents(s2) <- s2$SingleR.Luoma2
out <- lapply(broad.clusters, function(x) DE_volcano(s2, ident.1 ="UC_chang" , assay = "RNA", group.by= "orig.ident", subset.ident= x, FCcutoff = 0.8, xlim = c(-3, 3),
        title = paste('Cluster ',x, ': UC/CPI-colitis', sep=""),
        legendPosition= "none"))

# plot the DE
plist <- lapply(1:7, function(x) {out[[x]]$plot + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + labs(subtitle = paste("Cell Type: ",broad.clusters[x],sep=""),title=NULL)} )
p <- grid.arrange(grobs= plist, ncol=3,
        bottom= textGrob(label= expression(paste("Log"[2]*" fold change")),
                gp = gpar(col = "black", fontsize = 20)), 
        left=textGrob(label= expression(paste("-Log"[10]*italic(P))), rot=90, 
                gp = gpar(col = "black", fontsize = 20)),
        top = textGrob(label= expression(bold("UC/CPI-Colitis Differential Expression")), just = "center",
                gp = gpar(col = "black", fontsize = 24)))

ggsave(paste("figures/UC_chang_021321/uccpirec2_DE_banno2_comparison.pdf" ,sep=""),p, width= 8*4, height= 8*3, units= "in")

# compute fgseaRes
fRes <- lapply(1:7, function(x) get_fgseaRes(msigdbr.gs= hallmarks, seurat4.markers= out[[x]]$markers))


# heatmap visualization of the leadingEdge genes in the GSEA analysis for TNF_signaling
library(pheatmap)
getLeadingEdgeGenes <- function(index) {
        fgseaResMain <- fRes[[index]]$result[fRes[[index]]$result$pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]
        if(nrow(fgseaResMain)==0) {
                return(0)
        } else {fgseaResMain[, leadingEdge := mapIdsList(
                                             x=org.Hs.eg.db, 
                                             keys=leadingEdge,
                                             keytype="ENTREZID", 
                                             column="SYMBOL")]
                return(fgseaResMain$leadingEdge[[1]]) 
        }
        
}
leadingEdgeGenes <- lapply(1:7, getLeadingEdgeGenes)
names(leadingEdgeGenes) <- broad.clusters

leadingEdgeGenes
nLEG <- lapply(leadingEdgeGenes, length) %>% unlist()
labels <- rep(broad.clusters, nLEG)
genes <- leadingEdgeGenes %>% unlist()
tab <- table(labels,genes)
if(colnames(tab)[1] == "0") tab <- tab[,-1] # remove the 0 if present
# rearrange rows
anno <- c("IgA plasma B","B", "CD4 T", "CD8 T", "ILCs","Myeloid","Mast")
tab.ordered <- tab[c(2,3,1,5,4,7,6),]
# save figure
pdf("figures/UCCPI_DE_GSEA_031221/fgsea_pheatmap_TNF.pdf",width=8,height=3)
pheatmap(tab.ordered, cluster_row=F, cluster_cols=T,fontsize= 8, angle_col=c(90),
        legend=F)
dev.off()



##############################################
# fGSEA per broad cell type annotation (controls comparison)
##############################################
library(fgsea)
library(data.table)
library(ggplot2)

s <- readRDS("saved_objects/UC_chang_021321/uccpirec-allCD45.rds")

# pulling the hallmarks geneset from the msigdbr
library(msigdbr)
hallmarks <- msigdbr(species = "Homo sapiens", category = "H")

# getting the genes
# subset UC colitis and CPI-colitis
s2 <- s[,s$colitis3=="UC-Control" | s$colitis3=="Control"]
broad.clusters <- s$SingleR.Luoma2 %>% table() %>% names()

DefaultAssay(s2) <- "RNA"
Idents(s2) <- s2$SingleR.Luoma2
out <- lapply(broad.clusters, function(x) DE_volcano(s2, ident.1 ="UC_chang" , assay = "RNA", group.by= "orig.ident", subset.ident= x, FCcutoff = 0.8, xlim = c(-3, 3),
        title = paste('Cluster ',x, ': UC/CPI-colitis', sep=""),
        legendPosition= "none"))

# plot the DE
plist <- lapply(1:7, function(x) {out[[x]]$plot + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + labs(subtitle = paste("Cell Type: ",broad.clusters[x],sep=""),title=NULL)} )
p <- grid.arrange(grobs= plist, ncol=3,
        bottom= textGrob(label= expression(paste("Log"[2]*" fold change")),
                gp = gpar(col = "black", fontsize = 20)), 
        left=textGrob(label= expression(paste("-Log"[10]*italic(P))), rot=90, 
                gp = gpar(col = "black", fontsize = 20)),
        top = textGrob(label= expression(bold("UC/CPI-Controls Differential Expression")), just = "center",
                gp = gpar(col = "black", fontsize = 24)))

ggsave(paste("figures/UCCPI_DE_GSEA_031221/uccpirec2_DE_banno2_comparison_of_controls.pdf" ,sep=""),p, width= 8*4, height= 8*3, units= "in")

# compute fgseaRes
fRes <- lapply(1:7, function(x) get_fgseaRes(msigdbr.gs= hallmarks, seurat4.markers= out[[x]]$markers))



# rough plotting function for all cell types
fgsea_plot_all <- function(anno.index) {
        topPathwaysUp <- fRes[[anno.index]]$result %>% 
                           subset(ES > 0 & padj < 0.05) %>%
                           arrange(pval) %>%
                           top_n(10, pval) %>%
                           pull(pathway)
        topPathwaysDown <- fRes[[anno.index]]$result %>% 
                           subset(ES < 0 & padj < 0.05) %>%
                           arrange(pval) %>%
                           top_n(10,pval) %>%
                           pull(pathway)
        topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
        pdf(paste("figures/UCCPI_DE_GSEA_031221/uccpirec2_fgsea_controls_celltype",anno.index, ".pdf",sep=""),width=10,height=5)
        plotGseaTable(fRes[[anno.index]]$pathways[topPathways], fRes[[anno.index]]$ranks, fRes[[anno.index]]$result, gseaParam=0.5)
        dev.off()
}

lapply(1:7,fgsea_plot_all)













#########################

