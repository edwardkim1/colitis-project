############################################################
## Featureplots for figures
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
plan(multicore, workers = 20)
options(future.globals.maxSize = 20000 * 1024^2)
options(ggrepel.max.overlaps = Inf)
key.colitis.colors = c("#a3d2ca", "#ffcda3", "#5eaaa8","#056676","#db6400")

############################## 
# 2 - Source file 
##############################
source("scripts/Stat111functions.R")
source("scripts/new_clustering_functions.R")

############################## 
# 3 - UC+CPI all CD45+
##############################
load('data/ssuo_CD45/seurat.object.new.RData')
DimPlot(Tcell, label=FALSE)
ggsave("figures/FIGURE_REDO_031321/CPI_orig_umap.pdf",width= 7, height= 6, units= "in")

s <- readRDS("saved_objects/UC_chang_021321/uccpirec-allCD45.rds")

#s$colitis3 <- ifelse(s$colitis3== "Colitis", "CPI-Colitis", s$colitis3)
#s$colitis3 <- ifelse(s$colitis3== "Control", "CPI-Control", s$colitis3)
#s$colitis3 <- ifelse(s$colitis3== "No-Colitis", "CPI-No-Colitis", s$colitis3)
#s$colitis3 <- factor(s$colitis3, levels= c("CPI-Control", "UC-Control" , "CPI-No-Colitis", "CPI-Colitis", "UC"))
#s$orig.ident <- ifelse(s$orig.ident== "UC_chang", "UC_Boland", s$orig.ident)
# s$SingleR.Luoma2 <- as.character(s$SingleR.Luoma2)
# s$SingleR.Luoma2 <- ifelse(s$SingleR.Luoma2 == "lgA.plasma.B", "Plasma.B",s$SingleR.Luoma2 )
# s$SingleR.Luoma2 <- ifelse(s$SingleR.Luoma2 == "lgA.plasma.B", "Plasma.B",s$SingleR.Luoma2 )
# s$SingleR.Luoma2 <- factor(s$SingleR.Luoma2)
#saveRDS(s, "saved_objects/UC_chang_021321/uccpirec-allCD45.rds")

# DE_heatmap
out <- DE_heatmap(s)
out$p
ggsave("figures/FIGURE_REDO_031321/UCCPI_DE_avgexp.pdf", width=4.5,height=8, units="in")


# distribution by sample
patients <- factor(s$patient.id, levels= c(str_glue("C{1:8}"), str_glue("NC{1:6}"),str_glue("CT{1:6}"),str_glue("C{9:16}"),str_glue("U{1:7}") ))
data.frame(cluster = s$seurat_clusters, patient= patients) %>%
        ggplot() + geom_bar(
                mapping = aes(x= patient, fill= cluster),
                position = "fill"
        ) + labs(y="proportion")+ theme_classic() + coord_flip()
ggsave("figures/FIGURE_REDO_031321/UCCPI_cluster_prop_per_sample.pdf",width= 4, height= 6, units= "in")

# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/FIGURE_REDO_031321/UCCPI_VlnPlot_percentmito_by_cluster.pdf", width= 4, height= 4, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/FIGURE_REDO_031321/UCCPI_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 4, height= 4, units= "in")


# group and split umap by dataset
DimPlot(s, reduction= "umap", label=FALSE, group.by="colitis3", ncol=1, cols=key.colitis.colors)
ggsave("figures/FIGURE_REDO_031321/UCCPI_UMAP_groupby_dataset.pdf", width= 8, height= 7, units= "in")

DimPlot(s, reduction= "umap", label=FALSE, split.by="colitis3", ncol=3)
ggsave("figures/FIGURE_REDO_031321/UCCPI_UMAP_splitby_dataset.pdf", width= 5*3, height= 5*2, units= "in")


# DA figures
p <- plot.propbar.errors(annotations = s$SingleR.Luoma2,
        sample.id = s$sample.id,
        conditions = s$colitis3,
        clusters= names(table(s$SingleR.Luoma2))[c(2,3,7,1,6,5,4)],
        group.names= c("CPI-Control", "UC-Control" , "CPI-No-Colitis", "CPI-Colitis", "UC"))
p + scale_fill_manual(values=key.colitis.colors)
ggsave("figures/FIGURE_REDO_031321/uccpi_CD45_barplot.pdf", width= 10, height= 5, units= "in")

#### Wilcoxon Test (Unpaired) ####
wilcoxon.p.values <- function(clusters, df.prop) {
        p.values <- unlist(lapply(clusters,function(x) {wilcox.test(proportion ~ status, data = subset(df.prop, cluster==x), paired = FALSE)$p.value}))
        names(p.values) <- unique(df.prop$cluster)
        return(p.values)
}
# get p-values using wilcoxon test
df.s <- wilcoxon.df(s$SingleR.Luoma2, s$sample.id, s$colitis3)
s.pvals <- wilcoxon.p.values(names(table(s$SingleR.Luoma2)),subset(df.s, status=="Colitis" | status=="UC"))
s.pvals <- wilcoxon.p.values(names(table(s$SingleR.Luoma2)),subset(df.s, status=="Control" | status=="UC-Control"))

# DA analysis with qlglm
s2 <- subset(s, colitis == "Colitis" | colitis == "UC")
s2$sample.id <- droplevels(s2$sample.id)
s2$UCoverCPI <- s2$colitis == "UC"
design <- function(y) model.matrix(~factor(UCoverCPI),y)
out <- DA_analysis(s2, design, title = "FIGURE_REDO_031321/uccpiCD45_DA", annotations=s2$SingleR.Luoma2)
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
out <- DA_analysis(s2, design, title = "FIGURE_REDO_031321/uccpiB_controls_DA", annotations=s2$SingleR.Luoma2)
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)



# violin plots!
Idents(s) <- s$colitis3
pdf(paste("figures/FIGURE_REDO_031321/violinplot_cxcl.pdf", sep=""), width=4*2, height=4*2)
StackedVlnPlot(object = s[,s$colitis3 != "CPI-No-Colitis"], features = c("CXCL2","CXCL8"), split.by = "colitis3", color.use = key.colitis.colors[c(1,2,4,5)])
dev.off()


############################## 
# 4 - UC+CPI B cell subset
##############################
s <- readRDS("saved_objects/UC_chang_021321/uccpirec-Bcell_subset.rds")

# s$colitis3 <- s$colitis
# s$colitis3 <- ifelse(s$colitis3 == "Control" & s$orig.ident== "UC_Boland", "UC-Control" ,s$colitis3)
# s$colitis3 <- ifelse(s$colitis3== "Colitis", "CPI-Colitis", s$colitis3)
# s$colitis3 <- ifelse(s$colitis3== "Control", "CPI-Control", s$colitis3)
# s$colitis3 <- ifelse(s$colitis3== "No-Colitis", "CPI-No-Colitis", s$colitis3)
# s$colitis3 <- factor(s$colitis3, levels= c("CPI-Control", "UC-Control" , "CPI-No-Colitis", "CPI-Colitis", "UC"))
# s$orig.ident <- ifelse(s$orig.ident== "UC_chang", "UC_Boland", s$orig.ident)
# saveRDS(s, "saved_objects/UC_chang_021321/uccpirec-Bcell_subset.rds")

s <- RunUMAP(s, return.model=TRUE,dims=1:25)
p <- DimPlot(s, reduction="umap", label=T, raster=F)
ggsave("figures/FIGURE_REDO_031321/UCCPIB_umap25pcs.pdf" ,width= 6, height= 6)

bcell.anno <- c("0: IgA Plasma B", 
        "1: MHCII+ Memory B",
        "2: IgA Plasma B",
        "3: Naive B",
        "4: IgA Plasma B",
        "5: Mixed B",
        "6: B-T Doublet",
        "7: IgG Plasma B",
        "8: IgA Plasma B",
        "9: Cycling B",
        "10: CD21-low B")
s$bcell.anno <- s$seurat_clusters
levels(s$bcell.anno) <- bcell.anno

bcell.anno2 <- c("IgA Plasma B", 
        "Memory B",
        "IgA Plasma B",
        "Naive B",
        "IgA Plasma B",
        "Mixed B",
        "B-T Doublet",
        "IgG Plasma B",
        "IgA Plasma B",
        "Cycling B",
        "CD21-low B")
s$bcell.anno2 <- s$seurat_clusters
levels(s$bcell.anno2) <- bcell.anno2


p <- DimPlot(s, reduction="umap", group.by="bcell.anno2", label=F, raster=F)
ggsave("figures/FIGURE_REDO_031321/UCCPIB_bcellanno.pdf" ,width= 6, height= 6)

DimPlot(s, reduction= "umap", label=FALSE, group.by="colitis3", ncol=1, cols=key.colitis.colors)
ggsave("figures/FIGURE_REDO_031321/UCCPIB_groupby_colitis3.pdf" ,width= 7, height= 6)


# important IgG-IgA ratio figure



# B cells vizualization feature plots
DefaultAssay(s) <- "RNA"
FeaturePlot(s, c("IGHA1","IGHA2", "IGHD", "IGHG1", "IGHG2", "CD79A"),ncol=3)
ggsave("figures/FIGURE_REDO_031321/uccpiB2_IGfeatures.pdf", width= 16, height= 8, units= "in")
# Feature plots from the UC Chang study
DefaultAssay(s) <- "RNA"
FeaturePlot(s, c("PRDM1","XBP1", "CD79A","IGHM", "CD19", "CD3E"),ncol=3)
ggsave("figures/FIGURE_REDO_031321/uccpiB2_Bcellfeatures.pdf", width= 16, height= 8, units= "in")

#violin plot
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "CR2", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/FIGURE_REDO_031321/UCCPIB_VlnPlot_CD21orCR2.pdf", width= 4, height= 4, units= "in")
p <- VlnPlot(s,features = "CXCR4", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/FIGURE_REDO_031321/UCCPIB_VlnPlot_CXCR4.pdf", width= 4, height= 4, units= "in")

p <- VlnPlot(s,features = "CD5", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/FIGURE_REDO_031321/UCCPIB_VlnPlot_CD5.pdf", width= 4, height= 4, units= "in")

p <- VlnPlot(s,features = "AICDA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/FIGURE_REDO_031321/UCCPIB_VlnPlot_AID.pdf", width= 4, height= 4, units= "in")

# CD22 (inhibitory receptor), CD32, CCR7, FCGR2A/B,
p <- VlnPlot(s,features = "CCR7", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/FIGURE_REDO_031321/UCCPIB_VlnPlot_CCR7.pdf", width= 4, height= 4, units= "in")
p <- VlnPlot(s,features = "CD22", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/FIGURE_REDO_031321/UCCPIB_VlnPlot_CD22.pdf", width= 4, height= 4, units= "in")

p <- VlnPlot(s,features = "FCGR2A", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/FIGURE_REDO_031321/UCCPIB_VlnPlot_FCGR2A.pdf", width= 4, height= 4, units= "in")
p <- VlnPlot(s,features = "FCGR2B", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/FIGURE_REDO_031321/UCCPIB_VlnPlot_FCGR2B.pdf", width= 4, height= 4, units= "in")



# DE_heatmap
DefaultAssay(s) <- "integrated"
out <- DE_heatmap(s)
out$p
ggsave("figures/FIGURE_REDO_031321/UCCPIB_DE_avgexp.pdf", width=4,height=7, units="in")




# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/FIGURE_REDO_031321/UCCPIB_VlnPlot_percentmito_by_cluster.pdf", width= 4, height= 4, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/FIGURE_REDO_031321/UCCPIB_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 4, height= 4, units= "in")

# distribution by sample
patients <- factor(s$patient.id, levels= c(str_glue("C{1:8}"), str_glue("NC{1:6}"),str_glue("CT{1:6}"),str_glue("C{9:16}"),str_glue("U{1:7}") ))
data.frame(cluster = s$seurat_clusters, patient= patients) %>%
        ggplot() + geom_bar(
                mapping = aes(x= patient, fill= cluster),
                position = "fill"
        ) + labs(y="proportion")+ theme_classic() + coord_flip()
ggsave("figures/FIGURE_REDO_031321/UCCPIB_cluster_prop_per_sample.pdf",width= 4, height= 6, units= "in")

# DA figures
#s$colitis3 <- ifelse(s$colitis3== "Colitis", "CPI-Colitis", s$colitis3)
#s$colitis3 <- ifelse(s$colitis3== "Control", "CPI-Control", s$colitis3)
#s$colitis3 <- ifelse(s$colitis3== "No-Colitis", "CPI-No-Colitis", s$colitis3)
#s$colitis3 <- factor(s$colitis3, levels= c("CPI-Control", "UC-Control" , "CPI-No-Colitis", "CPI-Colitis", "UC"))
#s$orig.ident <- ifelse(s$orig.ident== "UC_chang", "UC_Boland", s$orig.ident)
#saveRDS(s, "saved_objects/UC_chang_021321/uccpirec-Bcell_subset.rds")

p <- plot.propbar.errors(annotations = s$bcell.anno2,
        sample.id = s$patient.id,
        conditions = s$colitis3,
        clusters= levels(s$bcell.anno2),
        group.names= c("CPI-Control", "UC-Control" , "CPI-No-Colitis", "CPI-Colitis", "UC"))
p + scale_fill_manual(values=key.colitis.colors) +
	theme(axis.text.x = element_text(angle = 0, vjust = 1)) + 
	scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
ggsave("figures/FIGURE_REDO_031321/uccpi_Bcells_barplot.pdf", width= 10, height= 5, units= "in")

#### Wilcoxon Test (Unpaired) ####
wilcoxon.p.values <- function(clusters, df.prop) {
        p.values <- unlist(lapply(clusters,function(x) {wilcox.test(proportion ~ status, data = subset(df.prop, cluster==x), paired = FALSE)$p.value}))
        names(p.values) <- unique(df.prop$cluster)
        return(p.values)
}
# get p-values using wilcoxon test
df.s <- wilcoxon.df(s$SingleR.Luoma2, s$sample.id, s$colitis3)
s.pvals <- wilcoxon.p.values(names(table(s$SingleR.Luoma2)),subset(df.s, status=="Colitis" | status=="UC"))
s.pvals <- wilcoxon.p.values(names(table(s$SingleR.Luoma2)),subset(df.s, status=="Control" | status=="UC-Control"))

# DA analysis with qlglm
s2 <- subset(s, colitis3 == "CPI-Colitis" | colitis3 == "UC")
s2$sample.id <- droplevels(s2$sample.id)
s2$UCoverCPI <- s2$colitis3 == "UC"
design <- function(y) model.matrix(~factor(UCoverCPI),y)
out <- DA_analysis(s2, design, title = "FIGURE_REDO_031321/uccpiB_DA", annotations=s2$bcell.anno)
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
out <- DA_analysis(s2, design, title = "FIGURE_REDO_031321/uccpiCD45_controls_DA", annotations=s2$SingleR.Luoma2)
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)


# CADM1 violin plots
DefaultAssay(s) <- "RNA"
p <- VlnPlot(s,features = "CADM1", pt.size= 0, group.by="colitis3") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank()) + scale_fill_manual(values=key.colitis.colors)
ggsave("figures/FIGURE_REDO_031321/UCCPIB_VlnPlot_CADM1_by_condition.pdf", width= 4, height= 4, units= "in")
p <- VlnPlot(s,features = "CADM1", pt.size= 0, group.by="bcell.anno") + theme(axis.text.x = element_text(angle = 45),legend.position = "none", axis.title.x=element_blank()) 
ggsave("figures/FIGURE_REDO_031321/UCCPIB_VlnPlot_CADM1_by_celltype.pdf", width= 4, height= 4, units= "in")




############################## 
# 5 - UC+CPI T cell subset
##############################
s <- readRDS("saved_objects/UC_chang_013021/uccpirec2.rds")

# violin plot of IGKC per cluster
DefaultAssay(s) <- "RNA"
p <- VlnPlot(s,features = "rna_IGKC", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 0),legend.position = "none", axis.title.x=element_blank())
ggsave(paste("figures/FIGURE_REDO_031321/UCCPIT_VlnPlot_IGKC_by_cluster.pdf", sep=""), width= 4, height= 4, units= "in")

# distribution by sample
patients <- factor(s$patient.id, levels= c(str_glue("C{1:8}"), str_glue("NC{1:6}"),str_glue("CT{1:6}"),str_glue("C{9:16}"),str_glue("U{1:7}") ))
data.frame(cluster = s$seurat_clusters, patient= patients) %>%
        ggplot() + geom_bar(
                mapping = aes(x= patient, fill= cluster),
                position = "fill"
        ) + labs(y="proportion")+ theme_classic() + coord_flip()
ggsave("figures/FIGURE_REDO_031321/UCCPIT_cluster_prop_per_sample.pdf",width= 4, height= 6, units= "in")

# featureplot of CD4 CD8 split
DefaultAssay(sT) <- "RNA"
FeaturePlot(sT, c("CD4","CD8A", "CD8B","CD3E"),ncol=2)
ggsave("figures/FIGURE_REDO_031321/UCCPIT_CD4xCD8_split_featureplot.pdf", width= 10, height= 10, units= "in")

###################################
# 6 - CPI PD-1 mono+combo T cell subset
###################################
s <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
p <- DE_heatmap(s)
p$p
ggsave("figures/PD1_031221/onlyTcells1550_DE_avgexp.pdf")

###################################
# 7 - CDCPI (Martin) T cells
###################################
s.a <- readRDS("saved_objects/cdcpi3_anchor.rds")
# differential gene analysis by cluster split by dataset
x <- subset(s.a, colitis=="CPI_colitis")
out1 <- DE_heatmap(x)
ggsave("figures/FIGURE_REDO_031321/cdcpi3DE_CPI_colitis.pdf",out1$p, width=4,height=7, units="in")

x <- subset(s.a, colitis=="Control")
out2 <- DE_heatmap(x)
ggsave("figures/FIGURE_REDO_031321/cdcpi3DE_CPI_control.pdf",out2$p, width=4,height=7, units="in")

x <- subset(s.a, colitis == "CD_inflamed")
out3 <- DE_heatmap(x)
ggsave("figures/FIGURE_REDO_031321/cdcpi3DE_CD_inf.pdf",out3$p, width=4,height=7, units="in")

x <- subset(s.a, colitis == "CD_uninflamed")
out4 <- DE_heatmap(x)
ggsave("figures/FIGURE_REDO_031321/cdcpi3DE_CD_uninf.pdf",out4$p, width=4,height=7, units="in")

# proportion by dataset per cluster
s.a$orig.ident <- ifelse(s.a$orig.ident == "CD3_Tcell", "CPI_luoma", s.a$orig.ident)
data.frame(dataset = s.a$orig.ident, cluster= s.a$seurat_clusters) %>%
	ggplot() + geom_bar(
		mapping = aes(x= cluster, fill= dataset),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/FIGURE_REDO_031321/cdcpi3_dataset_proportion.pdf", width= 5, height= 4, units= "in")




############################## 
# 2 - Implement Dream 
##############################
library('variancePartition')
library('edgeR')
library('BiocParallel')



out <- get_DEG_dream(s,"Myeloid")

ggsave(paste("figures/FIGURE_REDO_031321/myeloid_DEgenes.pdf", sep=""), p, width= 3.5, height= 5.5, units= "in")
saveRDS(out, "saved_objects/UCCPI_DE_GSEA_031221/myeloid_DEgenes.rds")

# subset comparison conditions
out <- get_DEG_dream(s,"Myeloid")
s <- readRDS("saved_objects/UC_chang_021321/uccpirec-allCD45.rds")
ggsave(paste("figures/FIGURE_REDO_031321/plasmaBcells_DEgenes.pdf", sep=""), p, width= 3.5, height= 5.5, units= "in")


s.markers <- topTable(out$key.markers, coef='DiseaseUC', number=10000)

markers.up <- s.markers %>%
					subset(logFC > 1 & adj.P.Val < 0.05) %>%
	                   arrange(desc(logFC)) %>%
	                   arrange(P.Value) %>%
	                   arrange(adj.P.Val) %>%
	                   head(20) %>%
	                   rownames()

s.markers.c <- topTable(out$control.markers, coef='DiseaseUC-Control', number=10000)
	# get significant markers for removal
	markers.up.c <- s.markers.c %>% 
	                   subset(logFC > 1 & adj.P.Val < 0.05) %>%
	                   arrange(desc(logFC)) %>%
	                   arrange(P.Value) %>%
	                   arrange(adj.P.Val) %>% 
	                   rownames()
	markers.down.c <- s.markers.c %>% 
	                   subset(logFC < -1 & adj.P.Val < 0.05) %>%
	                   arrange(logFC) %>%
	                   arrange(P.Value) %>%
	                   arrange(adj.P.Val) %>% 
	                   rownames()
	# filter regular with control DE
	markers.up <- s.markers %>% subset(logFC > 1 & adj.P.Val < 0.05)
	markers.up <- markers.up[rownames(markers.up) %ni% markers.up.c,] %>%
	                   arrange(desc(logFC)) %>%
	                   arrange(P.Value) %>%
	                   arrange(adj.P.Val) %>%
	                   head(20) %>%
	                   rownames()
	markers.down <- s.markers %>% subset(logFC < -1 & adj.P.Val < 0.05)
	markers.down <- markers.down[rownames(markers.down) %ni% markers.down.c,] %>%
	                   arrange(logFC) %>%
	                   arrange(P.Value) %>%
	                   arrange(adj.P.Val) %>%
	                   head(20) %>%
	                   rownames()
	top.markers.filtered <- c(markers.up, rev(markers.down))



# attempt at contrasts with myeloid

get_DEG_dream_contrasts <- function(s, celltype) {
	s3 <- s[,s$SingleR.Luoma2==celltype]

	# extract metadata
	metadata <- data.frame(Individual = s3$patient.id, Condition= s3$colitis3, Celltype= s3$SingleR.Luoma2)

	# remove genes with fewer than 10% of cells expressing it
	x <- s3[["RNA"]]@counts 
	nnz_row <- tabulate(x@i + 1)
	keep <- rownames(x)[nnz_row> ncol(s3)*0.1 ]
	print(paste("Keeping",length(keep),"number of features for analysis."))
	# exactly
	geneExpr = DGEList( x[keep,] )
	geneExpr = calcNormFactors( geneExpr )
	# Specify parallel processing parameters
	# this is used implicitly by dream() to run in parallel
	param = SnowParam(20, "SOCK", progressbar=TRUE)
	register(param)
	# The variable to be tested must be a fixed effect
	form <- ~ 0 + Condition + (1|Individual) 

	# estimate weights using linear mixed model of dream
	print("Estimating weights...")
	vobjDream = voomWithDreamWeights( geneExpr, form, metadata )
	L1 = c(1, -1,0, -1, 1)
	print("Done.")
	
	print("Applying contrasts...")
	fitmm = dream( vobjDream, form, metadata, L=L1 )
	print("Done.")
	
	print("Selecting Markers...")
	s.markers <- topTable(fitmm, coef="L1", number=10000)
	markers.up <- s.markers %>% 
	                   subset(logFC > 0.8 & adj.P.Val < 0.05) %>%
	                   arrange(desc(logFC)) %>%
	                   arrange(P.Value) %>%
	                   arrange(adj.P.Val) %>% 
	                   head(30) %>%
	                   rownames()
	markers.down <- s.markers %>% 
	                   subset(logFC < -0.8 & adj.P.Val < 0.05) %>%
	                   arrange(logFC) %>%
	                   arrange(P.Value) %>%
	                   arrange(adj.P.Val) %>% 
	                   head(30) %>%
	                   rownames()
	top.markers <- c(markers.up, rev(markers.down))
	print("Done.")
	#Idents(s3) <- s3$colitis3
	#DefaultAssay(s) <- "RNA"
	#cluster.averages <- AverageExpression(s3, return.seurat = TRUE)
	# visualize
	#print("Saving plot...")
	#p <- DoHeatmap(cluster.averages,features = c(top.markers), size = 2, draw.lines = FALSE, raster=F)
	#print("Done. Fin")
 	return(list(fitmm=fitmm))
}

out2 <- get_DEG_dream_contrasts(s,"Myeloid")

ggsave(paste("figures/FIGURE_REDO_031321/myeloidcells_DEgenes_contrast.pdf", sep=""), out2$p, width= 3.5, height= 5.5, units= "in")



# executing DE for all broad cell types
celltypes <- levels(s$SingleR.Luoma2)
attempt1 <- lapply(celltypes[-5],function(x) get_DEG_dream(s,x))
saveRDS(attempt1, "saved_objects/UCCPI_DE_GSEA_031221/allbroadcelltypes_DEgenes.rds")
ggsave(paste("figures/FIGURE_REDO_031321/CD4T_DEgenes.pdf", sep=""), attempt1[[1]]$plot, width= 3.5, height= 5.5, units= "in")
ggsave(paste("figures/FIGURE_REDO_031321/IgAplasmaB_DEgenes.pdf", sep=""), attempt1[[2]]$plot, width= 3.5, height= 5.5, units= "in")
ggsave(paste("figures/FIGURE_REDO_031321/CD8T_DEgenes.pdf", sep=""), attempt1[[3]]$plot, width= 3.5, height= 5.5, units= "in")
ggsave(paste("figures/FIGURE_REDO_031321/B_DEgenes.pdf", sep=""), attempt1[[4]]$plot, width= 3.5, height= 5.5, units= "in")

ggsave(paste("figures/FIGURE_REDO_031321/Mast_DEgenes.pdf", sep=""), attempt1[[5]]$plot, width= 3.5, height= 5.5, units= "in")
ggsave(paste("figures/FIGURE_REDO_031321/ILCs_DEgenes.pdf", sep=""), attempt1[[6]]$plot, width= 3.5, height= 5.5, units= "in")


# contrasts for all broad cell types
attempt2 <- lapply(celltypes[-5],function(x) get_DEG_dream_contrasts(s,x))


cd4 <- get_DEG_dream_contrasts(s,celltypes[1])


s.markers <- topTable(cd4, coef="L1", number=10000)
	markers.up <- s.markers %>% 
	                   subset(logFC > 0.8 & adj.P.Val < 0.05) %>%
	                   arrange(desc(logFC)) %>%
	                   arrange(P.Value) %>%
	                   arrange(adj.P.Val) %>% 
	                   head(30) %>%
	                   rownames()
	markers.down <- s.markers %>% 
	                   subset(logFC < -0.8 & adj.P.Val < 0.05) %>%
	                   arrange(logFC) %>%
	                   arrange(P.Value) %>%
	                   arrange(adj.P.Val) %>% 
	                   head(30) %>%
	                   rownames()
	top.markers <- c(markers.up, rev(markers.down))



s3 <- s[,s$SingleR.Luoma2==celltypes[1]]
Idents(s3) <- s3$colitis3
DefaultAssay(s) <- "RNA"
cluster.averages <- AverageExpression(s3, return.seurat = TRUE)

print("Saving plot...")
p <- DoHeatmap(cluster.averages,features = c(top.markers), size = 2, draw.lines = FALSE, raster=F)
print("Done. Fin")
ggsave(paste("figures/FIGURE_REDO_031321/CD4T_DEgenes_contrast.pdf", sep=""), p, width= 3.5, height= 5.5, units= "in")
