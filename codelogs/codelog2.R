#### Clustering redo with regress out mito genes ####
dirname <- dir("../",pattern="GSM")
cells <- mclapply(paste(dirname,"postQC.rds", sep="_"),readRDS)
#### Merging cells excluding patient 5 and 16
s2.wo5 <- merge(cells[[3]], c(cells[[4]], cells[[5]], cells[[6]], cells[[7]], cells[[8]], cells[[9]], cells[[10]], cells[[11]], cells[[12]], cells[[13]], cells[[14]], cells[[15]], cells[[16]], cells[[17]], cells[[18]], cells[[19]], cells[[20]]))
s2 <- s2.wo5[['RNA']]@counts
# Filter based on mitochondrial reads
genename=as.matrix(rownames(s2))
rownames(genename)=rownames(s2)
colnames(genename)='Symbol'
lib.sizes = Matrix::colSums(s2)
mito.gene=rownames(genename)[grepl(x=rownames(genename), pattern = "^MT-")]
mt.counts = s2[mito.gene,]
mt.fraction = Matrix::colSums(mt.counts)/lib.sizes
mt.p = pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
mt.lim = min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 1e-2)])
s2.qc1 <- s2[,which(mt.fraction < mt.lim)]
#s2[["percent.mt"]] <- PercentageFeatureSet(s2, pattern = "^MT-")
#VlnPlot(s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(c(median(mt.fraction),mad(mt.fraction),mt.lim), "s2_wo5_mito_stats.rds")
# Filter the genes the appear infrequently
s2.qc2 <- CreateSeuratObject(s2.qc1,min.cells=10)
s2.qc2$orig.ident <- s2.wo5$orig.ident
# Filter cells that have low feature count
s2_postQC <- subset(s2.qc2, subset = nFeature_RNA > 200)
# Define Functions, Lists, and Operators
'%ni%' = Negate('%in%')
all.genes <- rownames(s2_postQC)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TR.V")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
# Normalized
s2_postQC <- NormalizeData(s2_postQC, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable genes
s2_postQC <- FindVariableFeatures(s2_postQC, selection.method = "vst", nfeatures = 2000)
hvf <- VariableFeatures(s2_postQC)
VariableFeatures(s2_postQC) <- hvf[hvf %ni% c(trvgenes,igvgenes)] 
# Adding patient ID and 
p.id <- ifelse(grepl(colnames(s2_postQC),pattern="_1$|_2$"),"Patient 6","")
p.id <- ifelse(grepl(colnames(s2_postQC),pattern="_3$|_4$"),"Patient 7",p.id)
p.id <- ifelse(grepl(colnames(s2_postQC),pattern="_5$|_6$"),"Patient 8",p.id)
p.id <- ifelse(grepl(colnames(s2_postQC),pattern="_7$|_8$"),"Patient 10",p.id)
p.id <- ifelse(grepl(colnames(s2_postQC),pattern="_9$|_10$"),"Patient 11",p.id)
p.id <- ifelse(grepl(colnames(s2_postQC),pattern="_11$|_12$"),"Patient 12",p.id)
p.id <- ifelse(grepl(colnames(s2_postQC),pattern="_13$|_14$"),"Patient 13",p.id)
p.id <- ifelse(grepl(colnames(s2_postQC),pattern="_15$|_16$"),"Patient 14",p.id)
p.id <- ifelse(grepl(colnames(s2_postQC),pattern="_17$|_18$"),"Patient 15",p.id)
s2_postQC$patient.id <- p.id
s2_postQC[["percent.mt"]] <- PercentageFeatureSet(s2_postQC, pattern = "^MT-")
## Scale to mean 0 var 1
options(future.globals.maxSize = 11000 * 1024^2)
s2_postQC <- ScaleData(s2_postQC, features = all.genes,vars.to.regress=c("patient.id","percent.mt"))
## PCA
s2_postQC <- RunPCA(s2_postQC, features = VariableFeatures(s2_postQC))
npc <- 20 # 20 PCS
s2_postQC <- FindNeighbors(s2_postQC, dims = 1:npc, k.param = 100)
s2_postQC <- FindClusters(s2_postQC, resolution = 0.4, graph.name ="RNA_snn")
## UMAP ##
s2_postQC <- RunUMAP(s2_postQC, min.dist =0.50, graph="RNA_snn")
pdf("s2v2_postQC_UMAP_3.pdf", width = 10, height = 10)
DimPlot(s2_postQC, reduction = "umap", label=T)
dev.off()
## Differential Gene Analysis with parallelization ##
library(future)
plan("multiprocess", workers = 40)
s2_postQC.markers <- FindAllMarkers(s2_postQC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
top10 <- s2_postQC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

### Label Inflammed and Uninflammed ###
inflammed <-c("_1$","_3$","_6$","_7$","_10$","_12$","_14$","_16$","_18$")
toMatch <- paste(inflammed,collapse="|")
s2_postQC$status <- ifelse(grepl(colnames(s2_postQC), pattern=toMatch),"inflamed","not inflamed")


DimPlot(s2_postQC, reduction = "umap", label=F, group_by="status")
pdf("s2v2_postQC_UMAP_features_2.pdf", width = 20, height = 20)
FeaturePlot(s2_postQC, features=c('CD14','CD68','C1QA','NKG7'))
dev.off()


#### Subset CD45+ cells ####
s3 <- subset(s2_postQC,idents= c(0,1,3,10,17,19))
cd.counts <- Matrix::colSums(s3[['RNA']]@counts[c('CD3E','CD3G','CD3D','CD4','CD8A','CD8B','CD28','ICOS','CD2','CD7'),])
s3_2 <- s3[,cd.counts>0]
s3_2 <- NormalizeData(s3_2, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable genes
s3_2 <- FindVariableFeatures(s3_2, selection.method = "vst", nfeatures = 2000)
hvf <- VariableFeatures(s3_2)
'%ni%' = Negate('%in%')
all.genes <- rownames(s3_2)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TR.V")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
# Normalized
VariableFeatures(s3_2) <- hvf[hvf %ni% c(trvgenes,igvgenes)] 
top10 <- head(VariableFeatures(s3_2), 10)
# Adding patient ID
p.id <- ifelse(grepl(colnames(s3_2),pattern="_1$|_2$"),"Patient 6","")
p.id <- ifelse(grepl(colnames(s3_2),pattern="_3$|_4$"),"Patient 7",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_5$|_6$"),"Patient 8",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_7$|_8$"),"Patient 10",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_9$|_10$"),"Patient 11",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_11$|_12$"),"Patient 12",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_13$|_14$"),"Patient 13",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_15$|_16$"),"Patient 14",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_17$|_18$"),"Patient 15",p.id)
s3_2$patient.id <- p.id
s3_2[["percent.mt"]] <- PercentageFeatureSet(s3_2, pattern = "^MT-")
## Scale to mean 0 var 1 and regress out patient variance
s3_2 <- ScaleData(s3_2, features = all.genes, vars.to.regress=c("patient.id","percent.mt"))
s3_2 <- RunPCA(s3_2, features = VariableFeatures(s3_2))
npc <- 30 # 20 PCS
s3_2 <- FindNeighbors(s3_2, dims = 1:npc, k.param = 50)
s3_2 <- FindClusters(s3_2, resolution = 0.5, graph.name ="RNA_snn")
## UMAP ##
s3_2 <- RunUMAP(s3_2, min.dist =0.30, graph="RNA_snn")
DimPlot(s3_2, reduction = "umap", label=T)

### Label Inflammed and Uninflammed ###
inflammed <-c("_1$","_3$","_6$","_7$","_10$","_12$","_14$","_16$","_18$")
toMatch <- paste(inflammed,collapse="|")
s3_2$status <- ifelse(grepl(colnames(s3_2), pattern=toMatch),"inflamed","not inflamed")
colnames(s3_2)[grepl(colnames(s3_2), pattern=toMatch)]

## Differential Gene Analysis with parallelization ##
library(future)
plan("multiprocess", workers = 40)
s3_2_regress.markers <- FindAllMarkers(s3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
top10 <- s3_2_regress.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


#### Subset T cells ####
s3 <- subset(s2_postQC,idents= c(1,2,3,5,6,8,10,12,15,17,18))
cd.counts <- Matrix::colSums(s3[['RNA']]@counts[c('CD3E','CD3G','CD3D','CD4','CD8A','CD8B','CD28','ICOS','CD2','CD7'),])
s3_2 <- s3[,cd.counts>0]
s3_2 <- NormalizeData(s3_2, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable genes
s3_2 <- FindVariableFeatures(s3_2, selection.method = "vst", nfeatures = 2000)
hvf <- VariableFeatures(s3_2)
'%ni%' = Negate('%in%')
all.genes <- rownames(s3_2)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TR.V")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
# Normalized
VariableFeatures(s3_2) <- hvf[hvf %ni% c(trvgenes,igvgenes)] 
top10 <- head(VariableFeatures(s3_2), 10)
# Adding patient ID
s3_2[["percent.mt"]] <- PercentageFeatureSet(s3_2, pattern = "^MT-")
## Scale to mean 0 var 1 and regress out patient variance
s3_2 <- ScaleData(s3_2, features = all.genes, vars.to.regress=c("patient.id","percent.mt"))
s3_2 <- RunPCA(s3_2, features = VariableFeatures(s3_2))
npc <- 15 # 20 PCS
s3_2 <- FindNeighbors(s3_2, dims = 1:npc, k.param = 100)
s3_2 <- FindClusters(s3_2, resolution = 0.5, graph.name ="RNA_snn")
## UMAP ##
s3_2 <- RunUMAP(s3_2, min.dist =0.50, graph="RNA_snn")
DimPlot(s3_2, reduction = "umap", label=T)

### Label Inflammed and Uninflammed ###
inflammed <-c("_1$","_3$","_6$","_7$","_10$","_12$","_14$","_16$","_18$")
toMatch <- paste(inflammed,collapse="|")
s3_2$status <- ifelse(grepl(colnames(s3_2), pattern=toMatch),"inflamed"," not inflamed")


#### Subsetting T cells to remove ILCs ####
s3_3 <- subset(s3_2,idents=c(0,1,2,3,4,5,8,9))
s3_3 <- NormalizeData(s3_3, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable genes
s3_3 <- FindVariableFeatures(s3_3, selection.method = "vst", nfeatures = 3000)
hvf <- VariableFeatures(s3_3)
'%ni%' = Negate('%in%')
all.genes <- rownames(s3_3)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TR.V")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
# Normalized
VariableFeatures(s3_3) <- hvf[hvf %ni% c(trvgenes,igvgenes)] 
top10 <- head(VariableFeatures(s3_3), 10)
# Adding patient ID
s3_3[["percent.mt"]] <- PercentageFeatureSet(s3_3, pattern = "^MT-")
## Scale to mean 0 var 1 and regress out patient variance
s3_3 <- ScaleData(s3_3, features = all.genes, vars.to.regress=c("patient.id","percent.mt"))

s3_3 <- cluster2(s3_3, npc=12, k.param=100, resolution=0.6, min.dist=0.3, parallel=T, n.cores=40)
s3_3 <- label.SingleR(s3_3)
df.prop.s3_3 <- wilcoxon.df(s3_3)

pdf("s3v3_barplot_sem_error_9.pdf", width = 10, height = 10)
df.prop.summary <- data_summary(df.prop.s3_3, varname="proportion", groupnames=c("cluster", "status"))
p <- ggplot(df.prop.summary, aes(x=cluster, y=proportion, fill=status)) + 
	geom_bar(stat="identity", position=position_dodge()) +
	geom_errorbar(aes(ymin=proportion-sem, ymax=proportion+sem), width=.2, position=position_dodge(.9)) +
	labs(y = "Fraction of total cells") + scale_fill_brewer(palette="Paired") + theme_minimal()
p
dev.off()

s3_3.pvals <- wilcoxon.p.values(0:9,df.prop.s3_3)





#### Subset CD8 and CD4 T cells ####
sCD8 <- subset(s3_3,subset= CD8A > 1 & CD4 < 1)
sCD8 <- cluster1(sCD8, 2000, parallel=T, n.cores=40)
#sCD8 <- cluster2(sCD8, npc=15, k.param=100, resolution=0.4,min.dist=0.3, parallel=T, n.cores=40) #v1
sCD8 <- cluster2(sCD8, npc=20, k.param=100, resolution=0.5,min.dist=0.3, parallel=T, n.cores=40) #v2
sCD8 <- label.SingleR(sCD8)
df.prop.sCD8 <- wilcoxon.df(sCD8)

pdf("sCD8_barplot_sem_error.pdf", width = 10, height = 10)
df.prop.summary <- data_summary(df.prop.sCD8, varname="proportion", groupnames=c("cluster", "status"))
p <- ggplot(df.prop.summary, aes(x=cluster, y=proportion, fill=status)) + 
	geom_bar(stat="identity", position=position_dodge()) +
	geom_errorbar(aes(ymin=proportion-sem, ymax=proportion+sem), width=.2, position=position_dodge(.9)) +
	labs(y = "Fraction of total cells") + scale_fill_brewer(palette="Paired") + theme_minimal()
p
dev.off()

sCD8.pvals <- wilcoxon.p.values(c(0,1,2,3),df.prop.sCD8)


sCD4 <- subset(s3_3,subset= CD4 > 1 & CD8A < 1)
sCD4 <- cluster1(sCD4, 2000, parallel=T, n.cores=40)
sCD4 <- cluster2(sCD4, npc=20, k.param=100, resolution=0.5,min.dist=0.3, parallel=T, n.cores=40)
sCD4 <- label.SingleR(sCD4)
df.prop.sCD4 <- wilcoxon.df(sCD4)

pdf("sCD4_barplot_sem_error.pdf", width = 10, height = 10)
df.prop.summary <- data_summary(df.prop.sCD4, varname="proportion", groupnames=c("cluster", "status"))
p <- ggplot(df.prop.summary, aes(x=cluster, y=proportion, fill=status)) + 
	geom_bar(stat="identity", position=position_dodge()) +
	geom_errorbar(aes(ymin=proportion-sem, ymax=proportion+sem), width=.2, position=position_dodge(.9)) +
	labs(y = "Fraction of total cells") + scale_fill_brewer(palette="Paired") + theme_minimal()
p
dev.off()

sCD4.pvals <- wilcoxon.p.values(c(0,1,2),df.prop.sCD4)


#### Save PDFs ####
pdf("sCD8_PCelbow_2000hvf.pdf", width = 10, height = 10)
ElbowPlot(sCD8, ndims=40)
dev.off()

pdf("s3v3_postQC_UMAP_9.pdf", width = 10, height = 10)
DimPlot(s3_3, reduction = "umap", label=T)
dev.off()

pdf("s3v3_postQC_UMAP_9_gbStatus.pdf", width = 10, height = 10)
DimPlot(s3_3, reduction = "umap", label=F, group.by="status")
dev.off()

pdf("s3v3_postQC_UMAP_9_gbPatient.pdf", width = 10, height = 10)
DimPlot(s3_3, reduction = "umap", label=F, group.by="patient.id")
dev.off()

pdf("sCD8_singler_Luoma.pdf", width = 10, height = 10)
DimPlot(sCD8, reduction = "umap", label=T, group.by="SingleR")
dev.off()

library(gridExtra)
pdf("sCD4_singler.pdf", width = 20, height = 10)
p1 <- DimPlot(sCD4, reduction = "umap", label=T, group.by="SingleR.Martin")
p2 <- DimPlot(sCD4, reduction = "umap", label=T, group.by="SingleR.Luoma")
grid.arrange(p1, p2, nrow = 1)
dev.off()

pdf("sCD8_singler.pdf", width = 20, height = 10)
p1 <- DimPlot(sCD8, reduction = "umap", label=T, group.by="SingleR.Martin")
p2 <- DimPlot(sCD8, reduction = "umap", label=T, group.by="SingleR")
grid.arrange(p1, p2, nrow = 1)
dev.off()

pdf("s3v3_singler.pdf", width = 20, height = 10)
p1 <- DimPlot(s3_3, reduction = "umap", label=T, group.by="SingleR.Martin")
p2 <- DimPlot(s3_3, reduction = "umap", label=T, group.by="SingleR.Luoma")
grid.arrange(p1, p2, nrow = 1)
dev.off()

pdf("sCD8_postQC_UMAP_2_gbStatus.pdf", width = 10, height = 10)
DimPlot(sCD8, reduction = "umap", label=F, group.by="status")
dev.off()

pdf("sCD4_postQC_UMAP_gbStatus.pdf", width = 10, height = 10)
DimPlot(sCD4, reduction = "umap", label=F, group.by="status")
dev.off()

pdf("sCD8_postQC_UMAP_2_gbPatient.pdf", width = 10, height = 10)
DimPlot(sCD8, reduction = "umap", label=F, group.by="patient.id")
dev.off()

pdf("sCD4_postQC_UMAP_2_gbPatient.pdf", width = 10, height = 10)
DimPlot(sCD4, reduction = "umap", label=F, group.by="patient.id")
dev.off()

#### functions ####
cluster1 <- function(seurat.object, n.var.features, parallel, n.cores) {
	if (parallel == TRUE) {
		require(future)
		plan("multiprocess", workers = n.cores)
	}
	seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
	## Highly variable genes
	seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = n.var.features)
	hvf <- VariableFeatures(seurat.object)
	'%ni%' = Negate('%in%')
	all.genes <- rownames(seurat.object)
	trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TR.V")]
	igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
	# Normalized
	VariableFeatures(seurat.object) <- hvf[hvf %ni% c(trvgenes,igvgenes)] 
	# Adding patient ID
	seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-")
	## Scale to mean 0 var 1 and regress out patient variance
	seurat.object <- ScaleData(seurat.object, features = all.genes, vars.to.regress=c("patient.id","percent.mt"))
	seurat.object
}

cluster2 <- function(seurat.object, npc, k.param, resolution, min.dist, parallel, n.cores) {
	if (parallel == TRUE) {
		require(future)
		plan("multiprocess", workers = n.cores)
	}
	seurat.object <- RunPCA(seurat.object, features = VariableFeatures(seurat.object))
	seurat.object <- FindNeighbors(seurat.object, dims = 1:npc, k.param = k.param)
	seurat.object <- FindClusters(seurat.object, resolution = resolution, graph.name ="RNA_snn")
	seurat.object <- RunUMAP(seurat.object, min.dist = min.dist, graph="RNA_snn")
	seurat.object
}

load.references <- function() {
	# Prepare Martin 2019 reference
	martin <- read.csv(file="martin_categories.csv", header=TRUE,sep=",")
	df.martin <- as.data.frame(martin)
	filtered.martin <- df.martin[which(df.martin$X != "28-Feb" & df.martin$X != "1-Mar"),] # need to remove the rows with "28-Feb" "1-Mar"
	rownames(filtered.martin) <- filtered.martin$X
	martin <- as.matrix(filtered.martin)
	martin <- martin[, -1]
	class(martin) <- "numeric"
	# Prepare Luoma reference
	load('../../CD3/colitis/seurat.object.RData')
	Tcell$names <- Tcell$seurat_clusters
	levels(Tcell$names) <- c("Type 1 cytokines Trm", "CD8 TRM", "IEL: CD8/gd", "Cytotoxic Trm", "LP Trm", "Type 3 cytokines Trm", "Treg", "Naive/CM", "TFH","Type 1 cytokines effector", "Cycling", "5", "MAIT")
}

label.SingleR <- function(seurat.object) {
 	sr.out.Luoma.cluster <- SingleR(method='cluster',
		sc_data=seurat.object[['RNA']]@data, 
		clusters= seurat.object$seurat_clusters,
		ref_data=as.matrix(Tcell[['RNA']]@data), 
		types=Tcell$names,
		numCores=40)
 	seurat.object$SingleR.Luoma <- seurat.object$seurat_clusters
	levels(seurat.object$SingleR.Luoma) <- sr.out.Luoma.cluster$labels

	sr.out.Martin.cluster <- SingleR(method='cluster',
		sc_data=seurat.object[['RNA']]@data, 
		clusters= seurat.object$seurat_clusters,
		ref_data=martin, 
		types=colnames(martin),
		numCores=40)
	seurat.object$SingleR.Martin <- seurat.object$seurat_clusters
	levels(seurat.object$SingleR.Martin) <- sr.out.Martin.cluster$labels

	return(seurat.object)
}


## wilcoxon functions ##
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE), sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}
wilcoxon.df <- function(s3_3) {
	viz.data <- prop.table(table(Idents(s3_3),s3_3$orig.ident), margin = 2)

	inflamed <- c(1,3,6,7,10,12,14,16,18)
	status <- rep("not inflamed", ncol(viz.data))
	status[inflamed] <- "inflamed"

	df.prop <- data.frame(cluster= rep(rownames(viz.data),ncol(viz.data)), proportion= as.vector(viz.data), status= rep(status, each=nrow(viz.data)))
	return(df.prop)
}

#### Wilcoxon Test (Paired) ####
#res <- wilcox.test(proportion ~ status, data = subset(df.prop, cluster==0), paired = TRUE)
wilcoxon.p.values <- function(clusters, df.prop) {
	p.values <- unlist(lapply(clusters,function(x) {wilcox.test(proportion ~ status, data = subset(df.prop, cluster==x), paired = TRUE)$p.value}))
	names(p.values) <- unique(df.prop$cluster)
	return(p.values)
}


## Differential Gene Analysis with parallelization ##
library(future)
plan("multiprocess", workers = 40)
sCD8.regress.markers <- FindAllMarkers(sCD8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
top10 <- sCD8.regress.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

pdf("s3v3_heatmap9.pdf", width = 30, height = 30)
DoHeatmap(s3_3, features = top10$gene) + NoLegend()
dev.off()


sCD4.regress.markers <- FindAllMarkers(sCD4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
top10 <- sCD4.regress.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


