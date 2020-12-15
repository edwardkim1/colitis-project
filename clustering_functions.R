#############################################################
## Clustering Functions
##
## add_inflamed_id
## add_patient_id
## filter_mito
## remove_emptyDroplets
## qc
## cluster1
## cluster2
## label_singleR
## save_figures
##
## edward_kim@college.harvard.edu - Dec. 20, 2019
############################################################

add_inflamed_id <- function(seurat.object) {
	inflammed <-c("_1$","_3$","_6$","_7$","_10$","_12$","_14$","_16$","_18$")
	toMatch <- paste(inflammed,collapse="|")
	seurat.object$status <- ifelse(grepl(colnames(seurat.object), pattern=toMatch),"inflamed","not inflamed")
	seurat.object
}

add_patient_id <- function(s2_postQC) {
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
	s2_postQC
}

filter_mito <- function(seurat.object, save = FALSE, filename = "") {
	genename=as.matrix(rownames(seurat.object))
	rownames(genename)=rownames(seurat.object)
	colnames(genename)='Symbol'
	lib.sizes = Matrix::colSums(seurat.object)
	mito.gene=rownames(genename)[grepl(x=rownames(genename), pattern = "^MT-")]
	mt.counts = seurat.object[mito.gene,]
	mt.fraction = Matrix::colSums(mt.counts)/lib.sizes
	mt.p = pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
	mt.lim = min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 1e-2)])
	filtered.object <- seurat.object[,which(mt.fraction < mt.lim)]
	if(save==TRUE) {
		saveRDS(c(median(mt.fraction),mad(mt.fraction),mt.lim), filename)
	}
	filtered.object
}


remove_emptyDroplets <- function(dirname) {
	sce <- read10xCounts(dirname)
	my.counts <- counts(sce)
	br.out <- barcodeRanks(my.counts)
	set.seed(100)
	e.out <- emptyDrops(my.counts)

	# OUPUTS # 
	out.dir <- paste("qc_files",dirname,sep="/")
	#save the e.out for the next step
	saveRDS(e.out,paste(out.dir,"e_out.rds",sep="_"))

	#save the plot of the knee
	out.filename = paste(out.dir,"brp.pdf",sep="_")
	pdf(out.filename, width = 7, height = 7)
	plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
	o <- order(br.out$rank)
	lines(br.out$rank[o], br.out$fitted[o], col="red")
	title(main=dirname)
	abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
	abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
	legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))
	dev.off()
}

qc <- function(dirname) {
	dirname
	e.out.filename <- paste(dirname,"e_out.rds",sep="_")
	e.out <- readRDS(e.out.filename)
	s1.data <- Read10X(data.dir = paste("..",dirname,sep="/"))

	## Quality Control (remove emptydroplets)
	is.cell <- e.out$FDR <= 0.01
	is.cell <- ifelse(is.na(is.cell),F,is.cell)
	s1.qc1 <- s1.data[,is.cell]

	## Quality Control (mito hypothesis testing)
	genename=as.matrix(rownames(s1.qc1))
	rownames(genename)=rownames(s1.qc1)
	colnames(genename)='Symbol'
	lib.sizes = Matrix::colSums(s1.qc1)
	mito.gene=rownames(genename)[grepl(x=rownames(genename), pattern = "^MT-")]
	mt.counts = s1.qc1[mito.gene,]
	mt.fraction = Matrix::colSums(mt.counts)/lib.sizes
	mt.p = pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
	mt.lim = min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 1e-2)])
	s1.qc2 <- s1.qc1[,which(mt.fraction < mt.lim)]
	saveRDS(c(median(mt.fraction),mad(mt.fraction),mt.lim), paste(dirname,"mito_stats.rds",sep="_"))

	## Create Seurat Object ##
	s1 <- CreateSeuratObject(counts = s1.qc2, project = dirname, min.cells = 0, min.features = 0)
	s1 <- NormalizeData(s1, normalization.method = "LogNormalize", scale.factor = 10000)
	s1 <- FindVariableFeatures(s1, selection.method = "vst", nfeatures = 3000) ## Highly variable genes
	all.genes <- rownames(s1) ## Scale to mean 0 var 1
	s1 <- ScaleData(s1, features = all.genes)
	s1 <- RunPCA(s1, features = VariableFeatures(object = s1)) ## PCA
	s1 <- FindNeighbors(s1, dims = 1:25) ## KNN + Louvian algorithm ##
	s1 <- FindClusters(s1, resolution = 0.4)
	s1 <- RunUMAP(s1, dims = 1:25) ## UMAP ##
	## Detecting Doublets: Doublet Finder ##
	## pK identification
	sweep.res.list <- paramSweep_v3(s1, PCs = 1:25)
	sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
	s1_pK <- find.pK(sweep.stats)
	pK_max <- s1_pK[which.max(s1_pK$BCmetric),]$pK
	pK_max <- as.numeric(levels(pK_max)[as.integer(pK_max)])
	## Output pK maximization using ggplot2 as pdf
	s1_pK$group <- numeric(length(s1_pK$ParamID))
	pdf(paste(dirname,"pK.pdf",sep="_"), width = 7, height = 7)
	p <- ggplot(data=s1_pK, aes(x=pK, y=BCmetric, group=0)) +
	  		geom_line()+ scale_x_discrete(breaks=seq(0,0.3,0.05)) +
	  		geom_point() + ggtitle(dirname) + theme_classic() 
	print(p)
	dev.off()
	## Calculation
	homotypic.prop <- modelHomotypic(s1@meta.data$seurat_clusters)   ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
	nExp_poi <- round(0.1*nrow(s1@meta.data))   ## Assuming 10% doublet formation rate - tailor for your dataset
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

	s1 <- doubletFinder_v3(s1, PCs = 1:25, pN = 0.25, pK= pK_max, nExp = nExp_poi, reuse.pANN = FALSE)
	first_pANN <- names(s1@meta.data)[grep("pANN", names(s1@meta.data))]
	s1 <- doubletFinder_v3(s1, PCs = 1:25, pN = 0.25, pK = pK_max, nExp = nExp_poi.adj, reuse.pANN = first_pANN)
	DFs <- names(s1@meta.data)[grep("DF.class", names(s1@meta.data))]

	s1@meta.data[,"DF_hi.lo"] <- s1[[DFs[1]]]
	s1@meta.data$DF_hi.lo[which(s1@meta.data$DF_hi.lo == "Doublet" & s1[[DFs[2]]] == "Singlet")] <- "Doublet_lo"
	s1@meta.data$DF_hi.lo[which(s1@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
	# Output predoublet removal
	saveRDS(s1,paste(dirname,"withdoublets.rds",sep="_"))
	s1 <- subset(s1, subset = DF_hi.lo != "Doublet_hi" ) #eliminating hi probability doublets
	# Output post doublet removal
	saveRDS(s1, paste(dirname,"postQC.rds", sep="_"))
}

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
	# Adding patient ID, inflamed ID, percent mito genes
	seurat.object <- add_patient_id(seurat.object)
	seurat.object <- add_inflamed_id(seurat.object)
	seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-")
	## Scale to mean 0 var 1 and regress out patient variance
	seurat.object <- ScaleData(seurat.object, features = all.genes, vars.to.regress=c("patient.id","percent.mt"))
	return(seurat.object)
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
	return(seurat.object)
}

load_references <- function() {
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

label_singleR <- function(seurat.object) {
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

save_figures <- function(seurat.object) {
	require(gridExtra)
	filename <- paste(seurat.object,".pdf",sep="")
	pdf(filename, width = 10, height = 10)
		p.elbow <- ElbowPlot(seurat.object, ndims=40)
		p.umap <- DimPlot(seurat.object, reduction = "umap", label=T)
		p1 <- DimPlot(seurat.object, reduction = "umap", label=T, group.by="SingleR.Martin")
		p2 <- DimPlot(seurat.object, reduction = "umap", label=T, group.by="SingleR")
		p.status <- DimPlot(seurat.object, reduction = "umap", label=F, group.by="status")
		p.patient <- DimPlot(seurat.object, reduction = "umap", label=F, group.by="patient.id")

		print(p.elbow)
		print(p.umap)
		print(p.status)
		print(p.patient)
		print(grid.arrange(p1, p2, nrow = 1))
	dev.off()
}

