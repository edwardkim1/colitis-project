#############################################################
## Updated Clustering Functions
##
## filter_stats # has to be done before scaling/normalizing
## filter_by_libsize # has to be done before scaling/normalizing
## filter_by_numfeats # has to be done before scaling/normalizing
## filter_by_mito# has to be done before scaling/normalizing
## find.pK.noPlot # 
## label_doublets # 
## cluster_umap #
##
## process_martin #
## add_inflamed_id #
## process_ssuo_Tcells #
## label_singleR #
##
## VariableFeaturePlot.Tcells # 
## VariableFeaturePlot.Tcells.SCT # 
## save_UMAP_figures #
##
## test_match_order #
## DE_heatmap
## DA_analysis
##
## filter_stats_CD # has to be done before scaling/normalizing
## qc_CD # 
## merge_CD # 
## save_figures_CD # 
## get_info_CD #
##
## edward_kim@college.harvard.edu - last edit: Oct. 7, 2020
############################################################
require(Seurat)

filter_stats_CD <- function(seurat.object, e.out, save = FALSE, filename = "") {
	## removing emptydroplets
	is.cell <- e.out$FDR <= 0.01
	is.cell <- ifelse(is.na(is.cell),F,is.cell)
	seurat.object <- seurat.object[,is.cell]
	# first, threshold by mitochondrial fraction
	# max threshold (20%)
	mt.fraction <- seurat.object$percent.mt
	mt.p = pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
	mt.lim1 = min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 1e-2)])
	mt.lim = min(20,mt.lim1)
	# thresholding on natural log of number of unique features; after removing dead cells
	# min threshold (log(200))
	lnf <- log(seurat.object$nFeature_RNA[which(mt.fraction < mt.lim)])
	lnf.p = pnorm(lnf, mean = median(lnf), sd = mad(lnf), lower.tail = TRUE)
	lnf.lim1 = max(lnf[which(p.adjust(lnf.p, method = "fdr") < 1e-2)])
	lnf.lim = max(log(200),lnf.lim1)
	# output
	filter.stats = list(is.cell = is.cell, mt.pre = mt.fraction, mt.post = mt.fraction[which(mt.fraction < mt.lim)], lnf.pre = lnf, lnf.post = lnf[which(lnf > lnf.lim)], cells.to.remove= !(log(seurat.object$nFeature_RNA) > lnf.lim & mt.fraction < mt.lim), mt.remove = sum(!(mt.fraction < mt.lim)) , mt.median = median(mt.fraction), mt.mad =mad(mt.fraction), mt.lim = mt.lim, lnf.remove = sum(!(lnf > lnf.lim)) , lnf.median = median(lnf), lnf.mad =mad(lnf), lnf.lim = lnf.lim, cells.to.remove.count = sum(!(log(seurat.object$nFeature_RNA) > lnf.lim & mt.fraction < mt.lim)))
	if(save==TRUE) {
		saveRDS(filter.stats, filename)
	}

	return(filter.stats)
}

qc_CD <- function(dirname, date) {
	temp <- Read10X(data.dir = paste("data/CD_martin/",dirname,sep=""))
	s <- CreateSeuratObject(counts = temp, project = dirname, min.cells = 0, min.features = 0) %>% PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")
	e.out <- readRDS(paste("saved_objects/CD_martin_qc/",dirname,"_e_out.rds",sep=""))
	stats <- filter_stats_CD(s, e.out= e.out, save=T, filename=paste("saved_objects/CD_martin_qc_", date, "/", dirname, "_filterstats.RDS", sep=""))
	# remove empty drops
	s <- s[,stats$is.cell]
	# remove low quality cells
	s$cells.to.remove <- stats$cells.to.remove
	s <- s[,!s$cells.to.remove]
	s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
	## Highly variable genes
	hvf <- VariableFeatures(s)
	'%ni%' = Negate('%in%')
	all.genes <- rownames(s)
	trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
	igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
	VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
	# Do the rest
	nPCs <- 20
	s <- ScaleData(s, features = all.genes, vars.to.regress="percent.mt") %>% RunPCA(features=VariableFeatures(s)) %>% RunUMAP(dims = 1:nPCs) %>% FindNeighbors(dims = 1:nPCs) %>% FindClusters(resolution = 0.5)
	s <- label_doublets(s,nPCs=nPCs, save.pK.figure=T, filename= paste("figures/CD_martin_", date, "/" ,dirname,"_pK_ggplot.pdf", sep=""))
	s1 <- subset(s, subset = DF_hi.lo == "Doublet_hi")
	saveRDS(s1,paste("saved_objects/CD_martin_qc_", date, "/", dirname, "_2star.RDS", sep=""))
	s <- subset(s, subset = DF_hi.lo != "Doublet_hi")
	saveRDS(s,paste("saved_objects/CD_martin_qc_", date, "/", dirname, "_3.RDS", sep=""))
}

merge_CD <- function(s.object1,s.object2) {
	require(Seurat)
	s <- merge(s.object1, s.object2)
	# extract and store Variable Features
	s <- FindVariableFeatures(s,selection.method = "vst", nfeatures = 2000)
	hvf <- VariableFeatures(s)
	'%ni%' = Negate('%in%')
	all.genes <- rownames(s)
	trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
	igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
	VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
	# extract and store PCA
	temp <- Embeddings(s.object1, reduction="pca")
	temp2 <- Embeddings(s.object2, reduction="pca")
	s[["pca"]]<- CreateDimReducObject(embeddings = rbind(temp,temp2), key = "PC_", assay = DefaultAssay(s.object1))	
	# extract and store UMAP
	temp <- Embeddings(s.object1, reduction="umap")
	temp2 <- Embeddings(s.object2, reduction="umap")
	s[["umap"]]<- CreateDimReducObject(embeddings = rbind(temp,temp2), key = "UMAP_", assay = DefaultAssay(s.object1))
	return(s)
}

save_figures_CD <- function(dirname, date) {
	## import relevant objects
	temp <- Read10X(data.dir = paste("data/CD_martin/",dirname,sep=""))
	s <- CreateSeuratObject(counts = temp, project = dirname, min.cells = 0, min.features = 0) %>% PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")
	stats <- readRDS(paste("saved_objects/CD_martin_qc_", date, "/", dirname, "_filterstats.RDS", sep=""))
	s1 <- s[,stats$is.cell]
	s1$cells.to.remove <- stats$cells.to.remove
	s2 <- readRDS(paste("saved_objects/CD_martin_qc_", date, "/", dirname, "_2star.RDS", sep=""))
	s3 <- readRDS(paste("saved_objects/CD_martin_qc_", date, "/", dirname, "_3.RDS", sep=""))
	# Mito% Pre/Post Violin Plot
	p <- ggplot(data=data.frame(x= c(rep("pre mito%", length(stats$mt.pre)),rep("post mito%", length(stats$mt.post))), percent.mito = c(stats$mt.pre,stats$mt.post))) + geom_violin(aes(x= x, y=percent.mito))
	ggsave(paste("figures/CD_martin_", date, "/",dirname,"_p101519_CD3_MITOviolin.pdf", sep=""))
	# Log(No. of Features) Pre/Post Violin Plot
	p <- ggplot(data=data.frame(x= c(rep("pre lnf", length(stats$lnf.pre)),rep("post lnf", length(stats$lnf.post))), log.num.feats = c(stats$lnf.pre,stats$lnf.post))) + geom_violin(aes(x= x, y=log.num.feats))
	ggsave(paste("figures/CD_martin_", date, "/",dirname,"_LNFviolin.pdf", sep=""))

	# Log.no.features vs. mito% Figure
	p <- s1 %>% FeatureScatter(feature1="nFeature_RNA", feature2="percent.mt", group.by="cells.to.remove")
	ggsave(paste("figures/CD_martin_", date, "/",dirname,"_LNFxMITO.pdf", sep=""))

	# Mito% Model Fitting Figure
	p <- plot_kde(stats$mt.pre, kernel = "gaussian", bw=0.5, lab.x = "mitochondrial gene percentage (per cell)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$mt.median,stats$mt.mad)}, 
	support=seq(0,100,0.1))
	p + geom_vline(xintercept= stats$mt.lim, color = "red") +
	geom_vline(xintercept= stats$mt.median-3*stats$mt.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$mt.median+3*stats$mt.mad, linetype = "dotted")
	ggsave(paste("figures/CD_martin_", date, "/",dirname,"_MITO-KDE-Normal.pdf", sep=""))
	# Log.no.features Model Fitting Figure
	p <- plot_kde(stats$lnf.pre, kernel = "gaussian", bw=0.01, lab.x = "Log(No. of Features)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$lnf.median,stats$lnf.mad)})
	p + geom_vline(xintercept= stats$lnf.lim, color = "red") +
	geom_vline(xintercept= stats$lnf.median-3*stats$lnf.mad, linetype = "dotted") +
	geom_vline(xintercept= stats$lnf.median+3*stats$lnf.mad, linetype = "dotted")
	ggsave(paste("figures/CD_martin_", date, "/",dirname,"_LNF-KDE-Normal.pdf", sep=""))
	# Variable Features Save & Figure
	VariableFeatures(s3)%>%saveRDS(paste("saved_objects/CD_martin_qc_", date, "/",dirname,"_hvf.RDS", sep=""))
	p <- VariableFeaturePlot.Tcells(s3) %>% LabelPoints(points=head(VariableFeatures(s3),10), repel= TRUE)
	ggsave(paste("figures/CD_martin_", date, "/",dirname,"_VariableFeatures.pdf", sep=""))
	# PCA Elbow Plot Figure
	p <- ElbowPlot(s3, ndims= 30)
	ggsave(paste("figures/CD_martin_", date, "/",dirname,"_elbowplot.pdf", sep=""))
	# Clustering and Doublet Detection Figure
	s4 <- merge_CD(s3, s2)
	p <- DimPlot(s4)
	q <- DimPlot(s4, reduction = "umap", group.by="DF_hi.lo")+ scale_colour_manual(values=c("red","yellow","gray"))
	out <- p + q
	ggsave(paste("figures/CD_martin_", date, "/",dirname,"_UMAP_doublets.pdf", sep=""), width= 12, height= 6, units= "in")
}

get_info_CD <- function(dirname,date) {
	stats <- readRDS(paste("saved_objects/CD_martin_qc_", date, "/", dirname, "_filterstats.RDS", sep=""))
	s2 <- readRDS(paste("saved_objects/CD_martin_qc_", date, "/", dirname, "_2star.RDS", sep=""))
	s3 <- readRDS(paste("saved_objects/CD_martin_qc_", date, "/", dirname, "_3.RDS", sep=""))
	x <- data.frame(
	"Orig.count" = length(stats$is.cell),
	"is.cell.fraction" = mean(stats$is.cell),
	"is.cell.counts" = sum(stats$is.cell),
	"Mito.U.lnf.removed" = stats$cells.to.remove.count,
	"Mito.upper.thres" = stats$mt.lim,
	"Feats.lower.thres" = exp(stats$lnf.lim),
	"Doublet_hi" = dim(s2)[2],
	"Doublet_lo" = table(s3$DF_hi.lo)[["Doublet_lo"]],
	"Singlet" = table(s3$DF_hi.lo)[["Singlet"]],
	"Final_count" = dim(s3)[2]
	)
	rownames(x) <- dirname
	return(x)
}

filter_stats <- function(seurat.object, save = FALSE, filename = "") {
	#thresholding number of unique features
	lnf <- log(seurat.object$nFeature_RNA)
	lnf.p = pnorm(lnf, mean = median(lnf), sd = mad(lnf), lower.tail = TRUE)
	lnf.lim = max(lnf[which(p.adjust(lnf.p, method = "fdr") < 1e-2)])
	#thresholding by mitochondrial fraction
	mt.fraction <- seurat.object$percent.mt
	mt.p = pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
	mt.lim = min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 1e-2)])
	# output
	filter.stats = list(mt.pre = mt.fraction, mt.post = mt.fraction[which(mt.fraction < mt.lim)], lnf.pre = lnf, lnf.post = lnf[which(lnf > lnf.lim)], cells.to.remove= !(lnf > lnf.lim & mt.fraction < mt.lim), mt.remove = sum(!(mt.fraction < mt.lim)) , mt.median = median(mt.fraction), mt.mad =mad(mt.fraction), mt.lim = mt.lim, lnf.remove = sum(!(lnf > lnf.lim)) , lnf.median = median(lnf), lnf.mad =mad(lnf), lnf.lim = lnf.lim, total.remove = sum(!(lnf > lnf.lim & mt.fraction < mt.lim)))
	if(save==TRUE) {
		saveRDS(filter.stats, filename)
	}

	return(filter.stats)
}

filter_by_libsize <- function(seurat.object, save = FALSE, filename = "") {
	lib.sizes = Matrix::colSums(seurat.object)
	lls <- log(lib.sizes)
	lls.p = pnorm(lls, mean = median(lls), sd = mad(lls), lower.tail = TRUE)
	lls.lim = max(lls[which(p.adjust(lls.p, method = "fdr") < 1e-2)])

	filtered.object <- seurat.object[,which(lls > lls.lim)]
	new.object = list(filtered.counts.mat = filtered.object, pre.lls = lls, post.lls = lls[which(lls > lls.lim)], total.removed = sum(!(lls > lls.lim)) , median = median(lls), mad =mad(lls), lls.lim = lls.lim)
	if(save==TRUE) {
		saveRDS(new.object, filename)
	}

	return(new.object)
}

filter_by_numfeats <- function(seurat.object, save = FALSE, filename = "") {
	n.feats = diff(seurat.object@p)
	lnf <- log(n.feats)
	lnf.p = pnorm(lnf, mean = median(lnf), sd = mad(lnf), lower.tail = TRUE)
	lnf.lim = max(lnf[which(p.adjust(lnf.p, method = "fdr") < 1e-2)])

	filtered.object <- seurat.object[,which(lnf > lnf.lim)]
	new.object = list(filtered.counts.mat = filtered.object, pre.lnf = lnf, post.lnf = lnf[which(lnf > lnf.lim)], total.removed = sum(!(lnf > lnf.lim)) , median = median(lnf), mad =mad(lnf), lnf.lim = lnf.lim)
	if(save==TRUE) {
		saveRDS(new.object, filename)
	}

	return(new.object)
}

filter_by_mito <- function(seurat.object, save = FALSE, filename = "") {
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
	new.object = list(filtered.counts.mat = filtered.object, pre.mt.fraction = mt.fraction, post.mt.fraction= mt.fraction[which(mt.fraction < mt.lim)], total.removed = sum(!(mt.fraction < mt.lim)) , median = median(mt.fraction), mad =mad(mt.fraction), mt.lim = mt.lim)
	if(save==TRUE) {
		saveRDS(new.object, filename)
	}

	return(new.object)
}

find.pK.noPlot <- function(sweep.stats) {
    "%ni%" <- Negate("%in%")
    if ("AUC" %ni% colnames(sweep.stats) == TRUE) {
        bc.mvn <- as.data.frame(matrix(0L, nrow = length(unique(sweep.stats$pK)), ncol = 5))
        colnames(bc.mvn) <- c("ParamID", "pK", "MeanBC", "VarBC", "BCmetric")
        bc.mvn$pK <- unique(sweep.stats$pK)
        bc.mvn$ParamID <- 1:nrow(bc.mvn)
        x <- 0
        for (i in unique(bc.mvn$pK)) {
            x <- x + 1
            ind <- which(sweep.stats$pK == i)
            bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
            bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
            bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind,"BCreal"])^2)
        }
        return(bc.mvn)
    }
    if ("AUC" %in% colnames(sweep.stats) == TRUE) {
        bc.mvn <- as.data.frame(matrix(0L, nrow = length(unique(sweep.stats$pK)), ncol = 6))
        colnames(bc.mvn) <- c("ParamID", "pK", "MeanAUC", "MeanBC",
            "VarBC", "BCmetric")
        bc.mvn$pK <- unique(sweep.stats$pK)
        bc.mvn$ParamID <- 1:nrow(bc.mvn)
        x <- 0
        for (i in unique(bc.mvn$pK)) {
            x <- x + 1
            ind <- which(sweep.stats$pK == i)
            bc.mvn$MeanAUC[x] <- mean(sweep.stats[ind, "AUC"])
            bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
            bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
            bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind,"BCreal"])^2)
        }
        return(bc.mvn)
    }
}

label_doublets <- function(seurat.object, nPCs= 30, doublet.formation.rate= 0.1, save.pK.figure = F, filename = "") {
	## Create Seurat Object ##
	require(DoubletFinder)
	require(dplyr)
	## pK identification
	seurat.object_pK <- paramSweep_v3(s, PCs = 1:nPCs) %>% summarizeSweep(GT = FALSE) %>% find.pK.noPlot()
	pK_max <- seurat.object_pK[which.max(seurat.object_pK$BCmetric),]$pK
	pK_max <- as.numeric(levels(pK_max)[as.integer(pK_max)])
	## Output pK maximization using ggplot2 as pdf
	seurat.object_pK$group <- numeric(length(seurat.object_pK$ParamID))
	if(save.pK.figure==TRUE) {
		p <- ggplot(data=seurat.object_pK, aes(x=pK, y=BCmetric, group=0)) +
	  		geom_line()+ scale_x_discrete(breaks=seq(0,0.3,0.05)) +
	  		geom_point() + theme_classic() 
		ggsave(filename)
	}
	## Computation of parameters
	homotypic.prop <- modelHomotypic(seurat.object@meta.data$seurat_clusters)
	nExp_poi <- round(doublet.formation.rate *nrow(seurat.object@meta.data)) #tailor doublet.formation.rate to your dataset
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	## Run Doublet Finder algorithm
	seurat.object <- doubletFinder_v3(seurat.object, PCs = 1:nPCs, pN = 0.25, pK= pK_max, nExp = nExp_poi, reuse.pANN = FALSE)
	first_pANN <- names(seurat.object@meta.data)[grep("pANN", names(seurat.object@meta.data))]
	seurat.object <- doubletFinder_v3(seurat.object, PCs = 1:nPCs, pN = 0.25, pK = pK_max, nExp = nExp_poi.adj, reuse.pANN = first_pANN)
	DFs <- names(seurat.object@meta.data)[grep("DF.class", names(seurat.object@meta.data))]
	## Saving DF status in the seurat.object
	seurat.object@meta.data[,"DF_hi.lo"] <- seurat.object[[DFs[1]]]
	seurat.object@meta.data$DF_hi.lo[which(seurat.object@meta.data$DF_hi.lo == "Doublet" & seurat.object[[DFs[2]]] == "Singlet")] <- "Doublet_lo"
	seurat.object@meta.data$DF_hi.lo[which(seurat.object@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
	## Returning the labeled object
	return(seurat.object)
}

VariableFeaturePlot.Tcells <- function(object, cols = c("black", "red"), pt.size = 1, log = NULL, selection.method = NULL, assay = NULL) {
    if (length(x = cols) != 2) {
        stop("'cols' must be of length 2")
    }
    hvf.info <- HVFInfo(object = object, assay = assay, selection.method = selection.method, status = TRUE)
    ### The following code was added by Edward Kim 6-30-2020 ###
    # Removes TCR variable genes and IG variable genes 
    '%ni%' = Negate('%in%')
	all.genes <- rownames(hvf.info) # do after importing counts matrix
	trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
	igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
	hvf.info <- hvf.info[all.genes %ni% c(trvgenes,igvgenes),]
	### End ###
    var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) + 1]
    hvf.info <- hvf.info[, c(1, 3)]
    axis.labels <- switch(EXPR = colnames(x = hvf.info)[2], variance.standarized = c("Average Expression",
        "Standardized Variance"), dispersion.scaled = c("Average Expression", "Dispersion"), residual_variance = c("Geometric Mean of Expression",
        "Residual Variance"))
    log <- log %||% (any(c("variance.standardized", "residual_variance") %in% colnames(x = hvf.info)))
    plot <- SingleCorPlot(data = hvf.info, col.by = var.status, pt.size = pt.size)
    plot <- plot + labs(title = NULL, x = axis.labels[1], y = axis.labels[2]) + scale_color_manual(labels = paste(c("Non-variable", "Variable"), "count:", table(var.status)), values = cols)
	if (log) {plot <- plot + scale_x_log10()}
    return(plot)
}
environment(VariableFeaturePlot.Tcells) <- environment(ElbowPlot)

VariableFeaturePlot.Tcells.SCT <- function(object, cols = c("black", "red"), pt.size = 1, log = NULL, selection.method = NULL, assay = NULL) {
    if (length(x = cols) != 2) {
        stop("'cols' must be of length 2")
    }
    hvf.info <- HVFInfo(object = object, assay = assay, selection.method = selection.method, status = TRUE)
    ### The following code was added by Edward Kim 6-30-2020 ###
    # Removes TCR variable genes and IG variable genes 
    hvf.info$variable <- rownames(hvf.info) %in% VariableFeatures(object)
    '%ni%' = Negate('%in%')
	all.genes <- rownames(hvf.info) # do after importing counts matrix
	trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
	igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
	hvf.info <- hvf.info[all.genes %ni% c(trvgenes,igvgenes),]
	### End ###
    var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) + 1]
    hvf.info <- hvf.info[, c(1, 3)]
    axis.labels <- switch(EXPR = colnames(x = hvf.info)[2], variance.standarized = c("Average Expression",
        "Standardized Variance"), dispersion.scaled = c("Average Expression", "Dispersion"), residual_variance = c("Geometric Mean of Expression",
        "Residual Variance"))
    log <- log %||% (any(c("variance.standardized", "residual_variance") %in% colnames(x = hvf.info)))
    plot <- SingleCorPlot(data = hvf.info, col.by = var.status, pt.size = pt.size)
    plot <- plot + labs(title = NULL, x = axis.labels[1], y = axis.labels[2]) + scale_color_manual(labels = paste(c("Non-variable", "Variable"), "count:", table(var.status)), values = cols)
	if (log) {plot <- plot + scale_x_log10()}
    return(plot)
}
environment(VariableFeaturePlot.Tcells.SCT) <- environment(ElbowPlot)

cluster_umap <- function(seurat.object, npc, k.param, resolution, min.dist) {
	require(future)
	plan("multicore", workers = 10)
	seurat.object <- FindNeighbors(seurat.object, dims = 1:npc, k.param = k.param) %>% FindClusters(resolution = resolution, graph.name ="RNA_snn") %>% RunUMAP(min.dist = min.dist, graph="RNA_snn")
	return(seurat.object)
}

process_martin <- function(martin) {
	# Prepare Martin 2019 reference
	df.martin <- as.data.frame(martin)
	filtered.martin <- df.martin[which(df.martin$X != "28-Feb" & df.martin$X != "1-Mar"),] # need to remove the rows with "28-Feb" "1-Mar"
	rownames(filtered.martin) <- filtered.martin$X
	martin <- as.matrix(filtered.martin)
	martin <- martin[, -1]
	class(martin) <- "numeric"
	return(martin)
}

process_ssuo_Tcells <- function(Tcell) {
	# Prepare Luoma reference
	Tcell$names <- Tcell$seurat_clusters
	levels(Tcell$names) <- c("Type 1 cytokines Trm", "CD8 TRM", "IEL: CD8/gd", "Cytotoxic effector", "LP Trm", "Type 3 cytokines Trm", "Treg", "Naive/CM", "TFH","Th1 effector", "Cycling", "5", "MAIT")
	return(Tcell)
}

label_singleR <- function(seurat.object, save.diagnostics=F, save.name=NULL) {
 	require(BiocParallel)
 	pred.CPI <- SingleR(method='cluster',
		test=seurat.object[['RNA']]@data, 
		clusters= seurat.object$seurat_clusters,
		ref=as.matrix(Tcell[['RNA']]@data), 
		labels=Tcell$names,
		de.method="wilcox",
		BPPARAM=MulticoreParam())
 	seurat.object$SingleR.Luoma <- seurat.object$seurat_clusters
	levels(seurat.object$SingleR.Luoma) <- pred.CPI$labels

	pred.Martin <- SingleR(method='cluster',
		test=seurat.object[['RNA']]@data, 
		clusters= seurat.object$seurat_clusters,
		ref=martin, 
		labels=colnames(martin),
		BPPARAM=MulticoreParam())
	seurat.object$SingleR.Martin <- seurat.object$seurat_clusters
	levels(seurat.object$SingleR.Martin) <- pred.Martin$labels

	# ref <- MonacoImmuneData()
	# pred.Monaco <- SingleR(method="cluster",
	# 	test=seurat.object[['RNA']]@data,
	# 	clusters= seurat.object$seurat_clusters, 
	# 	ref=ref, 
	# 	labels=ref$label.fine, 
	# 	BPPARAM=MulticoreParam())
	# seurat.object$SingleR.Monaco <- seurat.object$seurat_clusters
	# levels(seurat.object$SingleR.Monaco) <- pred.Monaco$labels

	if(save.diagnostics) {
		saveRDS(list(pred.CPI, pred.Martin),save.name)
	}

	return(seurat.object)
}

#add_labels() <- function(seurat.object, labels) {


#}

add_inflamed_id <- function(seurat.object) {
	inflammed <-c("_1$","_3$","_6$","_7$","_10$","_12$","_14$","_16$","_18$")
	toMatch <- paste(inflammed,collapse="|")
	seurat.object$status <- ifelse(grepl(colnames(seurat.object), pattern=toMatch),"CD_inflamed","CD_non-inflamed")
	seurat.object
}


save_UMAP_figures <- function(seurat.object) {
	require(gridExtra)
	filename <- paste(seurat.object,".pdf",sep="")
	pdf(filename, width = 10, height = 10)
		p.umap <- DimPlot(seurat.object, reduction = "umap", label=T)
		p1 <- DimPlot(seurat.object, reduction = "umap", label=T, group.by="SingleR.Martin")
		p2 <- DimPlot(seurat.object, reduction = "umap", label=T, group.by="SingleR.Luoma")
		p.status <- DimPlot(seurat.object, reduction = "umap", label=F, group.by="status")
		p.patient <- DimPlot(seurat.object, reduction = "umap", label=F, group.by="patient.id")

		print(p.elbow)
		print(p.umap)
		print(p.status)
		print(p.patient)
		print(grid.arrange(p1, p2, nrow = 1))
	dev.off()
}


test_match_order <- function(x,y) {
	if(isTRUE(all.equal(x,y))) print('Perfect match in same order')
	if(!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y)))) print('Perfect match in wrong order')
	if(!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y)))) print('No match')
}

DE_heatmap <- function(seurat.object, filename, n.cores=10) {
	require(Seurat)
	require(future)
	plan(multicore, workers= n.cores)
	s_regress.markers <- FindAllMarkers(seurat.object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
	top5 <- s_regress.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
	# change spaces in cluster names to underscores
	orig.levels <- levels(seurat.object)
	Idents(seurat.object) <- gsub(pattern = " ", replacement = "_", x = Idents(seurat.object))
	orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
	levels(seurat.object) <- orig.levels
	# store cluster averages in a Seurat object
	cluster.averages <- AverageExpression(seurat.object, return.seurat = TRUE)
	# visualize
	p <- DoHeatmap(cluster.averages,features = top5$gene, size = 3, draw.lines = FALSE, raster=F)
	ggsave(filename)
}

DA_analysis <- function(x, des = function(y) model.matrix(~factor(colitis2),y), title) {
	require(edgeR)
	abundances <- table(x$seurat_clusters, x$sample.id) %>% unclass()
	extra.info <- x@meta.data[match(colnames(abundances),x$sample.id),]
	y.ab <- DGEList(abundances,samples=extra.info)
	# filter just for procedural reasons; no clusters should actually be removed
	#keep <- filterByExpr(y.ab, group=y.ab$samples$colitis2) 
	#y.ab <- y.ab[keep,]
	# incorporate design matrix
	design <- des(y.ab$samples)
	y.ab <- estimateDisp(y.ab, design, trend="none")
	fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
	res <- glmQLFTest(fit.ab, coef=ncol(design))
	# figures
	pdf(paste("figures/",title,"_BVC.pdf",sep=""))
	plotBCV(y.ab, cex=1)
	dev.off()
	pdf(paste("figures/",title, "_QLDisp.pdf",sep=""))
	plotQLDisp(fit.ab, cex=1)
	dev.off()
	# output
	return(list(y.ab = y.ab, fit.ab = fit.ab,res = res))
}


matrix.sort.no.diag <- function(matrix) {
 nfile=matrix(0,ncol=max(dim(matrix)),nrow=max(dim(matrix)))
 nfile[,1:ncol(matrix)]=as.matrix(matrix)
 rownames(nfile)=rownames(matrix)
 addname=make.names(rep('A',max(dim(matrix))-ncol(matrix)), unique=T)
 colnames(nfile)=c(colnames(matrix),addname)
 ffile=matrix.sort(nfile)
 '%ni%' <- Negate('%in%')
 ncolnames=colnames(ffile)
 pos=which(ncolnames %ni% addname)
 matrix=ffile[,pos]
}

matrix.sort.no.diag.T <- function(matrix) {
 nfile=matrix(0,ncol=max(dim(matrix)),nrow=max(dim(matrix)))
 nfile[1:nrow(matrix),]=as.matrix(matrix)
 colnames(nfile)=colnames(matrix)
 addname=make.names(rep('A',max(dim(matrix))-nrow(matrix)), unique=T)
 rownames(nfile)=c(rownames(matrix),addname)
 ffile=matrix.sort(nfile)
 '%ni%' <- Negate('%in%')
 ncolnames=rownames(ffile)
 pos=which(ncolnames %ni% addname)
 matrix=ffile[pos,]
}

matrix.sort <- function(matrix) {
 if (nrow(matrix) != ncol(matrix)) stop("Not diagonal")
 if(is.null(rownames(matrix))) rownames(matrix) <- 1:nrow(matrix)

 row.max <- apply(matrix,1,which.max)
 if(all(table(row.max) != 1)) stop("Ties cannot be resolved")

 matrix[names(sort(row.max)),]
}


