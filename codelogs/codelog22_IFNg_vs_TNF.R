
# following function was gratiously provided by Shengbao Suo
ActivityScore <- function(SeuratObject,GeneSet,nCores=2) {
	# parameter:
	# SeuratObject: Seurat object.
	# GeneSet: should only contain one gene set and the class of this variable should be a list), the following is an example, this gene set contain 10 genes and the name of this gene set is geneSet1:
	# $geneSet1
	# [1] "Gene13" "Gene18" "Gene3" "Gene19" "Gene1" "Gene11" "Gene8" "Gene10"
	# [9] "Gene9" "Gene20"
	# nCores: how many core would be used, default is 2

	# This function returns a Seurat object, and the activity score of query geneset has been included in the meta.data of this Seurat object.
	require(AUCell)
	require(Seurat)
	data <- GetAssayData(object = SeuratObject, slot = "data")
	meta.data <- SeuratObject@meta.data
	data <- data[,rownames(meta.data)]
	aucellRankings1 <- AUCell_buildRankings(data, nCores=nCores, plotStats=FALSE)
	regulonAUC1 <- AUCell_calcAUC(GeneSet, aucellRankings1, nCores=nCores, aucMaxRank = ceiling(0.05 * nrow(aucellRankings1)))
	regulonMatix <- getAUC(regulonAUC1)
	meta.data$activity <- as.numeric(regulonMatix[names(GeneSet),])
	SeuratObject@meta.data <- meta.data
	return(SeuratObject)
}


# # example run
# pathways <- fgsea::gmtPathways('hallmark.gene.sets.all.v7.0.symbols.gmt')
# tumor=readRDS('seurat.object.rds')
# geneset=pathways['INTERFERON_GAMMA_RESPONSE']
# tumor=ActivityScore(tumor,geneset,nCores=20)

library(msigdbr)
hallmarks <- msigdbr(species = "Homo sapiens", category = "H")
hallmarks_list <- split(x = hallmarks$gene_symbol, f = hallmarks$gs_name)
geneset1 <- hallmarks_list["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]
geneset2 <- hallmarks_list["HALLMARK_INTERFERON_GAMMA_RESPONSE"]
s <- readRDS("saved_objects/UC_chang_021321/uccpirec-allCD45.rds")
s <- ActivityScore(s, geneset1, nCores=20)
s2 <- ActivityScore(s, geneset2, nCores=20)

# visualization

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}


# plot for all cells
df <- data.frame(activity = s$activity, sample = s$sample.id, condition=s$colitis3)
df.summary <- summarySE(df, measurevar="activity", groupvars=c("sample","condition"))
df.meanofmeans.TNF <- summarySE(df.summary, measurevar="activity", groupvars=c("condition"))
df.meanofmeans.TNF$pathway <- "TNFA_SIGNALING_VIA_NFKB"

df <- data.frame(activity = s2$activity, sample = s2$sample.id, condition=s2$colitis3)
df.summary <- summarySE(df, measurevar="activity", groupvars=c("sample","condition"))
df.meanofmeans.IFNg <- summarySE(df.summary, measurevar="activity", groupvars=c("condition"))
df.meanofmeans.IFNg$pathway <- "INTERFERON_GAMMA_RESPONSE"

df.combined <- rbind(df.meanofmeans.TNF,df.meanofmeans.IFNg)
df.combined$condition <- factor(df.combined$condition, levels= c("CPI-Control", "UC-Control" , "CPI-No-Colitis", "CPI-Colitis", "UC"))
# plot
key.colitis.colors = c("#a3d2ca", "#ffcda3", "#5eaaa8","#056676","#db6400")
p <- ggplot(df.combined, aes(x=pathway, y=activity, fill=condition)) + 
	        geom_bar(stat="identity", position=position_dodge()) +
	        geom_errorbar(aes(ymin=activity-ci, ymax=activity+ci), width=.2, position=position_dodge(.9)) +
	        labs(y = "AUCell activity score") + theme_classic() + scale_fill_manual(values=key.colitis.colors) +
	        theme(axis.title.x = element_blank())
ggsave("figures/FIGURE_REDO_031321/IFNg-TNF_activity_score_allcells.pdf", width= 7, height= 4, units= "in")

# TNF: split across celltypes
df <- data.frame(activity = s$activity, sample = s$sample.id, condition=s$colitis3, celltype = s$SingleR.Luoma2)
df.summary <- summarySE(df, measurevar="activity", groupvars=c("sample","condition","celltype"))
df.meanofmeans.TNF <- summarySE(df.summary, measurevar="activity", groupvars=c("condition","celltype"))
df.meanofmeans.TNF$condition <- factor(df.meanofmeans.TNF$condition, levels= c("CPI-Control", "UC-Control" , "CPI-No-Colitis", "CPI-Colitis", "UC"))
df.meanofmeans.TNF$celltype <- factor(df.meanofmeans.TNF$celltype, levels=c("CD4.T" ,"CD8.T" ,"Plasma.B", "B", "Myeloid","Mast", "ILCs"))

key.colitis.colors = c("#a3d2ca", "#ffcda3", "#5eaaa8","#056676","#db6400")
p <- ggplot(df.meanofmeans.TNF, aes(x=celltype, y=activity, fill=condition)) + 
	        geom_bar(stat="identity", position=position_dodge()) +
	        geom_errorbar(aes(ymin=activity-ci, ymax=activity+ci), width=.2, position=position_dodge(.9)) +
	        labs(y = "AUCell activity score") + theme_classic() + scale_fill_manual(values=key.colitis.colors) +
	        theme(axis.title.x = element_blank())
ggsave("figures/FIGURE_REDO_031321/activity_score_TNF_across_allcells.pdf", width= 7, height= 4, units= "in")


# IFNg: split across celltypes
df <- data.frame(activity = s2$activity, sample = s2$sample.id, condition=s2$colitis3, celltype = s2$SingleR.Luoma2)
df.summary <- summarySE(df, measurevar="activity", groupvars=c("sample","condition","celltype"))
df.meanofmeans.IFNg <- summarySE(df.summary, measurevar="activity", groupvars=c("condition","celltype"))
df.meanofmeans.IFNg$condition <- factor(df.meanofmeans.IFNg$condition, levels= c("CPI-Control", "UC-Control" , "CPI-No-Colitis", "CPI-Colitis", "UC"))
df.meanofmeans.IFNg$celltype <- factor(df.meanofmeans.IFNg$celltype, levels=c("CD4.T" ,"CD8.T" ,"Plasma.B", "B", "Myeloid","Mast", "ILCs"))

key.colitis.colors = c("#a3d2ca", "#ffcda3", "#5eaaa8","#056676","#db6400")
p <- ggplot(df.meanofmeans.IFNg, aes(x=celltype, y=activity, fill=condition)) + 
	        geom_bar(stat="identity", position=position_dodge()) +
	        geom_errorbar(aes(ymin=activity-ci, ymax=activity+ci), width=.2, position=position_dodge(.9)) +
	        labs(y = "AUCell activity score") + theme_classic() + scale_fill_manual(values=key.colitis.colors) +
	        theme(axis.title.x = element_blank())
ggsave("figures/FIGURE_REDO_031321/activity_score_IFNg_across_allcells.pdf", width= 7, height= 4, units= "in")


# calculate p-values: two sample t-test
celltypes <- c("CD4.T" ,"CD8.T" ,"Plasma.B", "B", "Myeloid","Mast", "ILCs")
t.test.conditions <- function(cell.type, df = df.summary) {
	group1 <- df %>% 
		filter(celltype==cell.type) %>%
		filter(condition=="CPI-Colitis")

	group2 <- df %>% 
		filter(celltype==cell.type) %>%
		filter(condition=="UC")

	test.result <- t.test(group1$activity,group2$activity)
	test.result$p.value
}
t.test.controls <- function(cell.type, df = df.summary) {
	group1 <- df %>% 
		filter(celltype==cell.type) %>%
		filter(condition=="CPI-Control")

	group2 <- df %>% 
		filter(celltype==cell.type) %>%
		filter(condition=="UC-Control")

	test.result <- t.test(group1$activity,group2$activity)
	test.result$p.value
}

df <- data.frame(activity = s$activity, sample = s$sample.id, condition=s$colitis3, celltype = s$SingleR.Luoma2)
df.summary <- summarySE(df, measurevar="activity", groupvars=c("sample","condition","celltype"))
p.values.TNF <- c(unlist(lapply(celltypes,t.test.conditions)),
unlist(lapply(celltypes,t.test.controls)))
p.values.TNF <- p.adjust(p.values.TNF, method = "bonferroni", n = length(p.values.TNF))
df <- data.frame(activity = s2$activity, sample = s2$sample.id, condition=s2$colitis3, celltype = s2$SingleR.Luoma2)
df.summary <- summarySE(df, measurevar="activity", groupvars=c("sample","condition","celltype"))
p.values.IFNg <- c(unlist(lapply(celltypes,t.test.conditions)),
unlist(lapply(celltypes,t.test.controls)))
p.values.IFNg <- p.adjust(p.values.IFNg, method = "bonferroni", n = length(p.values.IFNg))
p.values <- c(p.values.TNF, p.values.IFNg)
matrix(p.values,ncol=7)

# pvalues for all cells comparison
t.test.colitis <- function(df = df.summary) {
	group1 <- df %>% 
		filter(condition=="CPI-Colitis")

	group2 <- df %>% 
		filter(condition=="UC")

	test.result <- t.test(group1$activity,group2$activity)
	test.result$p.value
}
t.test.ctrl <- function(df = df.summary) {
	group1 <- df %>% 
		filter(condition=="CPI-Control")

	group2 <- df %>% 
		filter(condition=="UC-Control")

	test.result <- t.test(group1$activity,group2$activity)
	test.result$p.value
}

df <- data.frame(activity = s$activity, sample = s$sample.id, condition=s$colitis3)
df.summary.TNF <- summarySE(df, measurevar="activity", groupvars=c("sample","condition"))
df <- data.frame(activity = s2$activity, sample = s2$sample.id, condition=s2$colitis3)
df.summary.IFNg <- summarySE(df, measurevar="activity", groupvars=c("sample","condition"))
p <- c(t.test.colitis(df.summary.TNF),
t.test.ctrl(df.summary.TNF),
t.test.colitis(df.summary.IFNg),
t.test.ctrl(df.summary.IFNg))
adjp <- p.adjust(p, method = "bonferroni", n = length(p))
