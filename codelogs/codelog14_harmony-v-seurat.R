############################################################
## Integration comparison for cdcpi3_harmony vs cdcpi3_anchor
## edward_kim@college.harvard.edu - October 2020
#############################################################

##############################
# 0 - Load librairies
##############################
library(tidyverse)
library(Seurat)
library(future)
library(future.apply)
library(harmony)
#library(scran)
library(pheatmap)

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

#############################################
# 3 - Fixing sample.id labels
#############################################
s.h <- readRDS("saved_objects/cdcpi3_harmony.rds")
s.a <- readRDS("saved_objects/cdcpi3_anchor.rds")
s1 <- readRDS("saved_objects/CD_martin_qc_100720/martin_naive_CD3_15100.rds")
s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")

tab <- table(paste("after", s.h$seurat_clusters[s.h$orig.ident=="CD_martin"]),
    paste("before", s1$seurat_clusters))

#### colname exploration
#str_sub(colnames(s.h),-6) %>% table()
str_sub(colnames(s.a),-6) %>% table()
#str_sub(colnames(s2),-4) %>% table()
#str_sub(names(x),-6) %>% table()

#### adding sample ID for the Harmony integration
x <- s2$sampletype
names(x) <- str_glue("{names(x)}_2")
s.h$sample.id.backup <- s.h$sample.id
s.h$sample.id <- x
s.h$sample.id <- ifelse(s.h$orig.ident == "CD_martin", as.character(s.h$sample.id.backup), s.h$sample.id)
# check: table(s.h$sample.id, useNA="ifany")

#### adding CPI sample ID for the Seurat integration
x <- s2$sampletype
names(x) <- str_glue("{names(x)}_1")
s.a$sample.id.backup <- s.a$sample.id
s.a$sample.id <- x
# check if correctly added: table(s.a$sample.id[s.h$orig.ident=="CD3_martin"],useNA="ifany")
#### adding CD sample ID correctly
x <- s1$sample.id
table(x) # you should have CD__# . If not do the following:
# x <- str_glue("CD_{str_sub(x,-3)}")
names(x) <- str_glue("{names(x)}_2")
s.a$sample.id.backup <- x
table(s.a$sample.id.backup, useNA="ifany")
s.a$sample.id <- ifelse(s.a$orig.ident == "CD_martin", as.character(s.a$sample.id.backup), s.a$sample.id)
# check: table(s.a$sample.id, useNA="ifany")

saveRDS(s.h,"saved_objects/cdcpi3_harmony.rds")
saveRDS(s.a,"saved_objects/cdcpi3_anchor.rds")


############################################
# 4 - Heatmap comparing Seurat and Harmony
############################################
s.h <- readRDS("saved_objects/cdcpi3_harmony.rds")
s.a <- readRDS("saved_objects/cdcpi3_anchor.rds")
# buffer storage with just the clusters:
a <- s.h$seurat_clusters
b <- s.a$seurat_clusters
levels(a) <- paste("harmony", levels(a))
levels(b) <- paste("seurat", levels(b))
# investigation:
names(a[s.h$orig.ident=="CD_martin"])
names(a[s.h$orig.ident=="CD3_Tcell"])
names(b[s.a$orig.ident=="CD_martin"])
names(b[s.a$orig.ident=="CD3_Tcell"])
# swapping indices in a to match the those in b:
names(a) <- ifelse(s.h$orig.ident=="CD_martin",
	str_glue("{str_sub(names(a),1,-2)}2"), 
	str_glue("{str_sub(names(a),1,-2)}1"))
# check match and order appropiately
test_match_order(names(sort(a)),names(b)) # found that they are in different orders
temp.df <- merge(data.frame(cell=names(a),a=a),data.frame(cell=names(b),b=b), by="cell")
# producing the clustering correlation figure
tab <- table(temp.df$a,temp.df$b)
heat <- pheatmap(matrix.sort.no.diag.T(log10(tab+1)), cluster_row=FALSE, cluster_col=FALSE, main="Integration Comparison", silent=TRUE)

pdf("figures/cdcpi3_harmony/cdcpi3_harmonyvanchor_heatmapv3.pdf")
gridExtra::grid.arrange(heat[[4]])
dev.off()
# after comparing with and without doing the swap step, it turns out that index swap was not necessary. In addition, the match order step was not necessary either.


#######################################
# Metric 1a: DA of CD inf vs uninf
#######################################
s.h <- readRDS("saved_objects/cdcpi3_harmony.rds")
s.a <- readRDS("saved_objects/cdcpi3_anchor.rds")
s1 <- readRDS("saved_objects/CD_martin_qc_100720/martin_naive_CD3_15100.rds")
# remove cluster 8 in s.a
s.a <- subset(s.a, idents='8',invert=T)
s.a$seurat_clusters <- droplevels(s.a$seurat_clusters)

# DA analysis: CD inflamed vs CD non-inflamed
# CD only dataset
library(edgeR)
x <- subset(s1, orig.ident == "CD_martin")
x$sample.id <- droplevels(x$sample.id)
x$CD_inf <- x$colitis == "CD_inflamed"
design <- function(y) model.matrix(~factor(CD_inf),y)
out <- DA_analysis(x, design, title = "martin_naive_CD315100_CD_inf_vs_uninf")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

# Seurat v3 anchors
library(edgeR)
x <- subset(s.a, orig.ident == "CD_martin")
x$CD_inf <- x$colitis == "CD_inflamed"
design <- function(y) model.matrix(~factor(CD_inf),y)
out <- DA_analysis(x, design, title = "cdcpi3_anchor_CD_inf_vs_uninf")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

# Harmony
library(edgeR)
x <- subset(s.h, orig.ident == "CD_martin")
x$sample.id <- droplevels(x$sample.id)
x$CD_inf <- x$colitis == "CD_inflamed"
design <- function(y) model.matrix(~factor(CD_inf),y)
out <- DA_analysis(x, design, title = "cdcpi3h_CD_inf_vs_uninf")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

#############################################
# Metric 1b: DA of CD inf vs uninf w/o CD_P16
#############################################
s.h <- readRDS("saved_objects/cdcpi3_harmony.rds")
s.a <- readRDS("saved_objects/cdcpi3_anchor.rds")
s1 <- readRDS("saved_objects/CD_martin_qc_100720/martin_naive_CD3_15100.rds")

# DA analysis: CD inflamed vs CD non-inflamed
# CD only dataset
library(edgeR)
x <- subset(s1, orig.ident == "CD_martin" & patient.id != "P16")
x$sample.id <- droplevels(x$sample.id)
x$CD_inf <- x$colitis == "CD_inflamed"
design <- function(y) model.matrix(~factor(CD_inf),y)
out <- DA_analysis(x, design, title = "martin_naive_CD315100_wo16_CD_inf_vs_uninf")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

# Seurat v3 anchors
library(edgeR)
x <- subset(s.a, orig.ident == "CD_martin" & patient.id != "CD_P16")
x$CD_inf <- x$colitis == "CD_inflamed"
design <- function(y) model.matrix(~factor(CD_inf),y)
out <- DA_analysis(x, design, title = "cdcpi3_anchor_wo16_CD_inf_vs_uninf")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

# Harmony
library(edgeR)
x <- subset(s.h, orig.ident == "CD_martin" & patient.id != "CD_P16")
x$sample.id <- droplevels(x$sample.id)
x$CD_inf <- x$colitis == "CD_inflamed"
design <- function(y) model.matrix(~factor(CD_inf),y)
out <- DA_analysis(x, design, title = "cdcpi3h_wo16_CD_inf_vs_uninf")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

#############################################
# Metric 1c: DA of CD inf vs uninf w/o CD_P16
#############################################
# CD only dataset
library(edgeR)
x <- subset(s1, orig.ident == "CD_martin" & patient.id != "P16" & patient.id != "P5")
x$sample.id <- droplevels(x$sample.id)
x$CD_inf <- x$colitis == "CD_inflamed"
design <- function(y) model.matrix(~factor(CD_inf),y)
out <- DA_analysis(x, design, title = "martin_naive_CD315100_wo5n16_CD_inf_vs_uninf")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

# Seurat v3 anchors
library(edgeR)
x <- subset(s.a, orig.ident == "CD_martin" & patient.id != "CD_P16" & patient.id != "CD_P5")

x$CD_inf <- x$colitis == "CD_inflamed"
design <- function(y) model.matrix(~factor(CD_inf),y)
out <- DA_analysis(x, design, title = "cdcpi3_anchor_wo5n16_CD_inf_vs_uninf")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

# Harmony
library(edgeR)
x <- subset(s.h, orig.ident == "CD_martin" & patient.id != "CD_P16" & patient.id != "CD_P5")
x$sample.id <- droplevels(x$sample.id)
x$CD_inf <- x$colitis == "CD_inflamed"
design <- function(y) model.matrix(~factor(CD_inf),y)
out <- DA_analysis(x, design, title = "cdcpi3h_wo5n16_CD_inf_vs_uninf")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)


######################################################
# Metric 2: TFH, Treg, Naive Proportion (w/o cluster8)
######################################################
# integrated seurat v3 vs CD dataset (order: Treg, TFH, Naive)
eval_metric2(a=s.a$seurat_clusters[s.a$orig.ident=="CD_martin"],
	b=s1$seurat_clusters,
	pre.ttn = c("2","4","6"),
	post.ttn = c("6","9","2"))
# integrated seurat v3 vs CPI dataset (order: Treg, TFH, Naive)
eval_metric2(a=s.a$seurat_clusters[s.a$orig.ident=="CD3_Tcell"],
	b=s2$seurat_clusters,
	pre.ttn = c("6","8","7"),
	post.ttn = c("6","9","2"))
# integrated harmony vs CD dataset (order: Treg, TFH, Naive)
eval_metric2(a=s.h$seurat_clusters[s.h$orig.ident=="CD_martin"],
	b=s1$seurat_clusters,
	pre.ttn = c("2","4","6"),
	post.ttn = c("4","2"))
# integrated harmony vs CPI dataset (order: Treg, TFH, Naive)
eval_metric2(a=s.h$seurat_clusters[s.h$orig.ident=="CD3_Tcell"],
	b=s2$seurat_clusters,
	pre.ttn = c("6","8","7"),
	post.ttn = c("4","2"))

######################################################
# Metric 3: Check that there are 2x CD vs CPI
######################################################
# harmony
x <- subset(s.h, idents="4")
table(x$orig.ident)/ncol(x)

# seurat v3 anchor
x <- subset(s.a, idents="6")
table(x$orig.ident)/ncol(x)


############################################################
# Metric 4: Check that there are far more gd cells from CPI
############################################################
# harmony
x <- subset(s.h, idents="8")
table(x$orig.ident)/ncol(x)

# seurat v3 anchor
x <- subset(s.a, idents="5")
table(x$orig.ident)/ncol(x)


