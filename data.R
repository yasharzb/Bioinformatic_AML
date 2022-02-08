### Required Libraries

library(GEOquery)
library(limma)
library(Biobase)
library(reshape2)
library(plyr)
library(ggplot2)
library(stringr)
### Loading Data

gset <- getGEO("GSE48558", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("0000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXX3XXX3XXXXX",
               "XXXXXXXXXXXXXXXXXX5X1XXX3X3225X1XX11XX11X5X1X5X1X4",
               "XXX4XXX4XXXXXXXXXXXXXXXXXXXXXXXXXXXXX3333333001000",
               "55555551222231111111")
gsms <- gsub("", ",", gsms)
gsms <- gsub("0", "Test", gsms)
gsms <- gsub("1", "Normal_T_Cells", gsms)
gsms <- gsub("2", "Normal_Monocytes", gsms)
gsms <- gsub("3", "Normal_Granubcytes", gsms)
gsms <- gsub("4", "Normal_CD34_HSPC", gsms)
gsms <- gsub("5", "Normal_B_Cells", gsms)
sml <- strsplit(gsms, split=",")[[1]]
sml <- sml[-1]

# filter out excluded samples (marked as "X")

sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]
### Analysis

gene_exp <- exprs(gset)
boxplot(gene_exp)

### PCA

evaluate_pca <- function(g_expr, sml) {
  pca <- prcomp(g_expr)
  summary(pca)
  plot(pca)
  plot(pca$x[,1:2])
  ggplot(data.frame(pca$rotation[,1:3], group = sml), aes(PC1 , PC2 , color = group)) + geom_point()
}

evaluate_pca(gene_exp, sml)

gene_exp.scaled <- t(scale(t(gene_exp), scale = FALSE))
evaluate_pca(gene_exp.scaled, sml)
### Correlation Heatmap

heatmap(cor(gene_exp), labRow = sml, labCol = sml)

### Difference

pca <- prcomp(gene_exp.scaled)
pca_rot <- data.frame(pca$rotation[,1:3], group = sml)
sml_near <- sml
sml_near[pca_rot$PC2 > 0.13 & pca_rot$PC1 < -0.025 & pca_rot$group == "Test"] <- "Near"
sml_near <- factor(sml_near)
gset$group <- sml_near

design <- model.matrix(~ group + 0, gset)
colnames(design) <- levels(sml_near)

### Model fitting

fit <- lmFit(gset, design)

cont <- makeContrasts(Near - Normal_CD34_HSPC, levels = design)
fit_cont <- contrasts.fit(fit, cont)
fit_cont <- eBayes(fit_cont, 0.01)

top_table <- topTable(fit_cont, adjust = "fdr", sort.by = "B", number = Inf)

### Result

summary_table <- subset(top_table , select = c("Gene.symbol" , "Gene.title", "adj.P.Val"  , "logFC"))

exp_high <- subset(summary_table, logFC > 1 & adj.P.Val < 0.05)
high_genes <-unique( as.character(strsplit2((exp_high$Gene.symbol),"///")))
exp_low <- subset(summary_table, logFC < -1& adj.P.Val < 0.05)
low_genes <- unique(as.character(strsplit2( (exp_low$Gene.symbol),"///")))

print("High exp:")
high_genes

print("Low exp:")
low_genes
