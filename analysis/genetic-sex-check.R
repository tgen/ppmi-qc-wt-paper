# DESCRIPTION ----------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: March 21st, 2019

#Purpose: Check genetic sex assignment in PPMI dataset using genetic sex-specific genes

#Imports:
# PPMI TPM table

#Exports:
# t-SNE analysis and plots

## load libraries  -------------------------------------------------------------------------
library(tidyverse) #basis for data manipulation
library(ggpubr) #adds some visuals to ggplot, ggaggregate makes publication panels
library(Rtsne) #runs t-SNE analysis
# library(ggfortify)
# library(DESeq2)
library(FactoMineR) #PCA function needed for t-SNE input
library(factoextra) #additional PCA functions
# library(corrplot)
# library(ggrepel)
library(devtools) #load R scripts from github

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes

myPalette <- getPalette(3, "complementary")[c(3,1)]
dev.off()

##Import expression table -------------------------------------------------------------------

exprTable <- read_tsv("data/archive/featureCountsB38_deseq2Normalized.tsv")

##Import metadata table and filter out pools ------------------------------------------------

meta <- read_csv("data/megaMetaTable.csv")

meta <- meta %>%
  filter(QCflagIR3 != "remove") %>%
  filter(is.na(PoolAssign))
dim(meta)


##Subset to genes of interest
geneIDs <- c("ENSG00000229807.11" = "XIST",
             "ENSG00000129824.15" = "RPS4Y1",
             "ENSG00000280969.1" = "RPS4Y2",
             "ENSG00000012817.15" = "KDM5D",
             "ENSG00000067048.16" = "DDX3Y",
             "ENSG00000114374.12" = "USP9Y")


expr.sexGenes <- exprTable[exprTable$ens_id %in% names(geneIDs),]


expr.sexGenes.t <- t(column_to_rownames(expr.sexGenes, var = "ens_id")) #transpose
expr.sexGenes.t.meta <- inner_join(meta[,c("HudAlphaSampleName", "SEX")], rownames_to_column(as.data.frame(expr.sexGenes.t), var = "HudAlphaSampleName")) #join with metadata

expr.sexGenes.t.meta <- expr.sexGenes.t.meta %>% plyr::rename(geneIDs)
head(expr.sexGenes.t.meta)



##Run PCA ---------------------------------------------------------------------------------
pc <- expr.sexGenes.t.meta[,-c(1:2)]

res.pca <- PCA(pc, graph = FALSE)
eig.val <- get_eigenvalue(res.pca)

var <- get_pca_var(res.pca)

#ok now to look at individual contributions
ind <- get_pca_ind(res.pca)

fviz_pca_biplot(res.pca, 
                # Individuals
                geom.ind = "point",
                fill.ind = expr.sexGenes.t.meta$SEX, col.ind = "black",
                pointshape = 21,
                pointsize = 2,
                palette = myPalette,
                addEllipses = TRUE,
                ellipse.level = 0.5,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                
                legend.title = list(fill = "Reported Sex", color = "Contrib",
                                    alpha = "Contrib"),
                
                ggtheme = theme_bw()
)

##Run t-SNE -------------------------------------------------------------------
tsne <- Rtsne(pc, dims = 3, perplexity = 1300, verbose = TRUE, max_iter = 5000, check_duplicates = FALSE)


tsne.out <- tibble(Sample = rownames(expr.sexGenes.t.meta),
                   Reported_Sex = expr.sexGenes.t.meta$SEX,
                   x = tsne$Y[,1],
                   y = tsne$Y[,2],
                   z = tsne$Y[,3],
                   flag = FALSE)

#identify swaps and flag
swap1 <- subset(tsne.out, Reported_Sex == "Male" & x > 0)
swap2 <- subset(tsne.out, Reported_Sex == "Female" & x < 0)
swappedIDs <- c(swap1$Sample, swap2$Sample)

tsne.out$flag[tsne.out$Sample %in% swappedIDs] <- "TRUE"
tsne.out$flag <- factor(tsne.out$flag, levels = c("TRUE", "FALSE"))
tsne.out <- tsne.out %>% arrange(desc(flag))

write_csv(tsne.out, "data/sexCheck_tsne_output.csv")


