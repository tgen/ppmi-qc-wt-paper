# DESCRIPTION -------------------------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: August 29th, 2019

#Purpose: assess technical variation in the pools

#Imports:
# PPMI count table
# metadata table

#Exports:
# genes that display inter-plate variability

## load libraries  ----------------------------------------------------------------------------------------------
library(tidyverse) #42; base for data/table manipulation
library(edgeR) #filtering by expression level
library(limma) #
library(FactoMineR) #PCA plots
library(factoextra) #PCA plots
library(devtools) #load R scripts from github

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes

toolkitPath <- "~/Dropbox/repos/toolkit_ehutchins/" #for custom scripts and gene annotation files
source(paste(toolkitPath, "R-functions/createColorPalettes.R", sep = "")) #color palettes

duoPalette <- getPalette(21, "multi")[c(1,3)]
names(duoPalette) <- NULL
dev.off()

## load count table  ----------------------------------------------------------------------------------------------
countTable <- read_tsv("data/PPMI_ALL_geneCountB38.tsv")
head(countTable[,c(1:10)])

## biotypes and gene names ----------------------------------------------------------------------------------------
tmpG <- read_tsv("data/GRCh38_GENCODE29_geneInfo.txt",
                 skip = 1)
genes.anno <- tmpG[,c("gene_id", "gene_name", "gene_type")]
genes.anno <- unique(genes.anno)

## remove rRNA and globin ----------------------------------------------------------------------------------------------
globinGenes <- c("CYGB", "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ", "MB")
rRNAGenes <- as.character(genes.anno$gene_name[genes.anno$gene_type %in% c("rRNA", "Mt_rRNA")])

globinGenes.ens <- genes.anno$gene_id[genes.anno$gene_name %in% globinGenes]
rRNAGenes.ens <- genes.anno$gene_id[genes.anno$gene_name %in% rRNAGenes]

countTable.filt <- countTable[!countTable$Ensembl_ID %in% c(globinGenes.ens, rRNAGenes.ens),]
dim(countTable.filt)


## normalize expression, calculate PCA ----------------------------------------------------------------------------------------------
meta <- read_csv("data/megaMetaTable.csv")
dim(meta)

meta <- meta[meta$QCflagIR3 == "pass",] #filter out fails
dim(meta) #4756

meta.Pools <- meta[meta$PoolAssign %in% c("PD_POOL", "HC_POOL"),] #filter out pools
dim(meta.Pools) #106

poolList <- meta.Pools$HudAlphaSampleName
length(poolList)

sampleData <- meta.Pools[,c("HudAlphaSampleName", "PoolAssign")]

counts.pool <- column_to_rownames(countTable.filt, var = "Ensembl_ID")[ ,sampleData$HudAlphaSampleName]

identical(names(counts.pool), meta.Pools$HudAlphaSampleName) #double checking

meta.Pools$diversity <- colSums(counts.pool > 0) #add diversity


meta.Pools.PD <- meta[meta$PoolAssign %in% c("PD_POOL"),] #filter out PD pools
dim(meta.Pools.PD) #53

poolList.PD <- meta.Pools.PD$HudAlphaSampleName
length(poolList.PD)

counts.pool.PD <- column_to_rownames(countTable.filt, var = "Ensembl_ID")[ ,meta.Pools.PD$HudAlphaSampleName]


meta.Pools.HC <- meta[meta$PoolAssign %in% c("HC_POOL"),] #filter out HC pools
dim(meta.Pools.HC) #53

poolList.HC <- meta.Pools.HC$HudAlphaSampleName
length(poolList.HC)

counts.pool.HC <- column_to_rownames(countTable.filt, var = "Ensembl_ID")[ ,meta.Pools.HC$HudAlphaSampleName]


## logCPM normalization ----------------------------------------------------------------------------------------------

keep <- filterByExpr(counts.pool)
sum(keep) #20652 genes

dge <- DGEList(counts.pool[keep,])
dge  <- calcNormFactors(dge)

logCPM <- cpm(dge, log = TRUE, prior.count = 3)

keep.pd <- filterByExpr(counts.pool.PD)
sum(keep.pd) #20309 genes

dge.PD <- DGEList(counts.pool.PD[keep.pd,])
dge.PD  <- calcNormFactors(dge.PD)

logCPM.PD <- cpm(dge.PD, log = TRUE, prior.count = 3)

keep.HC <- filterByExpr(counts.pool.HC)
sum(keep.HC) #20632 genes

dge.HC <- DGEList(counts.pool.HC[keep.HC,])
dge.HC  <- calcNormFactors(dge.HC)

logCPM.HC <- cpm(dge.HC, log = TRUE, prior.count = 3)

## pca plot ----------------------------------------------------------------------------------------------
pca.PD <- prcomp(t(logCPM.PD)) # perform a PCA on the data in assay(x) for the selected genes
percentVar <- pca.PD$sdev^2 / sum( pca.PD$sdev^2 ) # the contribution to the total variance for each component
sum(round(percentVar * 100, digits = 3)[c(1:10)])

var.PD <- get_pca_var(pca.PD)

fviz_eig(pca.PD, addlabels = TRUE, ylim = c(0, 100))



fviz_contrib(pca.PD, choice = "var", axes = 1:10, top = 100) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("gene_id")
ggsave(paste(outdir, "PD_POOL_geneContributions.png", sep = ""))


contrib.PD <- rownames_to_column(data.frame(var.PD$contrib[,c(1:10)]), var = "gene_id")
contrib.PD$sum.contrib <- rowSums(contrib.PD[,-1]) #contribution sums for the first 10 PCs
contrib.PD <- right_join(genes.anno, contrib.PD)
write_csv(contrib.PD, "PD_POOL_PCA_gene_contributions.csv", sep = "")


fviz_pca_biplot(pca.PD, 
                # Individuals
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2,

                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                
                ggtheme = theme_bw(base_size = 18)
)



pca.HC <- prcomp(t(logCPM.HC)) # perform a PCA on the data in assay(x) for the selected genes
percentVar <- pca.HC$sdev^2 / sum( pca.HC$sdev^2 ) # the contribution to the total variance for each component
sum(round(percentVar * 100, digits = 3)[c(1:10)])

var.HC <- get_pca_var(pca.HC)

fviz_eig(pca.HC, addlabels = TRUE, ylim = c(0, 100))



fviz_contrib(pca.HC, choice = "var", axes = 1:10, top = 100) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("gene_id")

contrib.HC <- rownames_to_column(data.frame(var.HC$contrib[,c(1:10)]), var = "gene_id")
contrib.HC$sum.contrib <- rowSums(contrib.HC[,-1]) #contribution sums for the first 10 PCs
contrib.HC <- right_join(genes.anno, contrib.HC)
write_csv(contrib.HC, "HC_POOL_PCA_gene_contributions.csv", sep = "")



fviz_pca_biplot(pca.HC, 
                # Individuals
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2,

                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                
                ggtheme = theme_bw(base_size = 18)
)

## genes that make sum of 0.1 percent in first 10 PCs --------------------------------------------------------------------------------

perc0.1 <- c(contrib.HC$gene_id[contrib.HC$sum.contrib > 0.1], contrib.PD$gene_id[contrib.PD$sum.contrib > 0.1]) 
perc0.1 <- unique(perc0.1)
write_csv(as.data.frame(perc0.1), "filteredGenes.csv")

## output filtered count table ------------------------------------------------------------------------------------------------------------
countTable.filt.contrib <- countTable.filt[!countTable.filt$Ensembl_ID %in% perc0.1,]
dim(countTable.filt.contrib) #56041 genes
write_csv(countTable.filt.contrib, "data/rawCounts_poolPCAContributeRemoved.csv")


# filter out genes ------------------------------------------------------------------------------------------------------------
calcCPM <- function(df){
  keep <- filterByExpr(df)
  sum(keep)
  
  dge <- DGEList(df[keep,])
  dge  <- calcNormFactors(dge)
  
  cpm.df <- cpm(dge)
  
  return(cpm.df)
  
}


counts.pool.filt <- counts.pool[!row.names(counts.pool) %in% perc0.1,]

counts.pool.filt.PD <- counts.pool.filt[ ,meta.Pools.PD$HudAlphaSampleName]
counts.pool.filt.HC <- counts.pool.filt[ ,meta.Pools.HC$HudAlphaSampleName]


CPM.PD.before <- calcCPM(counts.pool.PD)
CPM.HC.before <- calcCPM(counts.pool.HC)
CPM.PD.after <- calcCPM(counts.pool.filt.PD)
CPM.HC.after <- calcCPM(counts.pool.filt.HC)


## overall spearman's correlation before and after ----------------------------------------------------------------------------------------------
getMeanCor <- function(x){
  corMat <- cor(x)
  corMat[apply(corMat, 2, function(x) x == 1)] <- NA
  tmp <- colMeans(corMat, na.rm = TRUE)
  meanCor <- mean(tmp)
  return(meanCor)
}

cor.HC.before <- getMeanCor(logCPM.HC)
cor.HC.before
cor.PD.before <- getMeanCor(logCPM.PD)
cor.PD.before
cor.HC.after <- getMeanCor(CPM.HC.after)
cor.HC.after
cor.PD.after <- getMeanCor(CPM.PD.after)
cor.PD.after

df <- data.frame("pool" = c("HC", "PD"),
           "spearmans_before" = c(cor.HC.before, cor.PD.before),
           "spearmans_after" = c(cor.HC.after, cor.PD.after)
           )
kable(df)
write_csv(df, paste(outdir, "spearmans_beforeAfter.csv", sep = ""))

df.p <- ggtexttable(df, rows = NULL, theme = ttheme("lBlackWhite", base_size = 14))
df.p

