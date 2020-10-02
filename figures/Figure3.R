# DESCRIPTION -------------------------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: April 28th, 2020

#Purpose: generate figure 3

#Imports:
# metadata
# normalized counts
# ensembl ID to gene name table
# variance partition output
# t-SNE output

#Exports:
# figure 3

## load libraries  ----------------------------------------------------------------------------------------------
library(tidyverse) #42; base for data/table manipulation
library(ggpubr) #adds some visuals to ggplot, ggaggregate makes publication panels
library(variancePartition) #package for examining variance
library(edgeR) #for count normalization
library(devtools) #load R scripts from github

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes

#load gene annotations ----------------------------------------------------------------------------------------------
tmpG <- read_tsv("data/GRCh38_GENCODE29_geneInfo.txt",
                 skip = 1)
genes.anno <- tmpG[,c("gene_id", "gene_name", "gene_type")]
genes.anno <- unique(genes.anno)

## sample metadata  ----------------------------------------------------------------------------------------------
meta <- read_csv("data/megaMetaTable.csv")
dim(meta)

myPalette <- getPalette(21, "multi")
dev.off()
names(myPalette) <- NULL

diseasePalette <- myPalette[c(1:7)]
names(diseasePalette) <- c("Healthy Control", "Genetic Unaffected",
                           "Idiopathic PD", "Genetic PD",
                           "Prodromal", "SWEDD", "Other ND")

enrollPalette <- diseasePalette[-c(2,4,7)]
enrollPalette <- append(enrollPalette, values = c(myPalette[c(13:14)], "red4", "darkblue"))
names(enrollPalette) <- c("Healthy Control", "De novo PD",
                        "Prodromal", "SWEDD",
                        "Genetic Registry", "Genetic Cohort",
                        "HC Pool", "PD Pool")

complementary <- getPalette(4, "complementary")
dev.off()

# PCA plot ----------------------------------------------------------------------------------------------
meta$SEX[meta$DIAGNOSIS == "HCPOOL"] <- "HC Pool"
meta$SEX[meta$DIAGNOSIS == "PDPOOL"] <- "PD Pool"
meta$SEX[meta$QCflagIR3 == "fail"] <- "failed QC"

meta$SEX <- factor(meta$SEX, levels = c("Female", "Male", "HC Pool", "PD Pool", "failed QC"))

sample.PCA <- ggplot(meta, aes(x = PC1, y = PC2, color = SEX)) +
  geom_point(alpha = 0.6) +
  theme_bw(base_size = 26) +
  xlab("PC1: 25.4% variance") +
  ylab("PC2: 16.7% variance") +
  scale_color_manual(values = c(complementary[c(4,1,3)], "#222222", "gray70")) +
  theme(legend.position = "bottom", legend.title = element_blank())
sample.PCA


#sex check --------------------------------------------------------------------------------------------
tsne.out <- read_csv("data/sexCheck_tsne_output.csv")


tsne_p <- ggplot(tsne.out) +
  geom_point(aes(x = x, y = y, fill = Reported_Sex),
             shape = 21,
             size = 2,
             color = "black") +
  scale_fill_manual(values = c("#FC4E07", "#00AFBB")) +
  theme_bw(base_size = 26) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  xlab("t-SNE dimension 1") +
  ylab("t-SNE dimension 2")
tsne_p

# variance partition --------------------------------------------------------------------------------------------
colPal <- getPalette(14, "multi")
dev.off()
names(colPal) <- NULL

colPal <- append(colPal, values = "lightgrey") #residuals light grey

vp <- read_csv("data/variancePartition_output.csv")
names(colPal) <- c(names(vp)[-c(1,17:18)])

vpPlot <- plotVarPart(vp[,-c(1, 17,18)], col = colPal) +
  theme_bw(base_size = 32) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank())
vpPlot

#pool PCA --------------------------------------------------------------------------------------------
#color palette
poolPalette <- c("#EE9402", "#222222")
names(poolPalette) <- c("HC Pool", "PD Pool")

#load count table, filter out globin and rRNA genes
countTable <- read_tsv("data/PPMI_ALL_geneCountB38.tsv")
head(countTable[,c(1:10)])

globinGenes <- c("CYGB", "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ", "MB")
rRNAGenes <- as.character(genes.anno$gene_name[genes.anno$gene_type %in% c("rRNA", "Mt_rRNA")])

globinGenes.ens <- genes.anno$gene_id[genes.anno$gene_name %in% globinGenes]
rRNAGenes.ens <- genes.anno$gene_id[genes.anno$gene_name %in% rRNAGenes]

countTable.filt <- countTable[!countTable$Ensembl_ID %in% c(globinGenes.ens, rRNAGenes.ens),]
dim(countTable.filt)


meta.pass <- meta[meta$QCflagIR3 == "pass",] #filter out fails
dim(meta.pass) #4756

meta.Pools <- meta.pass[meta.pass$PoolAssign %in% c("PD_POOL", "HC_POOL"),] #filter to pools
meta.Pools$PoolAssign <- gsub("_POOL", " Pool", meta.Pools$PoolAssign)
dim(meta.Pools) #106

poolList <- meta.Pools$HudAlphaSampleName
length(poolList)

sampleData <- meta.Pools[,c("HudAlphaSampleName", "PoolAssign")]

counts.pool <- column_to_rownames(countTable.filt, var = "Ensembl_ID")[ ,sampleData$HudAlphaSampleName]

keep <- filterByExpr(counts.pool)
sum(keep) #20652 genes

dge <- DGEList(counts.pool[keep,])
dge  <- calcNormFactors(dge)

logCPM <- cpm(dge, log = TRUE, prior.count = 3)

pca <- prcomp(t(logCPM)) # perform a PCA on the data in assay(x) for the selected genes
percentVar <- pca$sdev^2 / sum( pca$sdev^2 ) # the contribution to the total variance for each component

# assemble the data for the plot
d <- data.frame(HudAlphaSampleName = row.names(pca$x),
                PC1=pca$x[,1],
                PC2=pca$x[,2])
d <- left_join(d, meta.Pools[,c("HudAlphaSampleName", "PoolAssign", "Plate")])

poolPCA <- ggplot(data=d, aes_string(x = "PC1", y = "PC2", color = "PoolAssign", label = "Plate")) +
  geom_point(size = 3) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100),"% variance")) +
  scale_colour_manual(values = poolPalette) +
  theme_bw(base_size = 26) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  geom_label() +
  stat_ellipse()
poolPCA

# final figure --------------------------------------------------------------------------------------------
ggarrange(ggarrange(tsne_p, sample.PCA, poolPCA, 
                    labels = c("a", "b", "c"),
                    font.label = list(size = 28),
                    ncol = 3, nrow = 1),
          vpPlot,
          labels = c("", "d"),
          font.label = list(size = 28),
          ncol = 1, nrow = 2,
          heights = c(3,4))

ggsave("Figure3.png", width = 24, height = 16, dpi = 600)


