# DESCRIPTION -------------------------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: September 10th, 2019

#Purpose: variance partition

#Imports:
# metadata table
# filtered count table after inter-plate-variability.R ("data/rawCounts_poolPCAContributeRemoved.csv")

#Exports:
# 

## load libraries  ----------------------------------------------------------------------------------------------
library(tidyverse) #42; base for data/table manipulation
library(tidyr) #for pivot tables
library(variancePartition) #package for examining variance
library(doParallel) #use multiple cores
library(edgeR) #expression normalization
library(limma) #more expression normalization
library(corrplot) #correlation plot matrix
library(beepr) #beep alerts when job is complete
library(ggpubr) #paper level visuals
library(seqsetvis) #ggplot2 venn
library(devtools) #load R scripts from github

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes


# function to run VP ---------------------------------------------------------------------------------------------
runVP <- function(exprObj, form, info, formCor, colPal, outName){
  
  #normalize with limma/voom ----------------------------------------------------------------------------------
  # optional step to run analysis in parallel on multicore machines
  # Here use 16 threads
  # This is strongly recommended since the analysis
  # can be computationally intensive
  cl <- makeCluster(16)
  registerDoParallel(cl)
  
  # identify genes that pass expression cutoff
  geneCounts <- column_to_rownames(exprObj, var = "Ensembl_ID")
  isexpr <- rowSums(cpm(geneCounts)>0.5) >= 0.25 * ncol(geneCounts)
  
  # create data structure with only expressed genes
  gExpr <- DGEList(counts=geneCounts[isexpr,])
  
  # Perform TMM normalization
  gExpr <- calcNormFactors(gExpr)
  
  # Specify variables to be included in the voom() estimates of
  # uncertainty.
  # Recommend including variables with a small number of categories
  # that explain a substantial amount of variation
  #leaving NULL for the time being
  #design <- model.matrix( ~ info)
  
  
  # Estimate precision weights for each gene and sample
  # This models uncertainty in expression measurements
  vobjGenes <- voom(gExpr)
  
  
  # variancePartition seamlessly deals with the result of voom()
  # by default, it seamlessly models the precision weights
  # This can be turned off with useWeights=FALSE
  varPart <- fitExtractVarPartModel( vobjGenes, form, info )
  
  #export as rds
  saveRDS(varPart, file = "data/varPart.rds")
  
  # sort variables (i.e. columns) by median fraction
  # of variance explained
  vp <- sortCols( varPart )
  
  # Figure 1a
  # Bar plot of variance fractions for the first 10 genes
  plotPercentBars(vp[1:20,], col = colPal) +
    theme_bw(base_size = 18)
  ggsave(paste(outdir, "varianceExplained_first20_", outName, ".png", sep = ""))
  
  #write results to table
  vrf <- vp
  vrf <- left_join(rownames_to_column(vrf, var = "gene_id"), genes.anno)
  write_csv(vrf, "variancePartition_output.csv")

  
  stopCluster(cl)
  
  return(vp)
  
}




## sample metadata: participant data  ----------------------------------------------------------------------------------------------

meta.noPools <- meta[!meta$PoolAssign %in% c("PD_POOL", "HC_POOL"),] #filter out pools
dim(meta.noPools) #4650

#add in med use info
medUseYN <- read_csv("data/Use_of_PD_Medication.csv")
medUseYN$PATNO_VISIT <- paste(medUseYN$PATNO, medUseYN$EVENT_ID, sep = "_")
meta.noPools <- left_join(meta.noPools, medUseYN[,c("PATNO_VISIT", "PDMEDYN")])
dim(meta.noPools) #4650

## import count table and filter out variable genes ----------
countTable.filtered <- read_csv("data/rawCounts_poolPCAContributeRemoved.csv")

countTable.filt.noPools <- countTable.filtered[,c("Ensembl_ID", meta.noPools$HudAlphaSampleName)]
dim(countTable.filt.noPools) #56041 genes

## factors for VP ----------------------------------------------------------------------------------------------
info <- data.frame(row.names = meta.noPools$HudAlphaSampleName,
                   ID = meta.noPools$HudAlphaSampleName,
                   
                   Individual = as.character(meta.noPools$PATNO),
                   
                   
                   Plate = meta.noPools$Plate,
                   Position = as.character(meta.noPools$POSITION),
                   
                   Disease = meta.noPools$Disease_Status,
                   
                   Usable_bases = meta.noPools$PCT_USABLE_BASES,
                   Site = meta.noPools$SITE,
                   Visit = meta.noPools$CLINICAL_EVENT,
                   Sex = meta.noPools$SEX,
                   Ethnicity = meta.noPools$Hispanic_or_Latino,
                   Age = meta.noPools$age_at_consent,
                   Lymphocytes = as.numeric(meta.noPools$`Lymphocytes (%)`),
                   Neutrophils = as.numeric(meta.noPools$`Neutrophils (%)`),
                   Genetics = meta.noPools$Genetic_Status,
                   PD_Medication = as.character(meta.noPools$PDMEDYN)
)

colPal <- getPalette(ncol(info) - 1, "multi")
dev.off()
names(colPal) <- NULL

colPal <- append(colPal, values = "lightgrey") #residuals light grey
names(colPal) <- c(names(info[c(2, 6, 3, 12:13, 7, 4:5, 8:11, 14:15)]), "Residuals")

# Define formula
form <- ~ (1|Plate) + (1|Position) + (1|Site) + (1|Visit) + (1|Sex) + (1|Ethnicity) + (1|Individual) + Age + Lymphocytes + Neutrophils + Usable_bases + (1|Disease) + (1|Genetics) + (1|PD_Medication)

# formula for correlation plot
formCor <- ~ Plate + Position + Site + Visit + Sex + Ethnicity + Individual + Age + Lymphocytes + Neutrophils + Usable_bases + Disease + Genetics + PD_Medication

#run variance partition  ----------------------------------------------------------------------------------------------
subject.vp <- runVP(countTable.filt.noPools, form, info, formCor, colPal)
beep(sound = 1)

