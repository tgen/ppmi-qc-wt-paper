# DESCRIPTION -------------------------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: April 28th, 2020

#Purpose: generate figure 2 for technical paper

#Imports:
# megaTable #metadata with plate info

#Exports:
# figure 2 for scientific data paper

## load libraries  ----------------------------------------------------------------------------------------------
library(tidyverse) #42; base for data/table manipulation
library(ggpubr) #adds some visuals to ggplot, ggaggregate makes publication panels
library(devtools) #load R scripts from github

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes

myPalette <- getPalette(21, "multi")
dev.off()
names(myPalette) <- NULL

## sample metadata  ----------------------------------------------------------------------------------------------
meta <- read_csv("data/megaMetaTable.csv")
meta <- subset(meta, QCflagIR3 == "pass")
names(meta)

meta$freq <- 1 #add frequency of 1 for histograms

meta.noPool <- meta[!meta$PoolAssign %in% c("PD_POOL", "HC_POOL"),]
meta.noPool <- meta.noPool[!is.na(meta.noPool$PATNO),]

# site plot  ----------------------------------------------------------------------------------------------
site.completecases <- meta.noPool[!is.na(meta.noPool$SITE),]
site.completecases$SITE <- as.character(site.completecases$SITE)

site.completecases$enrollArm <- site.completecases$Study_Arm
site.completecases$enrollArm <- gsub("GENUN|GENPD", "Genetic Cohort", site.completecases$enrollArm)
site.completecases$enrollArm <- gsub("REGUN|REGPD", "Genetic Registry", site.completecases$enrollArm)
site.completecases$enrollArm <- gsub("PRODROMA", "Prodromal", site.completecases$enrollArm)
site.completecases$enrollArm <- gsub("PD", "De novo PD", site.completecases$enrollArm)
site.completecases$enrollArm <- gsub("HC", "Healthy Control", site.completecases$enrollArm)

tmp <- site.completecases[,c("SITE", "Study_Arm", "Disease_Status", "enrollArm")] %>%
  group_by(SITE, enrollArm) %>%
  tally() %>%
  pivot_wider(names_from = enrollArm, values_from = n)

tmp <- as.data.frame(tmp)
tmp[is.na(tmp)] <- 0
tmp$total <- rowSums(tmp[,-1]) #total number assigned

siteSummary <- as.data.frame(column_to_rownames(tmp, var = "SITE"))

#convert to percentages
sitePercentages <- siteSummary
sitePercentages[, -c(ncol(sitePercentages))] <- sweep(sitePercentages[, -c(ncol(sitePercentages))], 1, sitePercentages[, ncol(sitePercentages)], "/")
sitePercentages <- sitePercentages[,c(1:(ncol(sitePercentages) - 1))] * 100
sitePercentages <- rownames_to_column(sitePercentages, var = "SITE")
rowSums(sitePercentages[,-1]) #sanity check

sitePercentages.longer <- sitePercentages %>% 
  pivot_longer(-SITE, names_to = "Study_Arm", values_to = "percent")

sitePercentages.longer$Study_Arm <- factor(sitePercentages.longer$Study_Arm,
                                           levels = c("Healthy Control", "De novo PD",
                                                      "Genetic Registry", "Genetic Cohort",
                                                      "Prodromal", "SWEDD"))
sitePalette <- diseasePalette[-c(2,4,7)]
sitePalette <- append(sitePalette, values = myPalette[c(13:14)])
names(sitePalette) <- c("Healthy Control", "De novo PD",
                        "Prodromal", "SWEDD",
                        "Genetic Registry", "Genetic Cohort")


sitePercent <- ggplot(sitePercentages.longer, aes(x = SITE, y = percent, fill = Study_Arm)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sitePalette) + 
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        #axis.text.y = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank()) +
  labs(x = "Site", y = "Enrollment \nGroup by Site (%)")
sitePercent

# insert size  ----------------------------------------------------------------------------------------------
platePalette <- colorRampPalette(c("#ffffff", myPalette[3]))(length(unique(meta$Plate)))
insertSize <- ggplot(meta.noPool, aes(x = Plate, y = MEAN_INSERT_SIZE, fill = Plate)) + 
  geom_boxplot() +
  scale_fill_manual(values = platePalette) + 
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        #axis.text.y = element_text(size = 14),
        legend.position = "none") +
  labs(x = "Plate ID", y = "Insert Size")
insertSize


#visit  ----------------------------------------------------------------------------------------------
labels <- c(BL = "Baseline", V02 = "6 Month Visit",
            V04 = "12 Month Visit", V06 = "24 Month Visit",
            V08 = "36 Month Visit")
visitPalette <- getPalette(12, "multi")[c(8:12)]
dev.off()
names(visitPalette) <- names(labels)
visit <- ggplot(meta.noPool, aes(x = Plate, y = freq, fill = VISIT)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ VISIT, nrow = 3, labeller = labeller(VISIT = labels)) +
  scale_fill_manual(values = visitPalette) + 
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        #axis.text.y = element_text(size = 14),
        legend.position = "none") +
  labs(x = "Plate ID", y = "Samples Per Plate")
visit

# study enrollment  ----------------------------------------------------------------------------------------------
meta.noPool$enrollArm <- meta.noPool$Study_Arm
meta.noPool$enrollArm <- gsub("GENUN|GENPD", "Genetic Cohort", meta.noPool$Study_Arm)
meta.noPool$enrollArm <- gsub("REGUN|REGPD", "Genetic Registry", meta.noPool$enrollArm)
meta.noPool$enrollArm <- gsub("PRODROMA", "Prodromal", meta.noPool$enrollArm)
meta.noPool$enrollArm <- gsub("PD", "De novo PD", meta.noPool$enrollArm)
meta.noPool$enrollArm <- gsub("HC", "Healthy Control", meta.noPool$enrollArm)
meta.noPool$enrollArm <- factor(meta.noPool$enrollArm, levels = c("Healthy Control", "De novo PD",
                                                                  "Genetic Registry", "Genetic Cohort",
                                                                  "Prodromal", "SWEDD"))
study <- ggplot(meta.noPool, aes(x = Plate, y = freq, fill = enrollArm)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ enrollArm, nrow = 3) +
  scale_fill_manual(values = sitePalette) + 
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        #axis.text.y = element_text(size = 14),
        legend.position = "none") +
  labs(x = "Plate ID", y = "Samples Per Plate")
study

#put figure together
## make figure
ggarrange(sitePercent,
          insertSize,
          study,
          visit,
          labels = c("a", "b", "c", "d"),
          font.label = list(size = 22),
          ncol = 2, nrow = 2)
ggsave("scientificData_paper/finalFigures/Figure2.png", width = 24, height = 12, dpi = 600)
ggsave("scientificData_paper/finalFigures/Figure2.eps", width = 24, height = 12, dpi = 600)
ggsave("scientificData_paper/finalFigures/Figure2.pdf", width = 24, height = 12, dpi = 600)
