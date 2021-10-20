library(tidyverse)
if (!require("pacman")) install.packages("pacman")

pacman::p_load_gh("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", "ropensci/rnaturalearthhires")
pacman::p_load("cowplot", "ggrepel", "ggspatial", "maps", "paletteer", "patchwork", "rgdal", "rnaturalearth", "sf", "tidyverse", "reshape2", "MCMC.OTU", "pairwiseAdonis", "RColorBrewer", "Redmonder", "flextable", "lubridate", "officer", "adegenet", "dendextend", "gdata", "ggdendro", "hierfstat", "Imap", "kableExtra", "poppr", "reshape2", "StAMPP", "vcfR", "vegan", "boa", "measurements", "magick", "rgeos")

#```{r, Analysis of Molecular Variance}
#reading in bcf file
pastVcf = read.vcfR("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/data_files/transcriptome_mapped/pastNoClones.bcf", verbose = TRUE)
#convert to genlight files for poppr
pastGenlightPopulation = vcfR2genlight(pastVcf, n.cores = 1)

#taking metadata file, without technical replicates & clones, reads in population data for each sample
popData = read.csv("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/data_files/transcriptome_mapped/poritesastreoidesMetaData_clonesremoved_trans.csv") %>% select("sample" = tube_id, "pop" = region)

#setting up amova
strata(pastGenlightPopulation) = data.frame(popData)
setPop(pastGenlightPopulation) = ~pop
#Runs AMOVA looking at samples by region
amova <- poppr.amova(pastGenlightPopulation, ~pop)
amova

set.seed(694)
amovasignif <- randtest(amova, nrepet = 99)
amovasignif$names

amovasignif$obs

amovasignif$pvalue

amovaPerc = paste(round(amova$componentsofcovariance$`%`[1], 2), "%",sep="")
amovaP = amovasignif$pvalue[3]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#```{r, PCoA with IBS}
pastMa = as.matrix(read.table("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/data_files/transcriptome_mapped/pastNoClones.ibsMat"))
pastMds = cmdscale(pastMa, eig = TRUE, x.ret = TRUE)

# Determine percent variation captured on each axis
# Calculate the eigenvalues so later we can figure out % variation shown on each Principal Coordinate
pastPcoaVar = round(pastMds$eig/sum(pastMds$eig)*100, 1)
head(pastPcoaVar)
# Format data to plot
pastPcoaValues = pastMds$points
head(pastPcoaValues)
# 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
pastI2P = read.csv("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/data_files/transcriptome_mapped/poritesastreoidesMetaData_clonesremoved_trans.csv") %>% select("sample" = tube_id, "pop" = region)

row.names(pastI2P) = pastI2P[,1]
pastPcoaValues=cbind(pastI2P, pastPcoaValues)
pastPcoaValues =as.data.frame(pastPcoaValues, sample = rownames(pastPcoaValues))
colnames(pastPcoaValues)[c(3,4)] = c("PCo1", "PCo2")
head(pastPcoaValues)

pastPCoA = merge(pastPcoaValues, aggregate(cbind(mean.x=PCo1,mean.y=PCo2)~pop, pastPcoaValues, mean), by="pop")

pastPCoA$pop = factor(pastPCoA$pop)
pastPCoA$pop = factor(pastPCoA$pop, levels(pastPCoA$pop)[c(4, 2, 5, 1, 3)])

flPal = paletteer_d("rcartocolor::Sunset")[c(7, 6, 4, 3, 1)]

pastPcoaPlotA = ggplot(pastPCoA, aes(x = PCo1, y = PCo2, color = pop, fill = pop)) +
  geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
  geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
  #ellipses
  stat_ellipse(data = pastPCoA, type = "t", geom = "polygon", alpha = 0.1) +
  #individual's indicated by small transparent circles
  geom_point(aes(x = PCo1, y = PCo2), size = 3, alpha = 0.3, show.legend = FALSE) +
  #population centroids indicated by large circles
  geom_point(aes(x = mean.x, y = mean.y), size = 5, color = "black", shape = 21) +
  annotate(geom = "text", x = 0.1, y = -0.25, label = bquote("AMOVA:"~.(amovaPerc)*","~italic(p)~" = "~.(amovaP))) +
  scale_fill_manual(values = flPal, name = "Region") +
  scale_color_manual(values = flPal, guide = NULL) +
  xlab(paste ("PCo 1 (", pastPcoaVar[1],"%)", sep = "")) + #Prints percent variation explained by first axis
  ylab(paste ("PCo 2 (", pastPcoaVar[2],"%)", sep = "")) + #Prints percent variation explained by second axis
  guides(shape = guide_legend(order = 2), linetype = guide_legend(override.aes = list(linetype = c(1,2), alpha = 1, color = "black", fill = NA), order = 3), fill = guide_legend(override.aes = list(shape = 22, size = 5, color = NA, alpha = NA), order = 1))+
  theme_bw()

pastPcoaPlot = pastPcoaPlotA +
  theme(axis.title.x = element_text(color = "black", size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "right",
        panel.border = element_rect(color = "black", size = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# pastPcoaPlot

ggsave("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/figures/transcriptome_mapped/pcoaPlot_trans.png", plot = pastPcoaPlot, height = 3.5, width = 7, units = "in", dpi = 300)
ggsave("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/figures/transcriptome_mapped/pcoaPlot_trans.pdf", plot = pastPcoaPlot, height = 3.5, width = 7, units = "in", dpi = 300)
ggsave("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/figures/transcriptome_mapped/pcoaPlot_trans.tiff", plot = pastPcoaPlot, height = 3.5, width = 7, units = "in", dpi = 300)
