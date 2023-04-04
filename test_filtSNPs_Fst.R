############No clones dendrogram?#############################
#{r, Dendrogram to identify clones}
#packages: tidyverse, ggdendro, paletteer
cloneBams = read.csv("data_files/poritesastreoidesMetaData_clonesremoved.csv")
cloneMa = as.matrix(read.table("data_files/pastNoClones.ibsMat"))
dimnames(cloneMa) = list(cloneBams[,1],cloneBams[,1])
clonesHc = hclust(as.dist(cloneMa),"ave")
clonePops = cloneBams$lineage
NocloneDend = cloneMa %>% as.dist() %>% hclust(.,"ave") %>% 
  as.dendrogram()
cloneDData = NocloneDend %>% dendro_data()

#Making the branches hang shorter so we can easily see clonal groups
cloneDData$segments$yend2 = cloneDData$segments$yend
for(i in 1:nrow(cloneDData$segments)) {
  if (cloneDData$segments$yend2[i] == 0) {
    cloneDData$segments$yend2[i] = (cloneDData$segments$y[i] - 0.01)}}

cloneDendPoints = cloneDData$labels
cloneDendPoints$pop = clonePops[order.dendrogram(NocloneDend)]
rownames(cloneDendPoints) = cloneDendPoints$label

# Making points at the leaves to place symbols for populations
point = as.vector(NA)
for(i in 1:nrow(cloneDData$segments)) {
  if (cloneDData$segments$yend[i] == 0) {
    point[i] = cloneDData$segments$y[i] - 0.01
  } else {
    point[i] = NA}}

cloneDendPoints$y = point[!is.na(point)]

kColPal3 = c("#AFDE62", "#FF8C8D", "mediumpurple3")

cloneDendPoints$pop = factor(cloneDendPoints$pop)
#cloneDendPoints$pop = factor(cloneDendPoints$pop,levels(cloneDendPoints$pop)[c(4,3,5,1,2)])

NocloneDendA = ggplot() +
  geom_segment(data = segment(cloneDData), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = cloneDendPoints, aes(x = x, y = y, fill = pop), size = 4, stroke = 0.25, shape = 24) +
  scale_fill_manual(values = kColPal3, name = "Lineage") +
#  geom_hline(yintercept = 0.165, color = "red", lty = 5, size = 0.75) + # creating a dashed line to indicate a clonal distance threshold
#  geom_text(data = subset(cloneDendPoints, subset = label %in% techReps), aes(x = (x + 0.53), y = (y - .0055), label = "*"), angle = 90, size = 10) + # spacing technical replicates further from leaf
  #  geom_text(data = subset(cloneDendPoints, subset = !label %in% techReps), aes(x = x, y = (y - .010), label = label), angle = 90) +
  labs(y = "Genetic distance (1 - IBS)") +
  theme_classic()

NocloneDend = NocloneDendA + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.line.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(size = 16, color = "black", angle = 90),
  axis.text.y = element_text(size = 14, color = "black"),
  axis.line.y = element_line(),
  axis.ticks.y = element_line(),
  panel.grid = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  plot.background = element_blank(),
  legend.key = element_blank(),
  legend.background = element_blank(),
  legend.title = element_text(size = 18),
  legend.text = element_text(size = 16),
  legend.position = "bottom")

NocloneDend

ggsave("figures/NocloneDend_GENOME.png", plot = NocloneDend, height = 5, width = 12, units = "in", dpi = 300)

###########Adding neutral structure in as factor in dbRDA??#####

##############can't get figure to look right
pastVcf = read.vcfR("data_files/pastFiltSnpsNoAdmixFINAL.vcf.recode.vcf", verbose = FALSE)
pastGenlightPopulation = vcfR2genlight(pastVcf, n.cores = 1)
locNames(pastGenlightPopulation) = paste(pastVcf@fix[,1],pastVcf@fix[,2],sep="_")

popData = read.csv("data_files/poritesastreoidesMetaData_clonesremoved_NOADMIX.csv") %>% select("sample" = tube_id, "pop" = lineage)

strata(pastGenlightPopulation) = data.frame(popData)
setPop(pastGenlightPopulation) = ~pop

pastGenlightPopulation$pop = factor(pastGenlightPopulation$pop)
levels(pastGenlightPopulation$pop)
pastGenlightPopulation$pop = factor(pastGenlightPopulation$pop,
                                    levels(pastGenlightPopulation$pop)[c(2, 1, 3)])

set.seed(694)

#99 permutations
sefl.fst <- stamppFst(pastGenlightPopulation, nboots = 99, percent = 95, nclusters = 2)
sefl.fst$Fsts

sefl.fst$Pvalues

#Generating heat map of pairwise Fst values
#reordering my samples to stay formatted properly for the matrix, and still have them go north to south
pop.order = c("Green", "Pink", "Purple")

#reads in Fst matrix
snpFstMa <- as.matrix(sefl.fst$Fsts)

#rebuilding the matrix based on order of populations
upperTriangle(snpFstMa, byrow = TRUE) <- lowerTriangle(snpFstMa)
snpFstMa <- snpFstMa[,pop.order] %>% .[pop.order,]
snpFstMa[upper.tri(snpFstMa)] <- NA
snpFstMa <- as.data.frame(snpFstMa)

snpFstMa$Pop = factor(row.names(snpFstMa))

snpQMa <- as.matrix(sefl.fst$Pvalues)
upperTriangle(snpQMa, byrow=TRUE) <- lowerTriangle(snpQMa)
snpQMa <- snpQMa[,pop.order] %>%
  .[pop.order,]
snpQMa[upper.tri(snpQMa)] <- NA
snpQMa <- as.data.frame(snpQMa)
snpQMa$Pop = factor(row.names(snpQMa), levels = unique(pop.order))

snpFstMa$Pop = factor(row.names(snpFstMa), levels = unique(pop.order))
snpFst = melt(snpFstMa, id.vars = "Pop", value.name = "Fst", variable.name = "Pop2", na.rm = FALSE)
snpFst$Fst = round(snpFst$Fst, 3)
snpFst = snpFst %>% mutate(Fst = replace(Fst, Fst < 0, 0))

snpQ = melt(snpQMa, id.vars = "Pop", value.name = "Pval", variable.name = "Pop2", na.rm = FALSE)
snpQ$Qval = p.adjust(snpQ$Pval, method = "BH")

snpFst$region = snpFst$Pop
snpFst$region = factor(gsub("\\n.*", "", snpFst$region))
snpFst$region = factor(snpFst$region, levels = levels(snpFst$region)[c(2, 1, 3)])

snpFst$region2 = snpFst$Pop2
snpFst$region2 = factor(gsub("\\n.*", "", snpFst$region2))
snpFst$region2 = factor(snpFst$region2, levels = levels(snpFst$region2)[c(2, 1, 3)])

snpFst$Fst = sprintf('%.3f', snpFst$Fst)
snpFst$Fst = factor(gsub("\\NA", NA, snpFst$Fst))
snpFst$Fst = factor(gsub("\\.000", "", snpFst$Fst))
snpFst$Fst = factor(gsub("\\-", "", snpFst$Fst))

kColPal3 = c("#FF8C8D", "#AFDE62", "mediumpurple3")

snpHeatmapA = ggplot(data = snpFst, aes(Pop, Pop2, fill = as.numeric(as.character(Fst))))+
  geom_tile(color = "white") +
  geom_segment(data = snpFst, aes(x = 0.48, xend = -0.43, y = Pop, yend = Pop, color = region), size = 23) + #edits y-axis pop titles, x changes x-axis width of color bars into the right/Fst side, xend changes width on the left side
  geom_segment(data = snpFst, aes(x = Pop2, xend = Pop2, y = 0.45, yend = -0.7, color = region2), size = 77.5) + #size changes width of color boxes
  scale_color_manual(values = kColPal3[], guide = NULL) +
  scale_fill_gradient(low = "white", high = "red", limit = c(0, 0.07), space = "Lab", name = expression(paste(italic("F")[ST])), na.value = "white",  guide = "colourbar") +
  geom_text(data = snpFst, aes(Pop, Pop2, label = Fst), color = ifelse(snpQ$Qval <= 0.05,"black", "darkgrey"), size = ifelse(snpQ$Qval < 0.05, 6, 5), fontface = ifelse (snpQ$Qval < 0.05, "bold", "plain")) +
  guides(fill = guide_colorbar(barwidth = 12, barheight = 1, title.position = "top", title.hjust = 0.5)) +
  scale_y_discrete(position = "left", limits = (levels(snpFst$Pop2))[c(1:2)]) +
  scale_x_discrete(limits = rev(levels(snpFst$Pop))[c(1:3)]) +
  coord_cartesian(xlim = c(1, 2), ylim = c(1, 2), clip = "off") +
  theme_minimal()

snpHeatmap = snpHeatmapA + theme(
  axis.text.x = element_text(vjust = 1, size = 16, hjust = 0.5, color = "black"),
  axis.text.y = element_text(size = 16, color = "black"),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_blank(),
  plot.background = element_rect(fill = "white"),
  axis.ticks = element_blank(),
  legend.position = c(0.7, 0.87),
  legend.direction = "horizontal",
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  plot.title = element_text(size = 16))

snpHeatmap


###################outlier SNPs###################
#Identifying outlier SNPs, indicative of some type of selection taking place
#{r, outlier}
# packages: tidyverse
bayescan = read.table("data_files/pastFilt.baye_fst.txt",header=T) %>% mutate(loc = rownames(.), out.05 = ifelse(qval < 0.05, 1, 0), out.1 = ifelse(qval < 0.1, 1, 0))
bayescan[bayescan[, 3]<=0.0001, 3] = 0.0001

bayescanPlotA = ggplot(data = bayescan, aes(x = log10(qval), y = fst, color = as.factor(out.05), alpha = as.factor(out.05))) +  
  geom_point(size = 1) +
#  geom_vline(xintercept = log10(0.05), linetype = 2, color = "purple") +
  xlab(expression(log[10]*"("*italic("q")*"-value)")) +
  ylab(expression(italic("F")[ST])) +
  scale_x_reverse() +
  scale_color_manual(values = c("grey45", "purple", "pink")) +
  scale_alpha_manual(values = c(0.25, 0.25, 0.5)) +
  theme_bw()

bayescanPlot = bayescanPlotA +
  theme(axis.title.x = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        axis.ticks.x = element_line(color = "black"),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 10),
        axis.ticks.y = element_line(color = "black"),
        legend.position = "none",
        legend.key.size = unit(0.3, 'cm'),
        panel.border = element_rect(color = "black"),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

bayescanPlot

##no outliers


###############Depth distribution of lineages??#############
### Depth distribtuion###

depthAov = aov(depthm ~ cluster, data = subset(pcangsd, subset = pcangsd$cluster!="Admixed"))

summary(depthAov)
TukeyHSD(depthAov)

chisq.test(x = subset(pcangsd, subset = pcangsd$cluster!="Admixed")$depth, y = subset(pcangsd, subset = pcangsd$cluster!="Admixed")$cluster)

lineageViolinA = ggplot(data = subset(pcangsd, subset = pcangsd$cluster!="Admixed"), aes(x = cluster, y = depthm, fill = cluster, group = cluster)) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 30, ymax = Inf,  fill = "black", alpha = 0.10, color = NA) +
  geom_beeswarm(shape = 21, size = 2, cex = 1.75, alpha = 1) +
  geom_violin(adjust = 1, linewidth = 0, color = "black", alpha = 0.35, width = 0.9, trim = F, scale = "width") +
  geom_violin(adjust = 1, linewidth = 0.4, color = "black", alpha = 1, width = 0.9, trim = F, fill = NA, scale = "width") +
  geom_boxplot(width = 0.2, color = "black", outlier.colour = NA, linewidth = 0.6, alpha = 0.5) +
  scale_fill_discrete(type = kColPal, name = "Lineage") +
  scale_color_discrete(type = kColPal, name = "Lineage") +
  xlab("Lineage") +
  ylab("Depth (m)") +
  scale_y_reverse(breaks = seq(10, 50, 5)) +
  theme_classic()
# theme_bw()

lineageViolin = lineageViolinA + theme(
  # axis.title.x = element_text(color = "black", size = 12),
  axis.title = element_text(color = "black", size = 16),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  # axis.title.y = element_text(color = "black", size = 12),
  # axis.text.x = element_text(color = "black", size = 10),
  axis.text.y = element_text(color = "black", size = 14),
  # axis.ticks.y = element_blank(),
  legend.position = "none",
  legend.key.size = unit(0.3, 'cm'),
  # panel.border = element_rect(color = "black"),
  panel.background = element_blank(),
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

# lineageViolin

ggsave("../figures/poster/lineage.svg", plot = lineageViolin, height = 4  , width = 6, units = "in", dpi = 600)
