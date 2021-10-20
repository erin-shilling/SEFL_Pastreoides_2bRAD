library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggdendro)
cloneBams = read.csv("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/data_files/transcriptome_mapped/poritesastreoidesMetaData_trans.csv")

cloneMa = as.matrix(read.table("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/data_files/transcriptome_mapped/pastClones.ibsMat"))
dimnames(cloneMa) = list(cloneBams[,1],cloneBams[,1])
clonesHc = hclust(as.dist(cloneMa),"ave")
clonePops = cloneBams$region
cloneDend = cloneMa %>% as.dist() %>% hclust(.,"ave") %>% as.dendrogram()
cloneDData = cloneDend %>% dendro_data()

# Making the branches hang shorter so we can easily see clonal groups
cloneDData$segments$yend2 = cloneDData$segments$yend
for(i in 1:nrow(cloneDData$segments)) {
  if (cloneDData$segments$yend2[i] == 0) {
    cloneDData$segments$yend2[i] = (cloneDData$segments$y[i] - 0.01)}}

cloneDendPoints = cloneDData$labels
cloneDendPoints$pop = clonePops[order.dendrogram(cloneDend)]
rownames(cloneDendPoints) = cloneDendPoints$label

# Making points at the leaves to place symbols for populations
point = as.vector(NA)
for(i in 1:nrow(cloneDData$segments)) {
  if (cloneDData$segments$yend[i] == 0) {
    point[i] = cloneDData$segments$y[i] - 0.01
  } else {
    point[i] = NA}}

cloneDendPoints$y = point[!is.na(point)]

techReps = c("P028-1", "P028-2", "P028-3", "P046-1", "P046-2", "P046-3", "P077-1", "P077-2", "P077-3")

cloneDendPoints$pop = factor(cloneDendPoints$pop)
cloneDendPoints$pop = factor(cloneDendPoints$pop,levels(cloneDendPoints$pop)[c(5,3,6,2,4,1)])

cloneDendA = ggplot() +
  geom_segment(data = segment(cloneDData), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = cloneDendPoints, aes(x = x, y = y, fill = pop), size = 4, stroke = 0.25, shape = 24) +
  scale_fill_brewer(palette = "Dark2", name = "Population") +
  geom_hline(yintercept = 0.175, color = "red", lty = 5, size = 0.75) + # creating a dashed line to indicate a clonal distance threshold
  geom_text(data = subset(cloneDendPoints, subset = label %in% techReps), aes(x = x, y = (y - .015), label = label), angle = 90) + # spacing technical replicates further from leaf
  geom_text(data = subset(cloneDendPoints, subset = !label %in% techReps), aes(x = x, y = (y - .010), label = label), angle = 90) +
  labs(y = "Genetic distance (1 - IBS)") +
  theme_classic()

cloneDend = cloneDendA + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.line.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(size = 12, color = "black", angle = 90),
  axis.text.y = element_text(size = 10, color = "black"),
  axis.line.y = element_line(),
  axis.ticks.y = element_line(),
  panel.grid = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  plot.background = element_blank(),
  legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "bottom")

cloneDend

ggsave("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/figures/transcriptome_mapped/cloneDend_trans.png", plot = cloneDend, height = 8, width = 35, units = "in", dpi = 300)
ggsave("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/figures/transcriptome_mapped/cloneDend_trans.eps", plot = cloneDend, height = 8, width = 35, units = "in", dpi = 300)
