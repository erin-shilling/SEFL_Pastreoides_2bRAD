# packages: tidyverse, reshape2, RColorBrewer, paletteer
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(paletteer)
# read in sample metadata with clones, remove technical replicates
setwd("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/allie")
popDataClones = read.csv("rrcMetadataR_fix.csv") %>% select("sample" = sampleID, "pop" = Site)
#treat populations as a factor & order them North-->South
popDataClones$pop = factor(popDataClones$pop)
popDataClones$pop = factor(popDataClones$pop, levels = levels(popDataClones$pop)[c(3, 2, 5, 4, 1)])
#read in zoox reads data with clones incl, technical replicates removed
zoox = read.delim("zooxReads", header = FALSE, check.names = FALSE)
head(zoox)
zoox$V2[is.na(zoox$V2)] <- as.character(zoox$V1[is.na(zoox$V2)])
zoox$V1 = gsub("P.*", "chr", zoox$V1)
zoox$V2 = gsub(".trim.*", "", zoox$V2)
zoox = zoox %>% filter(zoox$V1 != "*")

zooxLst = split(zoox$V2, as.integer(gl(length(zoox$V2), 20, length(zoox$V2))))

zooxMaps = NULL

for(i in zooxLst){
  zooxMaps = rbind(zooxMaps, data.frame(t(i)))
}

colnames(zooxMaps) = c("sample", zoox$V1[c(2:20)])

for(i in c(2:20)){
  zooxMaps[,i] = as.numeric(zooxMaps[,i])
}
str(zooxMaps)
zooxMaps$Symbiodinium = rowSums(zooxMaps[2:6])
zooxMaps$Breviolum = rowSums(zooxMaps[7:10])
zooxMaps$Cladocopium = rowSums(zooxMaps[11:16])
zooxMaps$Durusdinium = rowSums(zooxMaps[17:20])

zooxMaps = zooxMaps[,c(1, 21:24)]
zooxProp = zooxMaps

zooxProp$sum = apply(zooxProp[, c(2:length(zooxProp[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})

#turn them into proportions
zooxProp = cbind(zooxProp$sample, (zooxProp[, c(2:(ncol(zooxProp)-1))]
                                   / zooxProp$sum))
colnames(zooxProp)[1] = "sample"
head(zooxProp)
apply(zooxProp[, c(2:(ncol(zooxProp)))], 1, function(x) {
  sum(x, na.rm = T)
})

#bind it to metadata
dfZoox = popDataClones %>% left_join(zooxProp)
dfZoox$pop = as.factor(dfZoox$pop)
levels(dfZoox$pop)

# setting up for plotting
dfZoox = dfZoox[order(dfZoox$pop),]
sampleCounts = plyr::count(dfZoox, c('pop'))
meltedList = reshape2::melt(lapply(sampleCounts$freq,function(x){c(1:x)}))
dfZoox$barPlotOrder = meltedList$value
dfZoox = dfZoox[c(1,ncol(dfZoox),2:(ncol(dfZoox)-1))]
zDat = melt(dfZoox, id.vars = c("sample", "pop", "barPlotOrder"), variable.name = "Symbiont", value.name = "Fraction")

# plotting
colPalZoox = c("#D0D3D4", "#ADDC91", "#6A9FA1", "#F1A7DC")
flPal = paletteer_d("rcartocolor::Sunset")[c(7, 6, 4, 3, 1)]

popAnno = data.frame(x1 = c(0.5, 0.5, 0.5, 0.5, 0.5), x2 = c(14.5, 18.5, 14.5, 23.5, 26.5),
                     y1 = -0.065, y2 = -0.065, pop = c("NorthECA", "MiddleECA", "SouthECA", "Sand", "Looe"))

popAnno$pop = factor(popAnno$pop)
popAnno$pop = factor(popAnno$pop, levels = levels(popAnno$pop)[c(3, 2, 5, 4, 1)])

dfZoox = zDat %>% left_join(popAnno, by = "pop")

zooxSNPA = ggplot(data = dfZoox, aes(x = barPlotOrder, y = Fraction, fill = Symbiont, order = barPlotOrder)) +
  geom_bar(stat = "identity", position = "stack", colour = "grey25", width = 1, size = 0.2) +
  xlab("Population") +
  scale_x_discrete(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(-.001, -0.001)) +
  scale_color_manual(values = flPal) +
  geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2, color = pop), size = 7) +
  scale_fill_manual(values = colPalZoox, name = "Symbiodiniaceae genus") +
  coord_cartesian(ylim = c(-.01,1.01), clip = "off") +
  facet_grid(~ pop, drop = TRUE, space = "free", scales = "free", switch = "both") +
  guides(color = "none") +
  theme_bw()

zooxSNP = zooxSNPA + theme(plot.title = element_text(),
                           panel.grid = element_blank(),
                           panel.background = element_rect(fill = "gray25", colour = "grey25"),
                           panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
                           panel.spacing.x = grid:::unit(0.05, "lines"),
                           panel.spacing.y = grid:::unit(0.05, "lines"),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title = element_blank(),
                           strip.background.x = element_blank(),
                           strip.background.y = element_blank(),
                           strip.text = element_text(size = 10),
                           strip.text.y.left = element_text(size = 10, angle = 90),
                           strip.text.x.bottom = element_text(vjust = -.05, color = "black"),
                           legend.key.size = unit(0.75, "line"),
                           legend.title = element_text(size = 10),
                           legend.text = element_text(size = 8),
                           legend.key = element_blank(),
                           legend.position = "bottom")

# zooxSNP figure with clones
zooxSNP
ggsave("zooxs.png", plot = zooxSNP, width = 20, height = 8, units = "cm", dpi = 300)

