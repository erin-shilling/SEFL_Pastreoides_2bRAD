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
colPalZoox= c("chocolate4", "burlywood1","olivedrab","magenta4")
#colPalZoox = c("#D0D3D4", "#ADDC91", "#6A9FA1", "#F1A7DC")
#flPal = paletteer_d("rcartocolor::Sunset")[c(7, 6, 4, 3, 1)]
flPal = paletteer_d("LaCroixColoR::PeachPear")[c(1,2,3,4,5)]
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

###### now by resistance level)
# read in sample metadata
setwd("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/allie")
popDataClones = read.csv("rrcMetadataR_fix.csv") %>% select("sample" = sampleID, "pop" = Resistance)

#treat populations as a factor & order them North-->South
popDataClones$pop = factor(popDataClones$pop)
popDataClones$pop = factor(popDataClones$pop, levels = levels(popDataClones$pop)[c(2, 1, 3)])
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
#dfZoox$pop = as.factor(dfZoox$pop)
dfZoox$pop <- factor(dfZoox$pop, levels=c("High", "Medium", "Low"))
levels(dfZoox$pop)
# setting up for plotting
dfZoox = dfZoox[order(dfZoox$pop),]
sampleCounts = plyr::count(dfZoox, c('pop'))
meltedList = reshape2::melt(lapply(sampleCounts$freq,function(x){c(1:x)}))
dfZoox$barPlotOrder = meltedList$value
dfZoox = dfZoox[c(1,ncol(dfZoox),2:(ncol(dfZoox)-1))]
zDat = melt(dfZoox, id.vars = c("sample", "pop", "barPlotOrder"), variable.name = "Symbiont", value.name = "Fraction")
# colors
names(colPalZoox) = levels(zDat$Symbiont)
# plotting
colPalZoox= c("chocolate4", "burlywood1","olivedrab","magenta4")
#colPalZoox = c("#247EA3", "#FFBF46", "#6A9FA1", "Purple3")
flPal = paletteer_d("rcartocolor::Sunset")[c(7, 6, 4)]


popAnno = data.frame(x1 = c(0.5, 0.5, 0.5), x2 = c(31.5, 32.5, 32.5),
                     y1 = -0.065, y2 = -0.065, pop = c("High", "Medium", "Low"))
popAnno$pop = factor(popAnno$pop)
popAnno$pop = factor(popAnno$pop, levels = levels(popAnno$pop)[c(4, 2, 5, 1, 3)])
dfZoox = zDat %>% left_join(popAnno, by = "pop")
zooxSNPA = ggplot(data = dfZoox, aes(x = barPlotOrder, y = Fraction, fill = Symbiont, order = barPlotOrder)) +
  geom_bar(stat = "identity", position = "stack", colour = "grey25", width = 1, size = 0.2) +
  xlab("Population") +
  scale_x_discrete(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(-.001, -0.001)) +
  scale_color_manual(values = c("palegreen2", "khaki", "indianred")) +
  geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2, color = pop), size = 7) +
  scale_fill_manual(values = colPalZoox, name = "Symbiodiniaceae genus") +
  coord_cartesian(ylim = c(-.01,1.01), clip = "off") +
  facet_grid(~ pop, drop = TRUE, space = "free", scales = "free", switch = "both") +
  guides(color = "none") +
  theme_bw()
zooxSNPdis = zooxSNPA + theme(plot.title = element_text(),
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
zooxSNPdis
ggsave("zooxsdisease.png", plot = zooxSNPdis, width = 20, height = 8, units = "cm", dpi = 300)

###############Site Map###########
Map and shapefile data
Here we are setting up shape files and base maps as well as sample meta data.
```{r, map data}
library(tidyverse)
library(lubridate)
library(measurements)
library(sf)
library(rnaturalearth)
library(paletteer)
library(ggspatial)
#packages: tidyverse, lubridate, measurements, sf, rnaturalearth
setwd("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/allie")
ofavSamples = read.csv("RRC_LATLONG/updatedRRC_LatLong.csv", header = TRUE) %>% select(sample = sampleID, !sampleID)

ofavSamples$Site = factor(ofavSamples$Site)
ofavSamples$Site = factor(ofavSamples$Site, levels = levels(ofavSamples$Site)[c(3, 2, 5, 1, 4)])
ofavSamples$Resistance = factor(ofavSamples$Resistance, levels = levels(ofavSamples$Resistance)[c(3, 2, 1)])
#ofavSamples$collection_date = mdy(pastSamples$collection_date) %>% format("%d %b %Y")
#ofavSamples$depthM = conv_unit(pastSamples$collection_depth_ft, from = "ft", to = "m") %>% round(1)

ofavSites = ofavSamples %>% group_by(Site)%>% summarize(latDD = first(latDD), longDD = first(lonDD), n = n()) %>% droplevels()

states = st_as_sf(ne_states(country = c("United States of America")))
countries = st_as_sf(ne_countries(country = c("Cuba", "The Bahamas")))
florida = read_sf("shp/flCountiesLo.shp") %>% st_transform(crs = 4326)
frt = read_sf("shp/flReefs.shp") %>% st_transform(crs = 4326)

countyNames = st_as_sf(maps::map("county", plot = FALSE, fill = TRUE)) %>% st_transform(crs = 4326) %>%
filter(grepl("florida", ID)) %>% filter(grepl("broward|miami-dade", ID))
countyNames$ID = c("Broward \nCounty", "Miami- \nDade \nCounty")

popLabs = ofavSamples %>% group_by(Site)%>% summarize(latDD = first(latDD), longDD = first(lonDD))
levels(popLabs$Site) = c("NorthECA", "MiddleECA", "SouthECA", "Looe", "Sand")
  
##### build hi-res polygon
floridaMap = ggplot() +
  geom_sf(data = states, fill = "white", size = 0.25) +
  geom_sf(data = countries, fill = "white", size = 0.25) +
  geom_rect(aes(xmin = -81.9, xmax = -79.2, ymin = 24, ymax = 26), color = paletteer_d("vapoRwave::vapoRwave")[8], fill = NA) +
  coord_sf(xlim = c(-87, -77), ylim = c(23, 31)) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_minimal(), height = unit(1.2, "cm"), pad_x = unit(-0.25, "cm")) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(color = "black", size = 0.75, fill = NA),
        plot.background = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
floridaMap


#### final Map
library(magick)
library(cowplot)
flPal = paletteer_d("LaCroixColoR::PeachPear")[c(1,2,3,4,5)]

latLab = paste(format(abs(seq(24, 26.4, by = 0.2)), nsmall = 1), "ºN", sep = "")
longLab = paste(format(abs(seq(-82, -79, by = 0.2)), nsmall = 1), "ºW", sep = "")

siteMap = ggplot() +
  geom_sf(data = florida, fill = "white", color = "gray40", size = 0.25) +
  geom_sf(data = countries, fill = "white", color = "gray40", size = 0.25) +
  geom_sf(data = frt, fill = "gray80", size = 0) +
  geom_point(data = popLabs, aes(x = longDD, y = latDD, fill = Site), shape = 21, size = 4) +
  geom_sf_text(data = countyNames, aes(label = ID), nudge_x = c(0.1, 0, 0.1),
               nudge_y = c(0, -.025, -0.1)) +
  scale_fill_manual(values = flPal, name = "Population:") +
  coord_sf(xlim = c(-82.0, -79.0), ylim = c(24.0, 26.4)) +
  scale_x_continuous(breaks = c(seq(-82.0, -79.0, by = 0.2)), labels = longLab) +
  scale_y_continuous(breaks = c(seq(24.0, 26.4, by = 0.2)), labels = latLab) +
  annotation_scale(location = "br") +
  guides(fill = guide_legend(ncol = 3)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(color = "black", size = 0.75, fill = NA),
        plot.background = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")

#add a picture of O. faveolata
ofavPic = image_read("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/allie/Ofavpic1.jpg") %>%
  image_border("black", "10x10")

ofavMap = ggdraw() +
  draw_plot(siteMap) +
  draw_plot(floridaMap, x = 0.73, y = 0.655, width = 0.25, height = 0.25)+
  draw_image(ofavPic, x = 0.63, y = 0.21, width = 0.336, height = 0.336)
  
ofavMap
ggsave("ofavMap_test.png", plot = ofavMap, width = 12, height = 13, units = "cm", dpi = 300)
ggsave("C:/Users/erin_/Documents/GitHub/SEFL_Pastreoides_2bRAD/allie/ofavMapfix.png", plot = ofavMap, width = 12, height = 13, units = "cm", dpi = 300)