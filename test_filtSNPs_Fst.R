
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
