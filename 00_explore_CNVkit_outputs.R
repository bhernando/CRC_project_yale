################################################################################
## This script aims to explore the otuputs of CNVkit
## Started: 22/10/25
## Updated: 22/10/25
##
## We have CRC samples from a TP53-KO mice model, and profiles have been derived
## with CNVkit. Here I just explore profiles derived by CNVkit.
##
################################################################################

### LIBRARIES
rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
# library(QDNAseq)
library(GenomicRanges)
library(foreach)
library(doMC)
library(YAPSA)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(lemon)
library(RColorBrewer)
library(igraph)
library(qgraph)
library(lsa)
library(dplyr)
library(mclust)
library(reshape2)

### PLOTTING THEME
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6))


### PATHS
BASE=dirname(this.path())
INPUT_DIR=file.path(BASE,"CNVkit_outputs")
OUTPUT_DIR=file.path(BASE,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(BASE,"plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)


### LOAD INPUT DATA

## RAW OUTPUT: log2 per probe
# Example of 1 sample
CNVkit = fread(file.path(INPUT_DIR, "/raw_output/MD7040c.dupmarked.1_19_XY.txt"))

CNVkit$length=CNVkit$end-CNVkit$start
CNVkit$offRegions=ifelse(CNVkit$gene=="Antitarget", TRUE, FALSE)
quantile(CNVkit$length[CNVkit$offRegions==T])
quantile(CNVkit$length[CNVkit$offRegions==F])
#       0%    25%    50%    75%   100%
# Off: 9374  16462  35194 132547 224986
# On:  120   172     201  232     399

dtLength=as.data.frame(melt(CNVkit[,8:9]))

png(file=paste0(PLOTS_DIR, "/BoxPlot_BinResolution.png"), width = 80/25.4, height = 80/25.4, units = "in", res = 300)
p = ggplot(dtLength, aes(x=offRegions, y=value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(pch = 21, color = "grey20", fill = "grey20", position = position_jitter(0.2), alpha = 0.4) +
    xlab("Off-target regions") + ylab("Length") +
    labs(title="", subtitle = "")+
    theme(text = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"))
print(p)
dev.off()

# Bin window is variable. The resolution is superhigh in on-target regions, which
# can add noise. The resolution in off-target regions goes from 9-220kb 


## ABSCN OUTPUT: log2 per probe
# Example of 1 sample
CNVkit = fread(file.path(INPUT_DIR, "/abs_output/MD7040c.dupmarked.1_19_XY.purityadj.call.cns"))

