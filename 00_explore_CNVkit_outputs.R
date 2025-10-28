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


### LOAD INPUT DATA & PROCESS 

#### RAW OUTPUT: log2 per probe ####
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


#### ABSCN OUTPUT: log2 and abscn per copy number segment ####
# The output provides log2 values (with 95%CI) and also absolute copy number values
# The absCN are generated taking into account the purity and ploidy estimates by PureCN
# CNVkit has been run using a normal sample as control, so noise (and thus false segments) may be actually reduced a lot

## IMP! All samples are near-diploid (ploidies ~ 2)
## This is also the case for most cases in Din et al Clin Cancer Res, 2018 (Supp Table 6) => 

## Copy number profiles 
# Extract information
files=dir(paste0(INPUT_DIR,"/abs_output"),full=T,pattern="purityadj.call") #25 samples

CNVkit=as.data.frame(do.call(rbind,lapply(files, function(thisFile){
    cnv = fread(thisFile)
    cnv$sample = gsub(".dupmarked.1_19_XY.purityadj.call.cns","",basename(thisFile))
    return(cnv)
})))


# Construct segTables
CNVprofiles=CNVkit[,c(1:3,8,12)]
CNVprofiles$length=CNVprofiles$end-CNVprofiles$start



## Explore copy number profiles 
# I see segments with negative values... most of them are in sex chromosomes 
# that can be that they are not present at all, but it's true that fitting
# sex chromosomes is more challenging and are normally excluded 
CNVprofiles=CNVprofiles[!CNVprofiles$chromosome%in%c("X","Y"),]


# Distribution of cnvalues and segment size
p1 = ggplot(CNVprofiles[CNVprofiles$cn!=2,], aes(x=length, y=cn)) +
    geom_point(pch = 21, color = "grey20", fill = "grey20") +
    xlab("Segment size") + ylab("Segment copy number") +
    labs(title="TP53-ko colitis-CRC model (n=25)", subtitle = "segVals != 2")+
    scale_x_continuous(trans='log2') +
    theme(text = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"))

# long segments at 0 or negative cnvals may be not real... cells cannot survive => we normally use a threshold of 10Mb as maximum with 0 copies
# this can be because the ploidy is not well estimated? and should be higher?
# negative values can be noise or also because the ploidy should be higher
p2 = ggplot(CNVprofiles[CNVprofiles$cn<2,], aes(x=length, y=cn)) +
    geom_point(pch = 21, color = "grey20", fill = "grey20") +
    xlab("Segment size") + ylab("Segment copy number") +
    labs(title="TP53-ko colitis-CRC model (n=25)", subtitle = "segVals < 2")+
    scale_x_continuous(trans='log2') +
    theme(text = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"))

p3 = ggplot(CNVprofiles[CNVprofiles$cn>2,], aes(x=length, y=cn)) +
    geom_point(pch = 21, color = "grey20", fill = "grey20") +
    xlab("Segment size") + ylab("Segment copy number") +
    labs(title="TP53-ko colitis-CRC model (n=25)", subtitle = "segVals > 2")+
    scale_x_continuous(trans='log2') +
    theme(text = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"))

png(file=paste0(PLOTS_DIR, "/ScatterPlot_segVal_vs_segSize.png"), width = 80/25.4, height = 160/25.4, units = "in", res = 300)
p=ggpubr::ggarrange(p1,p2,p3,ncol=1,hjust=T)
print(p)
dev.off()


# Distribution per sample
png(file=paste0(PLOTS_DIR, "/BoxPlot_segVals_per_sample.png"), width = 160/25.4, height = 80/25.4, units = "in", res = 300)
p=ggplot(CNVprofiles[CNVprofiles$cn!=2,], aes(x = sample, y = cn)) + 
    # geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = sample), pch = 21, position = position_jitter(0.2), alpha = 1) +
    xlab("") + ylab("Segment copy number") +
    scale_y_continuous(trans='log10') +
    labs(title="TP53-ko colitis-CRC model (n=25)",
         subtitle = "segVals != 2")+
    theme(text = element_text(size = 7),
          axis.text = element_text(size = 6),
          axis.text.x = element_text(angle = 90, hjust=0.5, vjust=0.5),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"),
          legend.position = "none")
print(p)
dev.off()

# Distribution per sample
png(file=paste0(PLOTS_DIR, "/BoxPlot_segSize_per_sample.png"), width = 160/25.4, height = 80/25.4, units = "in", res = 300)
p=ggplot(CNVprofiles[CNVprofiles$cn!=2,], aes(x = sample, y = length)) + 
    # geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = sample), pch = 21, position = position_jitter(0.2), alpha = 1) +
    xlab("") + ylab("Segment size") +
    scale_y_continuous(trans='log10') +
    labs(title="TP53-ko colitis-CRC model (n=25)",
         subtitle = "segVals != 2")+
    theme(text = element_text(size = 7),
          axis.text = element_text(size = 6),
          axis.text.x = element_text(angle = 90, hjust=0.5, vjust=0.5),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"),
          legend.position = "none")
print(p)
dev.off()



## Number of copy number alterations per sample
dtCNAs = as.data.frame(table(CNVprofiles$sample[CNVprofiles$cn!=2]))
colnames(dtCNAs)=c("sample","nCNAs")
median(dtCNAs$nCNAs) #70
sd(dtCNAs$nCNAs) #47.55302
mean(dtCNAs$nCNAs) #76.96

dtCNAs$cohort="TP53ko mice"

png(file=paste0(PLOTS_DIR, "/BoxPlot_nCNAs.png"), width = 80/25.4, height = 100/25.4, units = "in", res = 300)
p=ggplot(dtCNAs, aes(x = cohort, y = nCNAs)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = sample), pch = 21, position = position_jitter(0.2), alpha = 1) +
    xlab("") + ylab("Number of CNAs (segVal !=2)") +
    scale_y_continuous(trans='log10') +
    labs(title="TP53-ko colitis-CRC model (n=25)",
         subtitle = "")+
    theme(text = element_text(size = 7),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"))
print(p)
dev.off()

## Fraction of genome altered per sample (mouse genome size ~ 2.7Gb)
wgii <- CNVprofiles %>% mutate(length=end-start) %>% filter(cn!=2) %>% group_by(sample) %>% summarise(wgii=sum(length/2.7e9))
wgii$cohort="TP53ko mice"

png(file=paste0(PLOTS_DIR, "/BoxPlot_FGA.png"), width = 80/25.4, height = 100/25.4, units = "in", res = 300)
p=ggplot(wgii, aes(x = cohort, y = wgii)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = sample), pch = 21, position = position_jitter(0.2), alpha = 1) +
    xlab("") + ylab("Fraction of genome altered") +
    scale_y_continuous(trans='log10') +
    labs(title="TP53-ko colitis-CRC model (n=25)",
         subtitle = "")+
    theme(text = element_text(size = 7),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"))
print(p)
dev.off()

# what is the CNAs in human data?
Din.data=fread(file.path(dirname(BASE), "extdata/Din_2018.txt"))
Din.data$Hypermutated=factor(Din.data$Hypermutated,levels = c("TRUE", "FALSE"), labels = c("Hypermutated", "No-Hypermutated"))

png(file=paste0(PLOTS_DIR, "/BoxPlot_nCNAs_Din2018.png"), width = 120/25.4, height = 100/25.4, units = "in", res = 300)
p=ggplot(Din.data, aes(x = TP53, y = nCNAs)) +
    facet_wrap(~Hypermutated, ncol = 2) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(size=2)+
    geom_point(color = "black", fill = "black", pch = 21, position = position_jitter(0.2), alpha = 0.8) +
    xlab("TP53 mutation") + ylab("Number of CNAs") +
    scale_y_continuous(trans='log10') +
    labs(title="Din et al 2018",
         subtitle = "")+
    theme(text = element_text(size = 7),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"))
print(p)
dev.off()


png(file=paste0(PLOTS_DIR, "/BoxPlot_FGA_Din2018.png"), width = 120/25.4, height = 100/25.4, units = "in", res = 300)
p=ggplot(Din.data, aes(x = TP53, y = FGA)) +
    facet_wrap(~Hypermutated, ncol = 2) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(size=2)+
    geom_point(color = "black", fill = "black", pch = 21, position = position_jitter(0.2), alpha = 0.8) +
    xlab("TP53 mutation") + ylab("Fraction of genome altered") +
    scale_y_continuous(trans='log10') +
    labs(title="Din et al 2018",
         subtitle = "")+
    theme(text = element_text(size = 7),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"))
print(p)
dev.off()


