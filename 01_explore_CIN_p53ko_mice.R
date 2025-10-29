################################################################################
## This script aims to analyse the copy number landscape of colitis-associated 
## colorectal cancer in a TP53-ko mice model
## Started: 22/10/25
## Updated: 29/10/25
##
## We have CRC samples from a TP53-KO mice model, and profiles have been derived
## with CNVkit. We want to explore CNAs in this mice model and to compare with
## human data already published. 
##
## The mice model was developed by Antonia Churchhouse
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
OUTPUT_DIR=file.path(BASE,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(BASE,"plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)


### LOAD INPUT DATA & PROCESS 
# Copy number profiles of the mice model
CNVprofiles=readRDS(file.path(OUTPUT_DIR,"00_segment_tables_antonia_p53KOmice.rds"))
# Load metadata from the cohort
sample_metadata=fread(file.path(BASE,"data/sample_metadata.txt"))
sample_pathology=fread(file.path(BASE,"data/sample_pathology_metadata.txt"))

# Copy number information from Din et al. 2018
Din.data=fread(file.path(dirname(BASE), "extdata/Din_2018.txt"))
Din.data$Hypermutated=factor(Din.data$Hypermutated,levels = c("TRUE", "FALSE"), labels = c("Hypermutated", "No-Hypermutated"))

# Copy number profiles from Al Bakir et al. 2025
# list with matrices per sample => information per 500kb bin
# the cn is relative => this 1=diploid, 0.5=single-copy loss, 1.5=single-copy gain and so on... we may need to transform this?
load(file.path(dirname(BASE),"extdata/AlBakir_Gut_2025/bin_locations_4401.Rdata")) 
AlBakir.data=c(readRDS(file.path(dirname(BASE), "extdata/AlBakir_Gut_2025/Colitis_Discovery_set_SS7GC.rds")),
               readRDS(file.path(dirname(BASE), "extdata/AlBakir_Gut_2025/Colitis_Validation_set_SS7GC.rds")))

AlBakir.profiles=as.data.frame(do.call(rbind,lapply(1:length(AlBakir.data),function(thisSamp){
    mat=as.data.frame(t(AlBakir.data[[thisSamp]]))
    mat=as.data.frame(cbind(bins_4401,mat[,c(1,5)]))
    mat$sample=names(AlBakir.data[thisSamp])
    return(mat)
})))


################################################################################
#### 1. EXPLORE CIN IN TP53-ko MICE MODEL ######################################

## Number of copy number alterations per sample
dtCNAs = as.data.frame(table(CNVprofiles$sample[CNVprofiles$cn!=2]))
colnames(dtCNAs)=c("sample","nCNAs")
median(dtCNAs$nCNAs) #70
sd(dtCNAs$nCNAs) #47.55302
mean(dtCNAs$nCNAs) #76.96

# Add pathology information
dtCNAs$hist_grade=sample_pathology$Histological_grade[match(dtCNAs$sample,sample_pathology$Sanger_label)]
dtCNAs$cohort="TP53ko mice"

png(file=paste0(PLOTS_DIR, "/BoxPlot_nCNAs_per_histology.png"), width = 80/25.4, height = 80/25.4, units = "in", res = 300)
p=ggplot(dtCNAs, aes(x = hist_grade, y = nCNAs)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = hist_grade), pch = 21, position = position_jitter(0.2), alpha = 1) +
    stat_compare_means(size=2)+
    xlab("Histology grade") + ylab("Number of CNAs (segVal !=2)") +
    scale_y_continuous(trans='log10') +
    labs(title="TP53-ko colitis-CRC model (n=25)",
         subtitle = "")+
    theme(text = element_text(size = 7),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"),
          legend.position = "none")
print(p)
dev.off()

## Fraction of genome altered per sample (mouse genome size ~ 2.7Gb)
wgii <- CNVprofiles %>% mutate(length=end-start) %>% filter(cn!=2) %>% group_by(sample) %>% summarise(wgii=sum(length/2.7e9))
# Add pathology information
wgii$hist_grade=sample_pathology$Histological_grade[match(wgii$sample,sample_pathology$Sanger_label)]
wgii$cohort="TP53ko mice"


png(file=paste0(PLOTS_DIR, "/BoxPlot_FGA_per_histology.png"), width = 80/25.4, height = 80/25.4, units = "in", res = 300)
p=ggplot(wgii, aes(x = hist_grade, y = wgii)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = hist_grade), pch = 21, position = position_jitter(0.2), alpha = 1) +
    stat_compare_means(size=2)+
    xlab("Histology grade") + ylab("Fraction of genome altered") +
    scale_y_continuous(trans='log10') +
    labs(title="TP53-ko colitis-CRC model (n=25)",
         subtitle = "")+
    theme(text = element_text(size = 7),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"),
          legend.position = "none")
print(p)
dev.off()


## Which genes are impacted? Compute number of copies per gene
gene.locations=fread(file.path(dirname(BASE),"extdata/mm10_r102_geneCoordinates.tsv"))
gene.locations=gene.locations[!gene.locations$chromosome%in%c("X","Y","MT"),]
gene.locations=gene.locations[!duplicated(gene.locations$gene_name),]

segs.gr = GRanges(seqnames = CNVprofiles$chromosome, IRanges(start = CNVprofiles$start, end = CNVprofiles$end))
genes.gr <- GRanges(seqnames = gene.locations$chromosome, IRanges(start = gene.locations$start, end = gene.locations$end))
genes.gr$gene_name <- gene.locations$gene_name
overlaps <- findOverlaps(genes.gr, segs.gr)

dfCN.genes = as.data.frame(do.call(rbind,lapply(1:length(gene.locations$gene_name), function(i){
    idx = overlaps@to[overlaps@from == i]
    cn = CNVprofiles[idx,]
    vals = c()
    for(s in unique(cn$sample)){
        val=mean(cn$cn[cn$sample==s])
        val=c(gene.locations$gene_name[i],s,val)
        vals=rbind(vals,val)
    }
    return(vals)
})))
colnames(dfCN.genes)=c("gene_name","sample","cnval")
dfCN.genes$cnval=as.numeric(dfCN.genes$cnval)
write.table(dfCN.genes, file.path(BASE,"01_cnValues_per_gene_antonia_p53KOmice.txt"), quote=F, sep="\t", col.names = T, row.names = F)

# Prevalence per gene
dfCN.genes=fread(file.path(BASE,"01_cnValues_per_gene_antonia_p53KOmice.txt"))
dfGoI=dfCN.genes[dfCN.genes$cnval!=2,]

# Gains
dfGoI.gains=dfGoI[dfGoI$cnval>2,]
dfGoI.gains=dfGoI.gains%>%group_by(gene_name)%>%summarise(n=n())
dfGoI.gains$freq=dfGoI.gains$n/length(unique(CNVprofiles$sample))

# Losses
dfGoI.losses=dfGoI[dfGoI$cnval<2,]
dfGoI.losses=dfGoI.losses%>%group_by(gene_name)%>%summarise(n=n())
dfGoI.losses$freq=dfGoI.losses$n/length(unique(CNVprofiles$sample))


################################################################################
#### 2. EXPLORE CIN IN HUMAN DATA FROM Din et al 2018 ##########################

# How is the CNA burden in human data?
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


