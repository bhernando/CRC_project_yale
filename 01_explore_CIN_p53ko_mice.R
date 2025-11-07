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
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))


### PATHS
BASE=dirname(this.path())
OUTPUT_DIR=file.path(BASE,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(BASE,"plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)


### FUNTIONS
# To unify boundaries & plot copy number profiles
plotHeatmap=function(cntable, bin.size, chr_sizes, genome, path, cluster, names, format, outdir, outname, dims=NULL){
    
    ## 1) Unify copy number segments across samples
    ## Unified matrix
    # Define parameters
    window=as.numeric(bin.size*2)
    
    # Unify segments
    source(paste0(path,"/uniqueEvents_functions.R")) #upload functions
    UNIFY.MAT<-unifySegments.multiple(profiles=cntable,window=window,bin.size=bin.size)
    
    # Scale chromosome sizes => Number of heatmap columns proportional to chr size
    unify=chrScale(UNIFY.MAT,window,chr_sizes)
    
    # Round copy number values
    unify[,4:ncol(unify)]=sapply(4:ncol(unify),function(x)round(unify[,x]))
    
    ## 2) Heatmap with copy number values
    # Set annotation
    ann=data.frame(unify$chromosome)
    colnames(ann)="chr"
    ann$chr=as.numeric(ann$chr)
    if(genome=="mmd10"){
        cols=c("#E046DC","#83A263","#64E399","#7E40E1","#E0AD4F","#E57155","#7BA9D9","#B56DDA","#D99CD8",
               "#878F8E","#DA9799","#C7E1E3","#CDE389","#6DC4DC","#E9D9B3","#85E658","#7C88DB","#DFC9E0",
               "#E55A9E") 
    }
    if(genome=="hg19"){
        cols=c("#E046DC","#83A263","#64E399","#7E40E1","#E0AD4F","#E57155","#7BA9D9","#B56DDA","#D99CD8",
               "#878F8E","#DA9799","#C7E1E3","#CDE389","#6DC4DC","#E9D9B3","#85E658","#7C88DB","#DFC9E0",
               "#E55A9E","#6CE8DD","#D7E149","#B0E3C0") 
    }
    cols=stats::setNames(cols,unique(ann$chr))
    colAnn=ComplexHeatmap::HeatmapAnnotation(chr=ann$chr,col=list(chr=cols))
    
    # Plot as discrete value
    unify=t(unify[,4:ncol(unify)])
    colMain=c("0"="#3182BD","1"="#9ECAE1","2"="#E5E5E5","3"="#FDCC8A","4"="#FC8D59","5"="#E34A33",
              "6"="#B30000","7"="#980043","8"="#DD1C77","9"="#DF65B0","10"="#C994C7","11"="#D4B9DA",
              "12"="#BEAEEA","13"="#cee393","14"="#9bb060","15"="#616e3c")
    unify[unify>14]=15
    unify[unify<0]=0
    
    #Dimensions
    if(is.null(dims)){heigth=10}
    if(!is.null(dims)){heigth=nrow(unify)/2}
    
    if(format == "png"){
        png(file=paste0(outdir, "/Heatmap_unify_profiles_",outname,".png"), width = 14, height = heigth, units = "in", res = 300)
        print(Heatmap(unify, col = colMain, column_title = paste0("Copy number profiles of ", length(unique(cntable$sample))," samples"), 
                      cluster_rows = cluster, cluster_columns = FALSE,
                      top_annotation = colAnn, show_column_names = FALSE, 
                      show_row_names = names, row_dend_width = unit(2, "cm")))
        dev.off()
    }
    
    if(format == "pdf"){
        pdf(file=paste0(outdir, "/Heatmap_unify_profiles_",outname,".pdf"), width = 14, height = heigth)
        print(Heatmap(unify, col = colMain, column_title = paste0("Copy number profiles of ", length(unique(cntable$sample))," samples"), 
                      cluster_rows = cluster, cluster_columns = FALSE,
                      top_annotation = colAnn, show_column_names = FALSE, 
                      show_row_names = cell_names, row_dend_width = unit(2, "cm")))
        dev.off()
    }
}

# Scale chrom sizes
chrScale=function(unify, window, chr_sizes){
    # Get chromosome sizes & split into fixed windows
    coordinates=as.data.frame(do.call(rbind,lapply(chr_sizes$chr, function(chr){
        start=seq(1,chr_sizes$length[chr_sizes$chr==chr],by=window)
        end=c(start[2:length(start)]-1,chr_sizes$length[chr_sizes$chr==chr])
        return(cbind(rep(chr,length(start)),cbind(start,end)))
    })))
    colnames(coordinates)=colnames(unify)[1:3]
    coordinates[,2:3]=sapply(2:3, function(x)as.numeric(coordinates[,x]))
    
    # Map unified segments with windows
    window.gr=GenomicRanges::GRanges(seqnames=coordinates$chromosome, IRanges(start=coordinates$start, end=coordinates$end))
    mat.gr=GenomicRanges::GRanges(seqnames=unify$chromosome, IRanges(start=unify$start, end=unify$end))
    overlaps=GenomicRanges::findOverlaps(window.gr, mat.gr)
    
    # Get copy number values
    df=matrix(data = NA, nrow = nrow(coordinates), ncol=ncol(unify)-3)
    colnames(df)=colnames(unify[4:ncol(unify)])
    for (i in 1:nrow(df)){
        idx = overlaps@to[overlaps@from == i]
        if(length(idx)!=0){
            df[i, ] <- colMeans(unify[idx, 4:ncol(unify)])
        }
    }
    df=cbind(coordinates,df)
    return(df[!is.na(df[,4]),])
}

# get copy number segments from bins
getSegments=function(bintable){
    segTable<-c()
    for(c in unique(bintable$chromosome))
    {
        sn<-bintable$segmented[bintable$chromosome==c]
        tab<-bintable[bintable$chromosome==c,]
        sn.rle<-rle(sn)
        starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
        ends <- cumsum(sn.rle$lengths)
        lapply(seq_len(length(sn.rle$lengths)), function(s) {
            from <- tab$start[starts[s]]
            to <- tab$end[ends[s]]
            segValue <- sn.rle$value[s]
            c(tab$chromosome[starts[s]], from, to, segValue)
        }) -> segtmp
        segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=TRUE),stringsAsFactors=FALSE)
        segTable<-rbind(segTable,segTableRaw)
    }
    colnames(segTable) <- c("chromosome", "start", "end", "segVal")
    segTable
}

### LOAD INPUT DATA & PROCESS 
# Copy number profiles of the mice model
CNVprofiles=readRDS(file.path(OUTPUT_DIR,"00_segment_tables_antonia_p53KOmice.rds"))
# Load metadata from the cohort
sample_metadata=fread(file.path(BASE,"data/sample_metadata.txt"))
sample_pathology=fread(file.path(BASE,"data/sample_pathology_metadata.txt"))
sample_plodies=fread(file.path(BASE,"CNVkit_outputs/all_purities.csv")) #all are near-diploid samples

# Copy number information from Din et al. 2018
Din.data=fread(file.path(dirname(BASE), "extdata/Din_2018.txt"))
Din.data$Hypermutated=factor(Din.data$Hypermutated,levels = c("TRUE", "FALSE"), labels = c("Hypermutated", "No-Hypermutated"))

# Copy number profiles from Al Bakir et al. 2025
# list with matrices per sample => information per 500kb bin
# the cn is relative => this 1=diploid, 0.5=single-copy loss, 1.5=single-copy gain and so on... we may need to transform this? --> no info of ploidy for samples
load(file.path(dirname(BASE),"extdata/AlBakir_Gut_2025/bin_locations_4401.Rdata")) 
AlBakir.data=c(readRDS(file.path(dirname(BASE), "extdata/AlBakir_Gut_2025/Colitis_Discovery_set_SS7GC.rds")),
               readRDS(file.path(dirname(BASE), "extdata/AlBakir_Gut_2025/Colitis_Validation_set_SS7GC.rds")))

AlBakir.profiles=as.data.frame(do.call(rbind,lapply(1:length(AlBakir.data),function(thisSamp){
    mat=as.data.frame(t(AlBakir.data[[thisSamp]]))
    mat=as.data.frame(cbind(bins_4401,mat[,c(1,5)]))
    mat$sample=names(AlBakir.data[thisSamp])
    return(mat)
})))

# Genes of interest
favourite_genes=fread(file.path(BASE,"data/gene_names_pathways.txt"))
# Chrom sizes 
mmd10_chr_sizes=fread(file.path(BASE,"data/mm10.chrom.sizes.txt"))


################################################################################
#### 1. EXPLORE CIN IN TP53-ko MICE MODEL ######################################

## Process copy number profiles
colnames(CNVprofiles)[4]="segVal"
CNVprofiles.ok=c()
for(s in unique(CNVprofiles$sample)){
    cnv=CNVprofiles[CNVprofiles$sample==s,]
    for(c in unique(cnv$chromosome)){
        cn=cnv[cnv$chromosome==c,]
        for(i in 1:(nrow(cn)-1)){
            cn$end[i]=ifelse(cn$end[i]==cn$start[i+1],cn$start[i+1]-1,cn$end[i])
        }
        CNVprofiles.ok=rbind(CNVprofiles.ok,cn)
    }
}
CNVprofiles=CNVprofiles.ok

# Plot heatmap with profiles
plotHeatmap(cntable = CNVprofiles, bin.size = 100e3, chr_sizes = mmd10_chr_sizes,
            genome = "mmd10", path = "C:/Users/bhernando/Desktop/CNIO/Projects/scCINSignatures/scripts",
            names = T, cluster = T, format = "png", outdir = PLOTS_DIR, outname = "p53KOmice")



## Number of copy number alterations per sample
dtCNAs = as.data.frame(table(CNVprofiles$sample[CNVprofiles$segVal!=2]))
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
wgii <- CNVprofiles %>% mutate(length=end-start) %>% filter(segVal!=2) %>% group_by(sample) %>% summarise(wgii=sum(length/2.7e9))
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
gene.locations=gene.locations[!duplicated(gene.locations$gene_name),]
GoI.locations=gene.locations[gene.locations$gene_name%in%favourite_genes$mouse_name,] #Ctnnb2 is not in mouse genome (only Ctnnb1)

segs.gr = GRanges(seqnames = CNVprofiles$chromosome, IRanges(start = CNVprofiles$start, end = CNVprofiles$end))
genes.gr <- GRanges(seqnames = GoI.locations$chromosome, IRanges(start = GoI.locations$start, end = GoI.locations$end))
genes.gr$gene_name <- GoI.locations$gene_name
overlaps <- findOverlaps(genes.gr, segs.gr)

dfCN.genes = as.data.frame(do.call(rbind,lapply(1:length(GoI.locations$gene_name), function(i){
    idx = overlaps@to[overlaps@from == i]
    cn = CNVprofiles[idx,]
    vals = c()
    for(s in unique(cn$sample)){
        val=mean(cn$segVal[cn$sample==s])
        val=c(GoI.locations$gene_name[i],s,val)
        vals=rbind(vals,val)
    }
    return(vals)
})))
colnames(dfCN.genes)=c("gene_name","sample","cnval")
dfCN.genes$cnval=as.numeric(dfCN.genes$cnval)
write.table(dfCN.genes, file.path(BASE,"01_cnValues_per_GoI_antonia_p53KOmice.txt"), quote=F, sep="\t", col.names = T, row.names = F)

# Prevalence per gene
dfCN.genes=fread(file.path(BASE,"01_cnValues_per_GoI_antonia_p53KOmice.txt"))

# All samples are diploid. So, I classify alterations as gains (3/4 copies), amps (>4 copies), losses (1 copy), deletions (0 copies)
dfGoI=dfCN.genes[dfCN.genes$cnval!=2,] #remove diploid 
dfGoI$status=ifelse(dfGoI$cnval>2 & dfGoI$cnval<5, "Gain", ifelse(dfGoI$cnval>=5, "Amplification",ifelse(dfGoI$cnval<=0, "Deletion", "Loss")))
dfGoI$status=factor(dfGoI$status, levels = c("Deletion","Loss","Gain","Amplification"))

dfGoI.stats=as.data.frame(dfGoI%>%group_by(gene_name,status, .drop=F)%>%summarise(n=n()))
dfGoI.stats$freq=dfGoI.stats$n/length(unique(CNVprofiles$sample))

palette=c("#799FCB","#AFC7D0", "#FEC9C9", "#F9665E")
names(palette)=c("Deletion","Loss","Gain","Amplification")

png(file=paste0(PLOTS_DIR, "/BarPlot_GoI_frequency_antonia_p53KOmice.png"), width = 250/25.4, height = 80/25.4, units = "in", res = 300)
p = ggplot(dfGoI.stats,aes(x=gene_name,y=n,fill=status))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = palette) +
    labs(subtitle = "TP53ko colitis-CRC mice model (n=25)",
         x = "Gene of interest", y = "Number of samples") +
    theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 6),
          legend.position = "bottom", legend.title = element_blank())
print(p)
dev.off()


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



################################################################################
#### 3. EXPLORE CIN IN HUMAN DATA FROM Al Bakir et al 2025 #####################

## Process copy number profiles
# We don't have info of ploidy...
# The cn is relative => this 1=diploid, 0.5=single-copy loss, 1.5=single-copy gain and so on...
# In addition, profiles are per bin (not per segment)

AlBakir.segTabs=c()
for(s in unique(AlBakir.profiles$sample)){
    bintable=AlBakir.profiles[AlBakir.profiles$sample==s,]
    segtable=getSegments(bintable)
    segtable$sample=s
    AlBakir.segTabs=rbind(AlBakir.segTabs,segtable)
}
AlBakir.segTabs.abs=AlBakir.segTabs
AlBakir.segTabs.abs$segVal=AlBakir.segTabs.abs$segVal*2
AlBakir.segTabs.abs=AlBakir.segTabs.abs[AlBakir.segTabs.abs$sample!="3010_OPSC12-A5.3010_OPSC12-A5",]

# Set everything very close to 2 to 2
AlBakir.segTabs.abs$segVal[AlBakir.segTabs.abs$segVal > (2-0.1) & AlBakir.segTabs.abs$segVal < (2+0.1)] = 2
AlBakir.segTabs.abs$segVal=round(AlBakir.segTabs.abs$segVal)

## Number of copy number alterations per sample
dtCNAs = as.data.frame(table(AlBakir.segTabs.abs$sample[AlBakir.segTabs.abs$segVal!=2]))
colnames(dtCNAs)=c("sample","nCNAs")
median(dtCNAs$nCNAs) #3
sd(dtCNAs$nCNAs) #5.601506
mean(dtCNAs$nCNAs) #5.374449


# Plot heatmap with profiles
# Remove samples with >1 copy number alterations
plotHeatmap(cntable = AlBakir.segTabs.abs[AlBakir.segTabs.abs$sample%in%dtCNAs$sample[dtCNAs$nCNAs>1],], 
            bin.size = 500e3, chr_sizes = CNpare:::chr_sizes,
            genome = "hg19", path = "C:/Users/bhernando/Desktop/CNIO/Projects/scCINSignatures/scripts",
            names = F, cluster = T, format = "png", outdir = PLOTS_DIR, outname = "AlBakir2025")


# How is the CNA burden in human data?
dtCNAs$cohort="AlBakir2025"
png(file=paste0(PLOTS_DIR, "/BoxPlot_nCNAs_AlBakir2025.png"), width = 80/25.4, height = 80/25.4, units = "in", res = 300)
p=ggplot(dtCNAs, aes(x = cohort, y = nCNAs)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(color = "black", fill = "black", pch = 21, position = position_jitter(0.2), alpha = 0.8) +
    xlab("") + ylab("Number of CNAs") +
    scale_y_continuous(trans='log10') +
    labs(title="Al Bakir et al 2025",
         subtitle = "")+
    theme(text = element_text(size = 7),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"))
print(p)
dev.off()


## Fraction of genome altered per sample (mouse genome size ~ 2.7Gb)
wgii <- AlBakir.segTabs.abs %>% mutate(length=end-start) %>% filter(segVal!=2) %>% group_by(sample) %>% summarise(wgii=sum(length/3e9))
wgii$cohort="AlBakir2025"

png(file=paste0(PLOTS_DIR, "/BoxPlot_FGA_AlBakir2025.png"), width = 80/25.4, height = 80/25.4, units = "in", res = 300)
p=ggplot(wgii, aes(x = cohort, y = wgii)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_point(color = "black", fill = "black", pch = 21, position = position_jitter(0.2), alpha = 1) +
    xlab("") + ylab("Fraction of genome altered") +
    scale_y_continuous(trans='log10') +
    labs(title="Al Bakir et al 2025",
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
library(biomaRt)
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 37)  # use GRCh37 as in the paper
# Query coordinates
GoI.locations <- getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
    filters = "hgnc_symbol",
    values = favourite_genes$human_name,
    mart = mart
)
colnames(GoI.locations)=c("gene_name","chromosome","start","end","strand")
GoI.locations=GoI.locations[GoI.locations$chromosome%in%c(1:22,"X","Y"),] #profiles from Al Bakir do not contain sex chromosomes

segs.gr = GRanges(seqnames = AlBakir.segTabs.abs$chromosome, IRanges(start = AlBakir.segTabs.abs$start, end = AlBakir.segTabs.abs$end))
genes.gr <- GRanges(seqnames = GoI.locations$chromosome, IRanges(start = GoI.locations$start, end = GoI.locations$end))
genes.gr$gene_name <- GoI.locations$gene_name
overlaps <- findOverlaps(genes.gr, segs.gr)

dfCN.genes = as.data.frame(do.call(rbind,lapply(1:length(GoI.locations$gene_name), function(i){
    idx = overlaps@to[overlaps@from == i]
    cn = AlBakir.segTabs.abs[idx,]
    vals = c()
    for(s in unique(cn$sample)){
        val=mean(cn$segVal[cn$sample==s])
        val=c(GoI.locations$gene_name[i],s,val)
        vals=rbind(vals,val)
    }
    return(vals)
})))
colnames(dfCN.genes)=c("gene_name","sample","cnval")
dfCN.genes$cnval=as.numeric(dfCN.genes$cnval)


# Prevalence per gene
# We are going to follow the same criteria for classifying alterations... but we don't know the ploidy (3/4 copies), amps (>4 copies), losses (1 copy), deletions (0 copies)
dfGoI=dfCN.genes[dfCN.genes$cnval!=2,] #remove no different from ploidy
dfGoI$status=ifelse(dfGoI$cnval>2 & dfGoI$cnval<5, "Gain", ifelse(dfGoI$cnval>=5, "Amplification",ifelse(dfGoI$cnval<=0, "Deletion", "Loss")))
dfGoI$status=factor(dfGoI$status, levels = c("Deletion","Loss","Gain","Amplification"))

dfGoI.stats=as.data.frame(dfGoI%>%group_by(gene_name,status, .drop=F)%>%summarise(n=n()))
dfGoI.stats$freq=dfGoI.stats$n/length(unique(CNVprofiles$sample))

palette=c("#799FCB","#AFC7D0", "#FEC9C9", "#F9665E")
names(palette)=c("Deletion","Loss","Gain","Amplification")

png(file=paste0(PLOTS_DIR, "/BarPlot_GoI_frequency_AlBakir2025.png"), width = 250/25.4, height = 80/25.4, units = "in", res = 300)
p = ggplot(dfGoI.stats,aes(x=gene_name,y=n,fill=status))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = palette) +
    labs(subtitle = "Al Bakir 2025 (n=227)",
         x = "Gene of interest", y = "Number of samples") +
    theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 6),
          legend.position = "bottom", legend.title = element_blank())
print(p)
dev.off()

