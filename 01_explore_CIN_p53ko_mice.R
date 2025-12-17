################################################################################
## This script aims to analyse the copy number landscape of colitis-associated 
## colorectal cancer in a TP53-ko mice model
## Started: 22/10/25
## Updated: 26/11/25
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
# cntable needs to be structured as chrm, start, end, segValue, sample
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
    
    # Plot as continuous value
    unify=t(unify[,4:ncol(unify)])
    unify[unify>log2(15/2)]=(log2(15/2)) #max value
    unify[unify<log2(0.2/2)]=log2(0.2/2) #min value
    colMain = circlize::colorRamp2(c(min(unify), -1, 0, +1, max(unify)), c("#3da9cd", "#a9d9e9", "white", "#ff7272", "#fb0000"))
   
    #Dimensions
    if(is.null(dims)){heigth=10}
    if(!is.null(dims)){heigth=nrow(unify)/2}
    
    if(format == "png"){
        png(file=paste0(outdir, "/Heatmap_unify_profiles_",outname,".png"), width = 14, height = heigth, units = "in", res = 300)
        print(Heatmap(unify, col = colMain, column_title = paste0("Copy number profiles of ", length(unique(cntable$sample))," samples"), 
                      cluster_rows = cluster, cluster_columns = FALSE,
                      top_annotation = colAnn, show_column_names = FALSE, 
                      heatmap_legend_param = list(title = "relCN"),
                      show_row_names = names, row_dend_width = unit(2, "cm")))
        dev.off()
    }
    
    if(format == "pdf"){
        pdf(file=paste0(outdir, "/Heatmap_unify_profiles_",outname,".pdf"), width = 14, height = heigth)
        print(Heatmap(unify, col = colMain, column_title = paste0("Copy number profiles of ", length(unique(cntable$sample))," samples"), 
                      cluster_rows = cluster, cluster_columns = FALSE,
                      top_annotation = colAnn, show_column_names = FALSE, 
                      heatmap_legend_param = list(title = "relCN"),
                      show_row_names = names, row_dend_width = unit(2, "cm")))
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
AlBakir.metadata=rbind(fread(file.path(dirname(BASE), "extdata/AlBakir_Gut_2025/Discovery_samples_used.csv")),
                       fread(file.path(dirname(BASE), "extdata/AlBakir_Gut_2025/Validation_samples_used.csv")))

AlBakir.profiles=as.data.frame(do.call(rbind,lapply(1:length(AlBakir.data),function(thisSamp){
    mat=as.data.frame(t(AlBakir.data[[thisSamp]]))
    mat=as.data.frame(cbind(bins_4401,mat[,c(1,5)]))
    mat$sample=names(AlBakir.data[thisSamp])
    return(mat)
})))

# Genes of interest
favourite_genes=fread(file.path(BASE,"data/gene_names_pathways.txt"))
interest_genes=fread(file.path(BASE,"data/genes_of_interest_Dietlein2020.txt"))

# Chrom sizes 
mmd10_chr_sizes=fread(file.path(BASE,"data/mm10.chrom.sizes.txt"))


################################################################################
#### 1. EXPLORE CIN IN TP53-ko MICE MODEL ######################################

## Process copy number profiles
# We have absolute copy number values, but we want to transform that to relative copy number values (log2 transformed)
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

# avoid artefacts
CNVprofiles$cn[CNVprofiles$cn<0]=0

# relative copy number values
CNVprofiles$segVal=log2(pmax(CNVprofiles$cn, 0.1) / 2) # divided by ploidy (all diploid) & avoid infinite values
CNVprofiles=CNVprofiles[,c(1:3,7,5,6)]

## IMP! for classifying gains and losses
# 2 -> 0 
# 4 -> +1
# 1 -> -1

# Plot heatmap with relative values
plotHeatmap(cntable = CNVprofiles, bin.size = 100e3, chr_sizes = mmd10_chr_sizes,
            genome = "mmd10", path = "C:/Users/bhernando/Desktop/CNIO/Projects/scCINSignatures/scripts",
            names = T, cluster = T, format = "pdf", outdir = PLOTS_DIR, outname = "p53KOmice")



## Number of copy number alterations per sample
dtCNAs = as.data.frame(table(CNVprofiles$sample[CNVprofiles$segVal!=0]))
colnames(dtCNAs)=c("sample","nCNAs")
median(dtCNAs$nCNAs) #70
sd(dtCNAs$nCNAs) #47.55302
mean(dtCNAs$nCNAs) #76.96

# Add pathology information
dtCNAs$hist_grade=sample_pathology$Histological_grade[match(dtCNAs$sample,sample_pathology$Sanger_label)]
# dtCNAs$hist_grade[dtCNAs$hist_grade=="Invasive - mucinous"]="Invasive"
dtCNAs = dtCNAs[dtCNAs$hist_grade!="?",]
dtCNAs$cohort="TP53ko mice"

pdf(file=paste0(PLOTS_DIR, "/BoxPlot_nCNAs_per_histology.pdf"), width = 60/25.4, height = 60/25.4)
p=ggplot(dtCNAs, aes(x = hist_grade, y = nCNAs)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_point(color="black", fill="black", pch = 21, alpha = 0.8) +
    ggrepel::geom_text_repel(
        aes(label = sample),
        size = 2,
        max.overlaps = Inf,
        min.segment.length = 0,     # ensures the line is always drawn
        segment.color = "grey40",
        segment.size = 0.3
    ) +
    stat_compare_means(size=2)+
    xlab("Histology grade") + ylab("Number of CNAs (relVal !=2)") +
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
wgii <- CNVprofiles %>% mutate(length=end-start) %>% filter(segVal!=0) %>% group_by(sample) %>% summarise(wgii=sum(length/2.7e9))
# Add pathology information
wgii$hist_grade=sample_pathology$Histological_grade[match(wgii$sample,sample_pathology$Sanger_label)]
wgii = wgii[wgii$hist_grade!="?",]
wgii$cohort="TP53ko mice"


pdf(file=paste0(PLOTS_DIR, "/BoxPlot_FGA_per_histology.pdf"), width = 60/25.4, height = 60/25.4)
p=ggplot(wgii, aes(x = hist_grade, y = wgii)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_point(color="black", fill="black", pch = 21, alpha = 0.8) +
    ggrepel::geom_text_repel(
        aes(label = sample),
        size = 2,
        max.overlaps = Inf,
        min.segment.length = 0,     # ensures the line is always drawn
        segment.color = "grey40",
        segment.size = 0.3
    ) +
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
gene.locations=gene.locations[gene.locations$gene_source%in%c("ensembl","ensembl_havana"),]

segs.gr = GRanges(seqnames = CNVprofiles$chromosome, IRanges(start = CNVprofiles$start, end = CNVprofiles$end))
genes.gr <- GRanges(seqnames = gene.locations$chromosome, IRanges(start = gene.locations$start, end = gene.locations$end))
genes.gr$gene_name <- gene.locations$gene_name
overlaps <- findOverlaps(genes.gr, segs.gr)

dfCN.genes = as.data.frame(do.call(rbind,lapply(1:length(gene.locations$gene_name), function(i){
    idx = overlaps@to[overlaps@from == i]
    cn = CNVprofiles[idx,]
    vals = c()
    for(s in unique(cn$sample)){
        val=mean(cn$segVal[cn$sample==s])
        val=c(gene.locations$gene_name[i],s,val)
        vals=rbind(vals,val)
    }
    return(vals)
})))
colnames(dfCN.genes)=c("gene_name","sample","cnval")
dfCN.genes$cnval=as.numeric(dfCN.genes$cnval)
write.table(dfCN.genes, file.path(OUTPUT_DIR,"01_cnValues_per_genes_antonia_p53KOmice.txt"), quote=F, sep="\t", col.names = T, row.names = F)

# Prevalence per gene
dfCN.genes=fread(file.path(OUTPUT_DIR,"01_cnValues_per_genes_antonia_p53KOmice.txt"))
dfCN.genes=dfCN.genes[dfCN.genes$cnval <= -1 | dfCN.genes$cnval >= 1,] #remove no significant changes 
dfGoI=dfCN.genes[dfCN.genes$gene_name%in%favourite_genes$gene,]

# I classify alterations as gains (> +1 copies), amps (> +4 copies), losses (< -1 copy), deletions (0 copies)
dfGoI$status=ifelse(dfGoI$cnval>=log2(3/2) & dfGoI$cnval<log2(6/2), "Gain", ifelse(dfGoI$cnval>=log2(6/2), "Amplification",ifelse(dfGoI$cnval<log2(0.5/2), "Deletion", "Loss")))
dfGoI$status=factor(dfGoI$status, levels = c("Deletion","Loss","Gain","Amplification"))
dfGoI=dfGoI[order(dfGoI$gene_name),]
dfGoI$gene_name=factor(dfGoI$gene_name, levels = unique(dfGoI$gene_name))

dfGoI.stats=as.data.frame(dfGoI%>%group_by(gene_name,status, .drop=F)%>%summarise(n=n()))
dfGoI.stats$freq=dfGoI.stats$n/length(unique(CNVprofiles$sample))

# Plot prevalence
palette=c("#3da9cd", "#a9d9e9", "#ff7272", "#fb0000")
names(palette)=c("Deletion","Loss","Gain","Amplification")

pdf(file=paste0(PLOTS_DIR, "/BarPlot_GoI_frequency_antonia_p53KOmice.pdf"), width = 250/25.4, height = 80/25.4)
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
write.table(dfGoI, paste0(OUTPUT_DIR,"/favourite_genes_CNAs_per_sample.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)



## Pathway analyses => enrichment over background genes
dfCN.genes=fread(file.path(OUTPUT_DIR,"01_cnValues_per_genes_antonia_p53KOmice.txt"))
dfCN.genes$hit=0
dfCN.genes$hit[dfCN.genes$cnval <= -1 | dfCN.genes$cnval >= 1]=1 #remove no significant changes 

matHits=dcast(dfCN.genes, gene_name ~ sample, value.var = "hit") #25902 genes
row.names(matHits)=matHits$gene_name
matHits=matHits[,-1]
matHits_n=rowSums(matHits)   #counts
matHits_freq=rowSums(matHits) / ncol(matHits)   #frequency 0–1

paths=unique(favourite_genes$pathway)
path_enrich=as.data.frame(do.call(rbind,lapply(paths, function(thisPath) {
    print(thisPath)
    genes_pw=favourite_genes$gene[favourite_genes$pathway==thisPath]
    genes_bg=names(matHits_freq)[!names(matHits_freq)%in%favourite_genes$gene] # no favourite genes in the background
    
    # More samples mutated than random? frequency over background genes
    pw_freqs=matHits_freq[genes_pw]
    pw_freqs=pw_freqs[!is.na(pw_freqs)]
    bg_freqs=matHits_freq[genes_bg]
    bg_freqs=bg_freqs[!is.na(bg_freqs)]
    wtest=wilcox.test(pw_freqs, bg_freqs, alternative="greater")
    
    # More pathway genes altered than no-pathway genes?
    is_mut=matHits_n > 0
    yes_pw=sum(is_mut[genes_pw],na.rm = T)
    no_pw=length(genes_pw) - yes_pw
    yes_bg=sum(is_mut[genes_bg],na.rm = T)
    no_bg=length(genes_bg) - yes_bg
    tab=matrix(c(yes_pw,no_pw,yes_bg,no_bg), nrow = 2)
    ftest=fisher.test(tab, alternative = "greater")
    
    out=c(thisPath,wtest$p.value,mean(pw_freqs),mean(bg_freqs),
          ftest$p.value,yes_pw,no_pw)
    return(out)
})))
colnames(path_enrich)=c("pathway","wilcoxon.pval","mean_freq_pw","mean_freq_bg",
                        "fisher.pval","yes_ngenes_pw","no_ngenes_pw")
path_enrich[,2:7]=sapply(2:7, function(x)as.numeric(path_enrich[,x]))
path_enrich$wilcoxon.padj <- p.adjust(path_enrich$wilcoxon.pval, method="BH")
path_enrich$fisher.padj <- p.adjust(path_enrich$fisher.pval, method="BH")


## Gene enrichment analyses => enrichment over other favourite genes
dfCN.genes=fread(file.path(OUTPUT_DIR,"01_cnValues_per_genes_antonia_p53KOmice.txt"))
dfCN.genes$hit=0
dfCN.genes$hit[dfCN.genes$cnval <= -1 | dfCN.genes$cnval >= 1]=1 #remove no significant changes 
# dfCN.genes=dfCN.genes[dfCN.genes$gene_name%in%favourite_genes$gene,]

matHits=dcast(dfCN.genes, gene_name ~ sample, value.var = "hit") #145 genes
row.names(matHits)=matHits$gene_name
matHits=matHits[,-1]
matHits_n=rowSums(matHits)   #counts
matHits_freq=rowSums(matHits) / ncol(matHits)   #frequency 0–1

genes_bg=names(matHits_freq)[!names(matHits_freq)%in%favourite_genes$gene]
genes=unique(favourite_genes$gene)
genes=genes[!genes%in%c("Sfrp3","Il13ra1","Rbm10","Amer1","Pcdha3","Mcf2")] #no altered in any sample
gene_enrich=as.data.frame(do.call(rbind,lapply(genes, function(thisGene) {
    print(thisGene)
    # genes_bg=genes[genes!=thisGene]
    
    # More samples mutated than random? frequency over background genes
    pw_freqs=matHits_freq[thisGene]
    pw_freqs=pw_freqs[!is.na(pw_freqs)]
    bg_freqs=matHits_freq[genes_bg]
    bg_freqs=bg_freqs[!is.na(bg_freqs)]
    wtest=wilcox.test(pw_freqs, bg_freqs, alternative="greater")
    
    out=c(thisGene,wtest$p.value,mean(pw_freqs),mean(bg_freqs))
    return(out)
})))
colnames(gene_enrich)=c("gene","wilcoxon.pval","freq_gene","mean_freq_bg")
gene_enrich[,2:4]=sapply(2:4, function(x)as.numeric(gene_enrich[,x]))
gene_enrich$wilcoxon.padj <- p.adjust(gene_enrich$wilcoxon.pval, method="BH")


################################################################################
#### 2. EXPLORE CIN IN HUMAN DATA FROM Din et al 2018 ##########################

# How is the CNA burden in human data?
pval=wilcox.test(Din.data$nCNAs[Din.data$Hypermutated=="Hypermutated"], 
                 Din.data$nCNAs[Din.data$Hypermutated!="Hypermutated"])$p.value

pdf(file=paste0(PLOTS_DIR, "/BoxPlot_nCNAs_Din2018.pdf"), width = 80/25.4, height = 60/25.4)
p=ggplot(Din.data, aes(x = TP53, y = nCNAs)) +
    facet_wrap(~Hypermutated, ncol = 2) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(size=2)+
    geom_point(color = "black", fill = "black", pch = 21, position = position_jitter(0.2), alpha = 0.8) +
    xlab("TP53 mutation") + ylab("Number of CNAs") +
    scale_y_continuous(trans='log10') +
    labs(title="Din et al 2018",
         subtitle = "Yes vs No hypermutated: Wilcoxon, p = 7.61e-4")+
    theme(text = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"))
print(p)
dev.off()

pval=wilcox.test(Din.data$FGA[Din.data$Hypermutated=="Hypermutated"], 
                 Din.data$FGA[Din.data$Hypermutated!="Hypermutated"])$p.value

pdf(file=paste0(PLOTS_DIR, "/BoxPlot_FGA_Din2018.pdf"), width = 80/25.4, height = 60/25.4)
p=ggplot(Din.data, aes(x = TP53, y = FGA)) +
    facet_wrap(~Hypermutated, ncol = 2) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(size=2)+
    geom_point(color = "black", fill = "black", pch = 21, position = position_jitter(0.2), alpha = 0.8) +
    xlab("TP53 mutation") + ylab("Fraction of genome altered") +
    scale_y_continuous(trans='log10') +
    labs(title="Din et al 2018",
         subtitle = "Yes vs No hypermutated: Wilcoxon, p = 5.41e-4")+
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

# relative copy number values
AlBakir.segTabs.abs$segVal=log2(AlBakir.segTabs.abs$segVal) # not divided by the ploidy because it's already relative change => it only needs log transformation
AlBakir.segTabs.abs=AlBakir.segTabs.abs[AlBakir.segTabs.abs$sample!="3010_OPSC12-A5.3010_OPSC12-A5",]

# Set everything very close to 0 to 0
AlBakir.segTabs.abs$segVal[AlBakir.segTabs.abs$segVal > log2(1.5/2) & AlBakir.segTabs.abs$segVal < log2(2.5/2)] = 0

# Number of copy number alterations per sample
dtCNAs = as.data.frame(table(AlBakir.segTabs.abs$sample[AlBakir.segTabs.abs$segVal!=0]))
colnames(dtCNAs)=c("sample","nCNAs")
median(dtCNAs$nCNAs) #3
sd(dtCNAs$nCNAs) #5.51601
mean(dtCNAs$nCNAs) #5.30531

# Plot heatmap with profiles
# Remove samples with >1 copy number alterations
plotHeatmap(cntable = AlBakir.segTabs.abs[AlBakir.segTabs.abs$sample%in%dtCNAs$sample[dtCNAs$nCNAs>1],], 
            bin.size = 500e3, chr_sizes = CNpare:::chr_sizes,
            genome = "hg19", path = "C:/Users/bhernando/Desktop/CNIO/Projects/scCINSignatures/scripts",
            names = F, cluster = T, format = "pdf", outdir = PLOTS_DIR, outname = "AlBakir2025")


## How is the CNA burden in human data? => split by Progressors and NoProgressors
dtCNAs$sample.id=sub(".*\\.", "", dtCNAs$sample)
dtCNAs$progression=AlBakir.metadata$Progression[match(dtCNAs$sample.id,AlBakir.metadata$Sample)]

pdf(file=paste0(PLOTS_DIR, "/BoxPlot_nCNAs_AlBakir2025.pdf"), width = 60/25.4, height = 60/25.4)
p=ggplot(dtCNAs[!is.na(dtCNAs$progression),], aes(x = progression, y = nCNAs)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(color = "black", fill = "black", pch = 21, position = position_jitter(0.2), alpha = 0.8) +
    stat_compare_means(size=2)+
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
wgii <- AlBakir.segTabs.abs %>% mutate(length=end-start) %>% filter(segVal!=0) %>% group_by(sample) %>% summarise(wgii=sum(length/3e9))
wgii$sample.id=sub(".*\\.", "", wgii$sample)
wgii$progression=AlBakir.metadata$Progression[match(wgii$sample.id,AlBakir.metadata$Sample)]


pdf(file=paste0(PLOTS_DIR, "/BoxPlot_FGA_AlBakir2025.pdf"), width = 60/25.4, height = 60/25.4)
p=ggplot(wgii[!is.na(wgii$progression),], aes(x = progression, y = wgii)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_point(color = "black", fill = "black", pch = 21, position = position_jitter(0.2), alpha = 1) +
    stat_compare_means(size=2)+
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
    values = favourite_genes$human_gene,
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


# Prevalence per gene => overall cohort
# We are going to follow the same criteria for classifying alterations... but we don't know the ploidy (3/4 copies), amps (>4 copies), losses (1 copy), deletions (0 copies)
# I classify alterations as gains (>+1 rel copies), amps (> +2 rel copies), losses (-1 rel copy), deletions (< -2 copies)
dfGoI=dfCN.genes[dfCN.genes$cnval<= log2(1/2) | dfCN.genes$cnval>=log2(3/2),] #remove subclonal copy number changes with <1-copy change
dfGoI$status=ifelse(dfGoI$cnval>=log2(3/2) & dfGoI$cnval<log2(6/2), "Gain", ifelse(dfGoI$cnval>=log2(6/2), "Amplification",ifelse(dfGoI$cnval<log2(0.5/2), "Deletion", "Loss")))
dfGoI$status=factor(dfGoI$status, levels = c("Deletion","Loss","Gain","Amplification"))
dfGoI=dfGoI[order(dfGoI$gene_name),]
dfGoI$gene_name=factor(dfGoI$gene_name, levels = unique(dfGoI$gene_name))

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


dfGoI$sample.id=sub(".*\\.", "", dfGoI$sample)
dfGoI$progression=AlBakir.metadata$Progression[match(dfGoI$sample.id,AlBakir.metadata$Sample)]
dfGoI=dfGoI[!is.na(dfGoI$progression),]
dfGoI.stats=as.data.frame(dfGoI%>%group_by(progression,gene_name,status, .drop=F)%>%summarise(n=n()))

palette=c("#799FCB","#AFC7D0", "#FEC9C9", "#F9665E")
names(palette)=c("Deletion","Loss","Gain","Amplification")

png(file=paste0(PLOTS_DIR, "/BarPlot_GoI_frequency_Progression_AlBakir2025.png"), width = 250/25.4, height = 120/25.4, units = "in", res = 300)
p = ggplot(dfGoI.stats,aes(x=gene_name,y=n,fill=status))+
    facet_wrap(~progression,ncol=1)+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = palette) +
    labs(subtitle = "Al Bakir 2025",
         x = "Gene of interest", y = "Number of samples") +
    theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 6),
          legend.position = "bottom", legend.title = element_blank())
print(p)
dev.off()
write.table(dfGoI, paste0(OUTPUT_DIR,"/favourite_genes_CNAs_per_sample_AlBakir2025.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
