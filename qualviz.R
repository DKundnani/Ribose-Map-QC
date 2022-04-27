#!/usr/bin/env Rscript

#Taking in arguments
library("optparse")
 
option_list = list(
  make_option(c("-q", "--QC"), type="character", default=NULL, 
              help="bed file in the format: chr\tstart\tstop\tid\tscore\tstrand \n Preprocess the bed file to filter the chormosomes you want to plot", metavar="filetype"),
  make_option(c("-l", "--libinfo"), type="character", default=NULL, 
              help="bed file in the format: chr\tstart\tstop\tid\tscore\tstrand \n Preprocess the bed file to filter the chormosomes you want to plot", metavar="filetype")
  #make_option(c("-b", "--bin_width"), type="integer", default=1e3, 
              #help="bin size/width,\n Recommended: 1e3 for yeast and 1e4 for human [default %default]", metavar="integer"),
  #make_option(c("-y", "--y_max"), type="integer", default=300, 
              #help="max limit of y axis for all chromosomes [default %default]", metavar="integer"),
  #make_option(c("-o", "--output_prefix"), type="character", default="out", 
              #help="output file name [default %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Checking for required input
writeLines("\n...Checking input...\n")
if (is.null(opt$QC)){
  print_help(opt_parser)
  stop("Please specify input quality file.n", call.=FALSE)
}

if (is.null(opt$libinfo)){
  print_help(opt_parser)
  stop("Please specify library info file.n", call.=FALSE)
}

writeLines(paste("Input file name:", opt$QC))
writeLines(paste("Library info:", opt$libinfo))
#writeLines(paste("Bin width used:", opt$bin_width))
#writeLines(paste("Maximum limit of yaxis:", opt$y_max))


#####Main Script
#1.Calling required Packages
writeLines("\n...Calling required Package(s)...\n")
library(ggplot2)
#library(gg.gap)
library(scales)
library(stringr) #wrapping Labels
library(tidyr)
library(ggpubr)

#2.Defining functions if any
writeLines("\n...Defining required functions...\n")

#3. Preprocessing input file
writeLines("...Preprocessing input file(s)...")
QC<-read.table(opt$QC,sep="\t",header=T)
libinfo<-read.table(opt$libinfo,sep="\t",header=T)

#4.Main code
writeLines("\n... Executing Main Code...\n")
ymax=round(max(QC$Counts_Nuclear/1e6))+1
perc=QC[,c("Library", "Counts_Nuclear", "Counts_Mito")]
perc=gather(perc, "Organelle", "counts",-Library)

#5.Plotting
#writeLines("...Plotting...")
writeLines("\n... plotting Raw counts...\n")

png("Raw Nuclear counts.png", width=1200, height=700)
ggplot(QC , aes(x=Library, y=Counts_Nuclear/1e6))+
  geom_bar(stat="identity", fill="steelblue")+
  scale_x_discrete(limits=QC$Library,labels = str_replace_all(QC$Library, "-", "\n"))+
  scale_y_continuous(expand=c(0,0), n.breaks=10, limits = c(0,ymax))+
  ggtitle("Ribonucleotides in Nucleus")+
  xlab("")+
  ylab("Raw ribonucleotide counts (x10^6)")+
  theme_classic()+
  theme(plot.title = element_text(size=30, face="bold",hjust = 0.5), 
        axis.text = element_text(size=20, colour = "black"), 
        axis.title = element_text(size=25))
dev.off()
png("Raw Mitochondria counts.png", width=1200, height=700)
ggplot(QC , aes(x=Library, y=Counts_Mito/1e6))+
  geom_bar(stat="identity", fill="steelblue")+
  scale_x_discrete(limits=QC$Library,labels = str_replace_all(QC$Library, "-", "\n"))+
  scale_y_continuous(expand=c(0,0), n.breaks=10, limits = c(0,ymax))+
  ggtitle("Ribonucleotides in Mitochondria")+
  xlab("")+
  ylab("Raw Mitochondria counts (x10^6)")+
  theme_classic()+
  theme(plot.title = element_text(size=30, face="bold",hjust = 0.5), 
        axis.text = element_text(size=20, colour = "black"), 
        axis.title = element_text(size=25))
dev.off()
png("Raw Nuclear counts_horizontal.png", width=1000, height=700)
ggplot(QC , aes(x=Library, y=Counts_Nuclear/1e6))+
  geom_bar(stat="identity", fill="steelblue")+
  coord_flip()+
  scale_x_discrete(limits=rev(QC$Library),labels= rev(QC$Library))+
  scale_y_continuous(expand=c(0,0), n.breaks=10, limits = c(0,ymax))+
  ggtitle("Ribonucleotides in Nucleus")+
  xlab("")+
  ylab("Raw ribonucleotide counts (x10^6)")+
  theme_classic()+
  theme(plot.title = element_text(size=30, face="bold",hjust = 0.5), 
        axis.text = element_text(size=20, colour = "black"), 
        axis.title = element_text(size=25))
dev.off()
png("Raw Mitochondria counts_horizontal.png", width=1000, height=700)
ggplot(QC , aes(x=Library, y=Counts_Mito/1e6))+
  geom_bar(stat="identity", fill="steelblue")+
  coord_flip()+
  scale_x_discrete(limits=rev(QC$Library),labels= rev(QC$Library))+
  scale_y_continuous(expand=c(0,0), n.breaks=10, limits = c(0,ymax))+
  ggtitle("Ribonucleotides in Mitochondria")+
  xlab("")+
  ylab("Raw Mitochondria counts (x10^6)")+
  theme_classic()+
  theme(plot.title = element_text(size=30, face="bold",hjust = 0.5), 
        axis.text = element_text(size=20, colour = "black"), 
        axis.title = element_text(size=25))
dev.off()




writeLines("\n... plotting Percentages...\n")
png("Nuclear and Mitochondria Percentage.png", width=1200, height=700)
ggplot(perc, aes(x=Library, y=counts, fill=Organelle))+
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(labels = c("Mitochondria", "Nucleus"),values=c("dodgerblue3","firebrick3"))+
  scale_x_discrete(limits=QC$Library,labels = str_replace_all(QC$Library, "-", "\n"))+
  scale_y_continuous(expand = c(0,0), labels = scales::percent_format())+
  ggtitle("Nuclear and Mitochondria Distribution")+
  ylab("Percentage")+
  xlab("")+
  theme_classic()+
  theme(plot.title = element_text(size=25, face="bold",hjust = 0.5), 
        axis.text = element_text(size=16, colour = "black"),axis.title = element_text(size=20),
        legend.text = element_text(colour="black", size=16), legend.title = element_text(colour="black", size=16))
dev.off()
png("Nuclear and Mitochondria Percentage_horizontal.png", width=1000, height=700)
ggplot(perc, aes(x=Library, y=counts, fill=Organelle))+
  geom_bar(stat="identity", position="fill")+
  coord_flip()+
  scale_fill_manual(labels = c("Mitochondria", "Nucleus"),values=c("dodgerblue3","firebrick3"))+
  scale_x_discrete(limits=rev(QC$Library),labels= rev(QC$Library))+
  scale_y_continuous(expand = c(0,0), labels = scales::percent_format())+
  ggtitle("Nuclear and Mitochondria Distribution")+
  ylab("Percentage")+
  xlab("")+
  theme_classic()+
  theme(plot.title = element_text(size=25, face="bold",hjust = 0.5), 
        axis.text = element_text(size=16, colour = "black"),axis.title = element_text(size=20),
        legend.text = element_text(colour="black", size=16), legend.title = element_text(colour="black", size=16))
dev.off()


writeLines("\n... plotting composition...\n")
png("Nuclear Grouped Composition.png", width=1200, height=700)
df=QC[,grep(pattern="Comp_nucl", x=colnames(QC))]
df$Group=libinfo$Group
data=gather(df,key="rNTP", value = "Percentage", -Group)
obj<-ggbarplot(data, x = "Group", y = "Percentage", fill="rNTP",
               add = c("mean_sd", "point"),
               #error.plot=c("crossbar"),
               add.params=list(width=0.4),
               palette = c("#D55E00","#0072B2", "#F0E442","#009E73"),
               position = position_dodge(0.75), 
               xlab="", ylab="Percentage", 
               width = 0.75,
               size = 1) +
  scale_y_continuous(expand = c(0,0),limits = c(0,100), breaks = seq(0,100,10)) +
  theme_classic(base_size = 30)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12))
ggpar(obj, legend=c("right"), legend.title = "")

dev.off()

png("Mitochondria Grouped Composition.png", width=1200, height=700)
df=QC[,grep(pattern="Comp_mito", x=colnames(QC))]
df$Group=libinfo$Group
data=gather(df,key="rNTP", value = "Percentage", -Group)
obj<-ggbarplot(data, x = "Group", y = "Percentage", fill="rNTP",
               add = c("mean_sd", "point"),
               #error.plot=c("crossbar"),
               add.params=list(width=0.4),
               palette = c("#D55E00","#0072B2", "#F0E442","#009E73"),
               position = position_dodge(0.75), 
               xlab="", ylab="Percentage", 
               width = 0.75,
               size = 1) +
  scale_y_continuous(expand = c(0,0),limits = c(0,100), breaks = seq(0,100,10)) +
  theme_classic(base_size = 30)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12))
ggpar(obj, legend=c("right"), legend.title = "")

dev.off()
png("Nuclear Grouped Normalized Composition.png", width=1200, height=700)
df=QC[,grep(pattern="Freq_nucl", x=colnames(QC))]
df$Group=libinfo$Group
data=gather(df,key="rNTP", value = "Percentage", -Group)
obj<-ggbarplot(data, x = "Group", y = "Percentage", fill="rNTP",
               add = c("mean_sd", "point"),
               #error.plot=c("crossbar"),
               add.params=list(width=0.4),
               palette = c("#D55E00","#0072B2", "#F0E442","#009E73"),
               position = position_dodge(0.75), 
               xlab="", ylab="Normalized Percentage", 
               width = 0.75,
               size = 1) +
  scale_y_continuous(expand = c(0,0),limits = c(0,100), breaks = seq(0,100,10)) +
  theme_classic(base_size = 30)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12))
ggpar(obj, legend=c("right"), legend.title = "")

dev.off()

png("Mitochondria Grouped  Normalized Composition.png", width=1200, height=700)
df=QC[,grep(pattern="Freq_mito", x=colnames(QC))]
df$Group=libinfo$Group
data=gather(df,key="rNTP", value = "Percentage", -Group)
obj<-ggbarplot(data, x = "Group", y = "Percentage", fill="rNTP",
               add = c("mean_sd", "point"),
               #error.plot=c("crossbar"),
               add.params=list(width=0.4),
               palette = c("#D55E00","#0072B2", "#F0E442","#009E73"),
               position = position_dodge(0.75), 
               xlab="", ylab="Normalized Percentage", 
               width = 0.75,
               size = 1) +
  scale_y_continuous(expand = c(0,0),limits = c(0,100), breaks = seq(0,100,10)) +
  theme_classic(base_size = 30)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12))
ggpar(obj, legend=c("right"), legend.title = "")

dev.off()




writeLines("\nJob Finished.\n")
writeLines(paste("Output files in current folder")