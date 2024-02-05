#!/usr/bin/env Rscript

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
if(!require(weatherData)) {install.packages("ggbreak", repos = "http://cran.us.r-project.org")}
library(ggbreak) 
library(cowplot)
library(grid)
library(gridExtra) 

#2.Defining functions if any
writeLines("\n...Defining required functions...\n")

#3. Preprocessing input file
writeLines("...Preprocessing input file(s)...")
QC<-read.table(opt$QC,sep="\t",header=T)
libinfo<-read.table(opt$libinfo,sep="\t",header=T)
rownames(QC)=QC$Library
QC=QC[libinfo$Library,]
QC$Library<-libinfo$Label

#4.Main code
writeLines("\n... Executing Main Code...\n")
ymax=round(max(QC$Counts_Nuclear/1e6))+1
perc=QC[,c("Library", "Counts_Nuclear", "Counts_Mito")]
groupedcounts=perc
groupedcounts$Group=libinfo$Group
groupedcounts$Nuclperc=groupedcounts$Counts_Nuclear/(groupedcounts$Counts_Nuclear+groupedcounts$Counts_Mito)
groupedcounts$Mitoperc=groupedcounts$Counts_Mito/(groupedcounts$Counts_Nuclear+groupedcounts$Counts_Mito)
perc=gather(perc, "Genome", "counts",-Library)


#5.Plotting
#writeLines("...Plotting...")
writeLines("\n... plotting Raw counts...\n")

png("Raw Nuclear counts.png", width=12, height=7,units= "in",  res=600)
ggplot(QC , aes(x=Library, y=Counts_Nuclear/1e6))+
  geom_bar(stat="identity", fill="steelblue")+
  scale_x_discrete(limits=QC$Library,labels = str_replace_all(QC$Library, "-", "\n"))+
  scale_y_continuous(expand=c(0,0), n.breaks=10, limits = c(0,ymax))+
  ggtitle("Ribonucleotides in Nucleus")+
  xlab("")+
  ylab("Raw ribonucleotide counts (x10^6)")+
  theme_classic()+
  theme(plot.title = element_text(size=30, face="bold",hjust = 0.5), 
        axis.text = element_text(size=16, colour = "black"), 
        axis.title = element_text(size=25))
dev.off()


png("Raw Mitochondria counts.png", width=12, height=7,units= "in",  res=600)
ggplot(QC , aes(x=Library, y=Counts_Mito/1e6))+
  geom_bar(stat="identity", fill="steelblue")+
  scale_x_discrete(limits=QC$Library,labels = str_replace_all(QC$Library, "-", "\n"))+
  scale_y_continuous(expand=c(0,0), n.breaks=10, limits = c(0,ymax))+
  ggtitle("Ribonucleotides in Mitochondria")+
  xlab("")+
  ylab("Raw Mitochondria counts (x10^6)")+
  theme_classic()+
  theme(plot.title = element_text(size=30, face="bold",hjust = 0.5), 
        axis.text = element_text(size=16, colour = "black"), 
        axis.title = element_text(size=25))
dev.off()

png("Raw Nuclear counts_horizontal.png", width=12, height=7,units= "in",  res=600)
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
        axis.text = element_text(size=16, colour = "black"), 
        axis.title = element_text(size=25))
dev.off()

png("Raw Nuclear counts_groupedhorizontal.png", width=18, height=7,units= "in",  res=600, bg="transparent")
ggplot(groupedcounts , aes(x=Group, y=Counts_Nuclear/1e6))+
  geom_boxplot(fill="steelblue")+
  #geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1))+
  coord_flip()+
  scale_x_discrete(limits=rev(unique(groupedcounts$Group)),labels= rev(unique(groupedcounts$Group)))+
  scale_y_continuous(expand=c(0,0), n.breaks=10, limits = c(0,ymax))+
  scale_y_break(c(2, 5), scales=1)+
  ggtitle("Ribonucleotides in Nucleus")+
  xlab("")+
  ylab("Raw ribonucleotide counts (x10^6)")+
  theme_bw(base_size=25)+
  #theme(plot.title = element_text(size=30, face="bold",hjust = 0.5), 
        #axis.text = element_text(size=12, colour = "black"), 
        #axis.title = element_text(size=25))
	#theme_classic(base_size=20)+
	theme(legend.title=element_blank(),
            legend.background = element_rect(fill="transparent"), #element_rect(fill = "transparent",colour = NA),
            legend.key=element_rect(colour="transparent"), #legend.key = element_rect(fill = "transparent"), 
            legend.key.size = unit(0.5, "cm"),
            axis.line = element_line(size=0.2),
            axis.title=element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank(),
            strip.background = element_blank(),
            plot.margin = unit(c(1, 1, 0.0, 0.0), "cm"), #t, r, b, l
            axis.text = element_text(color="black"),
            axis.text.y = element_text(size=25)
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
            )
dev.off()



png("Raw Mitochondria counts_horizontal.png", width=12, height=7,units= "in",  res=600)
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
        axis.text = element_text(size=16, colour = "black"), 
        axis.title = element_text(size=25))
dev.off()




writeLines("\n... plotting Percentages...\n")+
png("Nuclear and Mitochondria Percentage_legend.png",width=2, height=2,units= "in",  res=600, bg="transparent")

b=ggplot(perc, aes(x=Library, y=counts, fill=Genome))+
  geom_bar(stat="identity", position="fill")+
  #scale_fill_manual(labels = c("Mitochondria", "Nucleus"),values=c("dodgerblue3","firebrick3"))+
  scale_fill_manual(labels = c("Mitochondria", "Nucleus"),values=c("#93A1E5","#E16A86"))+
  scale_x_discrete(limits=QC$Library,labels = str_replace_all(QC$Library, "-", "\n"))+
  scale_y_continuous(expand = c(0,0), labels = scales::percent_format())+
  ggtitle("Nuclear and Mitochondria Distribution")+
  ylab("Percentage")+
  xlab("")+
  theme_classic()+
  theme(plot.title = element_text(size=30, face="bold", hjust = 0.5), 
        axis.text = element_text(size=16, colour = "black"),axis.title = element_text(size=25),
        legend.text = element_blank(), legend.title = element_blank())

legend<-cowplot::get_legend(b)
grid.newpage()
grid.draw(legend)
dev.off()


png("Nuclear and Mitochondria Percentage_horizontal.png", width=12, height=8,units= "in",  res=1000, )
ggplot(perc, aes(x=Library, y=counts, fill=Genome))+
  geom_bar(stat="identity", position="fill", width = 0.9)+
  #geom_text(nudge_y= -.01, color="black",size = 5,fontface="bold")+
  coord_flip()+
  scale_fill_manual(labels = c("Mitochondria", "Nucleus"),values=c("#93A1E5","#E16A86"))+#c("#517AC9","#C05D5D"))+
  scale_x_discrete(limits=rev(QC$Library),labels= rev(QC$Library))+
  scale_y_continuous(expand = c(0,0),labels = scales::percent_format())+
  ggtitle("Nuclear and Mitochondria Distribution")+
  ylab("Percentage")+
  xlab("")+
  theme_classic()+
  theme(plot.title = element_text(size=30, face="bold", hjust = 0.5), 
        axis.text = element_text(size=16, colour = "black"),axis.title = element_text(size=25),
        plot.margin = unit(c(1, 1, 0.0, 0.0), "cm"), #t, r, b, l
        legend.text = element_blank(), legend.title = element_blank()) + 
  theme(legend.position = "none")
dev.off()

png("Nuclear percent_groupedhorizontal.png", width=12, height=7,units= "in",  res=600)
ggplot(groupedcounts , aes(x=Group, y=Nuclperc*100))+
  geom_boxplot(fill="#E16A86")+
  geom_point()+
  coord_flip()+
  scale_x_discrete(limits=rev(unique(groupedcounts$Group)),labels= rev(unique(groupedcounts$Group)))+
  scale_y_continuous(expand=c(0,0), n.breaks=10, limits = c(0,120))+
  #scale_y_break(c(80, 90), scales=2)+
  ggtitle("Ribonucleotide percentage in Nucleus")+
  xlab("")+
  ylab("Percentage")+
  theme_classic()+
  theme(plot.title = element_text(size=30, face="bold",hjust = 0.5), 
        axis.text = element_text(size=12, colour = "black"), 
        axis.title = element_text(size=25))
dev.off()


writeLines("\n... plotting composition...\n")

for (comp in c("Comp_mito", "Comp_nucl", "Freq_mito", "Freq_nucl")) {
  png(paste(comp,".png",sep=""),width=length(unique(libinfo$Group))*2.1, height=6,units= "in",  res=600)
  df=QC[,grep(pattern=comp, x=colnames(QC))]
  df$Group=libinfo$Group
  data=gather(df,key="rNTP", value = "Percentage", -Group)
  obj<-ggbarplot(data, x = "Group", y = "Percentage", fill="rNTP",
                add = c("mean_se"),
                add.params=list(width=0.4,size=0.5),
                palette = c("#D55E00","#0072B2", "#F0E442","#009E73"),
                position = position_dodge(0.75), 
                xlab="", ylab="Percentage", 
                width = 0.75,
                size = 0.5) +
    geom_point(data,mapping=aes(Group,Percentage,color=rNTP), position=position_dodge(0.75),size=1)+guides(color = FALSE)+
    scale_color_manual(values=c("black","black","black","black"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,100), breaks = seq(0,100,10)) +
    theme_classic(base_size = 20)+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 12))
  p=ggpar(obj, legend=c("right"), legend.title = "")
  print(p)
  dev.off()

}

writeLines("\nJob Finished.\n")
writeLines(paste("Output files in current folder"))
