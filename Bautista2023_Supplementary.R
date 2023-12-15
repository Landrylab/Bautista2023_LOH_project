##Supplementary Figures Bautista_2023###

#In order to place the figures in the same folder as the data tables,
#you need to specify the name of this directory in the 'Set Directory' section below.

####Packages####
rm(list=ls())
#install.packages("agricolae")
library(agricolae)
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
library(BiocManager)
#install.packages("broom")
library(broom)
#install.packages("car",type="binary")
library(car)
#install.packages("cowplot")
library(cowplot)
#install.packages("data.table")
library(data.table)
#install.packages("dplyr")
library(dplyr)
#install.packages("drc")
library(drc)
#install.packages("egg")
library(egg)
#install.packages("factoextra")
library(factoextra)
#install.packages("foreign")
library(foreign)
#install.packages("flowCore")
library(flowCore)
#install.packages("flowViz")
library(flowViz)
#install.packages("gdata")
library(gdata)
#install.packages("ggimage")
library(ggimage)
#install.packages("ggsignif")
library(ggsignif)
#install.packages("ggtext")
library(ggtext)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("ggprism")
library(ggprism)
#install.packages("grid")
library(ggpubr)
#install.packages("ggthemes")
library(ggthemes)
#install.packages("glue")
library(glue)
#install.packages("grid")
library(grid)
#install.packages("gridExtra")
library(gridExtra)
#install.packages("growthcurver") 
library(growthcurver)
#install.packages("gtable")
library(gtable)
#install.packages("janitor")
library(janitor)
#install.packages("lme4",type="binary")
library(lme4)
#install.packages("magick")
library(magick)
#install.packages("magrittr")
library(magrittr)
#install.packages("Matrix") 
library(Matrix)
#install.packages("modelr")
library(modelr)
#install.packages("multcomp")
library(multcomp)
#install.packages("multcompView")
library(multcompView)
#install.packages("nlme") 
library(nlme)
#install.packages("pegas")
library(pegas)
#install.packages("png")
library(png)
#install.packages("quantreg", type="binary")
library(quantreg)
#install.packages("readr")
library(readr)
#install.packages("readxl")
library(readxl)
#install.packages("reshape2")
library(reshape2)
#install.packages("rstatix")
library(rstatix)
#install.packages("stats")
library(stats)
#install.packages("stringr")
library(stringr)
#install.packages("svglite")
library(svglite)
#install.packages("tidyr")
library(tidyr)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("UpSetR")
library(UpSetR)
#BiocManager::install("VariantAnnotation")
library(VariantAnnotation)
#install.packages("vcfR")
library(vcfR)
#################

####Set directory####
#Define the directory where you save all the data from Bautista_2023
setwd("")
#################

####Get Data####
Supp_coverage <- read_csv("1Supp_coverage.csv")
Supp_ploidy <- read.csv("2Supp_ploidy.csv",stringsAsFactors = F)
Supp_ploidy_heatmap <- read.csv("2Supp_ploidy_heatmap.csv")
Supp_genomics <- read_csv("3Supp_growth.csv")
Supp_data_aneuploidies<- read.csv("3Supp_data_aneuploidies.csv")
Supp_Table_counts_aneuploidies <- read.csv("3Supp_Table_counts_aneuploidies.csv")
Supp_Aneuploidy_chrom <- read_excel("3Supp_Aneuploidy_chrom.xlsx")
Supp_new_genome3c_new<-read.csv("3Supp_new_genome3c_new.csv")
Supp4_GO <- image_read_pdf("4Supp_GO.pdf")
Supp_raw_cov<-read.csv("5Supp_raw_cov.csv")
Supp_growth_plasmids<-read.csv("6Supp_growth_plasmids.csv")
Supp_growth_curves_plasmids<-read.csv("6Supp_growth_curves_plasmids.csv")
Supp_growth_plasmids2<-read.csv("6Supp_growth_plasmids2.csv")
Supp_growth_curves_plasmids2<-read.csv("6Supp_growth_curves_plasmids2.csv")
##########

############Figure Supplementary 1############
#Define some functions
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

Supp_coverage$library <- "Sample"

FigA <- Supp_coverage  %>% 
  ggplot(aes(y=Mean_Read_Depth, x=library))+
  theme_prism()+
  geom_boxplot(outlier.size=0.1,col="black",alpha=0.2)+
  stat_summary(fun=mean, geom="point", shape=20, size=3, col="red4",fill="red4") +
  annotate("text",
           y = 125,
           x = 1,
           label = "105.25",
           family = "", fontface = 1, size=3,
           col="red4")+
  ylab("Mean read depth")+
  scale_x_discrete(breaks=c("Sample"),labels=c("All the \n samples")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=11, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank())+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text= element_text(size=9)) 

Supp_coverage$Type <- as.factor(Supp_coverage$Type)
Supp_coverage$Strain <- as.factor(Supp_coverage$Strain)

amod <- aov(Mean_Read_Depth~Type*Strain, data=Supp_coverage)
summary(amod)

legend_title <- " "
FigA1 <- Supp_coverage %>% 
  ggplot(aes(x=Type,y=Mean_Read_Depth, col=Strain, group=interaction(Strain,Type)))+
  geom_boxplot(outlier.size=0.1,col="black",aes(fill=Strain),alpha=0.2)+ 
  theme_prism() + 
  theme_bw(base_size=24) +
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + 
  ylab("Mean read depth")+
  theme(legend.position="bottom",
        axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 25, hjust=1,size=10)) +
  scale_x_discrete(limit =c("Ancestor","Evolved_control", "Evolved_NQO"),
                   label= c("Ancestor","Evolved in control","Evolved in UV mimetic"))+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size=12), 
        legend.text = element_text(size=10))

leg1 <- get_legend(FigA1)

FigB <- Supp_coverage %>% 
  ggplot(aes(x=Type,y=Mean_Read_Depth, col=Strain, group=interaction(Strain,Type)))+
  geom_boxplot(outlier.size=0.1,col="black",aes(fill=Strain),alpha=0.2)+ 
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_prism() + 
  theme(axis.text.x = element_text()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Mean read depth")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(size=13),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 25, hjust=1,size=13)) +
  scale_x_discrete(limit =c("Ancestor","Evolved_control", "Evolved_NQO"),
                   label= c("Ancestor","Evolved\n in control","Evolved in \n UV mimetic"))+
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=10))+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=11, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank())+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text= element_text(size=9)) +
  annotate("text",
           y = c(450),
           x = c(1),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=3) 
#################
############Assemble and save Supplementary Figure 1############
FigA_label <- plot_grid(FigA, labels = c("a"),label_size=20)
FigA2_label <- plot_grid(FigB, labels = c("b"),label_size=20)

Fig1 <- plot_grid(FigA_label,FigA2_label,nrow = 1)
Fig1leg <- plot_grid(leg1,Fig1, nrow = 2, rel_heights = c(2,12))

ggsave (plot = Fig1leg, filename = "Supplementary_Fig1_low_quality.jpg", units = "cm", device = "jpg",width = 20, height =12, dpi = 300,bg = "white")
ggsave (plot = Fig1leg, filename = "Supplementary_Fig1.png", units = "cm", device = "png",width = 20, height =12, dpi = 1000,bg = "white")
ggsave (plot = Fig1leg, filename = "Supplementary_Fig1.jpg", units = "cm", device = "jpg",width = 20, height =12, dpi = 1000,bg = "white")
ggsave (plot = Fig1leg, filename = "Supplementary_Fig1.svg", units = "cm", device = "svg",width = 20, height =12, dpi = 1000,bg = "white")
ggsave (plot = Fig1leg, filename = "Supplementary_Fig1.pdf", units = "cm", device = "pdf",width = 20, height =12, dpi = 1000,bg = "white")
#################

############Figure Supplementary 2############
fdata <- Supp_ploidy

x<-fdata%>%filter(!(plate=="A"&evolved=="4n"))
x<-x%>%filter(!(plate=="C"&evolved=="4n"))
x<-x%>%filter(!(plate=="B"&evolved=="2n"&replicate=="Control"&other_info=="WT ploidy2"))
x<-x%>%filter(!(plate=="B"&evolved=="2n"&replicate=="Control"))
x<-x%>%filter(!(plate=="A"&evolved=="2n"&replicate=="Control"))
x<-x%>%filter(!(plate=="C"&evolved=="2n"&specie=="Spar haploid nat alfa"))
x<-x%>%filter(!(plate=="A"&evolved=="3n"))
x<-x%>%filter(!(plate=="C"&evolved=="3n"))
x<-x%>%filter(!(plate=="A"&evolved=="n"))
x<-x%>%filter(!(plate=="B"&evolved=="n"))
x<-x%>%filter(!(plate=="C"&evolved=="n"&specie=="Spar haploid hyg alfa"))
x <- filter()

#Graphs comparing evolved vs ancestor vs control
hybrid <- x %>% filter(specie=="Hybrid")
scer <- x %>% filter(specie=="Scer")
spar <- x %>% filter(specie=="Spar")

x1<- rbind(hybrid,scer,spar)
x2 <- filter(x1, replicate!="Control")
x2$replicate <- as.numeric(x2$replicate)

#Heatmap choosing the line complete 
AB1 <- filter(x2, specie=="Scer" & replicate=="11" & evolved=="Evolved_control")
AB2 <- filter(x2, specie=="Scer" & replicate=="11" & evolved=="Ancestor")
AB3 <- filter(x2, specie=="Scer" & replicate=="11" & evolved=="Evolved_NQO")
AB4 <- filter(x, evolved =="2n")
AB5 <- filter(x, evolved =="3n")
AB<- rbind(AB1,AB2,AB3,AB4,AB5)
AB$log_FL1A = log1p(AB$FL1_A)

AB<- rbind(AB3,AB5)

#Get counts      
newdat <- AB
newdat$log_FL1A = log1p(newdat$FL1_A)

#First I need to calculate the interval
newdat2 <- newdat %>% 
  group_by(plate,well) %>% 
  mutate(LMH = cut(log(FL1_H), breaks = 100, labels = FALSE))

newdat3 <-newdat2 %>% 
  group_by(plate,well,LMH,evolved,specie,replicate,other_info,same_Replicate) %>% dplyr:::summarise(n()) %>% ungroup()

colnames(newdat3)[9]  <- "counts"

newdat4 <- filter(newdat3, counts > 100)
newdat4$LMH <- as.numeric(newdat4$LMH)

#Do this for all the samples
x$log_FL1A = log1p(x$FL1_A)

newdat <- x %>% 
  group_by(specie,replicate,evolved,well,plate) %>% 
  mutate(LMH = cut(log_FL1A, breaks = 500, labels = FALSE))

newdat <-newdat %>% 
  group_by(specie,replicate,evolved,well,plate,LMH) %>% dplyr:::summarise(n()) %>% ungroup()

colnames(newdat)[7]  <- "counts"

x$same_rep = paste(x$specie, x$replicate, sep="_")
newdat$same_rep = paste(newdat$specie, newdat$replicate, sep="_")

con <- filter(newdat, replicate=="Control")
neworder <- c("n","2n","3n","4n")

iris2 <- arrange(transform(con,
                           evolved=factor(evolved,levels=neworder)),evolved)

legend_title <- "Cell count (Density)"

#Graphs control
con2<- iris2 %>% ggplot(aes(LMH, evolved))+
  geom_tile(aes(fill = counts)) + scale_fill_gradient(legend_title,low = 'white', high = 'darkblue',limits=c(0,800)) +
  ylab("Samples") +
  xlab("Fluorescence")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.line=element_blank())+
  theme(strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank()) +
  ylab("Controls")+
  xlab("DNA content (Fluorescence)")+
  theme(axis.text.y=element_text(size=5),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position="none")+
  theme(legend.title = element_text(angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=3)+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  xlim(300, 500)+
  geom_vline(xintercept = c(392, 420), linetype = "dashed", color = "black",size=0.5)

con3<- iris2 %>% ggplot(aes(LMH, evolved))+
  geom_tile(aes(fill = counts)) + scale_fill_gradient(legend_title,low = 'white', high = 'darkblue',limits=c(0,800)) +
  ylab("Samples") +
  xlab("Fluorescence")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.line=element_blank())+
  theme(axis.text = element_text(face="bold"))+
  theme(strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank()) +
  xlab("DNA content (Fluorescence)")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.position="none")+
  theme(legend.title = element_text(angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=3)+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  xlim(300, 500)+
  geom_vline(xintercept = c(392, 420), linetype = "dashed", color = "black",size=0.5)

newdat$replicate<-as.factor(newdat$replicate)

legend_title <- "Cell count (Density)"

cond.labs <- c(
  `Evolved_control` = "Evolved in control",
  `Evolved_NQO` = "Evolved in UV mimetic",
  `Ancestor` = "Ancestor",
  `1Scer` = "S. cerevisiae",
  `2Spar` = "S. paradoxus",
  `3Hybrid` = "Hybrid"
)

Supp_ploidy_heatmap<-Supp_ploidy_heatmap %>% filter(replicate!="Control")
Supp_ploidy_heatmap$specie <- factor(Supp_ploidy_heatmap$specie,
                                     labels = c("bolditalic(S.cerevisiae)",
                                                "bolditalic(S.paradoxus)",
                                                "bold(Hybrid)"))


Supp_ploidy_heatmap$evolved <- factor(Supp_ploidy_heatmap$evolved,
                                      labels = c("bold(Ancestor)",
                                                 "bold(Control)",
                                                 "bold(UV_mimetic)"))

Samples_heatmap_leg <- Supp_ploidy_heatmap %>% 
  ggplot(aes(LMH, replicate))+
  geom_tile(aes(fill = counts)) + scale_fill_gradient(legend_title,low = 'white', high = 'darkblue') +
  ylab("Samples") +
  xlab("Fluorescence")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.title = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.title = element_text(angle = 0, vjust= 0.3,hjust= 0.01),legend.title.align=0.5,
        legend.box.just = "center",
        legend.direction = "vertical",
        legend.position = "right")+
  theme(legend.title = element_text(angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=3)+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(plot.title = element_text(hjust = 0.8, size=12, face = "bold"))+
  facet_grid(evolved~specie, scales= "free", space = "free",labeller=(label_parsed))+
  xlim(300, 500)

legheatmap <- get_legend(Samples_heatmap_leg)

Samples_heatmap <- Supp_ploidy_heatmap %>% filter(replicate!="Control") %>% ggplot(aes(LMH, replicate))+#,label = interaction(specie,info,glycerol), position = "identity")) + 
  geom_tile(aes(fill = counts)) + scale_fill_gradient(legend_title,low = 'white', high = 'darkblue') +
  ylab("Samples") +
  xlab("Fluorescence")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
  )+
  theme(axis.text = element_text(face="bold"))+
  theme(strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank()) +
  xlab("DNA content (A.U. Fluorescence)")+
  theme(axis.text.y=element_text(size=5))+
  theme(legend.position="none")+
  theme(legend.title = element_text(angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=3)+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  scale_y_discrete(name ="Line", 
                   breaks=c("1","2","3","4","5","6",
                            "7","8","9","10","11",
                            "12","13","14","15","16",
                            "17","18","19","20","21",
                            "22","23","24","25","26",
                            "27","28","29","30"),
                   limits=c("1","2","3","4","5","6",
                            "7","8","9","10","11",
                            "12","13","14","15","16",
                            "17","18","19","20","21",
                            "22","23","24","25","26",
                            "27","28","29","30"),
                   label=c("1","2","3","4","5","6",
                           "7","8","9","10","11",
                           "12","13","14","15","16",
                           "17","18","19","20","21",
                           "22","23","24","25","26",
                           "27","28","29","30"))+
  facet_grid(evolved~specie, scales= "free", space = "free",labeller=label_parsed)+
  xlim(300, 500)+
  geom_vline(xintercept = c(392, 420), linetype = "dashed", color = "black",size=0.5)

Samples_heatmap2 <- Samples_heatmap +  theme(panel.spacing.x = unit(0.52, "cm"))

#################
############Assemble and save Supplementary Figure 2############
leg <- plot_grid(con2,con3,con3, ncol = 4,rel_widths = c(13.6,12.2,12,1.6))
Fig3SuppA <- plot_grid(leg,Samples_heatmap2, nrow = 2, rel_heights = c(2,10))
Fig3SuppAleg <- plot_grid(Fig3SuppA,legheatmap, nrow = 1, rel_widths  = c(10,2))

ggsave (plot = Fig3SuppAleg, filename = "Supplementary_Fig2_low_quality.jpg", units = "cm", device = "jpg",width = 25, height = 25, dpi = 300,bg = "white")
ggsave (plot = Fig3SuppAleg, filename = "Supplementary_Fig2.png", units = "cm", device = "png",width = 25, height =25, dpi = 1000,bg = "white")
ggsave (plot = Fig3SuppAleg, filename = "Supplementary_Fig2.jpg", units = "cm", device = "jpg",width = 25, height = 25, dpi = 1000,bg = "white")
ggsave (plot = Fig3SuppAleg, filename = "Supplementary_Fig2.svg", units = "cm", device = "svg",width = 25, height =25, dpi = 1000,bg = "white")
ggsave (plot = Fig3SuppAleg, filename = "Supplementary_Fig2.pdf", units = "cm", device = "pdf",width = 25, height =25, dpi = 1000,bg = "white")
#################

############Figure Supplementary 3############
genomics1 <- dplyr:::filter(Supp_genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics2<-genomics1 %>% dplyr:::filter(Evolved=="Evolved_NQO")
genomics1 <- dplyr:::filter(Supp_genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics3<-genomics1 %>% dplyr:::filter(Evolved=="Ancestor")

genomics4<- left_join(genomics2,genomics3, by=c("Specie"="Specie","condition"="condition","Replicate"="Replicate"))

genomics4$fitness_gain<- (genomics4$rval.x - genomics4$rval.y) / genomics4$rval.y
genomics4$percentage <- genomics4$fitness_gain * 100

Supp_Table_counts_aneuploidies_NQO <- Supp_Table_counts_aneuploidies %>% 
  dplyr:::filter(Type=="Evolved_NQO")

every_aneu_NQO<- right_join(Supp_Aneuploidy_chrom,Supp_Table_counts_aneuploidies_NQO, by=c("Specie"="Strain",
                                                                                       "Replicate"="Replicate"))

FITNESS_LOH <-left_join(every_aneu_NQO,genomics4,by=c("Specie"="Specie","Replicate"="Replicate"))

#Define some functions
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

NQO <- dplyr:::filter(FITNESS_LOH, Type=="Evolved_NQO")

#Get statistics
NQO <- dplyr:::filter(FITNESS_LOH, Type=="Evolved_NQO")
NQO_scer <- dplyr:::filter(NQO, Specie=="1Scer")
ggscatter(NQO_scer, x = "counts", y = "percentage",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n"))

NQO_spar <- dplyr:::filter(NQO, Specie=="2Spar")
ggscatter(NQO_spar, x = "counts", y = "percentage",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n"))

NQO_hyb <- dplyr:::filter(NQO, Specie=="3Hybrid")
ggscatter(NQO_hyb, x = "counts", y = "percentage",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n"))

legend_title <- " "

FigSupp3 <- FITNESS_LOH %>% filter(Type=="Evolved_NQO")%>%
  ggplot(aes(x=as.factor(counts),y = percentage))+
  stat_smooth(method="lm", alpha=0.1,aes(group=Specie,col=Specie))+
  geom_point(position=position_jitterdodge(seed=2, jitter.width=0.2), 
             pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  xlab("Number of aneuploidies / line")+ylab("Fitness gain (%)")+
  border()  +
  scale_fill_manual(legend_title,values=c("green4","dodgerblue1","#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus","Hybrid"))) +
  scale_color_manual(legend_title,values=c("green4","dodgerblue1","#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus","Hybrid"))) +
  theme(legend.position="none")+
  stat_smooth(method="lm", alpha=0.1)+
  theme_prism()+
  theme_bw(base_size=24) +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(-15, -15,-15,-15))+
  annotate("text",
           y = 190,x = 3,
           label = c("rs = -0.18, p > 0.05"),
           family = "", fontface = 3, size=3, col="darkgreen") +
  annotate("text",
           y = 180,x = 3,
           label = c("rs = -0.22, p > 0.05"),
           family = "", fontface = 3, size=3, col="dodgerblue3") +
  annotate("text",
           y = 170,x = 3,
           label = c("rs = -0.27, p > 0.05"),
           family = "", fontface = 3, size=3, col="lightpink4")+
  stat_smooth(method="lm", alpha=0.2) 
#################
############Assemble and save Supplementary Figure 3############
ggsave (plot = FigSupp3, filename = "Supplementary_Fig3_low_quality.jpg", units = "cm", device = "jpg",width = 15, height =10, dpi = 300,bg = "white")
ggsave (plot = FigSupp3, filename = "Supplementary_Fig3.png", units = "cm", device = "png",width = 15, height =10, dpi = 1000,bg = "white")
ggsave (plot = FigSupp3, filename = "Supplementary_Fig3.jpg", units = "cm", device = "jpg",width = 15, height =10, dpi = 1000,bg = "white")
ggsave (plot = FigSupp3, filename = "Supplementary_Fig3.svg", units = "cm", device = "svg",width = 15, height =10, dpi = 1000,bg = "white")
ggsave (plot = FigSupp3, filename = "Supplementary_Fig3.pdf", units = "cm", device = "pdf",width = 15, height =10, dpi = 1000,bg = "white")
#################

############Figure Supplementary 4############
#Special Packages of enrichment

#install.packages("AnnotationDbi")
library(AnnotationDbi)
#install.packages("AnnotationHub")
library(AnnotationHub)
#install.packages("biomaRt")
library(biomaRt)
#install.packages("biomartr")
library(biomartr) 
#install.packages("clusterProfiler")
library(clusterProfiler)
#install.packages("cola")
library(cola)
#install.packages("DO.db")
library(DO.db)
#install.packages("DOSE")
library(DOSE)
#install.packages("enrichplot")
library(enrichplot)
#install.packages("gage")
library(gage)
#install.packages("gageData")
library(gageData)
#install.packages("ggnewscale")
library(ggnewscale)
#install.packages("gsEasy")
library(gsEasy)
#install.packages("hu6800.db")
library(hu6800.db)
#install.packages("org.Sc.sgd.db")
library(org.Sc.sgd.db)
#install.packages("ReactomePA")
library(ReactomePA)
#install.packages("simplifyEnrichment")
library(simplifyEnrichment) 

m1 <- read.csv("4Supp_GO_scer_NQO_without_outliers.csv")
m1_miss <- filter(m1, Consequence=="missense_variant")
m1 <- m1_miss

m3 <- read.csv("4Supp_GO_spar_NQO_without_outliers.csv")
m3_miss <- filter(m3, Consequence=="missense_variant")
m3 <- m3_miss

m5 <- read.csv("4Supp_GO_hybrid_NQO_without_outliers_hybrid_concatenated.csv")
m5_miss <- filter(m5, Consequence=="missense_variant")
m5 <- m5_miss

m1a<- m1$ensembl_gene_id
m1a<-as.data.frame(m1a)
colnames(m1a) <- c("gene")
m1a<- m1$ensembl_gene_id
hyb_NQO <- enrichGO(m1a,"org.Sc.sgd.db", keyType = "ENSEMBL",
                    ont = "ALL",
                    qvalueCutoff = 0.5,pvalueCutoff=0.5)

x1 <- pairwise_termsim(hyb_NQO) 
set.seed(123)
A1<- emapplot(x1,color = "p.adjust",group_category = T, group_legend = T, cex_label_group = 1.5, node_label = "ONTOLOGY",nCluster = 5)+
  ggtitle("Scer 4-NQO")

m3a<- m3$ensembl_gene_id
m3a<-as.data.frame(m3a)
colnames(m3a) <- c("gene")
m3a<- m3$ensembl_gene_id
hyb_NQO <- enrichGO(m3a,"org.Sc.sgd.db", keyType = "ENSEMBL",
                    ont = "ALL",
                    qvalueCutoff = 0.5,pvalueCutoff=0.5)
x2 <- pairwise_termsim(hyb_NQO) 
set.seed(123)
A2<- emapplot(x2,color = "p.adjust",group_category = T, group_legend = T, cex_label_group = 1.5, node_label = "ONTOLOGY",nCluster = 5)+
  ggtitle("Spar 4-NQO")

m5a<- m5$ensembl_gene_id
m5a<-as.data.frame(m5a)
colnames(m5a) <- c("gene")
m5a<- m5$ensembl_gene_id
hyb_NQO <- enrichGO(m5a,"org.Sc.sgd.db", keyType = "ENSEMBL",
                    ont = "ALL",
                    qvalueCutoff = 0.5,pvalueCutoff=0.5)
x3 <- pairwise_termsim(hyb_NQO) 
set.seed(123)
A3<- emapplot(x3,color = "p.adjust",group_category = T, group_legend = T, cex_label_group = 1.5, node_label = "ONTOLOGY",nCluster = 5)+
  ggtitle("Hybrid 4-NQO")

network<- plot_grid(A1,A2,A3, nrow=1)
#################
############Assemble and save raw Supplementary Figure 4############
ggsave (plot =  network, filename = "Supplementary_Fig4_raw_low_quality.jpg", units = "cm", device = "jpg",width = 60, height =25, dpi = 400,bg = "white")
ggsave (plot =  network, filename = "Supplementary_Fig4_raw.png", units = "cm", device = "png",width = 60, height =25, dpi = 1000,bg = "white")
ggsave (plot =  network, filename = "Supplementary_Fig4_raw.jpg", units = "cm", device = "jpg",width = 60, height =25, dpi = 1000,bg = "white")
ggsave (plot =  network, filename = "Supplementary_Fig4_raw.svg", units = "cm", device = "svg",width = 60, height =25, dpi = 1000,bg = "white")
ggsave (plot =  network, filename = "Supplementary_Fig4_raw.pdf", units = "cm", device = "pdf",width = 60, height =25, dpi = 1000,bg = "white")
#################

#Calculate ratio of enrichment
x1b<-as.data.frame(x1)
x1b$Specie<-"Scer"
x2b<-as.data.frame(x2)
x2b$Specie<-"Spar"
x3b<-as.data.frame(x3)
x3b$Specie<-"Hybrid"

xb<-rbind(x1b,x2b,x3b)

fraction_a_number <- function(fraction) {
  eval(parse(text = fraction))
}
View(xb)

xb$GeneRatio_number <- sapply(xb$GeneRatio, fraction_a_number)
xb$BgRatio_number <- sapply(xb$BgRatio, fraction_a_number)
xb$enrichment_Ratio <-  xb$GeneRatio_number / xb$BgRatio_number

#So the most 2 enriched GO term (higher value of this division of ratios) per specie:
#S.cerevisiae: replication-born double-strand break repair via sister chromatid exchange and regulation of cell differentiation
#S.paradoxus: ER-associated misfolded protein catabolic process, endocytic vesicle
#Hybrid: trehalose metabolic process, ABC-type transporter activity (normally involved in efflux pumps poner)

#################
############Assemble and save Supplementary Figure 4############
#We added some aesthetic modifications so we import it
gpp <- rasterGrob(Supp4_GO, interpolate=TRUE)

FigSupp4<-plot_grid(gpp)

ggsave (plot = FigSupp4, filename = "Supplementary_Fig4_low_quality.jpg", units = "cm", device = "jpg",width = 20, height =40, dpi = 300,bg = "white")
ggsave (plot = FigSupp4, filename = "Supplementary_Fig4.png", units = "cm", device = "png",width = 20, height =40, dpi = 1000,bg = "white")
ggsave (plot = FigSupp4, filename = "Supplementary_Fig4.jpg", units = "cm", device = "jpg",width = 20, height =40, dpi = 300,bg = "white")
ggsave (plot = FigSupp4, filename = "Supplementary_Fig4.svg", units = "cm", device = "svg",width = 20, height =40, dpi = 1000,bg = "white")
ggsave (plot = FigSupp4, filename = "Supplementary_Fig4.pdf", units = "cm", device = "pdf",width = 20, height =40, dpi = 300,bg = "white")
#################

############Figure Supplementary 5############
##Panel A
Supp_new_genome3c_new2 <-  dplyr:::filter(Supp_new_genome3c_new, Strain=="3Hybrid")
Supp_new_genome3c_new2 <-  dplyr:::filter(Supp_new_genome3c_new2, Type=="Evolved_NQO")

exampleHybrid <- dplyr:::filter(Supp_new_genome3c_new2, Strain=="3Hybrid")
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=1)
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=21)
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=10)
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=25)
exampleHybrid<-filter(exampleHybrid,Type=="Evolved_NQO")

nq_snps_HYBRID <-  dplyr:::filter(exampleHybrid, Strain=="3Hybrid")
nq_snps_HYBRID<- nq_snps_HYBRID %>% mutate(Normalized_read=if_else(Normalized_read >1.5, 1.5,Normalized_read))
nq_snps_HYBRID<- nq_snps_HYBRID %>% mutate(Normalized_read=if_else(Normalized_read < -1.5, -1.5,Normalized_read))

nq_snps_HYBRID<- filter(nq_snps_HYBRID,Replicate==28)
nq_genome<- nq_snps_HYBRID

legend_title <- "Relative Read Depth"

FigSupp5A <-nq_genome %>% dplyr:::filter(Type=="Evolved_NQO") %>% 
  ggplot(aes(x = as.numeric(start2),y = as.factor(chrom4),fill=Normalized_read)) + 
  theme_prism() + 
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_tile(aes(fill =Normalized_read),height=0.65,size=0.2) +
  scale_fill_gradient2(legend_title,low = "darkslateblue", mid = "lavender", high = "orange3", midpoint = 0, limits = c(-2.9, 2.9),breaks= c(-2,-1,0,1,2))+
  scale_y_discrete(limits=c("32","31","30","29","28","27","26","25","24","23","22","21","20","19","18","17",
                            "16","15","14","13","12","11","10","9","8","7","6","5","4","3","2","1"),
                   labels=c("Spar XVI","Scer chrXVI","Spar chrXV","Scer chrXV",
                            "Spar chrXIV","Scer chrXIV","Spar chrXIII","Scer chrXIII",
                            "Spar chrXII","Scer chrXII","Spar chrXI","Scer chrXI",
                            "Spar chrX","Scer chrX","Spar chrIX","Scer chrIX",
                            "Spar chrVIII","Scer chrVIII","Spar chrVII","Scer chrVII",
                            "Spar chrVI","Scer chrVI","Spar chrV","Scer chrV",
                            "Spar chrIV","Scer chrIV","Spar chrIII","Scer chrIII",
                            "Spar chrII","Scer chrII","Spar chrI","Scer chrI"))+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),
                     labels=c(0,500,1000,1500))+
  facet_wrap(.~ Replicate, ncol=13)+
  xlab("Coordinate (kbp)")+
  ylab("Chromosome") +
  theme(strip.text = element_text(
    size = 20))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  theme(legend.text = element_text(size = 30))+
  theme(legend.background = element_rect(fill="snow2", 
                                         size=0.4, linetype="solid",colour ="grey"))+
  theme(legend.title = element_text(angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=0.1)+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10))+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text =element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(vjust = 0.5, hjust=1))

#Panel B
aneuploidy<-Supp_raw_cov
aneuploidy$Replicate<-as.numeric(aneuploidy$Replicate)

all_mean_genome_chrom<-aneuploidy %<>% mutate(chrom2 = ifelse(chrom=="utg351_pilon","Spar chrI",
                                                              ifelse(chrom=="utg1271_pilon","Spar chrII",
                                                                     ifelse(chrom=="utg584_pilon","Spar chrIII",
                                                                            ifelse(chrom=="utg988_pilon","Spar chrIV",
                                                                                   ifelse(chrom=="utg69_pilon","Spar chrV",
                                                                                          ifelse(chrom=="utg639_pilon","Spar chrVI",
                                                                                                 ifelse(chrom=="utg199_pilon","Spar chrVII",
                                                                                                        ifelse(chrom=="utg176_pilon","Spar chrVIII",
                                                                                                               ifelse(chrom=="utg1121_pilon","Spar chrIX",
                                                                                                                      ifelse(chrom=="utg245_pilon","Spar chrX",
                                                                                                                             ifelse(chrom=="utg298_pilon","Spar chrXI",
                                                                                                                                    ifelse(chrom=="utg110_pilon","Spar chrXII 1/2",
                                                                                                                                           ifelse(chrom=="utg11_pilon","Spar chrXII 2/2",
                                                                                                                                                  ifelse(chrom=="utg210_pilon","Spar chrXIII",
                                                                                                                                                         ifelse(chrom=="utg48_pilon","Spar chrXIV",
                                                                                                                                                                ifelse(chrom=="utg155_pilon","Spar chrXV 1/2",
                                                                                                                                                                       ifelse(chrom=="utg122_pilon","Spar chrXV 2/2",
                                                                                                                                                                              ifelse(chrom=="utg675_pilon","Spar chrXVI",
                                                                                                                                                                                     ifelse(chrom=="chrI","Scer chrI",
                                                                                                                                                                                            ifelse(chrom=="chrII","Scer chrII",
                                                                                                                                                                                                   ifelse(chrom=="chrIII","Scer chrIII",
                                                                                                                                                                                                          ifelse(chrom=="chrIV","Scer chrIV",
                                                                                                                                                                                                                 ifelse(chrom=="chrV","Scer chrV",
                                                                                                                                                                                                                        ifelse(chrom=="chrVI","Scer chrVI",
                                                                                                                                                                                                                               ifelse(chrom=="chrVII","Scer chrVII",
                                                                                                                                                                                                                                      ifelse(chrom=="chrVIII","Scer chrVIII",
                                                                                                                                                                                                                                             ifelse(chrom=="chrIX","Scer chrIX",
                                                                                                                                                                                                                                                    ifelse(chrom=="chrX","Scer chrX",
                                                                                                                                                                                                                                                           ifelse(chrom=="chrXI","Scer chrXI",
                                                                                                                                                                                                                                                                  ifelse(chrom=="chrXII","Scer chrXII",
                                                                                                                                                                                                                                                                         ifelse(chrom=="chrXIII","Scer chrXIII",
                                                                                                                                                                                                                                                                                ifelse(chrom=="chrXIV","Scer chrXIV",
                                                                                                                                                                                                                                                                                       ifelse(chrom=="chrXV","Scer chrXV",
                                                                                                                                                                                                                                                                                              ifelse(chrom=="chrXVI","Scer chrXVI",NA)))))))))))))))))))))))))))))))))))
all_pools_samples<- all_mean_genome_chrom 
all_pools_samples$log2 <- log2(all_pools_samples$mean_window)

all_pools_samples$gainloss <-ifelse(all_pools_samples$log2>7,"Gain",NA)
all_pools_samples$gainloss <-ifelse(all_pools_samples$log2<2,"Loss",all_pools_samples$gainloss)

all_pools_samples<-all_pools_samples %<>% mutate(chrom2 = ifelse(chrom=="utg351_pilon",17,
                                                                 ifelse(chrom=="utg1271_pilon",18,
                                                                        ifelse(chrom=="utg584_pilon",19,
                                                                               ifelse(chrom=="utg988_pilon",20,
                                                                                      ifelse(chrom=="utg69_pilon",21,
                                                                                             ifelse(chrom=="utg639_pilon",22,
                                                                                                    ifelse(chrom=="utg199_pilon",23,
                                                                                                           ifelse(chrom=="utg176_pilon",24,
                                                                                                                  ifelse(chrom=="utg1121_pilon",25,
                                                                                                                         ifelse(chrom=="utg245_pilon",26,
                                                                                                                                ifelse(chrom=="utg298_pilon",27,
                                                                                                                                       ifelse(chrom=="utg110_pilon",28,
                                                                                                                                              ifelse(chrom=="utg11_pilon",29,
                                                                                                                                                     ifelse(chrom=="utg210_pilon",30,
                                                                                                                                                            ifelse(chrom=="utg48_pilon",31,
                                                                                                                                                                   ifelse(chrom=="utg155_pilon",32,
                                                                                                                                                                          ifelse(chrom=="utg122_pilon",33,
                                                                                                                                                                                 ifelse(chrom=="utg675_pilon",34,
                                                                                                                                                                                        ifelse(chrom=="chrI",1,
                                                                                                                                                                                               ifelse(chrom=="chrII",2,
                                                                                                                                                                                                      ifelse(chrom=="chrIII",3,
                                                                                                                                                                                                             ifelse(chrom=="chrIV",4,
                                                                                                                                                                                                                    ifelse(chrom=="chrV",5,
                                                                                                                                                                                                                           ifelse(chrom=="chrVI",6,
                                                                                                                                                                                                                                  ifelse(chrom=="chrVII",7,
                                                                                                                                                                                                                                         ifelse(chrom=="chrVIII",8,
                                                                                                                                                                                                                                                ifelse(chrom=="chrIX",9,
                                                                                                                                                                                                                                                       ifelse(chrom=="chrX",10,
                                                                                                                                                                                                                                                              ifelse(chrom=="chrXI",11,
                                                                                                                                                                                                                                                                     ifelse(chrom=="chrXII",12,
                                                                                                                                                                                                                                                                            ifelse(chrom=="chrXIII",13,
                                                                                                                                                                                                                                                                                   ifelse(chrom=="chrXIV",14,
                                                                                                                                                                                                                                                                                          ifelse(chrom=="chrXV",15,
                                                                                                                                                                                                                                                                                                 ifelse(chrom=="chrXVI",16,NA)))))))))))))))))))))))))))))))))))
all_pools_samples1 <- filter(all_pools_samples, Unique_value=="Yes")
all_pools_samples1$Replicate <- as.numeric(all_pools_samples1$Replicate) 

NQ_hyb<-all_pools_samples%>% filter(Type=="Evolved_NQO")
NQ_hyb<-NQ_hyb %>% filter(Strain=="3Hybrid")
NQ_hyb28 <- NQ_hyb %>% filter(Replicate == 28)

hybnq <-NQ_hyb28

FigSupp5B <-hybnq %>%
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read \n depth (log2)")+
  ylim(0,10)+#+scale_x_discrete(breaks=c("a","b"),labels=c("2","3"))
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#Panel C
ploidy <- Supp_ploidy
x<-ploidy%>%filter(!(plate=="A"&evolved=="4n"))
x<-x%>%filter(!(plate=="C"&evolved=="4n"))
x<-x%>%filter(!(plate=="B"&evolved=="2n"&replicate=="Control"&other_info=="WT ploidy2"))
x<-x%>%filter(!(plate=="B"&evolved=="2n"&replicate=="Control"))
x<-x%>%filter(!(plate=="A"&evolved=="2n"&replicate=="Control"))
x<-x%>%filter(!(plate=="C"&evolved=="2n"&specie=="Spar haploid nat alfa"))
x<-x%>%filter(!(plate=="A"&evolved=="3n"))
x<-x%>%filter(!(plate=="C"&evolved=="3n"))
x<-x%>%filter(!(plate=="A"&evolved=="n"))
x<-x%>%filter(!(plate=="B"&evolved=="n"))
x<-x%>%filter(!(plate=="C"&evolved=="n"&specie=="Spar haploid hyg alfa"))

fdata <- x
newdat53n <- filter(newdat5, LMH>184)

triploides <- newdat53n %>% group_by(evolved,same_Replicate,specie) %>% dplyr:::summarise(n()) %>% ungroup()

colnames(newdat)[5]  <- "counts"
control3n<-filter(x, evolved=="3n")
control2n<-filter(x, evolved=="2n")
controles <- filter(x, evolved=="2n" | evolved=="3n")

fdata$log_FL1A = log1p(fdata$FL1_A)
control2n$log_FL1A = log1p(control2n$FL1_A)
control3n$log_FL1A = log1p(control3n$FL1_A)

legend_title <- " "

NQ_hyb<-fdata %>% filter(evolved=="Evolved_NQO")
NQ_hyb<-NQ_hyb %>% filter(specie=="3Hybrid")
NQ_hyb28 <- NQ_hyb %>% filter(replicate == 28)

tri28graph <-rbind(NQ_hyb28,control3n,control2n)

FigSupp5C<- tri28graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+
  ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Replicate 28"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(panel.spacing = unit(1, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#################
############Assemble and save Supplementary Figure 5############
FigA_label<- plot_grid(FigSupp5A,labels="a",label_size=25)
FigB_label<- plot_grid(FigSupp5B,labels="b",label_size=25)
FigC_label<- plot_grid(FigSupp5C,labels="c",label_size=25)

Figsupp5_left <- plot_grid(FigB_label,FigC_label,nrow=2)

Figsupp5_all <- plot_grid(FigA_label,Figsupp5_left,nrow=1, rel_widths = c(1,1.7))

ggsave (plot = Figsupp5_all, filename = "Supplementary_Fig5_low_quality.jpg", units = "cm", device = "jpg",width = 40, height =20, dpi = 300,bg = "white")
ggsave (plot = Figsupp5_all, filename = "Supplementary_Fig5.png", units = "cm", device = "png",width = 40, height =20, dpi = 1000,bg = "white")
ggsave (plot = Figsupp5_all, filename = "Supplementary_Fig5.jpg", units = "cm", device = "jpg",width = 40, height =20, dpi = 1000,bg = "white")
ggsave (plot = Figsupp5_all, filename = "Supplementary_Fig5.svg", units = "cm", device = "svg",width = 40, height =20, dpi = 1000,bg = "white")
ggsave (plot = Figsupp5_all, filename = "Supplementary_Fig5.pdf", units = "cm", device = "pdf",width = 40, height =20, dpi = 1000,bg = "white")
#################

############Figure Supplementary 6############
a1<- filter(Supp_growth_plasmids, Type=="WT + pmoby")
a2<- filter (Supp_growth_plasmids, Type=="Δpdr1 + pmoby")

all_data2<-rbind(a1,a2)

all_data2_day1<-all_data2%>% filter(Day==1)

fdata<-Supp_growth_curves_plasmids

fdata1 <- filter(fdata,Type!="WT")
fdata1 <- filter(fdata1,Type!="Δpdr1")

fdataapdr2<- mutate(fdata1, Type2=ifelse(grepl('^Control',Mutation), 'Control', "Mutation"))

fdataapdr2 <- dplyr::filter(fdataapdr2, !grepl('EMPTY', Type))

fdataapdr2 <- dplyr::filter(fdataapdr2,Day==2)

fdataapdr2$hour<-as.numeric(fdataapdr2$hour) 
fdataapdr2$hour_Rounded <- round(fdataapdr2$hour, digits = 1)

data_new3<- fdataapdr2%>%
  dplyr::group_by(Type,Specie,Mutation,hour_Rounded,Condition.y,Type2) %>%
  dplyr::summarise(mean_od = mean(od))

error_estandar <- fdataapdr2 %>%
  dplyr::group_by(Type,Specie,Mutation,hour_Rounded,Condition.y,Type2) %>%
  dplyr::summarise(error_estandar = sd(od) / sqrt(n()))

data_new4<-  full_join(data_new3,error_estandar,by=c("Type"="Type","Specie"="Specie","Mutation"="Mutation","hour_Rounded"="hour_Rounded","Condition.y"="Condition.y","Type2"="Type2"))

data_new4<- filter(data_new4,Mutation!="G282V")

data_new4<- data_new4 %<>% mutate(Mutation = ifelse(Mutation=="Control","WT",Mutation))

data_new4<- data_new4 %<>% mutate(Specie = ifelse(Specie=="BY4741 a","BY4741",Specie))

data_new4<- data_new4 %<>% filter(Specie!="BY4741 alpha")

data_new4<- data_new4 %<>% mutate(Type = ifelse(Type=="WT + pmoby","WT + pmoby_scer",Type))

data_new4<- data_new4 %<>% mutate(Type = ifelse(Type=="Δpdr1 + pmoby","Δpdr1 + pmoby_scer",Type))

new_names <- c("BY4741", "Scer", "Spar")
levels(data_new4$Specie) <- c("BY4741", "italic('S.cerevisiae')", "italic('S.paradoxus')")
levels(data_new4$Specie) <- c("italic('p=')*0.001", "italic('p=')*0.01", "italic('p=')*0.05")

data_bold<- data_new4
data_bold$Specie <- factor(data_bold$Specie,        # Change factor labels
                           labels = c("BY4741",
                                      "italic(S.cerevisiae)",
                                      "italic(S.paradoxus)"))
#Scer growth curves 
Figure2_v2<- data_bold%>% filter(Condition.y=="NQO_4μM") %>% 
  ggplot(aes(y = mean_od,x= hour_Rounded,ymin = mean_od - error_estandar, ymax = mean_od + error_estandar,group=Mutation)) +
  geom_ribbon(alpha = 0.5, aes(fill=Mutation)) +
  geom_line(aes(col=Mutation)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("OD (595 nm)")+xlab("Time (hours)")+
  theme(axis.text.x = element_text()) +
  theme_bw()+
  ylim(0.1,1.1)+
  facet_grid(Type~Specie,labeller=label_parsed)+
  theme(axis.title = element_text(size=20, face = "bold"),
    strip.text = element_text(face = "bold", size = 19),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.background = element_rect(fill = "white") )+
  scale_linetype_manual(values=c("dashed","solid"))+
  theme(legend.position="top",
        legend.direction = "horizontal")+
  scale_color_manual(legend_title,values = c("steelblue1",
                                "gold4","mediumaquamarine",
                                "violet","tan1"),name = " ",
                     breaks = c("G1041/2W","G280R","G280S","M308I","WT"),
                     labels = c("G1042W","G280R","G280S","M308I","WT"))+
  scale_fill_manual(legend_title,values = c("steelblue1",
                               "gold4","mediumaquamarine",
                               "violet","tan1"),name = " ",
                    breaks = c("G1041/2W","G280R","G280S","M308I","WT"),
                    labels = c("G1042W","G280R","G280S","M308I","WT"))+
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14))

all_data<-Supp_growth_plasmids

a1<- filter(all_data, Type=="WT + pmoby")
a2<- filter (all_data, Type=="Δpdr1 + pmoby")

all_data2<-rbind(a1,a2)

all_data2<-all_data2 %>% filter(Mutation!="G282V")

all_data2_day2<-all_data2%>% filter(Day==2)

all_data2_day2<- all_data2_day2 %<>% mutate(Type = ifelse(Type=="WT + pmoby","WT + pmoby_scer",Type))

all_data2_day2<- all_data2_day2 %<>% mutate(Type = ifelse(Type=="Δpdr1 + pmoby","Δpdr1 + pmoby_scer",Type))

check1_BY_day2<- all_data2_day2 %>% filter(Specie=="BY4741") %>% ggplot(aes(x = interaction(Specie,Type), y = Mutation, fill = aucexp)) +
  geom_tile() + 
  facet_grid(.~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "brown3", mid = "burlywood1", high = "green4", midpoint = 3,limits= c(0,36),breaks= c(1,9,18,27,36))+
  theme_bw()+
  ylab("Amino acid \n change")+
  theme(axis.title = element_text(size=14, face = "bold"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank())+
  scale_y_discrete(breaks = c("M308I","G280S","G280R","G1041/2W","Control"),
                   labels = c("M308I","G280S","G280R","G1042W","WT")) 

check1_Scer_day2<- all_data2_day2 %>% filter(Specie=="Scer") %>% ggplot(aes(x = interaction(Specie,Type), y = Mutation, fill = aucexp)) +
  geom_tile() + 
  facet_grid(.~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "brown3", mid = "burlywood1", high = "green4", midpoint = 5,limits= c(0,36),breaks= c(1,9,18,27,36))+
  theme_bw()+
  ylab("Amino acid \n change")+
  theme(axis.title = element_text(size=14, face = "bold"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank())+
  scale_y_discrete(breaks = c("M308I","G280S","G280R","G1041/2W","Control"),
                   labels = c("M308I","G280S","G280R","G1042W","WT")) 

check1_Spar_day2<- all_data2_day2 %>% filter(Specie=="Spar") %>% ggplot(aes(x = interaction(Specie,Type), y = Mutation, fill = aucexp)) +
  geom_tile() + 
  facet_grid(.~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "brown3", mid = "burlywood1", high = "green4", midpoint = 5,limits= c(0,36),breaks= c(1,9,18,27,36))+
  theme_bw()+
  ylab("Amino acid \n change")+
  theme(axis.title = element_text(size=14, face = "bold"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank())+
  scale_y_discrete(breaks = c("M308I","G280S","G280R","G1041/2W","Control"),
                   labels = c("M308I","G280S","G280R","G1042W","WT")) 

#Scer heatmap
combined_plot <- plot_grid(check1_BY_day2,check1_Scer_day2,check1_Spar_day2,ncol=3)

fdata<- Supp_growth_curves_plasmids2

fdata1 <- filter(fdata,Type!="WT")
fdata1 <- filter(fdata1,Type!="Δpdr1")

fdataapdr1 <- dplyr::filter(fdata1, !grepl('Δpdr1', Type))

fdataapdr2<- mutate(fdataapdr1, Type2=ifelse(grepl('^Control',Mutation), 'Control', "Mutation"))

fdataapdr2 <- dplyr::filter(fdataapdr2, !grepl('EMPTY', Type))

all_data<-Supp_growth_plasmids2

a1<- filter(all_data, Type=="WT + pmoby_spar")
a2<- filter (all_data, Type=="Δpdr1 + pmoby_spar")

all_data2<-rbind(a1,a2)

all_data2_day2<-all_data2%>% filter(Day==2)

check1_BY_day2<- all_data2_day2 %>% filter(Specie=="BY4741") %>% ggplot(aes(x = interaction(Specie,Type), y = Mutation, fill = aucexp)) +
  geom_tile() + 
  facet_grid(.~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "brown3", mid = "burlywood1", high = "green4", midpoint = 3,limits= c(0,36),breaks= c(1,9,18,27,36))+
  theme_bw()+
  ylab("Amino acid \n change")+
  theme(axis.title = element_text(size=14, face = "bold"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank())+
  scale_y_discrete(breaks = c("G282V", "G280S", "G280R", "Control"),
                   labels = c("G281V", "G279S", "G279R", "WT"))

all_data2_day2<-all_data2%>% filter(Day==2)

check1_Scer_day2<- all_data2_day2 %>% filter(Specie=="Scer") %>% ggplot(aes(x = interaction(Specie,Type), y = Mutation, fill = aucexp)) +
  geom_tile() + 
  facet_grid(.~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "brown3", mid = "burlywood1", high = "green4", midpoint = 5,limits= c(0,36),breaks= c(1,9,18,27,36))+
  theme_bw()+
  ylab("Amino acid \n change")+
  theme(axis.title = element_text(size=14, face = "bold"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank())+
  scale_y_discrete(breaks = c("G282V", "G280S", "G280R", "Control"),
                     labels = c("G281V", "G279S", "G279R", "WT"))

check1_Spar_day2<- all_data2_day2 %>% filter(Specie=="Spar") %>% ggplot(aes(x = interaction(Specie,Type), y = Mutation, fill = aucexp)) +
  geom_tile() + 
  facet_grid(.~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "brown3", mid = "burlywood1", high = "green4", midpoint = 5,limits= c(0,36),breaks= c(1,9,18,27,36))+
  theme_bw()+
  ylab("Amino acid \n change")+
  theme(axis.title = element_text(size=14, face = "bold"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank())+
  scale_y_discrete(breaks = c("G282V", "G280S", "G280R", "Control"),
                   labels = c("G281V", "G279S", "G279R", "WT"))

#Spar heatmap
combined_plot2 <- plot_grid(check1_BY_day2,check1_Scer_day2,check1_Spar_day2,ncol=3)

fdata<-Supp_growth_curves_plasmids2

fdata1 <- filter(fdata,Type!="WT")
fdata1 <- filter(fdata1,Type!="Δpdr1")

fdataapdr2<- mutate(fdata1, Type2=ifelse(grepl('^Control',Mutation), 'Control', "Mutation"))

fdataapdr2 <- dplyr::filter(fdataapdr2, !grepl('EMPTY', Type))
fdataapdr2 <- dplyr::filter(fdataapdr2, !grepl('CRISPR', Type))

fdataapdr2 <- dplyr::filter(fdataapdr2, !grepl('pmoby_scer', Type))

fdataapdr2 <- dplyr::filter(fdataapdr2,Day==2)

fdataapdr2$hour<-as.numeric(fdataapdr2$hour) 
fdataapdr2$hour_Rounded <- round(fdataapdr2$hour, digits = 1)

data_new3<- fdataapdr2%>%
  dplyr::group_by(Type,Specie,Mutation,hour_Rounded,Condition.y,Type2) %>%
  dplyr::summarise(mean_od = mean(od))

error_estandar <- fdataapdr2 %>%
  dplyr::group_by(Type,Specie,Mutation,hour_Rounded,Condition.y,Type2) %>%
  dplyr::summarise(error_estandar = sd(od) / sqrt(n()))

data_new4<-  full_join(data_new3,error_estandar,by=c("Type"="Type","Specie"="Specie","Mutation"="Mutation","hour_Rounded"="hour_Rounded","Condition.y"="Condition.y","Type2"="Type2"))

data_new4<- data_new4 %<>% mutate(Mutation = ifelse(Mutation=="Control","WT",Mutation))

data_new4 <- filter(data_new4,Specie!="White NQO")
data_new4 <- filter(data_new4,Specie!="White control")
data_bold<- data_new4
data_bold$Specie <- factor(data_bold$Specie,        
                           labels = c("BY4741",
                                      "italic(S.cerevisiae)",
                                      "italic(S.paradoxus)"))

data_bold<- filter(data_bold,hour_Rounded <41)

#Spar growth curves 
Figure2b_v2<- data_bold %>% filter(Condition.y=="NQO_4μM") %>% 
  ggplot(aes(y = mean_od,x= hour_Rounded,
             ymin = mean_od - error_estandar,
             ymax = mean_od + error_estandar,group=Mutation)) +
  geom_ribbon(alpha = 0.5, aes(fill=Mutation)) +
  geom_line(aes(col=Mutation)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("OD (595 nm)")+xlab("Time (hours)")+
  theme(axis.text.x = element_text()) +
  ylim(0.1,1.1)+
  theme_bw()+
  facet_grid(Type~Specie,labeller=label_parsed)+
  theme(axis.title = element_text(size=20, face = "bold"),
    strip.text = element_text(face = "bold", size = 19),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.background = element_rect(fill = "white") )+
  theme(legend.position="top",
        legend.direction = "horizontal")+
  scale_color_manual(legend_title,breaks = c("G280R", "G280S","G282V", "WT"),
                     labels = c("G279R", "G279S","G281V",  "WT"),
                     values = c("gold4","mediumaquamarine","purple","tan1",
                                "gold4","mediumaquamarine","purple","tan1"),name = " ")+
  scale_fill_manual(legend_title,breaks = c("G280R", "G280S","G282V", "WT"),
                    labels = c("G279R", "G279S","G281V",  "WT"),
                    values = c("gold4","mediumaquamarine","purple","tan1",
                               "gold4","mediumaquamarine","purple","tan1"),name = " ")+
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14))
#################
############Assemble and save Supplementary Figure 6############
Fig_supp<- plot_grid(Figure2_v2,Figure2b_v2,nrow=1)
Fig_suppB <- plot_grid(combined_plot,combined_plot2,nrow=2)
Fig_supp_A<- plot_grid(Fig_supp,labels="a",label_size=26)
Fig_suppB_B<- plot_grid(Fig_suppB,labels="b",label_size=26)
FigSupp6 <- plot_grid(Fig_supp_A,Fig_suppB_B,nrow=2)

ggsave (plot = FigSupp6, filename = "Supplementary_Fig6_low_quality.jpg", units = "cm", device = "jpg",width = 40, height =40, dpi = 300,bg = "white")
ggsave (plot = FigSupp6, filename = "Supplementary_Fig6.png", units = "cm", device = "png",width = 40, height =40, dpi = 1000,bg = "white")
ggsave (plot = FigSupp6, filename = "Supplementary_Fig6.jpg", units = "cm", device = "jpg",width = 40, height =40, dpi = 1000,bg = "white")
ggsave (plot = FigSupp6, filename = "Supplementary_Fig6.svg", units = "cm", device = "svg",width = 40, height =40, dpi = 700,bg = "white")
ggsave (plot = FigSupp6, filename = "Supplementary_Fig6.pdf", units = "cm", device = "pdf",width = 40, height =40, dpi = 700,bg = "white")
#################

############Figure Supplementary 7,8 and 9############
############Figure read coverage######## 
all_pools_samples <- Supp_raw_cov
all_pools_samples$log2 <- log2(all_pools_samples$mean_window)

all_pools_samples$gainloss <-ifelse(all_pools_samples$log2>7,"Gain",NA)
all_pools_samples$gainloss <-ifelse(all_pools_samples$log2<2,"Loss",all_pools_samples$gainloss)

all_pools_samples<-all_pools_samples %<>% mutate(chrom2 = ifelse(chrom=="utg351_pilon",17,
                                                                 ifelse(chrom=="utg1271_pilon",18,
                                                                        ifelse(chrom=="utg584_pilon",19,
                                                                               ifelse(chrom=="utg988_pilon",20,
                                                                                      ifelse(chrom=="utg69_pilon",21,
                                                                                             ifelse(chrom=="utg639_pilon",22,
                                                                                                    ifelse(chrom=="utg199_pilon",23,
                                                                                                           ifelse(chrom=="utg176_pilon",24,
                                                                                                                  ifelse(chrom=="utg1121_pilon",25,
                                                                                                                         ifelse(chrom=="utg245_pilon",26,
                                                                                                                                ifelse(chrom=="utg298_pilon",27,
                                                                                                                                       ifelse(chrom=="utg110_pilon",28,
                                                                                                                                              ifelse(chrom=="utg11_pilon",29,
                                                                                                                                                     ifelse(chrom=="utg210_pilon",30,
                                                                                                                                                            ifelse(chrom=="utg48_pilon",31,
                                                                                                                                                                   ifelse(chrom=="utg155_pilon",32,
                                                                                                                                                                          ifelse(chrom=="utg122_pilon",33,
                                                                                                                                                                                 ifelse(chrom=="utg675_pilon",34,
                                                                                                                                                                                        ifelse(chrom=="chrI",1,
                                                                                                                                                                                               ifelse(chrom=="chrII",2,
                                                                                                                                                                                                      ifelse(chrom=="chrIII",3,
                                                                                                                                                                                                             ifelse(chrom=="chrIV",4,
                                                                                                                                                                                                                    ifelse(chrom=="chrV",5,
                                                                                                                                                                                                                           ifelse(chrom=="chrVI",6,
                                                                                                                                                                                                                                  ifelse(chrom=="chrVII",7,
                                                                                                                                                                                                                                         ifelse(chrom=="chrVIII",8,
                                                                                                                                                                                                                                                ifelse(chrom=="chrIX",9,
                                                                                                                                                                                                                                                       ifelse(chrom=="chrX",10,
                                                                                                                                                                                                                                                              ifelse(chrom=="chrXI",11,
                                                                                                                                                                                                                                                                     ifelse(chrom=="chrXII",12,
                                                                                                                                                                                                                                                                            ifelse(chrom=="chrXIII",13,
                                                                                                                                                                                                                                                                                   ifelse(chrom=="chrXIV",14,
                                                                                                                                                                                                                                                                                          ifelse(chrom=="chrXV",15,
                                                                                                                                                                                                                                                                                                 ifelse(chrom=="chrXVI",16,NA)))))))))))))))))))))))))))))))))))










all_pools_samples1 <- filter(all_pools_samples, Unique_value=="Yes")
all_pools_samples1$Replicate <- as.numeric(all_pools_samples1$Replicate) 

#Scer - Ancestor
exampleScer <- filter(all_pools_samples1, Strain=="1Scer")
ancestor <- exampleScer %>% filter(Type=="Ancestor") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold.italic"))+
  ggtitle(expression(atop(paste(bolditalic("Saccharomyces cerevisiae"),bold(" - Ancestor")))))#+theme(strip.text.y = element_blank())

#Spar - Ancestor
exampleSpar <- filter(all_pools_samples1, Strain=="2Spar")
ancestor2 <- exampleSpar %>% filter(Type=="Ancestor") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Saccharomyces paradoxus")+ theme(plot.title = element_text(face = "italic"))+
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold.italic"))+
  ggtitle(expression(atop(paste(bolditalic("Saccharomyces paradoxus"),bold(" - Ancestor")))))#+theme(strip.text.y = element_blank())

#Hybrid - Ancestor
exampleHybrid <- filter(all_pools_samples1, Strain=="3Hybrid")
ancestor3 <- exampleHybrid %>% filter(Type=="Ancestor") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Hybrid - Ancestor")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold"))

M1<- plot_grid(ancestor,ancestor2,ancestor3,nrow=3)

#Scer - Evolved in Control
exampleScer <- filter(all_pools_samples1, Strain=="1Scer")
ancestor <- exampleScer %>% filter(Type=="Evolved_control") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold.italic"))+
  ggtitle(expression(atop(paste(bolditalic("Saccharomyces cerevisiae"),bold(" - Evolved in control")))))#+theme(strip.text.y = element_blank())

#Spar - Evolved in Control
exampleSpar <- filter(all_pools_samples1, Strain=="2Spar")
ancestor2 <- exampleSpar %>% filter(Type=="Evolved_control") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Saccharomyces paradoxus")+ theme(plot.title = element_text(face = "italic"))+
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold.italic"))+
  ggtitle(expression(atop(paste(bolditalic("Saccharomyces paradoxus"),bold(" - Evolved in control")))))#+theme(strip.text.y = element_blank())

#Hybrid - Evolved in Control
exampleHybrid <- filter(all_pools_samples1, Strain=="3Hybrid")
ancestor3 <- exampleHybrid %>% filter(Type=="Evolved_control") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Hybrid - Evolved in control")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold"))

M2<- plot_grid(ancestor,ancestor2,ancestor3,nrow=3)

#Scer - Evolved in UV mimetic
exampleScer <- filter(all_pools_samples1, Strain=="1Scer")
ancestor <- exampleScer %>% filter(Type=="Evolved_NQO") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold.italic"))+
  ggtitle(expression(atop(paste(bolditalic("Saccharomyces cerevisiae"),bold(" - Evolved in UV mimetic")))))#+theme(strip.text.y = element_blank())

#Spar - Evolved in UV mimetic
exampleSpar <- filter(all_pools_samples1, Strain=="2Spar")
ancestor2 <- exampleSpar %>% filter(Type=="Evolved_NQO") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Saccharomyces paradoxus")+ theme(plot.title = element_text(face = "italic"))+
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold.italic"))+
  ggtitle(expression(atop(paste(bolditalic("Saccharomyces paradoxus"),bold(" - Evolved in UV mimetic")))))#+theme(strip.text.y = element_blank())

#Hybrid - Evolved in UV mimetic
exampleHybrid <- filter(all_pools_samples1, Strain=="3Hybrid")
ancestor3 <- exampleHybrid %>% filter(Type=="Evolved_NQO") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Hybrid - Evolved in UV mimetic")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold"))

M3<- plot_grid(ancestor,ancestor2,ancestor3,nrow=3)
#################
############Assemble and save Supplementary Figure 7############
ggsave (plot = M1, filename = "Supplementary_Fig7_low_quality.jpg", units = "cm", device = "jpg",width = 70, height =90, dpi = 300,bg = "white")
ggsave (plot = M1, filename = "Supplementary_Fig7.png", units = "cm", device = "png",width = 70, height =90, dpi = 800,bg = "white")
ggsave (plot = M1, filename = "Supplementary_Fig7.jpg", units = "cm", device = "jpg",width = 70, height =90, dpi = 800,bg = "white")
ggsave (plot = M1, filename = "Supplementary_Fig7.svg", units = "cm", device = "svg",width = 70, height =90, dpi = 500,bg = "white")
ggsave (plot = M1, filename = "Supplementary_Fig7.pdf", units = "cm", device = "pdf",width = 70, height =90, dpi = 500,bg = "white")
#################
############Assemble and save Supplementary Figure 8############
ggsave (plot = M2, filename = "Supplementary_Fig8_low_quality.jpg", units = "cm", device = "jpg",width = 70, height =90, dpi = 300,bg = "white")
ggsave (plot = M2, filename = "Supplementary_Fig8.png", units = "cm", device = "png",width = 70, height =90, dpi = 800,bg = "white")
ggsave (plot = M2, filename = "Supplementary_Fig8.jpg", units = "cm", device = "jpg",width = 70, height =90, dpi = 800,bg = "white")
ggsave (plot = M2, filename = "Supplementary_Fig8.svg", units = "cm", device = "svg",width = 70, height =90, dpi = 500,bg = "white")
ggsave (plot = M2, filename = "Supplementary_Fig8.pdf", units = "cm", device = "pdf",width = 70, height =90, dpi = 500,bg = "white")
#################
############Assemble and save Supplementary Figure 9############
ggsave (plot = M3, filename = "Supplementary_Fig9_low_quality.jpg", units = "cm", device = "jpg",width = 70, height =90, dpi = 300,bg = "white")
ggsave (plot = M3, filename = "Supplementary_Fig9.png", units = "cm", device = "png",width = 70, height =90, dpi = 800,bg = "white")
ggsave (plot = M3, filename = "Supplementary_Fig9.jpg", units = "cm", device = "jpg",width = 70, height =90, dpi = 800,bg = "white")
ggsave (plot = M3, filename = "Supplementary_Fig9.svg", units = "cm", device = "svg",width = 70, height =90, dpi = 500,bg = "white")
ggsave (plot = M3, filename = "Supplementary_Fig9.pdf", units = "cm", device = "pdf",width = 70, height =90, dpi = 500,bg = "white")
#################

############Figure Supplementary 10############
Supp_raw_cov$Replicate<-as.numeric(Supp_raw_cov$Replicate)
all_mean_genome_chrom<-Supp_raw_cov%<>% mutate(chrom2 = ifelse(chrom=="utg351_pilon","Spar chrI",
                                                              ifelse(chrom=="utg1271_pilon","Spar chrII",
                                                                     ifelse(chrom=="utg584_pilon","Spar chrIII",
                                                                            ifelse(chrom=="utg988_pilon","Spar chrIV",
                                                                                   ifelse(chrom=="utg69_pilon","Spar chrV",
                                                                                          ifelse(chrom=="utg639_pilon","Spar chrVI",
                                                                                                 ifelse(chrom=="utg199_pilon","Spar chrVII",
                                                                                                        ifelse(chrom=="utg176_pilon","Spar chrVIII",
                                                                                                               ifelse(chrom=="utg1121_pilon","Spar chrIX",
                                                                                                                      ifelse(chrom=="utg245_pilon","Spar chrX",
                                                                                                                             ifelse(chrom=="utg298_pilon","Spar chrXI",
                                                                                                                                    ifelse(chrom=="utg110_pilon","Spar chrXII 1/2",
                                                                                                                                           ifelse(chrom=="utg11_pilon","Spar chrXII 2/2",
                                                                                                                                                  ifelse(chrom=="utg210_pilon","Spar chrXIII",
                                                                                                                                                         ifelse(chrom=="utg48_pilon","Spar chrXIV",
                                                                                                                                                                ifelse(chrom=="utg155_pilon","Spar chrXV 1/2",
                                                                                                                                                                       ifelse(chrom=="utg122_pilon","Spar chrXV 2/2",
                                                                                                                                                                              ifelse(chrom=="utg675_pilon","Spar chrXVI",
                                                                                                                                                                                     ifelse(chrom=="chrI","Scer chrI",
                                                                                                                                                                                            ifelse(chrom=="chrII","Scer chrII",
                                                                                                                                                                                                   ifelse(chrom=="chrIII","Scer chrIII",
                                                                                                                                                                                                          ifelse(chrom=="chrIV","Scer chrIV",
                                                                                                                                                                                                                 ifelse(chrom=="chrV","Scer chrV",
                                                                                                                                                                                                                        ifelse(chrom=="chrVI","Scer chrVI",
                                                                                                                                                                                                                               ifelse(chrom=="chrVII","Scer chrVII",
                                                                                                                                                                                                                                      ifelse(chrom=="chrVIII","Scer chrVIII",
                                                                                                                                                                                                                                             ifelse(chrom=="chrIX","Scer chrIX",
                                                                                                                                                                                                                                                    ifelse(chrom=="chrX","Scer chrX",
                                                                                                                                                                                                                                                           ifelse(chrom=="chrXI","Scer chrXI",
                                                                                                                                                                                                                                                                  ifelse(chrom=="chrXII","Scer chrXII",
                                                                                                                                                                                                                                                                         ifelse(chrom=="chrXIII","Scer chrXIII",
                                                                                                                                                                                                                                                                                ifelse(chrom=="chrXIV","Scer chrXIV",
                                                                                                                                                                                                                                                                                       ifelse(chrom=="chrXV","Scer chrXV",
                                                                                                                                                                                                                                                                                              ifelse(chrom=="chrXVI","Scer chrXVI",NA)))))))))))))))))))))))))))))))))))


all_pools_samples<- all_mean_genome_chrom 
all_pools_samples$log2 <- log2(all_pools_samples$mean_window)

all_pools_samples<-all_pools_samples %<>% mutate(chrom2 = ifelse(chrom=="utg351_pilon",17,
                                                                 ifelse(chrom=="utg1271_pilon",18,
                                                                        ifelse(chrom=="utg584_pilon",19,
                                                                               ifelse(chrom=="utg988_pilon",20,
                                                                                      ifelse(chrom=="utg69_pilon",21,
                                                                                             ifelse(chrom=="utg639_pilon",22,
                                                                                                    ifelse(chrom=="utg199_pilon",23,
                                                                                                           ifelse(chrom=="utg176_pilon",24,
                                                                                                                  ifelse(chrom=="utg1121_pilon",25,
                                                                                                                         ifelse(chrom=="utg245_pilon",26,
                                                                                                                                ifelse(chrom=="utg298_pilon",27,
                                                                                                                                       ifelse(chrom=="utg110_pilon",28,
                                                                                                                                              ifelse(chrom=="utg11_pilon",29,
                                                                                                                                                     ifelse(chrom=="utg210_pilon",30,
                                                                                                                                                            ifelse(chrom=="utg48_pilon",31,
                                                                                                                                                                   ifelse(chrom=="utg155_pilon",32,
                                                                                                                                                                          ifelse(chrom=="utg122_pilon",33,
                                                                                                                                                                                 ifelse(chrom=="utg675_pilon",34,
                                                                                                                                                                                        ifelse(chrom=="chrI",1,
                                                                                                                                                                                               ifelse(chrom=="chrII",2,
                                                                                                                                                                                                      ifelse(chrom=="chrIII",3,
                                                                                                                                                                                                             ifelse(chrom=="chrIV",4,
                                                                                                                                                                                                                    ifelse(chrom=="chrV",5,
                                                                                                                                                                                                                           ifelse(chrom=="chrVI",6,
                                                                                                                                                                                                                                  ifelse(chrom=="chrVII",7,
                                                                                                                                                                                                                                         ifelse(chrom=="chrVIII",8,
                                                                                                                                                                                                                                                ifelse(chrom=="chrIX",9,
                                                                                                                                                                                                                                                       ifelse(chrom=="chrX",10,
                                                                                                                                                                                                                                                              ifelse(chrom=="chrXI",11,
                                                                                                                                                                                                                                                                     ifelse(chrom=="chrXII",12,
                                                                                                                                                                                                                                                                            ifelse(chrom=="chrXIII",13,
                                                                                                                                                                                                                                                                                   ifelse(chrom=="chrXIV",14,
                                                                                                                                                                                                                                                                                          ifelse(chrom=="chrXV",15,
                                                                                                                                                                                                                                                                                                 ifelse(chrom=="chrXVI",16,NA)))))))))))))))))))))))))))))))))))

#Scer - Ancestor
Anc_scer<-all_pools_samples %>% filter(Type=="Ancestor")
Anc_scer<-Anc_scer %>% filter(Strain=="1Scer")
Anc_scer3 <- Anc_scer %>% filter(Replicate == 3) 
Anc_scer11 <- Anc_scer %>% filter(Replicate == 11) 
Anc_scer19 <- Anc_scer %>% filter(Replicate == 19) 

sceranc <-rbind(Anc_scer3,Anc_scer11,Anc_scer19)

sceranc_read <- sceranc %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#Scer - Evolved in Control
Con_scer<-all_pools_samples %>% filter(Type=="Evolved_control")
Con_scer<-Con_scer %>% filter(Strain=="1Scer")
Con_scer1 <- Con_scer %>% filter(Replicate == 1) 
Con_scer3 <- Con_scer %>% filter(Replicate == 11)

scercon <-rbind(Con_scer1,Con_scer3)

scercon_read <- scercon %>%  
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#Scer - Evolved in UV mimetic
NQ_scer<-all_pools_samples %>% filter(Type=="Evolved_NQO")
NQ_scer<-NQ_scer %>% filter(Strain=="1Scer")
NQ_scer5 <- NQ_scer %>% filter(Replicate == 5) 
NQ_scer11 <- NQ_scer %>% filter(Replicate == 11)
NQ_scer13 <- NQ_scer %>% filter(Replicate == 13) 
NQ_scer21 <- NQ_scer %>% filter(Replicate == 21)

scernq <-rbind(NQ_scer5,NQ_scer11,NQ_scer13,NQ_scer21)

scernq_read <-scernq %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,12)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#Spar - Ancestor
Anc_spar<-all_pools_samples %>% filter(Type=="Ancestor")
Anc_spar<-Anc_spar %>% filter(Strain=="2Spar")
Anc_spar1 <- Anc_spar %>% filter(Replicate == 1) 
Anc_spar2 <- Anc_spar %>% filter(Replicate == 2) 
Anc_spar4 <- Anc_spar %>% filter(Replicate == 4)

sparanc <-rbind(Anc_spar1,Anc_spar2,Anc_spar4)

sparanc_read <- sparanc %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#Spar - Evolved in UV mimetic
NQ_spar<-all_pools_samples %>% filter(Type=="Evolved_NQO")
NQ_spar<-NQ_spar %>% filter(Strain=="2Spar")
NQ_spar23 <- NQ_spar %>% filter(Replicate == 23)
NQ_spar24 <- NQ_spar %>% filter(Replicate == 24)

sparnq <-rbind(NQ_spar23,NQ_spar24)

sparnq_read <-sparnq %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#Hybrid - Ancestor
Anc_hyb<-all_pools_samples %>%filter(Type=="Ancestor")
Anc_hyb<-Anc_hyb %>% filter(Strain=="3Hybrid")
Anc_hyb23 <- Anc_hyb %>% filter(Replicate == 23) 
Anc_hyb24 <- Anc_hyb %>% filter(Replicate == 24) 
Anc_hyb27 <- Anc_hyb %>% filter(Replicate == 27)

hybanc <-rbind(Anc_hyb23,Anc_hyb24,Anc_hyb27)

hybanc_read <- hybanc  %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#Hybrid - Evolved in Control
Con_hyb<-all_pools_samples %>% filter(Type=="Evolved_control")
Con_hyb<-Con_hyb %>% filter(Strain=="3Hybrid")
Con_hyb20 <- Con_hyb %>% filter(Replicate == 20)
Con_hyb23 <- Con_hyb %>% filter(Replicate == 23) 
Con_hyb26 <- Con_hyb %>% filter(Replicate == 26) 
Con_hyb27 <- Con_hyb %>% filter(Replicate == 27)
Con_hyb28 <- Con_hyb %>% filter(Replicate == 28)

hybcon <-rbind(Con_hyb20,Con_hyb23,Con_hyb26,Con_hyb27,Con_hyb28)

hybcon_read <- hybcon %>%
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#Hybrid - Evolved in UV mimetic
NQ_hyb<-all_pools_samples%>% filter(Type=="Evolved_NQO")
NQ_hyb<-NQ_hyb %>% filter(Strain=="3Hybrid")
NQ_hyb26 <- NQ_hyb %>% filter(Replicate == 26) 
NQ_hyb27 <- NQ_hyb %>% filter(Replicate == 27) 
NQ_hyb28 <- NQ_hyb %>% filter(Replicate == 28)

hybnq <-rbind(NQ_hyb26,NQ_hyb27,NQ_hyb28)

hybnq_read <-hybnq %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#Ploidy graphs
ploidy <- Supp_ploidy

x<-ploidy%>%filter(!(plate=="A"&evolved=="4n"))
x<-x%>%filter(!(plate=="C"&evolved=="4n"))
x<-x%>%filter(!(plate=="B"&evolved=="2n"&replicate=="Control"&other_info=="WT ploidy2"))
x<-x%>%filter(!(plate=="B"&evolved=="2n"&replicate=="Control"))
x<-x%>%filter(!(plate=="A"&evolved=="2n"&replicate=="Control"))
x<-x%>%filter(!(plate=="C"&evolved=="2n"&specie=="Spar haploid nat alfa"))
x<-x%>%filter(!(plate=="A"&evolved=="3n"))
x<-x%>%filter(!(plate=="C"&evolved=="3n"))
x<-x%>%filter(!(plate=="A"&evolved=="n"))
x<-x%>%filter(!(plate=="B"&evolved=="n"))
x<-x%>%filter(!(plate=="C"&evolved=="n"&specie=="Spar haploid hyg alfa"))

fdata <- x
newdat53n <- filter(newdat5, LMH>184)

triploides <- newdat53n %>% group_by(evolved,same_Replicate,specie) %>% dplyr:::summarise(n()) %>% ungroup()

colnames(newdat)[5]  <- "counts"
control3n<-filter(x, evolved=="3n")
control2n<-filter(x, evolved=="2n")
controles <- filter(x, evolved=="2n" | evolved=="3n")

fdata$log_FL1A = log1p(fdata$FL1_A)
control2n$log_FL1A = log1p(control2n$FL1_A)
control3n$log_FL1A = log1p(control3n$FL1_A)

legend_title <- "Sample"

Anc_scer<-fdata %>% filter(evolved=="Ancestor")
Anc_scer<-Anc_scer %>% filter(specie=="1Scer")
Anc_scer3 <- Anc_scer %>% filter(replicate == 3) 
Anc_scer11 <- Anc_scer %>% filter(replicate == 11) 
Anc_scer19 <- Anc_scer %>% filter(replicate == 19)

Con_scer<-fdata %>% filter(evolved=="Evolved_control")
Con_scer<-Con_scer %>% filter(specie=="1Scer")
Con_scer1 <- Con_scer %>% filter(replicate == 1) 
Con_scer11 <- Con_scer %>% filter(replicate == 11)

NQ_scer<-fdata %>% filter(evolved=="Evolved_NQO")
NQ_scer<-NQ_scer %>% filter(specie=="1Scer")
NQ_scer5 <- NQ_scer %>% filter(replicate == 5) 
NQ_scer11 <- NQ_scer %>% filter(replicate == 11)
NQ_scer13 <- NQ_scer %>% filter(replicate == 13) 
NQ_scer21 <- NQ_scer %>% filter(replicate == 21)

Anc_spar<-fdata %>% filter(evolved=="Ancestor")
Anc_spar<-Anc_spar %>% filter(specie=="2Spar")
Anc_spar1 <- Anc_spar %>% filter(replicate == 1) 
Anc_spar2 <- Anc_spar %>% filter(replicate == 2) 
Anc_spar4 <- Anc_spar %>% filter(replicate == 4)

NQ_spar<-fdata %>% filter(evolved=="Evolved_NQO")
NQ_spar<-NQ_spar %>% filter(specie=="2Spar")
NQ_spar23 <- NQ_spar %>% filter(replicate == 23)
NQ_spar24 <- NQ_spar %>% filter(replicate == 24)

Anc_hyb<-fdata %>% filter(evolved=="Ancestor")
Anc_hyb<-Anc_hyb %>% filter(specie=="3Hybrid")
Anc_hyb19 <- Anc_hyb %>% filter(replicate == 19)
Anc_hyb23 <- Anc_hyb %>% filter(replicate == 23) 
Anc_hyb24 <- Anc_hyb %>% filter(replicate == 24) 
Anc_hyb27 <- Anc_hyb %>% filter(replicate == 27)

Con_hyb<-fdata %>% filter(evolved=="Evolved_control")
Con_hyb<-Con_hyb %>% filter(specie=="3Hybrid")
Con_hyb20 <- Con_hyb %>% filter(replicate == 20)
Con_hyb23 <- Con_hyb %>% filter(replicate == 23) 
Con_hyb26 <- Con_hyb %>% filter(replicate == 26) 
Con_hyb27 <- Con_hyb %>% filter(replicate == 27)
Con_hyb28 <- Con_hyb %>% filter(replicate == 28)

NQ_hyb<-fdata %>% filter(evolved=="Evolved_NQO")
NQ_hyb<-NQ_hyb %>% filter(specie=="3Hybrid")
NQ_hyb10 <- NQ_hyb %>% filter(replicate == 10) 
NQ_hyb25 <- NQ_hyb %>% filter(replicate == 25) 
NQ_hyb26 <- NQ_hyb %>% filter(replicate == 26) 
NQ_hyb27 <- NQ_hyb %>% filter(replicate == 27) 
NQ_hyb28 <- NQ_hyb %>% filter(replicate == 28)

tri1graph <-rbind(Anc_scer3,control3n,control2n)
tri2graph <-rbind(Anc_scer11,control3n,control2n)
tri3graph <-rbind(Anc_scer19,control3n,control2n)

tri4graph <-rbind(Con_scer1,control3n,control2n)
tri5graph <-rbind(Con_scer11,control3n,control2n)

tri6graph <-rbind(NQ_scer5,control3n,control2n)
tri7graph <-rbind(NQ_scer11,control3n,control2n)
tri8graph <-rbind(NQ_scer13,control3n,control2n)
tri9graph <-rbind(NQ_scer21,control3n,control2n)

tri10graph <-rbind(Anc_spar1,control3n,control2n)
tri11graph <-rbind(Anc_spar2,control3n,control2n)
tri12graph <-rbind(Anc_spar4,control3n,control2n)

tri13graph <-rbind(NQ_spar23,control3n,control2n)
tri14graph <-rbind(NQ_spar24,control3n,control2n)

tri15graph <-rbind(Anc_hyb19,control3n,control2n)
tri16graph <-rbind(Anc_hyb23,control3n,control2n)
tri17graph <-rbind(Anc_hyb24,control3n,control2n)
tri18graph <-rbind(Anc_hyb27,control3n,control2n)

tri19graph <-rbind(Con_hyb20,control3n,control2n)
tri20graph <-rbind(Con_hyb23,control3n,control2n)
tri21graph <-rbind(Con_hyb26,control3n,control2n)
tri22graph <-rbind(Con_hyb27,control3n,control2n)
tri23graph <-rbind(Con_hyb28,control3n,control2n)

tri24graph <-rbind(NQ_hyb10,control3n,control2n)
tri25graph <-rbind(NQ_hyb25,control3n,control2n)
tri26graph <-rbind(NQ_hyb26,control3n,control2n)
tri27graph <-rbind(NQ_hyb27,control3n,control2n)
tri28graph <-rbind(NQ_hyb28,control3n,control2n)

#Scer - Ancestor
tri1<- tri1graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+
  ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 3"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri2<- tri2graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 11"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

tri3<- tri3graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 19"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#Scer - Evolved in Control
tri4<- tri4graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 1"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri5<- tri5graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 11"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Scer - Evolved in UV mimetic
tri6<- tri6graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 5"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri7<- tri7graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 11"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tri8<- tri8graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 13"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri9<- tri9graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 21"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#Spar - Ancestor
tri10<- tri10graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 1"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri11<- tri11graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 2"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri12<- tri12graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 4"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Spar - Evolved in UV mimetic
tri13<- tri13graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 23"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tri14<- tri14graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 24"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#Hybrid - Ancestor
tri16<- tri16graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 23"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri17<- tri17graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 24"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri18<- tri18graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 27"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#Hybrid - Evolved in Control
tri19<- tri19graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 20"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tri20<- tri20graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 23"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri21<- tri21graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 26"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri22<- tri22graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 27"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tri23<- tri23graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 28"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#Hybrid - Evolved in UV mimetic
tri26<- tri26graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 26"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri27<- tri27graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 27"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri28<- tri28graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (A.U. Fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Line 28"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#################
############Assemble and save Supplementary Figure 10############
tri123<- plot_grid(tri1,tri2,tri3, nrow=3)
ggtitle(expression(atop(paste(bolditalic("Saccharomyces cerevisiae"),bold(" - Ancestor")))))#+theme(strip.text.y = element_blank())

title <- ggdraw() + draw_label(expression(atop(paste(bolditalic("Saccharomyces cerevisiae"),bold(" - Ancestor")))),size=15)
scer_ancestor<-plot_grid(sceranc_read,tri123,nrow=1)
scer_ancestor_title <-plot_grid(title,scer_ancestor,nrow=2,rel_heights=c(0.2,2))

tri45<- plot_grid(tri4,tri5, nrow=2)
title <- ggdraw() + draw_label(expression(atop(paste(bolditalic("Saccharomyces cerevisiae"),bold(" - Evolved in control")))),size=15)
scer_con<-plot_grid(scercon_read,tri45,nrow=1)
scer_con_title<-plot_grid(title,scer_con,nrow=2,rel_heights=c(0.2,2))

tri6789<- plot_grid(tri6,tri7,tri8,tri9, nrow=4)
title <- ggdraw() + draw_label(expression(atop(paste(bolditalic("Saccharomyces cerevisiae"),bold(" - Evolved in UV mimetic")))),size=15)
scer_nq<-plot_grid(scernq_read,tri6789,nrow=1)
scer_nq_title<-plot_grid(title,scer_nq,nrow=2,rel_heights=c(0.2, 2))

tri101112<- plot_grid(tri10,tri11,tri12, nrow=3)
title <- ggdraw() + draw_label(expression(atop(paste(bolditalic("Saccharomyces paradoxus"),bold(" - Ancestor")))),size=15)
spar_ancestor<-plot_grid(sparanc_read,tri101112,nrow=1)
spar_ancestor_title<-plot_grid(title,spar_ancestor,nrow=2,rel_heights=c(0.2, 2))

tri1314<- plot_grid(tri13,tri14, nrow=2)
title <- ggdraw() + draw_label(expression(atop(paste(bolditalic("Saccharomyces paradoxus"),bold(" - Evolved in UV mimetic")))),size=15)
spar_nq<-plot_grid(sparnq_read,tri1314,nrow=1)
spar_nq_title<-plot_grid(title,spar_nq,nrow=2,rel_heights=c(0.2, 2))

tri161718<- plot_grid(tri16,tri17,tri18, nrow=3)
title <- ggdraw() + draw_label(expression(paste(bold("Hybrid - Ancestor"))), fontface='bold',size=15)
hyb_ancestor<-plot_grid(hybanc_read,tri161718,nrow=1)
hyb_ancestor_title<-plot_grid(title,hyb_ancestor,nrow=2,rel_heights=c(0.2, 2))

tri1920212223<- plot_grid(tri19,tri20,tri21,tri22,tri23, nrow=5)
title <- ggdraw() + draw_label(expression(paste(bold("Hybrid - Evolved in control"))), fontface='bold',size=15)
hyb_con<-plot_grid(hybcon_read,tri1920212223,nrow=1)
hyb_con_title<-plot_grid(title,hyb_con,nrow=2,rel_heights=c(0.2, 2))

tri262728<- plot_grid(tri26,tri27,tri28, nrow=3)
title <- ggdraw() + draw_label(expression(paste(bold("Hybrid - Evolved in UV mimetic"))), fontface='bold',size=15)
hyb_nq<-plot_grid(hybnq_read,tri262728,nrow=1)
hyb_nq_title<-plot_grid(title,hyb_nq,nrow=2,rel_heights=c(0.2, 2))

A<-plot_grid(scer_ancestor_title,scer_nq_title,ncol=2)
z<- plot_grid(scer_con_title,spar_nq_title,ncol=1)
B<-plot_grid(z,spar_ancestor_title,ncol=2)
C<-plot_grid(hyb_ancestor_title,hyb_nq_title,ncol=1)
CD<-plot_grid(C,hyb_con_title,ncol=2)
BC <-plot_grid(A,B,CD,ncol=1,rel_heights=c(1.1,1.1,1.6))

ggsave (plot = BC, filename = "Supplementary_Fig10_low_quality.jpg", units = "cm", device = "jpg",width =70, height =96, dpi = 300,bg = "white")
ggsave (plot = BC, filename = "Supplementary_Fig10.png", units = "cm", device = "png",width =70, height =96, dpi = 800,bg = "white")
ggsave (plot = BC, filename = "Supplementary_Fig10.jpg", units = "cm", device = "jpg",width =70, height =96, dpi = 800,bg = "white")
ggsave (plot = BC, filename = "Supplementary_Fig10.svg", units = "cm", device = "svg",width =70, height =96, dpi = 500,bg = "white")
ggsave (plot = BC, filename = "Supplementary_Fig10.pdf", units = "cm", device = "pdf",width =70, height =96, dpi = 500,bg = "white")
#################