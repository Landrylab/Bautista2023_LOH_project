#Main Figures of the Manuscript Bautista_2024

##Figure 2###

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
#install.packages("flowCore")
library(flowCore)
#install.packages("flowViz")
library(flowViz)
#install.packages("foreign")
library(foreign)
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
#Define the directory where you save all the data from Bautista_2024
setwd("")
#################

####1.Get data####
genomics <- read_csv("2_growth.csv")
norm_aneu_heatmap <- read.csv("2A_data_aneuploidies_heatmap_fitness_gain.csv")
norm_aneu <- read.csv("2A_data_aneuploidies.csv")
B_Table_LOH_Line13<-read.csv("2B_Table_LOH_Line13.csv")
D_Table_LOH_Percentage<-read.csv("2D_Table_LOH_Percentage.csv")
E_Table_counts_LOH<-read.csv("2E_Table_counts_LOH.csv")
#################

######Figure 2A######
#Heatmap read depth
#Define some elements for the graph
cond.labs <- c(
  `Evolved_control` = "Evolved in control",
  `Evolved_NQO` = "Evolved in UV mimetic",
  `Ancestor` = "Ancestor"
)

make_labels <- function(labels) {
  result <- str_split(labels, "\\.")
  unlist(lapply(result, function(x) x[1]))
}

legend_title <- "Relative Read Depth"

# Generate separate panels for each genotypic backgrounds since Line numbers are ordered by fitness gain in 4-NQO
#HYBRID
g_heatmap<-norm_aneu_heatmap 
ghyb <- g_heatmap%>% dplyr:::filter(Strain=="3Hybrid") 
ghyb <- ghyb%>% dplyr:::filter(Type=="Evolved_NQO")
g <- ghyb %>% ggplot(aes(y = interaction(Replicate,Strain,fitness_gain_evolved_in_NQ0), x = as.factor(chrom4))) + facet_grid(Type~., labeller = as_labeller(cond.labs)) 
g <- g + geom_tile(aes(fill = log2_anc_norm_aneu), color = "white", lwd = .7) +theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ghybC  <- g + scale_fill_gradient2(legend_title,low = "darkslateblue", mid = "gray90", high = "#cc7000", midpoint = 0,limits= c(-2.9,2.9),breaks= c(-2,-1,0,1,2))+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
                            18,19,20,21,22,23,24,25,26,27,28,29,30,31,32),
                   labels=c("Scer chrI","Spar chrI","Scer chrII","Spar chrII",
                            "Scer chrIII","Spar chrIII","Scer chrIV","Spar chrIV",
                            "Scer chrV","Spar chrV","Scer chrVI","Spar chrVI",
                            "Scer chrVII","Spar chrVII","Scer chrVIII","Spar chrVIII",
                            "Scer chrIX","Spar chrIX","Scer chrX","Spar chrX",
                            "Scer chrXI","Spar chrXI","Scer chrXII","Spar chrXII",
                            "Scer chrXIII","Spar chrXIII","Scer chrXIV","Spar chrXIV",
                            "Scer chrXV","Spar chrXV",
                            "Scer chrXVI","Spar chrXVI"))+
  theme_bw()+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.title = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        strip.background = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Chromosome")+
  theme(axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=9))+
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  theme(legend.title = element_text(angle = 0),#, vjust= 0.2,hjust= 0.9),
        legend.title.align=3)+
  theme(legend.key.size = unit(0.4, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.6, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=7))+
  theme(axis.title = element_text(size=15, face = "bold"),
        axis.text = element_text(size=11),
        strip.text = element_text(face = "bold", size = 14),
        strip.background = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "top",
        legend.direction = "horizontal") #+ 
geom_tile(aes(fill =Normalized_read),height=0.65,size=0.2) 

legheatmap <- get_legend(ghybC)

ghybC  <- g + scale_y_discrete(labels = make_labels, name = "Replicate") +
  scale_fill_gradient2(legend_title,low = "darkslateblue", mid = "gray90", high = "#cc7000", midpoint = 0,limits= c(-2.9,2.9),breaks= c(-2,-1,0,1,2))+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
                            18,19,20,21,22,23,24,25,26,27,28,29,30,31,32),
                   labels=c("Scer chrI","Spar chrI","Scer chrII","Spar chrII",
                            "Scer chrIII","Spar chrIII","Scer chrIV","Spar chrIV",
                            "Scer chrV","Spar chrV","Scer chrVI","Spar chrVI",
                            "Scer chrVII","Spar chrVII","Scer chrVIII","Spar chrVIII",
                            "Scer chrIX","Spar chrIX","Scer chrX","Spar chrX",
                            "Scer chrXI","Spar chrXI","Scer chrXII","Spar chrXII",
                            "Scer chrXIII","Spar chrXIII","Scer chrXIV","Spar chrXIV",
                            "Scer chrXV","Spar chrXV",
                            "Scer chrXVI","Spar chrXVI"))+
  theme_bw()+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.title = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        strip.background = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+#+theme(strip.text.y = element_blank())
  xlab("Chromosome")+
  theme(legend.background = element_rect(fill="snow2", 
                                         size=0.4, linetype="solid",colour ="grey"))+
  theme(axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=9))+
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  theme(legend.title = element_text(angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=3)+
  theme (legend.position="none")

ghybC2  <- ghybC + geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = 9.5, ymax = 10.5),color = "firebrick4",fill = NA)+
  geom_rect(aes(xmin = 3.5, xmax = 4.5, ymin = 5.5, ymax = 6.5),color = "firebrick4",fill = NA)+
  geom_rect(aes(xmin = 27.5, xmax = 28.5, ymin = 3.5, ymax = 4.5),color = "blue4",fill = NA)+
  geom_rect(aes(xmin = 23.5, xmax = 24.5, ymin = 3.5, ymax = 4.5),color = "blue4",fill = NA)+
  #geom_rect(aes(xmin = 12.5, xmax = 13.5, ymin = 20.5, ymax = 21.5),color = "firebrick4",fill = NA)+
  #geom_rect(aes(xmin = 13.5, xmax = 14.5, ymin = 20.5, ymax = 21.5),color = "blue4",fill = NA)+
  geom_rect(aes(xmin = 14.5, xmax = 15.5, ymin = 21.5, ymax = 22.5),color = "blue4",fill = NA)


ghyb <- g_heatmap%>%  dplyr:::filter(Strain=="3Hybrid") 
ghyb <- ghyb%>%  dplyr:::filter(Type=="Evolved_control") 

g <- ghyb %>% ggplot(aes(y = interaction(Replicate,Strain), as.factor(chrom4))) + facet_grid(Type~., labeller = as_labeller(cond.labs)) 
g <- g + geom_tile(aes(fill = log2_anc_norm_aneu), color = "white", lwd = .7) +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ggtitle("Hybrid")
ghybA  <- g +
  scale_y_discrete(limits=c("14.3Hybrid","21.3Hybrid","24.3Hybrid","4.3Hybrid","26.3Hybrid",
                            "10.3Hybrid","11.3Hybrid","9.3Hybrid","2.3Hybrid","25.3Hybrid",
                            "5.3Hybrid","22.3Hybrid","30.3Hybrid","27.3Hybrid","13.3Hybrid",
                            "16.3Hybrid", "17.3Hybrid","6.3Hybrid","3.3Hybrid","7.3Hybrid",
                            "19.3Hybrid","1.3Hybrid","23.3Hybrid","29.3Hybrid","28.3Hybrid",
                            "15.3Hybrid","20.3Hybrid","12.3Hybrid","18.3Hybrid","8.3Hybrid"),
                   labels=c(14,21,24,4,26,10,11,9,2,25,5,22,30,27,13,
                            16,17,6,3,7,19,1,23,29,28,15,20,12,18,8)
  )+theme_bw()+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.title = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        strip.background = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+#+theme(strip.text.y = element_blank())
  theme (legend.position="none")+
  xlab("Chromosome")+
  theme(axis.text.y=element_text(size=8),
        axis.text.x=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  scale_fill_gradient2(low = "darkslateblue", mid = "gray90", high = "#cc7000", midpoint = 0,limits= c(-2.9,2.9),breaks= c(-2,-1,0,1,2))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


ghybA2  <- ghybA + geom_rect(aes(xmin = 21.5, xmax = 22.5, ymin = 23.5, ymax = 24.5),color = "firebrick4",fill = NA)


HYBRID <- plot_grid(ghybA2,ghybC2, nrow = 2,rel_heights = c(1,1.15))

#S.PARADOXUS
gspar <- g_heatmap%>%  dplyr:::filter(Strain=="2Spar") 
gspar <- gspar%>%  dplyr:::filter(Type=="Evolved_control") 
g <- gspar %>% ggplot(aes(y = interaction(Replicate,Strain,fitness_gain_evolved_in_NQ0), as.factor(chrom4))) + facet_grid(Type~., labeller = as_labeller(cond.labs)) 
g <- g + geom_tile(aes(fill = log2_anc_norm_aneu), color = "white", lwd = .7) +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ggtitle("Saccharomyces paradoxus")
gsparA <- g + scale_y_discrete(labels = make_labels, name = "Replicate") +
  theme_bw()+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.title = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_gradient2(low = "darkslateblue", mid = "gray90", high = "#cc7000", midpoint = 0,limits= c(-2.9,2.9),breaks= c(-2,-1,0,1,2))+
  theme (legend.position="none") +
  theme(strip.text.y = element_blank(),
        axis.text.y=element_text(size=8))+
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold.italic"))

gspar <- g_heatmap%>%  dplyr:::filter(Strain=="2Spar") 
gspar <- gspar%>%  dplyr:::filter(Type=="Evolved_NQO")  
g <- gspar %>% ggplot(aes(y = interaction(Replicate,Strain,fitness_gain_evolved_in_NQ0), as.factor(chrom4))) + facet_grid(Type~., labeller = as_labeller(cond.labs)) 
g <- g + geom_tile(aes(fill = log2_anc_norm_aneu), color = "white", lwd = .7) +theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
gsparC  <- g + scale_y_discrete(labels = make_labels, name = "Replicate") +
  scale_x_discrete(breaks=c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32),
                   labels=c("Spar chrI","Spar chrII",
                            "Spar chrIII","Spar chrIV",
                            "Spar chrV","Spar chrVI",
                            "Spar chrVII","Spar chrVIII",
                            "Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII",
                            "Spar chrXIII","Spar chrXIV",
                            "Spar chrXV","Spar chrXVI"))+
  theme_bw()+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.title = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme (legend.position="none")+
  xlab("Chromosome")+
  theme(strip.text.y = element_blank(),
        axis.text.y=element_text(size=8))+
  scale_fill_gradient2(low = "darkslateblue", mid = "gray90", high = "#cc7000", midpoint = 0,limits= c(-2.9,2.9),breaks= c(-2,-1,0,1,2))


gsparC2<- gsparC+geom_rect(aes(xmin = 12.5, xmax = 13.5, ymin = 3.5, ymax = 4.5),color = "firebrick4",fill = NA)+
  geom_rect(aes(xmin = 14.5, xmax = 15.5, ymin = 4.5, ymax = 5.5),color = "firebrick4",fill = NA)+
  geom_rect(aes(xmin = 14.5, xmax = 15.5, ymin = 5.5, ymax = 6.5),color = "firebrick4",fill = NA)+
  geom_rect(aes(xmin = 14.5, xmax = 15.5, ymin = 7.5, ymax = 8.5),color = "firebrick4",fill = NA)+
  geom_rect(aes(xmin = 12.5, xmax = 13.5, ymin = 21.5, ymax = 22.5),color = "firebrick4",fill = NA)+
  #geom_rect(aes(xmin = 12.5, xmax = 13.5, ymin = 20.5, ymax = 21.5),color = "firebrick4",fill = NA)+
  #geom_rect(aes(xmin = 13.5, xmax = 14.5, ymin = 20.5, ymax = 21.5),color = "blue4",fill = NA)+
  geom_rect(aes(xmin = 12.5, xmax = 13.5, ymin = 24.5, ymax = 25.5),color = "firebrick4",fill = NA)

SPAR <- plot_grid(gsparA,gsparC2, nrow = 2,rel_heights = c(1,1.15))

#S.CEREVISIAE
gscer <- g_heatmap%>%  dplyr:::filter(Strain=="1Scer") 
gscer <- gscer%>%  dplyr:::filter(Type=="Evolved_control") 

g <- gscer %>% ggplot(aes(y = interaction(Replicate,Strain,fitness_gain_evolved_in_NQ0), as.factor(chrom4))) + facet_grid(Type~., labeller = as_labeller(cond.labs)) 
g <- g + geom_tile(aes(fill = log2_anc_norm_aneu), color = "white", lwd = .7) +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ggtitle("Saccharomyces cerevisiae")
gscerA <- g +scale_y_discrete(labels = make_labels, name = "Line") +theme_bw()+
  theme(axis.text = element_text(face="bold"))+
  theme(strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=15,face = "bold"))+#+theme(strip.text.y = element_blank())
  scale_fill_gradient2(low = "darkslateblue", mid = "gray90", high = "#cc7000", midpoint = 0,limits= c(-2.9,2.9),breaks= c(-2,-1,0,1,2))+
  theme (legend.position="none") +
  theme(strip.text.y = element_blank())+
  ylab("Line")+
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold.italic"),
        axis.text.y=element_text(size=8))

gscer <- g_heatmap%>%  dplyr:::filter(Strain=="1Scer") 
gscer <- gscer%>%  dplyr:::filter(Type=="Evolved_NQO")

g <- gscer %>% ggplot(aes(y = interaction(Replicate,Strain,fitness_gain_evolved_in_NQ0), as.factor(chrom4))) + facet_grid(Type~., labeller = as_labeller(cond.labs)) 
g <- g + geom_tile(aes(fill = log2_anc_norm_aneu), color = "white", lwd = .7) +theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
gscerC  <- g +   scale_y_discrete(labels = make_labels, name = "Line") +
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold.italic"))+
  scale_x_discrete(breaks=c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31),
                   labels=c("Scer chrI","Scer chrII",
                            "Scer chrIII","Scer chrIV",
                            "Scer chrV","Scer chrVI",
                            "Scer chrVII","Scer chrVIII",
                            "Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII",
                            "Scer chrXIII","Scer chrXIV",
                            "Scer chrXV","Scer chrXVI"))+
  theme_bw()+
  theme(axis.text = element_text(face="bold"))+
  theme(strip.text = element_text(face = "bold", size = 16),
        strip.background = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_text(size=15,face = "bold"))+#+theme(strip.text.y = element_blank())
  theme (legend.position="none")+
  xlab("Line")+
  xlab("Chromosome")+
  theme(strip.text.y = element_blank(),
        axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=9),
        axis.title = element_text(size = 15,face="bold"))+
  xlab("Line")+
  scale_fill_gradient2(low = "darkslateblue", mid = "gray90", high = "#cc7000", midpoint = 0,limits= c(-2.9,2.9),breaks= c(-2,-1,0,1,2))


gscerC2<-gscerC+geom_rect(aes(xmin = 14.5, xmax = 15.5, ymin = 0.5, ymax = 1.5),color = "firebrick4",fill = NA)+
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = 4.5, ymax = 5.5),color = "firebrick4",fill = NA)+
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = 25.5, ymax = 26.5),color = "firebrick4",fill = NA)
  
SCER <- plot_grid(gscerA,gscerC2, nrow = 2,rel_heights = c(1,1.15))


#create new legend
data <- data.frame(
  x = c(1, 3),
  y = c("Chromosome gain","Chromosome loss"))

# Create the plot
leg2 <-ggplot(data, aes(x = x, y = y,col=as.factor(y))) +
  geom_point(shape=22,size=6)+
  theme(legend.direction ="horizontal")+
  scale_colour_manual(values=c("firebrick4","blue4"))+
  geom_rect(xmin = 0.5, xmax = 1.5, ymin = 4.5, ymax = 5.5,, color = "firebrick4", fill = "white") +
  geom_rect(xmin = 1.5, xmax = 2.5, ymin = 4.5, ymax = 5.5,, color = "blue4", fill = "white") +
  theme_bw()+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.title = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        strip.background = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Chromosome")+
  theme(axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=9))+
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  theme(legend.title = element_text(angle = 0),#, vjust= 0.2,hjust= 0.9),
        legend.title.align=3)+
  theme(legend.key.size = unit(0.4, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.6, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=7))+
  theme(axis.title = element_text(size=15, face = "bold"),
        axis.text = element_text(size=11),
        strip.text = element_text(face = "bold", size = 14),
        strip.background = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "top",
        legend.direction = "horizontal")+
  theme(legend.title=element_blank())


legheatmap2 <- get_legend(leg2)


# Combine the panels to create a merged plot
Fig2Afitness <- plot_grid(SCER,SPAR,HYBRID, nrow = 1,rel_widths = c(4,4,8))
Fig2A<- Fig2Afitness 

legheatmap3<-plot_grid(legheatmap,legheatmap2)
Fig2Afitness2 <- plot_grid(legheatmap3,Fig2Afitness, nrow=2,rel_heights =  c(0.4,5))
Fig2A<- Fig2Afitness2
#################

######Figure 2B######
###Example t-LOH on line 13
legend_title <- "Relative Read Depth"
Fig2Bb <-B_Table_LOH_Line13%>% ggplot(aes(x = as.numeric(start2),y = as.factor(chrom4),fill=Normalized_read)) + 
  theme_prism()+
  theme_bw(base_size=24) +
  geom_tile(aes(fill =Normalized_read),height=0.65,size=0.2) +
  scale_fill_gradient2(legend_title,low = "darkslateblue", mid = "lavender", high = "orange3", midpoint = 0, limits = c(-2.9, 2.9),breaks= c(-2,-1,0,1,2))+
  facet_wrap(.~ Replicate, ncol=6)+
  scale_y_discrete(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32),
                 labels=c("Scer chrI","Spar chrI",
                          "Scer chrII","Spar chrII",
                          "Scer chrII","Spar chrIII",
                          "Scer chrIV","Spar chrIV",
                          "Scer chrV","Spar chrV",
                          "Scer chrVI","Spar chrVI",
                          "Scer chrVII","Spar chrVII",
                          "Scer chrVIII","Spar chrVIII",
                          "Scer chrIX","Spar chrIX",
                          "Scer chrX","Spar chrX",
                          "Scer chrXI","Spar chrXI",
                          "Scer chrXII","Spar chrXII",
                          "Scer chrXIII","Spar chrXIII",
                          "Scer chrXIV","Spar chrXIV",
                          "Scer chrXV","Spar chrXV",
                          "Scer chrXVI","Spar chrXVI"),
                 limits=rev)+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),
                     labels=c(0,500,1000,1500))+
  xlab("Coordinate (kbp)")+
  ylab("Chromosome")+
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.6, 'cm'),
        axis.title = element_text(size=15, face = "bold"),
        axis.text = element_text(size=11),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(t = -0.2, r = -0.2, b = -0.2, l = -0.2),
        panel.spacing = margin(0, 0, 0, 0),
        strip.text = element_blank(),
        strip.background = element_blank())+
  theme(legend.position="none")

Fig2B <- plot_grid(legheatmap,Fig2Bb, nrow=2,rel_heights =  c(0.29,5))
#################

######Figure 2C######
###t-LOH scheme
img <-image_read_pdf("Fig2C.pdf")
gpp <- rasterGrob(img, interpolate = TRUE)
Fig2C <- plot_grid(gpp)
#################

######Figure 2D######
#Percentage of LOH across genomes
legend_title <- "% of hybrid lines \n bearing LOH"   
Fig2D <- D_Table_LOH_Percentage %>% dplyr:::filter(genome_parental=="Spar") %>% ggplot(aes(x = as.numeric(start2),y = as.factor(chrom4),fill=Percentage)) + 
  geom_tile(aes(fill =Percentage,height=0.65)) +
  scale_fill_gradient2(legend_title,low = "wheat", mid="indianred4",high="red3",midpoint = 50, limits = c(0, 100))+
  scale_y_discrete(breaks=c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32),
                 labels=c("chrI","chrII",
                          "chrIII","chrIV",
                          "chrV","chrVI",
                          "chrVII","chrVIII",
                          "chrIX","chrX",
                          "chrXI","chrXII",
                          "chrXIII","chrXIV",
                          "chrXV","chrXVI"),
                 limits=rev)+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),
                     labels=c(0,500,1000,1500))+
  xlab("Coordinate (kbp)")+
  ylab("Chromosome")+
  theme_bw(base_size=24) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.6, 'cm'),
        axis.title = element_text(size=15, face = "bold"),
        axis.text = element_text(size=11),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        strip.text = element_blank(),
        strip.background = element_blank())+
  theme(legend.spacing.x = unit(0.1, 'cm'),
        legend.text = element_text(margin = margin(t = 5)))


leg_percentage <- get_legend(Fig2D)

Fig2D <- D_Table_LOH_Percentage %>% dplyr:::filter(genome_parental=="Spar") %>% ggplot(aes(x = as.numeric(start2),y = as.factor(chrom4),fill=Percentage)) + 
  geom_tile(aes(fill =Percentage,height=0.65)) +
  scale_fill_gradient2(legend_title,low = "wheat", mid="indianred4",high="red3",midpoint = 50, limits = c(0, 100))+
  scale_y_discrete(breaks=c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32),
                   labels=c("chrI","chrII",
                            "chrIII","chrIV",
                            "chrV","chrVI",
                            "chrVII","chrVIII",
                            "chrIX","chrX",
                            "chrXI","chrXII",
                            "chrXIII","chrXIV",
                            "chrXV","chrXVI"),
                   limits=rev)+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),
                     labels=c(0,500,1000,1500))+
  xlab("Coordinate (kbp)")+
  ylab("Chromosome")+
  theme_bw(base_size=24) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.6, 'cm'),
        axis.title = element_text(size=15, face = "bold"),
        axis.text = element_text(size=11),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        strip.text = element_blank(),
        strip.background = element_blank())+
  theme(legend.position = "none")

Fig2Dleg <- plot_grid(leg_percentage,Fig2D, nrow=2,rel_heights =  c(0.45,5))
#################

######Figure 2E######
#Boxplot LOH control vs UV mimetic conditions

#Exclude triploids
control<-E_Table_counts_LOH %>% filter(Type=="Evolved_control")
control <- control %>% filter (Replicate!=20)
control <- control %>% filter (Replicate!=23)
control <- control %>% filter (Replicate!=26)
control <- control %>% filter (Replicate!=27)
control <- control %>% filter (Replicate!=28)


nqo<-E_Table_counts_LOH %>% filter(Type=="Evolved_NQO")
nqo <- nqo %>% filter (Replicate!=26)
nqo <- nqo %>% filter (Replicate!=27)
nqo <- nqo %>% filter (Replicate!=28)

E_Table_counts_LOH<- rbind(nqo, control)

wilcox.test(nqo$new_counts, control$new_counts)

Fig2E <- ggplot(E_Table_counts_LOH, aes(x = Type, y = new_counts))+
  theme_prism() + 
  geom_violin(colour="black",pch=21,size=1,alpha=0.5, fill="#FF9999",trim=TRUE)+#width=1.15)+
  theme(legend.position = "bottom",
      legend.direction = "horizontal")+
  geom_segment(aes(x=1, xend=2, y=3.2, yend=3.2)) +
  annotate("text",
           y = c(3.4),
           x = c(1.5),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=4) +
  ylab("Number of \n t-LOH / line")+
  #xlab(" ")+
  theme_bw(base_size=24) +
  theme(axis.title = element_text(size=8, face = "bold"),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_blank(),
    axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=8),
    panel.background = element_blank()) +theme(strip.text.x = element_blank())+
  theme(plot.title = element_text(size=15, hjust=0, face = "bold", color="black"),
        axis.text.x = element_text(size = 13, color="black",face = "bold"),
        axis.text.y = element_text(size = 13, color="black"),
        axis.title = element_text(size = 15,face="bold"),
        axis.title.x=element_blank())+
  scale_x_discrete("", labels=c("Evolved in \n control","  Evolved in \n UV mimetic")) 
#################

############Assemble and save Figure 2############
Fig2A2<- plot_grid(Fig2A,labels="a",label_size=20)
Fig2B2<- plot_grid(Fig2B,labels="b",label_size=20)
Fig2C2<- plot_grid(Fig2C,labels="c",label_size=20)
Fig2D2<- plot_grid(Fig2Dleg,labels="d",label_size=20)
Fig2E2<- plot_grid(Fig2E,labels="e",label_size=20)

Fig2DE <- plot_grid(Fig2D2,Fig2E2,nrow=2,rel_heights =c(2.5,1.4))
Fig2BC <- plot_grid(Fig2B2,Fig2C2,nrow=1)
Fig2BBCDE<- plot_grid(Fig2BC,Fig2DE,nrow=1, rel_widths  =c(2,1))
Fig2<- plot_grid(Fig2A2,Fig2BBCDE,nrow=2, rel_heights =c(2,2))


#Save the image in the previously set working directory
ggsave (plot = Fig2, filename = "Bautista2024_Figure2_low_quality.jpg", units = "cm", device = "jpg",width =33, height =35, dpi = 300,bg = "white")
ggsave (plot = Fig2, filename = "Bautista2024_Figure2.png", units = "cm", device = "png",width =33, height =35, dpi = 1000,bg = "white")
ggsave (plot = Fig2, filename = "Bautista2024_Figure2.jpg", units = "cm", device = "jpg",width =33, height =35, dpi = 1000,bg = "white")
ggsave (plot = Fig2, filename = "Bautista2024_Figure2.svg", units = "cm", device = "svg",width =33, height =35, dpi = 1000,bg = "white")
ggsave (plot = Fig2, filename = "Bautista2024_Figure2.pdf", units = "cm", device = "pdf",width =33, height =35, dpi = 1000,bg = "white")
#################
