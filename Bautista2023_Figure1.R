#Main Figures of the Manuscript Bautista_2023

##Figure 1###

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

####Define some functions####
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}
#################

#####1.Get data#####
img <- readPNG("Fig1A.png")
genomics <- read_csv("1_growth.csv")
expevol <- read_csv("1_expevol.csv")
#################

############Figure 1A############
gpp <- rasterGrob(img, interpolate=TRUE)
Fig1A<-plot_grid(gpp)
#################

##########Figure 1B############
#Filter by the growth in 4-NQO (UV mimetic)
AV<- dplyr:::filter(genomics, condition=="NQO")
AV<- dplyr:::filter(AV, Evolved=="Evolved_NQO")
AV<- dplyr:::filter(AV, day==2)

my_comparisons <- list( c("1Scer", "2Spar"), c("1Scer","3Hybrid"),
                        c("2Spar", "3Hybrid"))

Fig1B <- AV %>% 
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotypic background", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour)") + 
  scale_color_manual(values=c("green4", "dodgerblue1", "black"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  geom_segment(aes(x=1, xend=2, y=0.67, yend=0.67)) + 
  geom_segment(aes(x=2, xend=3, y=0.70, yend=0.70)) + 
  geom_segment(aes(x=1, xend=3, y=0.74, yend=0.74)) +
  annotate("text", x = 1, y = 0.78, size = 3,
           label = c("p < 0.0001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.69, 0.72, 0.76),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.01", "p < 0.001", "p < 0.0001"),
           family = "", fontface = 3, size=3) +  
  guides(color = guide_legend(title = "Genotype")) +
  theme_prism()+
  theme_bw(base_size=24) +
  theme(legend.position = "none",
        axis.title = element_text(size=11, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        panel.background = element_blank()) +theme(strip.text.x = element_blank())
##########

############Figure 1C############
fdataAREA21 <- dplyr:::filter(expevol,day==21)
fdataAREA21 <- dplyr:::filter(fdataAREA21,NQO=="Yes")

AV<- dplyr:::filter(genomics, condition=="NQO")
AV<- dplyr:::filter(AV, Evolved=="Evolved_NQO")
AV<- dplyr:::filter(AV, day==2)

AV$Replicate <-as.factor(AV$Replicate)
AV$Specie <-as.factor(AV$Specie)
fdataAREA21$rep <-as.factor(fdataAREA21$rep)
fdataAREA21$strain <-as.factor(fdataAREA21$strain)

joined <- full_join(fdataAREA21, AV, by=c("strain"="Specie","rep"="Replicate"))

legend_title <- " "

Fig1C <- joined %>% 
  ggplot(aes(x=rval.x,y = rval.y))+
  border()  +
  geom_point(pch=21, aes(fill=strain), col="black", size=2, alpha=0.5)+
  scale_fill_manual(legend_title,values=c("green4","dodgerblue1","#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus","Hybrid"))) +
  stat_smooth(method="lm", alpha=0.2) +
  ylim(0.2,0.8)+
  xlab("Growth rate (OD/hour) of entire populations") + ylab("Growth rate (OD/hour) of isolated clones") +
  theme_prism() + 
  theme(axis.title = element_text(size=14, face = "bold"),
        legend.position="top",
        legend.direction="horizontal",
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        panel.background = element_blank()) +theme(strip.text.x = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 5)))+
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14))

leg <- get_legend(Fig1C)

Fig1C <- joined %>% 
  ggplot(aes(x=rval.x,y = rval.y))+
  border()  +
  geom_point(pch=21, aes(fill=strain), col="black", size=2, alpha=0.5)+
  scale_fill_manual(legend_title,values=c("green4","dodgerblue1","#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus","Hybrid"))) +
  scale_color_manual(legend_title,values=c("green4","dodgerblue1","#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus","Hybrid"))) +
  annotate("text",
           y = 0.8,x = 0.25,
           label = c("rs = 0.76, p < 0.0001"),
           family = "", fontface = 3, size=3.5) + 
  annotate("text",
           y = 0.8,x = 0.55,
           label = c("rs = 0.37, p < 0.05"),
           family = "", fontface = 3, size=3, col="darkgreen") +
  annotate("text",
           y = 0.76,x = 0.55,
           label = c("rs = 0.43, p < 0.05"),
           family = "", fontface = 3, size=3, col="dodgerblue3") +
  annotate("text",
           y = 0.72,x = 0.55,
           label = c("rs = 0.63, p < 0.001"),
           family = "", fontface = 3, size=3, col="lightpink4")+
  stat_smooth(method="lm", alpha=0.2) +
  ylim(0.2,0.8)+
  xlab("Growth rate (OD/hour) of populations") + ylab("Growth rate (OD/hour) of isolated clones") +
  theme_prism()+
  theme_bw(base_size=24) +
  theme(legend.position = "none",
        axis.title = element_text(size=11, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        panel.background = element_blank()) +theme(strip.text.x = element_blank())
#################

############Assemble and save Figure 1############
Fig1A_label<- plot_grid(Fig1A,labels="a",label_size=20)
Fig1B_label<- plot_grid(Fig1B,labels="b",label_size=20)
Fig1C_label<- plot_grid(Fig1C,labels="c",label_size=20)

Fig1 <- plot_grid(Fig1B_label,Fig1C_label, nrow = 1, rel_widths = c(12,12))
Fig1leg <- plot_grid(Fig1,leg, nrow = 2, rel_heights = c(10,1))
Fig1_img <- plot_grid(Fig1A_label,Fig1leg, nrow = 2,rel_heights = c(1,1.1))

#Save the image in the previously set working directory
ggsave (plot = Fig1_img, filename = "Bautista2023_Figure1_low_quality.jpg", units = "cm", device = "jpg",width =28, height =25, dpi = 300,bg = "white")
ggsave (plot = Fig1_img, filename = "Bautista2023_Figure1.png", units = "cm", device = "png",width =28, height =25, dpi = 1000,bg = "white")
ggsave (plot = Fig1_img, filename = "Bautista2023_Figure1.jpg", units = "cm", device = "jpg",width =28, height =25, dpi = 1000,bg = "white")
ggsave (plot = Fig1_img, filename = "Bautista2023_Figure1.svg", units = "cm", device = "svg",width =28, height =25, dpi = 1000,bg = "white")
ggsave (plot = Fig1_img, filename = "Bautista2023_Figure1.pdf", units = "cm", device = "pdf",width =28, height =25, dpi = 1000,bg = "white")
#################