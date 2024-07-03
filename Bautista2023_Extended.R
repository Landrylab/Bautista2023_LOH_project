##Extended Figures Bautista_2024###

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
#BiocManager::install("flowCore")
library(flowCore)
#BiocManager::install("flowViz")
library(flowViz)
#install.packages("gdata")
library(gdata)
#install.packages("ggimage")
library(ggimage)
#install.packages("ggpattern")
library(ggpattern)
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
#BiocManager::install("ggtree")
library(ggtree)
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
#install.packages("pdftools")
library(pdftools)
#install.packages("pegas")
library(pegas)
#install.packages("plyr")  
library(plyr)
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
#BiocManager::install("trackViewer")
library(trackViewer)
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
#install.packages("dplyr")
library(dplyr)
#################

####Set directory####
#Define the directory where you save all the data from Bautista_2024
setwd("")
#################

####Get Data####
Extended_coverage <- read_csv("1Extended_coverage.csv")
Extended_genomics <- read_csv("1Extended_growth.csv")
Extended_ploidy <- read_csv("2Extended_ploidy.csv")
Extended_fold_change_ploidy <- read_csv("2Extended_fold_change_ploidy.csv")
Extended_LMH_ploidy <- read_csv("2Extended_LMH_ploidy.csv")
Extended_genome <- read_csv("3Extended_new_genome3c_new.csv")
Extended_LOH_window<-read_csv("3Extended_Table_LOH.csv")
Extended_LOH_chrom<-read_csv("4Extended_Table_LOH.csv")
Extended_counts_LOH <- read_csv("4Extended_Table_counts_LOH.csv")
Extended_SNPs<- read_csv("5Extended_SNP.csv")
Extended_GO<-read_csv("5Extended_GO_all_changes.csv")
Extended6_PDR1 <- readPNG("6Extended_PDR1.png")
Extended7_candida <-image_read_pdf("7Extended_candida.pdf")
Extended_growth<-read_csv("8Extended_growth.csv")
Extended_expression<-read_csv("8Extended_expression.csv")
Extended9_boxplot <- read_csv("9Extended_boxplot.csv")
Extended9_curves <- read_csv("9Extended_curves.csv")
Extended_expevol<-read_csv("10Extended_expevol.csv")
Extended10_Sanger <-image_read_pdf("10Extended_Sanger.pdf")
Extended10_LOH_table_parents <- read_excel("10Extended_LOH_parents_table.xlsx")
Extended10_LOH_parents_table_proportion<-read_excel("10Extended_LOH_parents_table_proportion.xlsx")
Extended10_LOH_parents_table_proportion_hybrid<-read_csv("10Extended_LOH_parents_table_proportion_hybrid.csv")
Extended10_LOH_line13_final<-read_excel("10Extended_LOH_line13_final.xlsx")
Extended10_LOH_line28_final<- read_excel("10Extended_LOH_line28_final.xlsx")
Extended10_LOH_line30_final<-read_excel("10Extended_LOH_line30_final.xlsx")
Extended_loh_table_day2<-read_csv("11Extended_loh_table_day2.csv")
Extended_fdata<-read_csv("11Extended_fdata.csv")
##########

############Figure Extended 1############
#NQO
AV<- dplyr::filter(Extended_genomics, condition=="NQO")
nq <- dplyr::filter(AV, Evolved=="Evolved_NQO")
con <- dplyr::filter(AV, Evolved=="Evolved_control")
anc <- dplyr::filter(AV, Evolved=="Ancestor")
AVB<-rbind(nq,con,anc)

#Define some functions
my_comparisons <- list( c("1.1Scer", "2.1Scer"), c("1.3Hybrid","2.3Hybrid"),
                        c("1.2Spar", "2.2Spar"))

cond.labs <- c(
  `Evolved_control` = "Evolved in control",
  `Evolved_NQO` = "Evolved in UV mimetic",
  `Ancestor` = "Ancestor"
)

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

my_comparisons <- list( c("1Scer", "2Spar"), c("1Scer","3Hybrid"),
                        c("2Spar", "3Hybrid"))

legend_title <- " "

Fig_ext1 <- AVB %>% dplyr::filter(day==2)%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in UV mimetic conditions") + 
  scale_color_manual(values=c("green4", "dodgerblue1", "black"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  theme_prism() + 
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'),
        legend.position="top",legend.text.align = 0) +
  guides(color = guide_legend(title = "Genotype")) +
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +ylim(0,0.8)+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=25,angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=0.1,
        legend.text = element_text(size=15),
        legend.direction="horizontal")+ylim(0,1)

legfitness <- get_legend(Fig_ext1)
daytwo <- AVB %>% dplyr::filter(day==2)

AVB1<-AVB %>%dplyr::filter(Evolved=="Evolved_NQO")
AVB1<-AVB1 %>%dplyr::filter(day==2)

AVB1$Specie<-as.factor(AVB1$Specie)
amod <- aov(rval ~ Specie, data = AVB1)
summary(amod)
inter.test1 <- glht(amod,  mcp(Specie = "Tukey"))
summary(inter.test1)
cld(inter.test1)

Fig_ext1_A <- daytwo %>% dplyr::filter(Evolved=="Evolved_NQO")%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in UV mimetic conditions") + 
  scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=0.67, yend=0.67)) + 
  geom_segment(aes(x=2, xend=3, y=0.70, yend=0.70)) + 
  geom_segment(aes(x=1, xend=3, y=0.74, yend=0.74)) +
  annotate("text", x = 1, y = 0.78, size = 3,
           label = c("p < 0.0001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.69, 0.72, 0.76),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.05", "p < 0.001", "p < 0.001"),
           family = "", fontface = 3, size=3)+
  guides(color = guide_legend(title = "Genotype")) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +ylim(0,0.8)+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ylim(0,1)

AVB1<-AVB %>%dplyr::filter(Evolved=="Evolved_control")
AVB1<-AVB1 %>%dplyr::filter(day==2)

AVB1$Specie<-as.factor(AVB1$Specie)
amod <- aov(rval ~ Specie, data = AVB1)
summary(amod)
inter.test1 <- glht(amod,  mcp(Specie = "Tukey"))
summary(inter.test1)
cld(inter.test1)

Fig_ext1_B <- daytwo %>% dplyr::filter(Evolved=="Evolved_control")%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in UV mimetic conditions") + scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                                                                                              labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=0.67, yend=0.67)) + 
  geom_segment(aes(x=2, xend=3, y=0.70, yend=0.70)) + 
  geom_segment(aes(x=1, xend=3, y=0.74, yend=0.74)) +
  annotate("text", x = 1, y = 0.78, size = 3,
           label = c("p < 0.0001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.69, 0.72, 0.76),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.0001", "p > 0.05", "p < 0.05"),
           family = "", fontface = 3, size=3) + 
  guides(color = guide_legend(title = "Genotype")) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +ylim(0,0.8)+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ylim(0,1)

AVB1<-AVB %>%dplyr::filter(Evolved=="Ancestor")
AVB1<-AVB1 %>%dplyr::filter(day==2)

AVB1$Specie<-as.factor(AVB1$Specie)
amod <- aov(rval ~ Specie, data = AVB1)
summary(amod)
inter.test1 <- glht(amod,  mcp(Specie = "Tukey"))
summary(inter.test1)
cld(inter.test1)

Fig_ext1_C <- daytwo %>% dplyr::filter(Evolved=="Ancestor")%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in UV mimetic conditions") + scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                                                                                              labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=0.67, yend=0.67)) + 
  geom_segment(aes(x=2, xend=3, y=0.70, yend=0.70)) + 
  geom_segment(aes(x=1, xend=3, y=0.74, yend=0.74)) +
  annotate("text", x = 1, y = 0.78, size = 3,
           label = c("p < 0.0001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.69, 0.72, 0.76),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.0001", "p < 0.0001", "p < 0.001"),
           family = "", fontface = 3, size=3) + 
  guides(color = guide_legend(title = "Genotype")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +ylim(0,0.8)+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ylim(0,1)

#CONTROL
AV<- dplyr::filter(Extended_genomics, condition=="Control")
nq <- dplyr::filter(AV, Evolved=="Evolved_NQO")
con <- dplyr::filter(AV, Evolved=="Evolved_control")
anc <- dplyr::filter(AV, Evolved=="Ancestor")
AVB<-rbind(nq,con,anc)

#Define some functions
my_comparisons <- list( c("1.1Scer", "2.1Scer"), c("1.3Hybrid","2.3Hybrid"),
                        c("1.2Spar", "2.2Spar"))

cond.labs <- c(
  `Evolved_control` = "Evolved in control",
  `Evolved_NQO` = "Evolved in UV mimetic",
  `Ancestor` = "Ancestor"
)

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

my_comparisons <- list( c("1Scer", "2Spar"), c("1Scer","3Hybrid"),
                        c("2Spar", "3Hybrid"))

legend_title <- " "

Fig_extD <- AVB %>% dplyr::filter(day==2)%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in control conditions") + 
  scale_color_manual(values=c("green4", "dodgerblue1", "black"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  theme_prism() + 
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'),
        legend.position="top",legend.text.align = 0) +
  guides(color = guide_legend(title = "Genotype")) +
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=25,angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=0.1,
        legend.text = element_text(size=15),
        legend.direction="horizontal")+ylim(0,1)

legfitness <- get_legend(Fig_ext1)

AVB1<-AVB %>%dplyr::filter(Evolved=="Evolved_NQO")
AVB1<-AVB1 %>%dplyr::filter(day==2)

AVB1$Specie<-as.factor(AVB1$Specie)
amod <- aov(rval ~ Specie, data = AVB1)
summary(amod)
inter.test1 <- glht(amod,  mcp(Specie = "Tukey"))
summary(inter.test1)
cld(inter.test1)

daytwo <- AVB %>% dplyr::filter(day==2)

Fig_ext1_D <- daytwo %>% dplyr::filter(Evolved=="Evolved_NQO")%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in control conditions") + 
  scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=0.77, yend=0.77)) + 
  geom_segment(aes(x=2, xend=3, y=0.80, yend=0.80)) + 
  geom_segment(aes(x=1, xend=3, y=0.84, yend=0.84)) +
  annotate("text", x = 1, y = 0.88, size = 3,
           label = c("p < 0.001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.79, 0.82, 0.86),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.01", "p < 0.01", "p > 0.05"),
           family = "", fontface = 3, size=3)+
  guides(color = guide_legend(title = "Genotype")) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ylim(0,1)

AVB1<-AVB %>%dplyr::filter(Evolved=="Evolved_control")
AVB1<-AVB1 %>%dplyr::filter(day==2)

AVB1$Specie<-as.factor(AVB1$Specie)
amod <- aov(rval ~ Specie, data = AVB1)
summary(amod)
inter.test1 <- glht(amod,  mcp(Specie = "Tukey"))
summary(inter.test1)
cld(inter.test1)

Fig_ext1_E <- daytwo %>% dplyr::filter(Evolved=="Evolved_control")%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in control conditions") + scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                                                                                              labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=0.87, yend=0.87)) + 
  geom_segment(aes(x=2, xend=3, y=0.90, yend=0.90)) + 
  geom_segment(aes(x=1, xend=3, y=0.94, yend=0.94)) +
  annotate("text", x = 1, y = 0.98, size = 3,
           label = c("p < 0.0001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.89, 0.92, 0.96),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.0001", "p > 0.05", "p < 0.01"),
           family = "", fontface = 3, size=3) + 
  guides(color = guide_legend(title = "Genotype")) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ylim(0,1)

AVB1<-AVB %>%dplyr::filter(Evolved=="Ancestor")
AVB1<-AVB1 %>%dplyr::filter(day==2)

AVB1$Specie<-as.factor(AVB1$Specie)
amod <- aov(rval ~ Specie, data = AVB1)
summary(amod)
inter.test1 <- glht(amod,  mcp(Specie = "Tukey"))
summary(inter.test1)
cld(inter.test1)

Fig_ext1_F <- daytwo %>% dplyr::filter(Evolved=="Ancestor")%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in control conditions") + scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                                                                                              labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=0.87, yend=0.87)) + 
  geom_segment(aes(x=2, xend=3, y=0.90, yend=0.90)) + 
  geom_segment(aes(x=1, xend=3, y=0.94, yend=0.94)) +
  annotate("text", x = 1, y = 0.98, size = 3,
           label = c("p < 0.001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.89, 0.92, 0.96),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.001", "p < 0.01", "p > 0.5"),
           family = "", fontface = 3, size=3) + 
  guides(color = guide_legend(title = "Genotype")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs))+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +ylim(0,1)
#################
############Assemble and save Extended Figure 1############
Fig_ext1<-plot_grid(Fig_ext1_C,Fig_ext1_B,Fig_ext1_A,
                    Fig_ext1_F,Fig_ext1_E,Fig_ext1_D,
                    nrow=2)
Fig_ext1leg <- plot_grid(legfitness, Fig_ext1, nrow=2,rel_heights = c(2,20))

ggsave (plot = Fig_ext1leg, filename = "Extended_Fig1_low_quality.jpg", units = "cm", device = "jpg",width =40, height =25, dpi = 300, bg = "white")
ggsave (plot = Fig_ext1leg, filename = "Extended_Fig1.png", units = "cm", device = "png",width =40, height =25, dpi = 1000, bg = "white")
ggsave (plot = Fig_ext1leg, filename = "Extended_Fig1.jpg", units = "cm", device = "jpg",width =40, height =25, dpi = 1000, bg = "white")
ggsave (plot = Fig_ext1leg, filename = "Extended_Fig1.svg", units = "cm", device = "svg",width =40, height =25, dpi = 1000, bg = "white")
ggsave (plot = Fig_ext1leg, filename = "Extended_Fig1.pdf", units = "cm", device = "pdf",width =40, height =25, dpi = 1000, bg = "white")
#################

############Figure Extended 2############
#Make some functions
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

foldchange3<- Extended_fold_change_ploidy 

#Test ANOVA for having p-value in top
amod1 <- aov(fold_change ~ evolved.x*specie, data = foldchange3)
summary(amod1)

#Test between control and UV mimetic
spar<- dplyr::filter(foldchange3,specie=="2Spar")
t.test(fold_change ~ evolved.x,data=spar)

scer<- dplyr::filter(foldchange3,specie=="1Scer")
t.test(fold_change ~ evolved.x,data=scer)

hyb<- dplyr::filter(foldchange3,specie=="3Hybrid")
t.test(fold_change ~ evolved.x,data=hyb)

#Bonferroni among genotypes
df_p_val1 <- foldchange3 %>%
  rstatix::group_by(specie) %>%
  rstatix::t_test(fold_change ~ evolved.x) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "specie", dodge = 0.8) # important for positioning!

df_p_val2 <- rstatix::t_test(foldchange3, fold_change ~ specie, 
                             p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position()

legend_title <- " "

Fig_Extended_2_A <- foldchange3 %>% ggplot(aes(x = factor(specie), y = fold_change))+ #fill=Strain, show.legend = F,col="black",alpha=0.2)) + 
  geom_boxplot_pattern(aes(pattern = evolved.x,pattern_fill=evolved.x,fill=specie),
                       outlier.shape = 21,
                       outlier.size = 2,
                       color = "black", pattern_fill = "grey100", pattern_angle = 45, pattern_density = 0.4, pattern_spacing = 0.03,alpha=0.25) +
  geom_segment(aes(x=0.8, xend=1.2, y=0.045, yend=0.045)) + 
  geom_segment(aes(x=1.8, xend=2.2, y=0.045, yend=0.045)) + 
  geom_segment(aes(x=2.8, xend=3.2, y=0.045, yend=0.045)) +
  annotate("text",
           y = c(0.049,0.049,0.049,0.088), #,0.064,0.074,0.084,
           x = c(1,2,3,0.9),#1.5,2.5,2,
           label = c("p > 0.05", "p < 0.05", "p > 0.05","p > 0.05"),
           family = "", fontface = 3, size=3) +
  guides(fill = FALSE, alpha = FALSE)+
  scale_pattern_manual(legend_title,values = c(Evolved_NQO = "crosshatch", Evolved_control = "none"),
                       labels = (c("Evolved in control",
                                   "Evolved in UV mimetic")))+
  ylab("Fold change of DNA content \n (A.U. Fluorescence)")+
  xlab("Genotype")+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid")))+
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank())+
  theme(legend.title = element_text(face = "italic"))+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12))+
  theme(plot.title = element_text(face = "italic", hjust = 0.5)) +
  theme(legend.position = "bottom",
        axis.title.x=element_blank())

#Fitness correlations
genomics1 <- dplyr:::filter(Extended_genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics2<-genomics1 %>% dplyr:::filter(Evolved=="Evolved_NQO")
genomics1 <- dplyr:::filter(Extended_genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics3<-genomics1 %>% dplyr:::filter(Evolved=="Ancestor")

genomics4<- left_join(genomics2,genomics3, by=c("Specie"="Specie","condition"="condition","Replicate"="Replicate"))

genomics4$fitness_gain<- (genomics4$rval.x - genomics4$rval.y) / genomics4$rval.y
genomics4$percentage <- genomics4$fitness_gain * 100

Extended_LMH_ploidy$replicate <- as.numeric(Extended_LMH_ploidy$replicate)
ploidy_fitness<- full_join(Extended_LMH_ploidy,genomics4,by=c("replicate"="Replicate","specie"="Specie","evolved"="Evolved.x"))

sc<-dplyr::filter(ploidy_fitness,specie=="1Scer")
sp<-dplyr::filter(ploidy_fitness,specie=="2Spar")
hy <-dplyr::filter(ploidy_fitness,specie=="3Hybrid")

ploidyfitnes2<-rbind(sc,sp,hy)

legend_title <- ""
Fig_Extended_2_B_left <-ploidyfitnes2 %>% dplyr::filter(evolved=="Evolved_NQO")%>%
  ggscatter(x = "LMH_median", y = "percentage",color ="specie",
            add = "reg.line",
            conf.int = TRUE,  # Add regressin line,
            add.params = list(color = "grey"),size=1.5,alpha=0.6)+
  scale_color_manual(legend_title,values=c("green4","dodgerblue1","#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus","Hybrid"))) +
  theme_bw()+
  annotate("text",
           y = 0,x = 172,
           label = c("rs = 0.029, p > 0.05"),
           family = "", fontface = 3, size=3.5) +
  ylim(0,200)+
  xlab("DNA content (A.U. Fluorescence)") + ylab("Increase in \n growth rate (%)") +
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        plot.title = element_text(face = "italic", hjust = 0.5),
        legend.position = "bottom")+
  xlim(170,185)

leg <- get_legend(Fig_Extended_2_B_left)

Fig_Extended_2_B_left <-ploidyfitnes2 %>% dplyr::filter(evolved=="Evolved_NQO")%>%
  ggscatter(x = "LMH_median", y = "percentage",color ="specie",
            add = "reg.line",
            conf.int = TRUE,  # Add regressin line,
            add.params = list(color = "grey"),size=1.5,alpha=0.6)+
  scale_color_manual(legend_title,values=c("green4","dodgerblue1","#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus","Hybrid"))) +theme_bw()+
  annotate("text",
           y = 5,x = 174,
           label = c("rs = 0.029, p > 0.05"),
           family = "", fontface = 3, size=3.5) +
  ylim(0,200)+
  xlab("DNA content (A.U. Fluorescence)") + ylab("Increase in \n growth rate (%)") +
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        plot.title = element_text(face = "italic", hjust = 0.5),
        legend.position = "none")+
  xlim(172,184)

ploidyfitnes2_nq <- ploidyfitnes2 %>% dplyr::filter(evolved=="Evolved_NQO")

cond.labs <- c(
  `Evolved_control` = "Evolved in control",
  `Evolved_NQO` = "Evolved in UV mimetic",
  `Ancestor` = "Ancestor",
  `1Scer` = "S. cerevisiae",
  `2Spar` = "S. paradoxus",
  `3Hybrid` = "Hybrid"
)

FigSuppB_scer <-ploidyfitnes2_nq %>% dplyr::filter(specie=="1Scer")%>%
  ggscatter(x = "LMH_median", y = "percentage",color ="specie",
            add = "reg.line",
            conf.int = TRUE,  # Add regressin line,
            add.params = list(color = "grey"),size=1.5,alpha=0.6)+
  scale_color_manual(legend_title,values=c("green4","dodgerblue1","#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus","Hybrid"))) +
  annotate("text",
           y = 5,x = 174,
           label = c("rs = -0.27, p > 0.05"),
           family = "", fontface = 3, size=3.5) +
  ylim(0,200)+
  xlab("DNA content (A.U. Fluorescence)") + ylab("Increase in \n growth rate (%)") +
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        plot.title = element_text(face = "italic", hjust = 0.5),
        legend.position = "none")+
  xlim(172,184)

FigSuppB_spar <-ploidyfitnes2_nq %>% dplyr::filter(specie=="2Spar")%>%
  ggscatter(x = "LMH_median",y = "percentage",color ="specie",
            add = "reg.line",
            conf.int = TRUE,  # Add regressin line,
            add.params = list(color = "grey"),size=1.5,alpha=0.6)+
  scale_color_manual(legend_title,values=c("dodgerblue1","#FF9999"), 
                     labels = toexpr(c("S. paradoxus","Hybrid"))) +
  annotate("text",
           y = 5,x = 174,
           label = c("rs = -0.05, p > 0.05"),
           family = "", fontface = 3, size=3.5) +
  ylim(0,200)+
  xlab("DNA content (A.U. Fluorescence)") + ylab("Increase in \n growth rate (%)") +
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        plot.title = element_text(face = "italic", hjust = 0.5),
        legend.position = "none")+
  xlim(172,184)

FigSuppB_hyb <-ploidyfitnes2_nq %>% dplyr::filter(specie=="3Hybrid")%>%
  ggscatter(x = "LMH_median", y = "percentage",color ="specie",
            add = "reg.line",
            conf.int = TRUE,  # Add regressin line,
            add.params = list(color = "grey"),size=1.5,alpha=0.6)+
  scale_color_manual(legend_title,values=c("#FF9999"), 
                     labels = toexpr(c("Hybrid"))) +
  annotate("text",
           y = 5,x = 174,
           label = c("rs = -0.084, p > 0.05"),
           family = "", fontface = 3, size=3.5) +
  ylim(0,200)+
  xlab("DNA content (A.U. Fluorescence)") + ylab("Increase in \n growth rate (%)") +
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        plot.title = element_text(face = "italic", hjust = 0.5),
        legend.position = "none")+
  xlim(172,184)
#################
############Assemble and save Extended Figure 2############
Fig_Extended_2_B_right <- plot_grid(FigSuppB_scer,FigSuppB_spar,FigSuppB_hyb,ncol=1)
Fig_Extended_Ploidyfitness <- plot_grid(Fig_Extended_2_B_left,Fig_Extended_2_B_right,ncol=2,rel_widths = c(2,2))
Fig_Extended_Ploidyfitness_leg <- plot_grid(leg,Fig_Extended_Ploidyfitness,ncol=1,rel_heights = c(2,20))

Fig_Extended_2A_label<- plot_grid(Fig_Extended_2_A,labels="a",label_size=15)
Fig_Extended_2B_label<- plot_grid(Fig_Extended_Ploidyfitness_leg,labels="b",label_size=15)

Fig_Extended2_ploidy<- plot_grid(Fig_Extended_2A_label,Fig_Extended_2B_label,nrow=2,rel_heights = c(1,1.5))

ggsave (plot = Fig_Extended2_ploidy, filename = "Extended_Fig2_low_quality.jpg", units = "cm", device = "jpg",width =22, height =22, dpi = 300, bg = "white")
ggsave (plot = Fig_Extended2_ploidy, filename = "Extended_Fig2.png", units = "cm", device = "png",width =22, height =22, dpi = 1000, bg = "white")
ggsave (plot = Fig_Extended2_ploidy, filename = "Extended_Fig2.jpg", units = "cm", device = "jpg",width =22, height =22, dpi = 1000, bg = "white")
ggsave (plot = Fig_Extended2_ploidy, filename = "Extended_Fig2.svg", units = "cm", device = "svg",width =22, height =22, dpi = 1000, bg = "white")
ggsave (plot = Fig_Extended2_ploidy, filename = "Extended_Fig2.pdf", units = "cm", device = "pdf",width =22, height =22, dpi = 1000, bg = "white")
#################

############Figure Extended 3############
hybrid_nq <-  dplyr:::filter(Extended_genome, Strain=="3Hybrid")
hybrid_nq <-  dplyr:::filter(hybrid_nq, Type=="Evolved_NQO")

exampleHybrid <- dplyr:::filter(hybrid_nq, Strain=="3Hybrid")
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=1)
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=21)
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=10)
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=25)

nq_snps_HYBRID <-  exampleHybrid

#Normalize max and min values
nq_snps_HYBRID<- nq_snps_HYBRID %>% mutate(Normalized_read=if_else(Normalized_read >1.5, 1.5,Normalized_read))
nq_snps_HYBRID<- nq_snps_HYBRID %>% mutate(Normalized_read=if_else(Normalized_read < -1.5, -1.5,Normalized_read))

nq_genome<- nq_snps_HYBRID
legend_title <- "Relative Read Depth"

Fig_ext3 <-nq_genome %>% dplyr:::filter(Type=="Evolved_NQO") %>% ggplot(aes(x = as.numeric(start2),y = as.factor(chrom4),fill=Normalized_read)) + 
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
  ggtitle("Hybrids evolved in UV mimetic conditions")+
  theme(strip.text = element_text(size = 20))+
  theme(plot.title = element_text(size = 40, face = "bold"))+
  theme(legend.background = element_rect(fill="snow2",size=0.4, linetype="solid",colour ="grey"))+
  theme(legend.key.size = unit(2, 'cm'),
        legend.key.height = unit(2, 'cm'),
        legend.key.width = unit(2, 'cm'),
        legend.title = element_text(size=25,angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=0.1,
        legend.text = element_text(size=25))
#################
############Assemble and save Extended Figure 3############
ggsave (plot = Fig_ext3, filename = "Extended_Fig3_low_quality.jpg", units = "cm", device = "jpg",width =95, height =50, dpi = 300, bg = "white")
ggsave (plot = Fig_ext3, filename = "Extended_Fig3.png", units = "cm", device = "png",width =95, height =50, dpi = 800, bg = "white")
ggsave (plot = Fig_ext3, filename = "Extended_Fig3.jpg", units = "cm", device = "jpg",width =95, height =50, dpi = 800, bg = "white")
ggsave (plot = Fig_ext3, filename = "Extended_Fig3.svg", units = "cm", device = "svg",width =95, height =50, dpi = 1000, bg = "white")
ggsave (plot = Fig_ext3, filename = "Extended_Fig3.pdf", units = "cm", device = "pdf",width =95, height =50, dpi = 1000, bg = "white")
#################

############Figure Extended 4############
#Figure Ext4_A
size_chrom <- Extended_genome %>% group_by(chrom3, Replicate,Strain,Type) %>% dplyr::summarise(Size_chrom = last(end))
size_chrom <- size_chrom %>% dplyr:::filter(Type!="Ancestor")

LOH_events_counts_correlation<- Extended_LOH_chrom %>% dplyr:::filter(LOH=="LOH") %>% group_by(Type,LOH,chrom3) %>%   
  dplyr:::summarise(n()) %>% ungroup()
colnames(LOH_events_counts_correlation)[4]  <- "counts"

LOH_events_counts_correlation$Strain <- "3Hybrid"
LOH_events_counts_correlation$Strain <- "3Hybrid"

corr_aneu5a<-LOH_events_counts_correlation  %<>% mutate(same_chrom = ifelse(chrom3=="Scer chrI",1,
                                                                            ifelse(chrom3=="Scer chrII",2,
                                                                                   ifelse(chrom3=="Scer chrIII",3,
                                                                                          ifelse(chrom3=="Scer chrIV",4,
                                                                                                 ifelse(chrom3=="Scer chrV",5,
                                                                                                        ifelse(chrom3=="Scer chrVI",6,
                                                                                                               ifelse(chrom3=="Scer chrVII",7,
                                                                                                                      ifelse(chrom3=="Scer chrVIII",8,
                                                                                                                             ifelse(chrom3=="Scer chrIX",9,
                                                                                                                                    ifelse(chrom3=="Scer chrX",10,
                                                                                                                                           ifelse(chrom3=="Scer chrXI",11,
                                                                                                                                                  ifelse(chrom3=="Scer chrXII",12,
                                                                                                                                                         ifelse(chrom3=="Scer chrXIII",13,
                                                                                                                                                                ifelse(chrom3=="Scer chrXIV",14,
                                                                                                                                                                       ifelse(chrom3=="Scer chrXV",15,
                                                                                                                                                                              ifelse(chrom3=="Scer chrXVI",16,
                                                                                                                                                                                     ifelse(chrom3=="Spar chrI",1,
                                                                                                                                                                                            ifelse(chrom3=="Spar chrII",2,
                                                                                                                                                                                                   ifelse(chrom3=="Spar chrIII",3,
                                                                                                                                                                                                          ifelse(chrom3=="Spar chrIV",4,
                                                                                                                                                                                                                 ifelse(chrom3=="Spar chrV",5,
                                                                                                                                                                                                                        ifelse(chrom3=="Spar chrVI",6,
                                                                                                                                                                                                                               ifelse(chrom3=="Spar chrVII",7,
                                                                                                                                                                                                                                      ifelse(chrom3=="Spar chrVIII",8,
                                                                                                                                                                                                                                             ifelse(chrom3=="Spar chrIX",9,
                                                                                                                                                                                                                                                    ifelse(chrom3=="Spar chrX",10,
                                                                                                                                                                                                                                                           ifelse(chrom3=="Spar chrXI",11,
                                                                                                                                                                                                                                                                  ifelse(chrom3=="Spar chrXII",12,
                                                                                                                                                                                                                                                                         ifelse(chrom3=="Spar chrXIII",13,
                                                                                                                                                                                                                                                                                ifelse(chrom3=="Spar chrXIV",14,
                                                                                                                                                                                                                                                                                       ifelse(chrom3=="Spar chrXV",15,
                                                                                                                                                                                                                                                                                              ifelse(chrom3=="Spar chrXVI",16,NA)))))))))))))))))))))))))))))))))
corr_aneu5b <- ddply(corr_aneu5a, .(same_chrom), summarise, sum_counts=sum(counts))

corr_aneu5b<-corr_aneu5b %<>% mutate(same_chrom2 = ifelse(same_chrom==3,"chrIII",
                                                          ifelse(same_chrom==7,"chrVII",
                                                                 ifelse(same_chrom==10,"chrX",
                                                                        ifelse(same_chrom==11,"chrXI",
                                                                               ifelse(same_chrom==12,"chrXII",
                                                                                      ifelse(same_chrom==14,"chrXIV",
                                                                                             ifelse(same_chrom==15,"chrXV",
                                                                                                    ifelse(same_chrom==16,"chrXVI",NA)))))))))
corr_aneu5b <- corr_aneu5b %>% mutate(Size_chrom2=ifelse(same_chrom==15,1048000,
                                                         ifelse(same_chrom==12,1003000,
                                                                ifelse(same_chrom==14,742000,
                                                                       ifelse(same_chrom==11,643000,
                                                                              ifelse(same_chrom==3,295000,
                                                                                     ifelse(same_chrom==7,1058000,
                                                                                            ifelse(same_chrom==10,699000,
                                                                                                   ifelse(same_chrom==16,908000,NA)))))))))
corr_aneu_5 <- corr_aneu5b
cond.labs <- c(
  `Evolved_control` = "Evolved in control",
  `Evolved_NQO` = "Evolved in UV mimetic",
  `Ancestor` = "Ancestor",
  `1Scer` = "Saccharomyces cerevisiae",
  `2Spar` = "Saccharomyces paradoxus",
  `3Hybrid` = "Hybrid"
)
corr_aneu_5$Strain <- "3Hybrid" 

ggscatter(corr_aneu_5, x = "Size_chrom2", y = "sum_counts",
          color = "black", shape = 21, size = 3,
          add = "reg.line",  
          add.params = list(color = "blue", fill = "lightgray"), 
          cor.coef = TRUE,
          cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n"))


Fig_ext4_A <-corr_aneu_5 %>%filter(sum_counts>0)%>%
  ggplot(aes(x=Size_chrom2,y = sum_counts))+
  stat_smooth(method="lm", alpha=0.1)+
  annotate("text",
           y = 8,x =300000,
           label = c("chrIII"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y =8,x =643000,
           label = c("chrXI"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y = -5,x =699000,
           label = c("chrX"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y = 16,x =742000,
           label = c("chrXIV"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y = 10,x =908000,
           label = c("chrXVI"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y = 46,x =1003000,
           label = c("chrXII"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y = 14,x =1058000,
           label = c("chrVII"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y = 34,x =1048000,
           label = c("chrXV"),
           family = "", fontface = 3, size=2.2)+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, aes(fill=Strain), col="black", show.legend = F, alpha=0.5,size=2)+
  xlab("Chromosome Size (kbp)")+ylab("Total number of LOH")+
  border()+
  scale_x_continuous(breaks=c(400000,600000,800000,1000000),
                     labels=c(400,600,800,1000))+
  scale_shape_manual("",values=c(1,19),
                     label=c("Evolved in control","Evolved in UV mimetic"))+theme(strip.text.x = element_blank())+
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme_prism()+
  theme_bw(base_size=24) +
  theme(legend.position="none",
        axis.title = element_text(size=10, face = "bold"),
        strip.background = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        panel.background = element_blank())+
  annotate("text",
           y = 40,x = 450000,
           label = c("rs = 0.81, p < 0.05"),
           family = "", fontface = 3, size=3.5) 

#Figure Ext4_B
genomics1 <- dplyr:::filter(Extended_genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics2<-genomics1 %>% dplyr:::filter(Evolved=="Evolved_NQO")
genomics1 <- dplyr:::filter(Extended_genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics3<-genomics1 %>% dplyr:::filter(Evolved=="Ancestor")

genomics4<- left_join(genomics2,genomics3, by=c("Specie"="Specie","condition"="condition","Replicate"="Replicate"))

genomics4$fitness_gain<- (genomics4$rval.x - genomics4$rval.y) / genomics4$rval.y
genomics4$percentage <- genomics4$fitness_gain * 100

newloh <-Extended_counts_LOH

newloh$Replicate <- as.numeric(newloh$Replicate)

FITNESS_LOH <-left_join(newloh,genomics4,by=c("Type"="Evolved.x","Replicate"="Replicate"))

#Some functions
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

FITNESS_LOH_nq<- FITNESS_LOH %>% dplyr::filter(Type=="Evolved_NQO")

FITNESS_LOH_nq<- FITNESS_LOH_nq %>% dplyr::filter(Replicate!=25)
FITNESS_LOH_nq<- FITNESS_LOH_nq %>% dplyr::filter(Replicate!=10)

FITNESS_LOH_nq_hyb<-dplyr::filter(FITNESS_LOH_nq,Specie=="3Hybrid")
FITNESS_LOH_nq_hyb$new_counts <- FITNESS_LOH_nq_hyb$counts /2

ggscatter(FITNESS_LOH_nq_hyb, x = "new_counts", y = "percentage",
          color = "black", shape = 21, size = 3,
          add = "reg.line",  
          add.params = list(color = "blue", fill = "lightgray"), 
          cor.coef = TRUE,
          cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n"))

Fig_ext4_B <- FITNESS_LOH_nq_hyb %>% dplyr::filter(Type=="Evolved_NQO")%>%
  ggplot(aes(x=new_counts,y = percentage))+
  stat_smooth(method="lm", alpha=0.1)+
  annotate("text",
           y = 180,x = 0.6,
           label = c("rs = 0.33, p > 0.05"),
           family = "", fontface = 3, size=3.5) + 
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, aes(fill=Type), col="black", show.legend = F, alpha=0.5,size=2)+
  xlab("Number of LOH / line")+ylab("Increase in growth rate (%)")+
  border()  +
  scale_fill_manual(values=c("#FF9999"), 
                    labels = toexpr(c("Hybrid"))) +
  scale_color_manual(values=c("#FF9999"), 
                     labels = toexpr(c("Hybrid"))) +
  theme(legend.position="none")+
  stat_smooth(method="lm", alpha=0.1)+
  theme_prism()+
  theme_bw(base_size=24) +
  theme(legend.position="none",
        axis.title = element_text(size=10, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        panel.background = element_blank()) +theme(strip.text.x = element_blank())
#################
############Assemble and save Extended Figure 4############
Fig4A_label<- plot_grid(Fig_ext4_A,labels="a",label_size=15)
Fig4B_label<- plot_grid(Fig_ext4_B,labels="b",label_size=15)

Fig4_supp_LOH<-plot_grid(Fig4A_label,Fig4B_label)

ggsave (plot = Fig4_supp_LOH, filename = "Extended_Fig4_low_quality.jpg", units = "cm", device = "jpg",width = 20, height =8, dpi = 300, bg = "white")
ggsave (plot = Fig4_supp_LOH, filename = "Extended_Fig4.png", units = "cm", device = "png",width = 20, height =8, dpi = 1000, bg = "white")
ggsave (plot = Fig4_supp_LOH, filename = "Extended_Fig4.jpg", units = "cm", device = "jpg",width = 20, height =8, dpi = 1000, bg = "white")
ggsave (plot = Fig4_supp_LOH, filename = "Extended_Fig4.svg", units = "cm", device = "svg",width = 20, height =8, dpi = 1000, bg = "white")
ggsave (plot = Fig4_supp_LOH, filename = "Extended_Fig4.pdf", units = "cm", device = "pdf",width = 20, height =8, dpi = 1000, bg = "white")
#################

############Figure Extended 5############
SNP_after_norm2<-Extended_SNPs %<>% mutate(chrom2 = ifelse(CHROM=="utg351_pilon","Spar chrI",
                                                             ifelse(CHROM=="utg1271_pilon","Spar chrII",
                                                                    ifelse(CHROM=="utg584_pilon","Spar chrIII",
                                                                           ifelse(CHROM=="utg988_pilon","Spar chrIV",
                                                                                  ifelse(CHROM=="utg69_pilon","Spar chrV",
                                                                                         ifelse(CHROM=="utg639_pilon","Spar chrVI",
                                                                                                ifelse(CHROM=="utg199_pilon","Spar chrVII",
                                                                                                       ifelse(CHROM=="utg176_pilon","Spar chrVIII",
                                                                                                              ifelse(CHROM=="utg1121_pilon","Spar chrIX",
                                                                                                                     ifelse(CHROM=="utg245_pilon","Spar chrX",
                                                                                                                            ifelse(CHROM=="utg298_pilon","Spar chrXI",
                                                                                                                                   ifelse(CHROM=="utg110_pilon","Spar chrXII 1/2",
                                                                                                                                          ifelse(CHROM=="utg11_pilon","Spar chrXII 2/2",
                                                                                                                                                 ifelse(CHROM=="utg210_pilon","Spar chrXIII",
                                                                                                                                                        ifelse(CHROM=="utg48_pilon","Spar chrXIV",
                                                                                                                                                               ifelse(CHROM=="utg155_pilon","Spar chrXV 1/2",
                                                                                                                                                                      ifelse(CHROM=="utg122_pilon","Spar chrXV 2/2",
                                                                                                                                                                             ifelse(CHROM=="utg675_pilon","Spar chrXVI",
                                                                                                                                                                                    ifelse(CHROM=="chrI","Scer chrI",
                                                                                                                                                                                           ifelse(CHROM=="chrII","Scer chrII",
                                                                                                                                                                                                  ifelse(CHROM=="chrIII","Scer chrIII",
                                                                                                                                                                                                         ifelse(CHROM=="chrIV","Scer chrIV",
                                                                                                                                                                                                                ifelse(CHROM=="chrV","Scer chrV",
                                                                                                                                                                                                                       ifelse(CHROM=="chrVI","Scer chrVI",
                                                                                                                                                                                                                              ifelse(CHROM=="chrVII","Scer chrVII",
                                                                                                                                                                                                                                     ifelse(CHROM=="chrVIII","Scer chrVIII",
                                                                                                                                                                                                                                            ifelse(CHROM=="chrIX","Scer chrIX",
                                                                                                                                                                                                                                                   ifelse(CHROM=="chrX","Scer chrX",
                                                                                                                                                                                                                                                          ifelse(CHROM=="chrXI","Scer chrXI",
                                                                                                                                                                                                                                                                 ifelse(CHROM=="chrXII","Scer chrXII",
                                                                                                                                                                                                                                                                        ifelse(CHROM=="chrXIII","Scer chrXIII",
                                                                                                                                                                                                                                                                               ifelse(CHROM=="chrXIV","Scer chrXIV",
                                                                                                                                                                                                                                                                                      ifelse(CHROM=="chrXV","Scer chrXV",
                                                                                                                                                                                                                                                                                             ifelse(CHROM=="chrXVI","Scer chrXVI",NA)))))))))))))))))))))))))))))))))))
SNP_after_norm2$sample <- paste(SNP_after_norm2$Strain,SNP_after_norm2$Replicate,SNP_after_norm2$condition, sep="")
SNP_after_norm2_nooutliers <- SNP_after_norm2 %>% dplyr::filter(sample!="1Scer10NQO")
SNP_after_norm2_nooutliers <- SNP_after_norm2_nooutliers %>% dplyr::filter(sample!="2Spar1NQO")
SNP_after_norm2_nooutliers <- SNP_after_norm2_nooutliers %>% dplyr::filter(sample!="2Spar3NQO")
SNP_after_norm2_nooutliers <- SNP_after_norm2_nooutliers %>% dplyr::filter(sample!="2Spar5NQO")
SNP_after_norm2_nooutliers <- SNP_after_norm2_nooutliers %>% dplyr::filter(sample!="2Spar7NQO")
SNP_after_norm2_nooutliers <- SNP_after_norm2_nooutliers %>% dplyr::filter(sample!="2Spar10NQO")

#Get some functions
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

all2<-SNP_after_norm2_nooutliers

all2$Type <- all2$condition

count_data <- all2 %>%
  group_by(Strain,Type, Replicate) %>%
  mutate(count = n())

rows.per.group  <- aggregate(rep(1, length(paste0(all2$Strain, all2$Type,all2$Replicate))),
                             by=list(all2$Strain, all2$Type,all2$Replicate), sum)

colnames(rows.per.group) <- c("Strain","Type","Replicate","Counts")

rows.per.group%>%
  group_by(Type,Strain)%>% 
  summarise(Mean=mean(Counts), Max=max(Counts), Min=min(Counts), Median=median(Counts), Std=sd(Counts))

p_meds <- rows.per.group%>%
  group_by(Type,Strain)%>% 
  summarise(Median=median(Counts))

#Remove some lines
all2_2_nohyb<-dplyr::filter(rows.per.group,Strain!="3Hybrid")
all2_2_hyb<-dplyr::filter(rows.per.group,Strain=="3Hybrid")
all2_2_hyb_con<-dplyr::filter(all2_2_hyb,Type!="NQO")
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb,Type=="NQO")
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb_nq,Replicate!=1)
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb_nq,Replicate!=10)
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb_nq,Replicate!=21)
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb_nq,Replicate!=25)

rows.per.groupB<-rbind(all2_2_nohyb,all2_2_hyb_con,all2_2_hyb_nq)

rows.per.group<-rows.per.groupB

p_meds <- rows.per.group%>%
  group_by(Type,Strain)%>% 
  summarise(Median=median(Counts))

rows.per.group_nqo<- dplyr::filter(rows.per.group, Type=="NQO")

p_meds <- rows.per.group_nqo%>%
  group_by(Type,Strain)%>% 
  dplyr::summarise(Median=median(Counts))

my_comparisons <- list( c("NQO.1Scer", "NQO.2Spar"), c("NQO.1Scer", "NQO.3Hybrid"), c("NQO.2Spar", "NQO.3Hybrid") )

legend_title <- " "
Fig5SuppA_1 <- rows.per.group %>% dplyr:::filter(Type=="NQO") %>%
  ggplot(aes(x=interaction(Type,Strain),y=Counts)) +
  geom_boxplot(aes(fill=Strain))+
  geom_point(colour="black",pch=21,size=2, aes(fill=Strain))+
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid")))+
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85") + ylab("Counts") +
  geom_label(data = p_meds, aes(x = interaction(Type,Strain), y = Median, label =  Median),fontface = "bold",
             size = 3.5)+
  ggtitle("All variants in UV mimetic")+
  ylab("Count")+xlab(" ")+
  scale_x_discrete(" ", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))+
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(size = 9,face = "bold"),
        strip.background = element_blank())+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12))+
  theme(plot.title = element_text(hjust = 0.5,face="bold"))+
  annotate("text",
           y = c(440),
           x = c(1),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=3)+
  theme(legend.position="none")+
  ylim(0,450)

rows.per.group_nqo<- dplyr::filter(rows.per.group, Type=="control")

p_meds <- rows.per.group_nqo%>%
  group_by(Type,Strain)%>% 
  dplyr::summarise(Median=median(Counts))

my_comparisons <- list( c("control.1Scer", "control.2Spar"), c("control.1Scer", "control.3Hybrid"), c("control.2Spar", "control.3Hybrid") )

legend_title<- " "
Fig5SuppA_2 <- rows.per.group %>% dplyr::filter(Type=="control") %>%
  ggplot(aes(x=interaction(Type,Strain),y=Counts)) +
  geom_boxplot(aes(fill=Strain))+
  geom_point(colour="black",pch=21,size=2, aes(fill=Strain))+
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid")))+
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85") + ylab("Counts") +
  geom_label(data = p_meds, aes(x = interaction(Type,Strain), y = Median, label =  Median),fontface = "bold",
             size = 3.5)+
  ggtitle("All variants in control")+
  ylab("Count")+xlab(" ")+
  scale_x_discrete(" ", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))+
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(size = 9,face = "bold"),
        strip.background = element_blank())+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5,face="bold"))+
  annotate("text",
           y = c(440),
           x = c(1),
           label = c("p < 0.001"),
           family = "", fontface = 3, size=3)+
  theme(legend.position = "top",
        legend.direction = "horizontal")+
  ylim(0,450)+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title.align=0.1,
        legend.text = element_text(size=15))

leg <- get_legend(Fig5SuppA_2)

Fig5SuppA_2 <- rows.per.group %>% dplyr::filter(Type=="control") %>%
  ggplot(aes(x=interaction(Type,Strain),y=Counts)) +
  geom_boxplot(aes(fill=Strain))+
  geom_point(colour="black",pch=21,size=2, aes(fill=Strain))+
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid")))+
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85") + ylab("Counts") +
  geom_label(data = p_meds, aes(x = interaction(Type,Strain), y = Median, label =  Median),fontface = "bold",
             size = 3.5)+
  ggtitle("All variants in control")+
  ylab("Count")+xlab(" ")+
  scale_x_discrete(" ", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))+
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(size = 9,face = "bold"),
        strip.background = element_blank())+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5,face="bold"))+
  annotate("text",
           y = c(440),
           x = c(1),
           label = c("p < 0.001"),
           family = "", fontface = 3, size=3)+
  theme(legend.position = "none")+
  ylim(0,450)

#Missense variants
all2_2 <- dplyr::filter(Extended_GO,Consequence=="missense_variant")
all2_2_nohyb<-dplyr::filter(all2_2,Strain!="3Hybrid")
all2_2_hyb<-dplyr::filter(all2_2,Strain=="3Hybrid")
all2_2_hyb_con<-dplyr::filter(all2_2_hyb,condition!="NQO")
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb,condition=="NQO")
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb_nq,Replicate!=25)
all2_2_b <-rbind(all2_2_nohyb,all2_2_hyb_con,all2_2_hyb_nq)
all2_2<-all2_2_b
rows.per.group  <- aggregate(rep(1, length(paste0(all2_2$Strain, all2_2$condition,all2_2$Replicate))),
                             by=list(all2_2$Strain, all2_2$condition,all2_2$Replicate), sum)
colnames(rows.per.group) <- c("Strain","Type","Replicate","Counts")
medians <- aggregate(Counts ~  Strain*Type,rows.per.group, median)
con <- rows.per.group %>% dplyr::filter(Type=="control")
con2 <- medians %>% dplyr::filter(Type=="control")
nq <- rows.per.group %>% filter(Type=="NQO")
nq2 <- medians %>% dplyr::filter(Type=="NQO")

#Get some functions
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

Fig5SuppB_1 <- con %>% dplyr::filter(Type=="control") %>%
  ggplot(aes(x=Strain,y=Counts)) +
  geom_boxplot(aes(fill=Strain))+
  geom_point(colour="black",pch=21,size=2, aes(fill=Strain))+
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85") + ylab("Counts") +
  geom_label(data = con2, aes(x = Strain, y = Counts, label =  Counts),fontface = "bold",
             size = 3.5)+
  theme_bw(base_size=15) + 
  ylab("Count")+xlab(" ")+
  ggtitle("Missense variants in control")+
  theme(legend.position = "none")+
  scale_x_discrete(" ", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))+
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid")))+
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(size = 9,face = "bold"),
        strip.background = element_blank())+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  theme(legend.position="none")+
  annotate("text",
           y = c(440),
           x = c(1),
           label = c("p < 0.001"),
           family = "", fontface = 3, size=3)+
  ylim(0,450)

Fig5SuppB_2 <- nq %>% dplyr::filter(Type=="NQO") %>%
  ggplot(aes(x=Strain,y=Counts)) +
  geom_boxplot(aes(fill=Strain))+
  geom_point(colour="black",pch=21,size=2, aes(fill=Strain))+
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85") + ylab("Counts") +
  geom_label(data = nq2, aes(x = Strain, y = Counts, label =  Counts),fontface = "bold",
             size = 3.5)+
  theme_bw(base_size=15) + 
  ylab("Count")+xlab(" ")+
  ggtitle("Missense variants in UV mimetic")+
  theme(legend.position = "none")+
  scale_x_discrete(" ", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))+
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid")))+
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(size = 9,face = "bold"),
        strip.background = element_blank())+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  annotate("text",
           y = c(440),
           x = c(1),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=3)+
  ylim(0,450)+
  theme(legend.position = "none")

#################
############Assemble and save Extended Figure 5############
Figure5a<-plot_grid(Fig5SuppA_2,Fig5SuppA_1,rel_widths = c(1,1))
Figure5b<-plot_grid(Fig5SuppB_1,Fig5SuppB_2,rel_widths = c(1,1))
Fig5B_label<- plot_grid(Figure5b,labels="b",label_size=15)
Fig5A_label<- plot_grid(Figure5a,labels="a",label_size=15)
Figure5_Extended_SNP <- plot_grid(Fig5A_label,Fig5B_label, nrow = 2)
Figure5_Extended_SNPleg <- plot_grid(leg, Figure5_Extended_SNP, nrow=2,rel_heights = c(2,20))

ggsave (plot = Figure5_Extended_SNPleg, filename = "Extended_Fig5_low_quality.jpg", units = "cm", device = "jpg",width = 22, height = 22, dpi = 300, bg = "white")
ggsave (plot = Figure5_Extended_SNPleg, filename = "Extended_Fig5.png", units = "cm", device = "png",width = 22, height = 22, dpi = 1000, bg = "white")
ggsave (plot = Figure5_Extended_SNPleg, filename = "Extended_Fig5.jpg", units = "cm", device = "jpg",width = 22, height = 22, dpi = 1000, bg = "white")
ggsave (plot = Figure5_Extended_SNPleg, filename = "Extended_Fig5.svg", units = "cm", device = "svg",width = 22, height = 22, dpi = 1000, bg = "white")
ggsave (plot = Figure5_Extended_SNPleg, filename = "Extended_Fig5.pdf", units = "cm", device = "pdf",width = 22, height = 22, dpi = 1000, bg = "white")
#################

############Figure Extended 6############
#Positions of mutations
SNP <- c(280,280, 282,282, 298,308,308, 308,308,484,516,691,729,762,820,867,1041,1041,1047,1045)

sample.gr <- GRanges("chr7", IRanges(SNP, width=1, names=paste0("", SNP)))
features <- GRanges("chr7", IRanges(c(0),
                                    width=c(1050),
                                    names=paste0("block", 1:1)))
sample.gr$color <- sample.int(6, length(SNP), replace=TRUE)

sample.gr$color <- (c("dodgerblue1", "#FF9999",
                      "green4", "#FF9999",
                      "#FF9999", "dodgerblue1",
                      "dodgerblue1", "#FF9999",
                      "#FF9999", "dodgerblue1",
                      "#FF9999", "dodgerblue1",
                      "green4", "dodgerblue1",
                      "#FF9999", "dodgerblue1",
                      "dodgerblue1", "green4",
                      "#FF9999","#FF9999"))

plot <- lolliplot(sample.gr, features)
#################
############Assemble and save Extended Figure 6############
#We added some aesthetic modifications so we import it
gpp <- rasterGrob(Extended6_PDR1, interpolate=TRUE)

Extended_Fig6<-plot_grid(gpp)

ggsave (plot = Extended_Fig6, filename = "Extended_Fig6_low_quality.jpg", units = "cm", device = "jpg",width =18, height =3, dpi = 300, bg = "white")
ggsave (plot = Extended_Fig6, filename = "Extended_Fig6.png", units = "cm", device = "png",width =18, height =3, dpi = 1000, bg = "white")
ggsave (plot = Extended_Fig6, filename = "Extended_Fig6.jpg", units = "cm", device = "jpg",width =18, height =3, dpi = 1000, bg = "white")
ggsave (plot = Extended_Fig6, filename = "Extended_Fig6.svg", units = "cm", device = "svg",width =18, height =3, dpi = 1000, bg = "white")
ggsave (plot = Extended_Fig6, filename = "Extended_Fig6.pdf", units = "cm", device = "pdf",width =18, height =3, dpi = 1000, bg = "white")
#################

############Figure Extended 7############
#Create tree
tree <- read.tree(text = "(Nakaseomyces\nglabratus:0.02,(Saccharomyces\nbayanus:0.2,(Saccharomyces\nkudriavzevii:0.1,(Saccharomyces\nmikatae:0.2,(Saccharomyces\ncerevisiae:0.1,Saccharomyces\nparadoxus:0.1))))C:10);")

ggtree_obj <- ggtree(tree)

a<- ggtree_obj + geom_tiplab(fontface=4) +hexpand(0.5,direction=1)
#################
############Assemble and save Extended Figure 7############
#We added some aesthetic modifications so we import it
gpp <- rasterGrob(Extended7_candida, width = unit(1, "npc"), height = unit(1, "npc"))

Fig_extended7<-plot_grid(gpp)

ggsave (plot = Fig_extended7, filename = "Extended_Fig7_low_quality.jpg", units = "cm", device = "jpg",width =45, height =45, dpi = 300, bg = "white")
ggsave (plot = Fig_extended7, filename = "Extended_Fig7.png", units = "cm", device = "png",width =45, height =45, dpi = 1000, bg = "white")
ggsave (plot = Fig_extended7, filename = "Extended_Fig7.jpg", units = "cm", device = "jpg",width =45, height =45, dpi = 1000, bg = "white")
ggsave (plot = Fig_extended7, filename = "Extended_Fig7.svg", units = "cm", device = "svg",width =45, height =45, dpi = 1000, bg = "white")
ggsave (plot = Fig_extended7, filename = "Extended_Fig7.pdf", units = "cm", device = "pdf",width =45, height =45, dpi = 1000, bg = "white")
#################

############Figure Extended 8############
#Extended 8A
Extended_8A<- Extended_growth%>% dplyr::filter(Specie=="BY4741")%>%
  ggplot(aes(x = hour_Rounded, y = mean_od, ymin = mean_od - error_estandar, ymax = mean_od + error_estandar, fill = Mutation)) +
  geom_ribbon(alpha = 0.5) +
  scale_shape_identity() +  
  geom_line(aes(col=as.factor(Mutation),alpha=0.4)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("OD (595 nm)")+xlab("Time (hours)")+
  theme(axis.text.x = element_text()) +
  theme_bw()+
  xlim(0,40)+
  theme(legend.position = "none",
        axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_color_manual(values = c("indianred1","tan1","gold4","mediumaquamarine",
                                "violet","steelblue1"),
                     breaks= c("Empty", "Wild type (WT)", "G280R","G280S","M308I","G1041/2W"),
                     labels = c("Empty", "WT", "G280R","G280S","M308I","G1042W"))+
  scale_fill_manual(values = c("indianred1","tan1","gold4","mediumaquamarine",
                               "violet","steelblue1"),
                    breaks= c("Empty", "Wild type (WT)", "G280R","G280S","M308I","G1041/2W"),
                    labels = c("Empty", "WT", "G280R","G280S","M308I","G1042W"))
#Extended 8B
#Statistics
Extended_expression$Mutation <- as.factor(Extended_expression$Mutation)
amod <- aov(`GRN-B-HLog`~Mutation, data=Extended_expression)
summary(amod)
inter.test1 <- glht(amod,  mcp(Mutation = "Tukey"))
summary(inter.test1)
cld(inter.test1)

orden_desired <- c("Empty", "WT", "G280R","G280S","M308I","G1041/2W")
Extended_expression$Mutation <- factor(Extended_expression$Mutation, levels = orden_desired)

legend_title <- " "

Extended_8B<-Extended_expression %>% ggplot(aes(x = Mutation,y=`GRN-B-HLog`,fill=Mutation),  group=well) +
  geom_violin(trim=FALSE,alpha=0.7) + xlab("Mutation")+ylab("Pdr5-GFP (U.A. Fluorescence)")+ 
  theme(panel.spacing = unit(0.72, "cm"))+
  scale_color_manual(legend_title,values = c("indianred1","tan1","gold4","mediumaquamarine",
                                "violet","steelblue1"),
                     breaks= c("Empty", "WT", "G280R","G280S","M308I","G1041/2W"),
                     labels = c("Empty", "WT", "G280R","G280S","M308I","G1042W"))+
  scale_fill_manual(legend_title,values = c("indianred1","tan1","gold4","mediumaquamarine",
                               "violet","steelblue1"),
                    breaks= c("Empty", "WT", "G280R","G280S","M308I","G1041/2W"),
                    labels = c("Empty", "WT", "G280R","G280S","M308I","G1042W"))+
  guides(fill = guide_legend(nrow = 1))+
  annotate("text",
           y = c(4),
           x = c(1),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=5) +
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        legend.position = "top")

leg <- get_legend(Extended_8B)

Extended_8B<- Extended_expression %>% ggplot(aes(x = Mutation,y=`GRN-B-HLog`,fill=Mutation),  group=well) +
  geom_violin(trim=FALSE,alpha=0.7) + 
  xlab("Amino acid change")+ylab("Pdr5-GFP (A.U. Fluorescence)")+ 
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="none")+
  scale_fill_manual(values = c("indianred1","tan1",
                               "gold4",
                               "mediumaquamarine",
                               "violet","steelblue1"))+
  scale_x_discrete(breaks= c("Empty", "WT", "G280R","G280S","M308I","G1041/2W"),
                   labels = c("Empty", "WT", "G280R","G280S","M308I","G1042W"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank())+
  annotate("text",
           y = c(4),
           x = c(1.3),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=4) 

#################
############Assemble and save Extended Figure 8############
FigA_label<- plot_grid(Extended_8A,labels="a",label_size=15)
FigB_label<- plot_grid(Extended_8B,labels="b",label_size=15)

Fig4AB<- plot_grid(FigA_label,FigB_label,nrow=1,rel_widths = c(1,1))

Fig4ABleg<- plot_grid(leg,Fig4AB,nrow=2,rel_heights = c(0.2,1))

ggsave (plot = Fig4ABleg, filename = "Extended_Fig8_low_quality.jpg", units = "cm", device = "jpg",width =22, height =12, dpi = 300, bg = "white")
ggsave (plot = Fig4ABleg, filename = "Extended_Fig8.png", units = "cm", device = "png",width =22, height =12, dpi = 1000, bg = "white")
ggsave (plot = Fig4ABleg, filename = "Extended_Fig8.jpg", units = "cm", device = "jpg",width =22, height =12, dpi = 1000, bg = "white")
ggsave (plot = Fig4ABleg, filename = "Extended_Fig8.svg", units = "cm", device = "svg",width =22, height =12, dpi = 1000, bg = "white")
ggsave (plot = Fig4ABleg, filename = "Extended_Fig8.pdf", units = "cm", device = "pdf",width =22, height =12, dpi = 1000, bg = "white")
#################

############Figure Extended 9############
#Figure A
all_arranged2<- Extended9_boxplot

all_arranged_BY<- all_arranged2 %>% filter(Specie=="BY")
Figure1E_BY <- all_arranged_BY %>% filter(Condition.y=="NQO_4μM")%>% 
  ggplot(aes(x=interaction(Type2,Specie), y=rval,group=interaction(Specie,Type2))) +
  geom_point(colour="black",pch=21,height = 0, size=2, aes(fill=as.factor(Type2)))+
  geom_boxplot(outlier.shape = NA,aes(col=Mutation,fill=Mutation,alpha=0.7, group=interaction(Specie,Type2)))+
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "Ancestral_homozygous" = "gray30", "Haploid"="antiquewhite1"))+
  scale_fill_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "Ancestral_homozygous" = "gray30", "Haploid"="antiquewhite1"))+
  ylab("Growth rate \n (OD/hour)") +
  xlab("")+
  theme_bw() +
  facet_grid(.~Specie)+
  theme(axis.title = element_text(size=20, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=12,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  ylim(0,0.7)+
  theme(axis.text.x = element_text(angle = 90))

#make pvalue
all_arranged_BY<- all_arranged2 %>% filter(Specie=="BY")
all_arranged_BY$Specie<-"S.cerevisiae - BY4743" 

all_arranged_BY$Type2 <- as.factor(all_arranged_BY$Type2)
amod <- aov(rval~Type2, data=all_arranged_BY)
summary(amod)
inter.test1 <- glht(amod,  mcp(Type2= "Tukey"))
summary(amod)
summary(inter.test1)
cld(inter.test1)

all_arranged_BY<- all_arranged2 %>% filter(Specie=="Hybrid")

all_arranged_BY$Type2 <- as.factor(all_arranged_BY$Type2)
amod <- aov(rval~Type2, data=all_arranged_BY)
summary(amod)
inter.test1 <- glht(amod,  mcp(Type2= "Tukey"))
summary(amod)
summary(inter.test1)
cld(inter.test1)

all_arranged_BY<- all_arranged2 %>% filter(Specie=="Spar")

all_arranged_BY$Type2 <- as.factor(all_arranged_BY$Type2)
amod <- aov(rval~Type2, data=all_arranged_BY)
summary(amod)
inter.test1 <- glht(amod,  mcp(Type2= "Tukey"))
summary(amod)
summary(inter.test1)
cld(inter.test1)

all_arranged_BY<- all_arranged2 %>% filter(Specie=="BY")
all_arranged_BY$Specie<-"S.cerevisiae - BY4743" 

legend_title<-""
Figure1E_BY <- all_arranged_BY %>% filter(Condition.y == "NQO_4μM") %>%
  ggplot(aes(x = Type2, y = rval)) +
  geom_boxplot(colour = "black", outlier.shape = NA, aes(colour = "black", fill = Type2)) +
  geom_point(colour = "black", pch = 21, size = 2, aes(fill = Type2)) +
  scale_fill_manual(legend_title,
                    values = c("gray30","mediumpurple3","#21908CFF"),
                    limits=c("Heterozygous",
                             "Homozygous",
                             "Ancestral_homozygous"),
                    breaks= c("Heterozygous",
                              "Homozygous",
                              "Ancestral_homozygous"),
                    labels= c("No mutation",
                              "Heterozygous",
                              "Homozygous")) +
  scale_x_discrete(limits=c("Ancestral_homozygous",
                            "Heterozygous",
                            "Homozygous"),
                   breaks= c("Ancestral_homozygous",
                             "Heterozygous",
                             "Homozygous"),
                   labels= c("No \n mutation",
                             "Heterozygous",
                             "Homozygous"))+
  ylab("Growth rate (OD/hour)") +theme_bw() +
  ylim(0,0.65)+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_text(size=14, face = "bold"),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  theme(strip.text = element_text(face = "italic"))+
  annotate("text",
           y = c(0.64),
           x = c(1),
           label = c("p < 0.001"),
           family = "", fontface = 3, size=5)+
  theme(legend.direction="horizontal",
        legend.key.width = unit(8, 'mm'),
        legend.key.height = unit(8, 'mm'),
        legend.text = element_text(size = 12))

leg <- get_legend(Figure1E_BY)

labels_custom <- function(variable, value) {
  return(lapply(value, function(x) {
    bquote(italic(S.~cerevisiae) ~ "- BY4743")
  }))
}

Figure1E_BY <- all_arranged_BY %>% filter(Condition.y == "NQO_4μM") %>%
  ggplot(aes(x = Type2, y = rval)) +
  geom_boxplot(colour = "black", outlier.shape = NA, aes(colour = "black", fill = Type2, alpha = 0.7)) +
  geom_point(colour = "black", pch = 21, size = 2, aes(fill = Type2)) +
  scale_fill_manual(
    values = c("mediumpurple3",
               "#21908CFF","gray30"),
    limits=c("Heterozygous",
             "Homozygous",
             "Ancestral_homozygous"),
    breaks= c("Heterozygous",
              "Homozygous",
              "Ancestral_homozygous")) +
  scale_x_discrete(limits=c("Ancestral_homozygous",
                            "Heterozygous",
                            "Homozygous"),
                   breaks= c("Ancestral_homozygous",
                             "Heterozygous",
                             "Homozygous"),
                   labels= c("No \n mutation",
                             "Heterozygous",
                             "Homozygous"))+
  ylab("Growth rate (OD/hour)") +theme_bw() +
  theme(legend.position = "none")+
  facet_grid(. ~ Specie, labeller = labels_custom)+
  ylim(0,0.65)+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_text(size=14, face = "bold"),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  theme(strip.text = element_text(face = "italic"))+
  annotate("text",
           y = c(0.64),
           x = c(1),
           label = c("p < 0.05"),
           family = "", fontface = 3, size=5)

all_arranged_BY2<- all_arranged2 %>% filter(Specie=="Hybrid")

outcome_labels <- c("No mutation", 
                    "Heterozygous /n Spar", 
                    "Heterozygous /n Scer",
                    "Homozygous")


Figure1E_BY2 <- all_arranged_BY2 %>% filter(Condition.y == "NQO_4μM") %>%
  ggplot(aes(x = Type2, y = rval)) +
  geom_boxplot(colour = "black", outlier.shape = NA, aes(colour = "black", fill = Type2, alpha = 0.7)) +
  geom_point(colour = "black", pch = 21, size = 2, aes(fill = Type2)) +
  scale_fill_manual(
    values = c("mediumpurple3", "mediumpurple3",
               "#21908CFF","gray30"),
    limits=c("Heterozygous_BY-M308I",
             "Heterozygous_Spar-M308I",
             "Homozygous",
             "Ancestral_homozygous"),
    breaks= c("Heterozygous_BY-M308I",
              "Heterozygous_Spar-M308I",
              "Homozygous",
              "Ancestral_homozygous")) +
  scale_x_discrete(limits=c("Ancestral_homozygous",
                            "Heterozygous_BY-M308I",
                            "Heterozygous_Spar-M308I",
                            "Homozygous"),
                   breaks= c("Ancestral_homozygous",
                             "Heterozygous_BY-M308I",
                             "Heterozygous_Spar-M308I",
                             "Homozygous"),
                   labels= c("No \n mutation",
                             "Heterozygous \n Scer-M308I",
                             "Heterozygous \n Spar-M307I",
                             "Homozygous"))+
  ylab("Growth rate (OD/hour)") +theme_bw() +
  theme(legend.position = "none")+
  facet_grid(. ~ Specie)+
  ylim(0,0.65)+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.ticks.y =element_blank(),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  annotate("text",
           y = c(0.64),
           x = c(1),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=5) 

all_arranged_BY3<- all_arranged2 %>% filter(Specie=="Spar")

labels_custom <- function(variable, value) {
  return(lapply(value, function(x) {
    bquote(italic(S.~paradoxus) ~ " ")
  }))
}

all_arranged_BY3$Specie<-"S.paradoxus" 
Figure1E_BY3 <- all_arranged_BY3 %>% filter(Condition.y == "NQO_4μM") %>%
  ggplot(aes(x = Type2, y = rval)) +
  geom_boxplot(colour = "black", outlier.shape = NA, aes(colour = "black", fill = Type2, alpha = 0.7)) +
  geom_point(colour = "black", pch = 21, size = 2, aes(fill = Type2)) +
  scale_fill_manual(
    values = c("mediumpurple3",
               "#21908CFF","gray30"),
    limits=c("Heterozygous",
             "Homozygous",
             "Ancestral_homozygous"),
    breaks= c("Heterozygous",
              "Homozygous",
              "Ancestral_homozygous")) +
  scale_x_discrete(limits=c("Ancestral_homozygous",
                            "Heterozygous",
                            "Homozygous"),
                   breaks= c("Ancestral_homozygous",
                             "Heterozygous",
                             "Homozygous"),
                   labels= c("No \n mutation",
                             "Heterozygous",
                             "Homozygous"))+
  ylab("Growth rate (OD/hour)") +theme_bw() +
  theme(legend.position = "none")+
  facet_grid(. ~ Specie, labeller = labels_custom)+
  ylim(0,0.65)+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.ticks.y =element_blank(),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  theme(strip.text = element_text(face = "italic"))+
  annotate("text",
           y = c(0.64),
           x = c(1),
           label = c("p < 0.001"),
           family = "", fontface = 3, size=5) 

#Arrange a
Figure1E_BY4<-plot_grid(Figure1E_BY,Figure1E_BY2,Figure1E_BY3,nrow=1,rel_widths = c(1.2,1.3,1))
FigureSuppAleg<-plot_grid(leg, Figure1E_BY4,nrow=2,rel_heights = c(1,7))

#Figure B
all_arranged2 <- Extended9_curves
all_arranged2<- all_arranged2%>% filter(Condition.y!="Control")
all_arranged_BY<- all_arranged2 %>% filter(Specie=="BY")

Figure1E_BY <- all_arranged_BY %>%  dplyr::filter(hour<21)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Type,Replicate,Mutation,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Mutation, alpha=0.4)) +
  ylab("OD \n (595 nm)")+xlab("Time (hours)")+
  theme_bw() +theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=7)) +
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "No mutation" = "gray30"),
                     labels = c("Heterozygous" = "Heterozygous PDR1 mutation","No mutation" = "No PDR1 mutation", "Homozygous" = "Homozygous PDR1 mutation")) +
  theme(axis.title = element_text(size=20, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=16,face = "bold"),
        axis.text.y = element_text(size=16,face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  theme(strip.text.x = element_text(size = 12),    
        strip.text.y = element_text(size = 12),      
        strip.background = element_rect(fill = "white"))+
  ylab("OD (595 nm)")+xlab("Time (hours)")

all_arranged_BY<- all_arranged2 %>% filter(Specie=="BY")
all_arranged_BY$Specie<-"S.cerevisiae - BY4743" 

labels_custom <- function(variable, value) {
  return(lapply(value, function(x) {
    bquote(italic(S.~cerevisiae) ~ "- BY4743")
  }))
}

Figure1E_BY <- all_arranged_BY %>%  dplyr::filter(hour<21)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Type2,Type,Replicate,Mutation,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Mutation, alpha=0.4)) +
  facet_grid(. ~ Specie, labeller = labels_custom)+
  ylab("OD \n (595 nm)")+xlab("Time (hours)")+
  theme_bw() +theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=7)) +
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "No mutation" = "gray30"),
                     labels = c("Heterozygous" = "Heterozygous PDR1 mutation","No mutation" = "No PDR1 mutation", "Homozygous" = "Homozygous PDR1 mutation")) +
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_text(size=14, face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  theme(strip.text = element_text(face = "italic"))+
  ylab("OD (595 nm)")+xlab("Time (hours)")

all_arranged_BY2<- all_arranged2 %>% filter(Specie=="Hybrid")

##add data for legend purposes
data_leg <- data.frame(
  Well = c("...1"),
  Day = c(1),
  Plate = c("Plate1"),
  Condition.x = c("Condition1"),
  time = c(0),
  od = c(NA),
  Specie = c("Hybrid"),
  Type = c("Type1"),
  Line = c("Line1"),
  Replicate = c(1),
  Mutation = c("Heterozygous"),
  Condition.y = c("Condition1"),
  Assay = c("Assay1"),
  hour = c(0),
  Type2 = c("Heterozygous")
)

all_arranged_BY2B <-all_arranged_BY2[, -c(1, 2)]

all_arranged_BY2B<-rbind(data_leg,all_arranged_BY2B)

legend_title <- ""
Figure1E_BY2 <- all_arranged_BY2B %>% dplyr::filter(hour<21)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Type2,Type,Replicate,Mutation,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Type2)) +
  scale_color_manual(legend_title,
                     values = c("mediumpurple3",
                                "darkblue",
                                "magenta3",
                                "#21908CFF",
                                "gray30"),
                     limits=c("Heterozygous",
                              "Heterozygous_BY-M308I",
                              "Heterozygous_Spar-M308I",
                              "Homozygous",
                              "Ancestral_homozygous"),
                     breaks= c("Heterozygous",
                               "Heterozygous_BY-M308I",
                               "Heterozygous_Spar-M308I",
                               "Homozygous",
                               "Ancestral_homozygous"),
                     labels= c("Heterozygous",
                               "Heterozygous Scer-M308I",
                               "Heterozygous Spar-M307I",
                               "Homozygous",
                               "No mutation")) +
  ylab("Growth rate (OD/hour)") +theme_bw() +
  theme(legend.direction =  "horizontal")+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_text(size=14, face = "bold"),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14),
        legend.key.width = unit(8, 'mm'),
        legend.key.height = unit(8, 'mm'),
        legend.text = element_text(size = 12))

leg <- get_legend(Figure1E_BY2)

Figure1E_BY2 <- all_arranged_BY2 %>% dplyr::filter(hour<21)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Type2,Type,Replicate,Mutation,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Type2, alpha=0.4)) +
  scale_color_manual(
    values = c("darkblue", "magenta3",
               "#21908CFF","gray30"),
    limits=c("Heterozygous_BY-M308I",
             "Heterozygous_Spar-M308I",
             "Homozygous",
             "Ancestral_homozygous"),
    breaks= c("Heterozygous_BY-M308I",
              "Heterozygous_Spar-M308I",
              "Homozygous",
              "Ancestral_homozygous")) +
  ylab("Growth rate (OD/hour)") +theme_bw() +
  facet_grid(. ~ Specie)+
  theme(legend.position = "none")+
  theme(axis.title = element_text(size=20, face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.text.y =element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  ylab("OD (595 nm)")+xlab("Time (hours)")

outcome_labels <- c("No mutation", 
                    "Heterozygous /n Spar", 
                    "Heterozygous /n Scer",
                    "Homozygous")

all_arranged_BY3<- all_arranged2 %>% filter(Specie=="Spar")

all_arranged_BY3$Specie<-"S.paradoxus" 

labels_custom <- function(variable, value) {
  return(lapply(value, function(x) {
    bquote(italic(S.~paradoxus) ~ " ")
  }))
}

Figure1E_BY3 <- all_arranged_BY3 %>%  dplyr::filter(hour<21)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Type2,Type,Replicate,Mutation,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Type2, alpha=0.4)) +
  facet_grid(. ~ Specie, labeller = labels_custom)+
  scale_color_manual(
    values = c("mediumpurple3",
               "#21908CFF","gray30"),
    limits=c("Heterozygous",
             "Homozygous",
             "Ancestral_homozygous"),
    breaks= c("Heterozygous",
              "Homozygous",
              "Ancestral_homozygous")) +
  ylab("Growth rate (OD/hour)") +theme_bw() +
  theme(legend.position = "none")+
  theme(axis.title = element_text(size=20, face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.text.y =element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  theme(strip.text = element_text(face = "italic"))+
  ylab("OD (595 nm)")+xlab("Time (hours)")

#arrange b
Figure1E_BY4<-plot_grid(Figure1E_BY,Figure1E_BY2,Figure1E_BY3,nrow=1,rel_widths = c(1.2,1.3,1))
FigureSuppBleg<-plot_grid(leg, Figure1E_BY4,nrow=2,rel_heights = c(1,7))

#################
############Assemble and save Extended Figure 9############
###Save panels 
FigA_label<- plot_grid(FigureSuppAleg,labels="a",label_size=32)
FigB_label<- plot_grid(FigureSuppBleg,labels="b",label_size=32)

FigureSuppleg<-plot_grid(FigA_label,FigB_label,nrow=2)

Fig_Extended_9<-FigureSuppleg

ggsave (plot = Fig_Extended_9, filename = "Extended_Fig9_low_quality.jpg", units = "cm", device = "jpg",width = 45, height =25, dpi = 300, bg = "white")
ggsave (plot = Fig_Extended_9, filename = "Extended_Fig9.png", units = "cm", device = "png",width = 45, height =25, dpi = 1000, bg = "white")
ggsave (plot = Fig_Extended_9, filename = "Extended_Fig9.jpg", units = "cm", device = "jpg",width = 45, height =25, dpi = 1000, bg = "white")
ggsave (plot = Fig_Extended_9, filename = "Extended_Fig9.svg", units = "cm", device = "svg",width = 45, height =25, dpi = 1000, bg = "white")
ggsave (plot = Fig_Extended_9, filename = "Extended_Fig9.pdf", units = "cm", device = "pdf",width = 45, height =25, dpi = 1000, bg = "white")
#################

############Figure Extended 10############
#Hybrid 13
Hybrid <- dplyr::filter(Extended_expevol, strain=="3Hybrid")
Hybrid_30 <- dplyr::filter(Hybrid, rep==13)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")

Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
scer25_1<- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

scer25_1<-scer25_1 + scale_color_grey() + theme_classic() 

loh_table <- Extended10_LOH_line13_final

loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- loh_table
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
scer25_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)") 

scer25_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() + 
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

data<-Extended10_LOH_parents_table_proportion_hybrid %>% dplyr::filter(Replicate==13)

scer25_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(scer25_1,labels="A",label_size=10)
Fig1B2<- plot_grid(scer25_2,labels="B",label_size=10)
Fig1C2<- plot_grid(scer25_3,labels="C",label_size=10)
Fig1D2<- plot_grid(scer25_4,labels="D",label_size=10)

Figure1<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("Hybrid Line 13", fontface='bold')
Figure1<-plot_grid(title,Figure1, ncol=1, rel_heights=c(0.1, 1))

#Hybrid 28
Hybrid <- dplyr::filter(Extended_expevol, strain=="3Hybrid")
Hybrid_30 <- dplyr::filter(Hybrid, rep==28)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")

Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
scer25_1<- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

scer25_1<-scer25_1 + scale_color_grey() + theme_classic() 

loh_table <- Extended10_LOH_line28_final

loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day

loh_tableB <- loh_table
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
scer25_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)") 

scer25_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() + 
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

data<-Extended10_LOH_parents_table_proportion_hybrid %>% dplyr::filter(Replicate==28)

scer25_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(scer25_1,labels="A",label_size=10)
Fig1B2<- plot_grid(scer25_2,labels="B",label_size=10)
Fig1C2<- plot_grid(scer25_3,labels="C",label_size=10)
Fig1D2<- plot_grid(scer25_4,labels="D",label_size=10)

Figure2<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("Hybrid Line 28", fontface='bold')
Figure2<-plot_grid(title,Figure2, ncol=1, rel_heights=c(0.1, 1))

#Hybrid 30
Hybrid <- dplyr::filter(Extended_expevol, strain=="3Hybrid")
Hybrid_30 <- dplyr::filter(Hybrid, rep==30)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")

Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
scer25_1<- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

scer25_1<-scer25_1 + scale_color_grey() + theme_classic() 

loh_table <- Extended10_LOH_line30_final

loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- loh_table
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
scer25_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)") 

scer25_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() + 
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

loh_table<-loh_table %<>% mutate(LOH_other_positions= ifelse(LOH_1=="Heterozygous","No",
                                                             ifelse(LOH_1=="Homozygous","Yes",NA)))

data<-Extended10_LOH_parents_table_proportion_hybrid %>% dplyr::filter(Replicate==30)

scer25_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(scer25_1,labels="A",label_size=10)
Fig1B2<- plot_grid(scer25_2,labels="B",label_size=10)
Fig1C2<- plot_grid(scer25_3,labels="C",label_size=10)
Fig1D2<- plot_grid(scer25_4,labels="D",label_size=10)

Figure3<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("Hybrid Line 30", fontface='bold')
Figure3<-plot_grid(title,Figure3, ncol=1, rel_heights=c(0.1, 1))

#Scer 25
Hybrid <- dplyr::filter(Extended_expevol, strain=="1Scer")
Hybrid_30 <- dplyr::filter(Hybrid, rep==25)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")

Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
scer25_1<- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

scer25_1<-scer25_1 + scale_color_grey() + theme_classic() 

loh_table<-Extended10_LOH_table_parents
loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- dplyr::filter(loh_table, Genotype=="Scer")
loh_tableB <- dplyr::filter(loh_tableB, Line==25)
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
scer25_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)") + facet_grid(.~Genotype)

scer25_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("#21908CFF", "gray30")) 

data<-dplyr::filter(Extended10_LOH_parents_table_proportion, Replicate==25)

scer25_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(scer25_1,labels="A",label_size=10)
Fig1B2<- plot_grid(scer25_2,labels="B",label_size=10)
Fig1C2<- plot_grid(scer25_3,labels="C",label_size=10)
Fig1D2<- plot_grid(scer25_4,labels="D",label_size=10)

Figure4<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("S. cerevisiae Line 25", fontface='bold')
Figure4<-plot_grid(title,Figure4, ncol=1, rel_heights=c(0.1, 1))

#Spar 22
Hybrid <- dplyr::filter(Extended_expevol, strain=="2Spar")
Hybrid_30 <- dplyr::filter(Hybrid, rep==22)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")
Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                  ifelse(day==4,20,
                                                                         ifelse(day==5,25,
                                                                                ifelse(day==6,30,
                                                                                       ifelse(day==7,35,
                                                                                              ifelse(day==8,40,
                                                                                                     ifelse(day==9,45,
                                                                                                            ifelse(day==10,50,
                                                                                                                   ifelse(day==12,55,
                                                                                                                          ifelse(day==13,65,
                                                                                                                                 ifelse(day==14,70,
                                                                                                                                        ifelse(day==15,75,
                                                                                                                                               ifelse(day==16,80,
                                                                                                                                                      ifelse(day==17,85,
                                                                                                                                                             ifelse(day==18,90,
                                                                                                                                                                    ifelse(day==19,95,
                                                                                                                                                                           ifelse(day==20,100,
                                                                                                                                                                                  ifelse(day==21,105,NA)))))))))))))))))))
spar22_1 <- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
    geom_line()+
    geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
    ylab("Growth rate (OD/hour)")
  
spar22_1<-spar22_1 + scale_color_grey() + theme_classic() 

loh_table<-Extended10_LOH_table_parents
loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- dplyr::filter(loh_table, Genotype=="Spar")
loh_tableB <- dplyr::filter(loh_tableB, Line==22)
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                          ifelse(day==3,15,
                                                                 ifelse(day==4,20,
                                                                        ifelse(day==5,25,
                                                                               ifelse(day==6,30,
                                                                                      ifelse(day==7,35,
                                                                                             ifelse(day==8,40,
                                                                                                    ifelse(day==9,45,
                                                                                                           ifelse(day==10,50,
                                                                                                                  ifelse(day==12,55,
                                                                                                                         ifelse(day==13,65,
                                                                                                                                ifelse(day==14,70,
                                                                                                                                       ifelse(day==15,75,
                                                                                                                                              ifelse(day==16,80,
                                                                                                                                                     ifelse(day==17,85,
                                                                                                                                                            ifelse(day==18,90,
                                                                                                                                                                   ifelse(day==19,95,
                                                                                                                                                                          ifelse(day==20,100,
                                                                                                                                                                                 ifelse(day==21,105,NA))))))))))))))))))))
spar22_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
    geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
    theme_classic() +
    ylab("Relative Frequency")+xlab("Time (Cycle)")
  
spar22_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
    geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
    theme_classic() + 
    ylab("Relative Frequency")+xlab("Time (Cycle)")+
    scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 
  
data<-dplyr::filter(Extended10_LOH_parents_table_proportion, Replicate==22)

spar22_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
    geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
    scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
    ylab("Relative Frequency")+xlab("Time (Cycles)") +
    scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 
  
Fig1A2<- plot_grid(spar22_1,labels="A",label_size=10)
Fig1B2<- plot_grid(spar22_2,labels="B",label_size=10)
Fig1C2<- plot_grid(spar22_3,labels="C",label_size=10)
Fig1D2<- plot_grid(spar22_4,labels="D",label_size=10)
  
Figure5<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("S. paradoxus Line 22", fontface='bold')
Figure5<-plot_grid(title,Figure5, ncol=1, rel_heights=c(0.1, 1))
  
#Scer 2
Hybrid <- dplyr::filter(Extended_expevol, strain=="1Scer")
Hybrid_30 <- dplyr::filter(Hybrid, rep==2)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")
Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                    ifelse(day==4,20,
                                                                           ifelse(day==5,25,
                                                                                  ifelse(day==6,30,
                                                                                         ifelse(day==7,35,
                                                                                                ifelse(day==8,40,
                                                                                                       ifelse(day==9,45,
                                                                                                              ifelse(day==10,50,
                                                                                                                     ifelse(day==12,55,
                                                                                                                            ifelse(day==13,65,
                                                                                                                                   ifelse(day==14,70,
                                                                                                                                          ifelse(day==15,75,
                                                                                                                                                 ifelse(day==16,80,
                                                                                                                                                        ifelse(day==17,85,
                                                                                                                                                               ifelse(day==18,90,
                                                                                                                                                                      ifelse(day==19,95,
                                                                                                                                                                             ifelse(day==20,100,
                                                                                                                                                                                    ifelse(day==21,105,NA)))))))))))))))))))
scer2_1 <- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
      geom_line()+
      geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
      ylab("Growth rate (OD/hour)")

scer2_1<-scer2_1+ scale_color_grey() + theme_classic() 

loh_table<-Extended10_LOH_table_parents
loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- dplyr::filter(loh_table, Genotype=="Scer")
loh_tableB <- dplyr::filter(loh_tableB, Line==2)
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                            ifelse(day==3,15,
                                                                   ifelse(day==4,20,
                                                                          ifelse(day==5,25,
                                                                                 ifelse(day==6,30,
                                                                                        ifelse(day==7,35,
                                                                                               ifelse(day==8,40,
                                                                                                      ifelse(day==9,45,
                                                                                                             ifelse(day==10,50,
                                                                                                                    ifelse(day==12,55,
                                                                                                                           ifelse(day==13,65,
                                                                                                                                  ifelse(day==14,70,
                                                                                                                                         ifelse(day==15,75,
                                                                                                                                                ifelse(day==16,80,
                                                                                                                                                       ifelse(day==17,85,
                                                                                                                                                              ifelse(day==18,90,
                                                                                                                                                                     ifelse(day==19,95,
                                                                                                                                                                            ifelse(day==20,100,
                                                                                                                                                                                   ifelse(day==21,105,NA))))))))))))))))))))
scer2_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
      geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
      theme_classic() +
      ylab("Relative Frequency")+xlab("Time (Cycle)")
    
    
scer2_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
      geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
      theme_classic() + 
      ylab("Relative Frequency")+xlab("Time (Cycle)")+
      scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 
    
data<-dplyr::filter(Extended10_LOH_parents_table_proportion, Replicate==2)
    
scer2_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
      geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
      scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
      ylab("Relative Frequency")+xlab("Time (Cycles)") +
      scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 
    
Fig1A2<- plot_grid(scer2_1,labels="A",label_size=10)
Fig1B2<- plot_grid(scer2_2,labels="B",label_size=10)
Fig1C2<- plot_grid(scer2_3,labels="C",label_size=10)
Fig1D2<- plot_grid(scer2_4,labels="D",label_size=10)
    
Figure6<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("S. cerevisiae Line 2", fontface='bold')
Figure6<-plot_grid(title,Figure6, ncol=1, rel_heights=c(0.1, 1))

#Spar 24
Hybrid <- dplyr::filter(Extended_expevol, strain=="2Spar")
Hybrid_30 <- dplyr::filter(Hybrid, rep==24)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")
Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                      ifelse(day==4,20,
                                                                             ifelse(day==5,25,
                                                                                    ifelse(day==6,30,
                                                                                           ifelse(day==7,35,
                                                                                                  ifelse(day==8,40,
                                                                                                         ifelse(day==9,45,
                                                                                                                ifelse(day==10,50,
                                                                                                                       ifelse(day==12,55,
                                                                                                                              ifelse(day==13,65,
                                                                                                                                     ifelse(day==14,70,
                                                                                                                                            ifelse(day==15,75,
                                                                                                                                                   ifelse(day==16,80,
                                                                                                                                                          ifelse(day==17,85,
                                                                                                                                                                 ifelse(day==18,90,
                                                                                                                                                                        ifelse(day==19,95,
                                                                                                                                                                               ifelse(day==20,100,
                                                                                                                                                                                      ifelse(day==21,105,NA)))))))))))))))))))
spar24_1 <- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
        geom_line()+
        geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
        ylab("Growth rate (OD/hour)")
      
spar24_1<-spar24_1 + scale_color_grey() + theme_classic() #+ ggtitle ("Line 30") + theme(plot.title = element_text(hjust = 0.5, size=20))

loh_table<-Extended10_LOH_table_parents
loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- dplyr::filter(loh_table, Genotype=="Spar")
loh_tableB <- dplyr::filter(loh_tableB, Line==24)
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                              ifelse(day==3,15,
                                                                     ifelse(day==4,20,
                                                                            ifelse(day==5,25,
                                                                                   ifelse(day==6,30,
                                                                                          ifelse(day==7,35,
                                                                                                 ifelse(day==8,40,
                                                                                                        ifelse(day==9,45,
                                                                                                               ifelse(day==10,50,
                                                                                                                      ifelse(day==12,55,
                                                                                                                             ifelse(day==13,65,
                                                                                                                                    ifelse(day==14,70,
                                                                                                                                           ifelse(day==15,75,
                                                                                                                                                  ifelse(day==16,80,
                                                                                                                                                         ifelse(day==17,85,
                                                                                                                                                                ifelse(day==18,90,
                                                                                                                                                                       ifelse(day==19,95,
                                                                                                                                                                              ifelse(day==20,100,
                                                                                                                                                                                     ifelse(day==21,105,NA))))))))))))))))))))
spar24_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
        geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
        theme_classic() +
        ylab("Relative Frequency")+xlab("Time (Cycle)")
      
spar24_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
        geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
        theme_classic() + 
        ylab("Relative Frequency")+xlab("Time (Cycle)")+
        scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30","red")) 
      
data<-dplyr::filter(Extended10_LOH_parents_table_proportion, Replicate==24)

spar24_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
        geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
        scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
        ylab("Relative Frequency")+xlab("Time (Cycles)") +
        scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 
      
Fig1A2<- plot_grid(spar24_1,labels="A",label_size=10)
Fig1B2<- plot_grid(spar24_2,labels="B",label_size=10)
Fig1C2<- plot_grid(spar24_3,labels="C",label_size=10)
Fig1D2<- plot_grid(spar24_4,labels="D",label_size=10)
      
Figure7<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("S. paradoxus Line 24", fontface='bold')
Figure7<-plot_grid(title,Figure7, ncol=1, rel_heights=c(0.1, 1))
      
#Spar 27
Hybrid <- dplyr::filter(Extended_expevol, strain=="2Spar")
Hybrid_30 <- dplyr::filter(Hybrid, rep==27)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")
Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                        ifelse(day==4,20,
                                                                               ifelse(day==5,25,
                                                                                      ifelse(day==6,30,
                                                                                             ifelse(day==7,35,
                                                                                                    ifelse(day==8,40,
                                                                                                           ifelse(day==9,45,
                                                                                                                  ifelse(day==10,50,
                                                                                                                         ifelse(day==12,55,
                                                                                                                                ifelse(day==13,65,
                                                                                                                                       ifelse(day==14,70,
                                                                                                                                              ifelse(day==15,75,
                                                                                                                                                     ifelse(day==16,80,
                                                                                                                                                            ifelse(day==17,85,
                                                                                                                                                                   ifelse(day==18,90,
                                                                                                                                                                          ifelse(day==19,95,
                                                                                                                                                                                 ifelse(day==20,100,
                                                                                                                                                                                        ifelse(day==21,105,NA)))))))))))))))))))
spar27_1 <- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
          geom_line()+
          geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
          ylab("Growth rate (OD/hour)")
        
spar27_1<-spar27_1 + scale_color_grey() + theme_classic()

loh_table<-Extended10_LOH_table_parents
loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- dplyr::filter(loh_table, Genotype=="Spar")
loh_tableB <- dplyr::filter(loh_tableB, Line==27)
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                                ifelse(day==3,15,
                                                                       ifelse(day==4,20,
                                                                              ifelse(day==5,25,
                                                                                     ifelse(day==6,30,
                                                                                            ifelse(day==7,35,
                                                                                                   ifelse(day==8,40,
                                                                                                          ifelse(day==9,45,
                                                                                                                 ifelse(day==10,50,
                                                                                                                        ifelse(day==12,55,
                                                                                                                               ifelse(day==13,65,
                                                                                                                                      ifelse(day==14,70,
                                                                                                                                             ifelse(day==15,75,
                                                                                                                                                    ifelse(day==16,80,
                                                                                                                                                           ifelse(day==17,85,
                                                                                                                                                                  ifelse(day==18,90,
                                                                                                                                                                         ifelse(day==19,95,
                                                                                                                                                                                ifelse(day==20,100,
                                                                                                                                                                                       ifelse(day==21,105,NA))))))))))))))))))))
spar27_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
          geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
          theme_classic() +
          ylab("Relative Frequency")+xlab("Time (Cycle)")
        
spar27_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
          geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
          theme_classic() + 
          ylab("Relative Frequency")+xlab("Time (Cycle)")+
          scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 
        
data<-dplyr::filter(Extended10_LOH_parents_table_proportion, Replicate==27)
        
spar27_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
          geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
          scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
          ylab("Relative Frequency")+xlab("Time (Cycles)") +
          scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 
        
Fig1A2<- plot_grid(spar27_1,labels="A",label_size=10)
Fig1B2<- plot_grid(spar27_2,labels="B",label_size=10)
Fig1C2<- plot_grid(spar27_3,labels="C",label_size=10)
Fig1D2<- plot_grid(spar27_4,labels="D",label_size=10)
        
Figure8<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("S. paradoxus Line 27", fontface='bold')
Figure8<-plot_grid(title,Figure8, ncol=1, rel_heights=c(0.1, 1))
        
#Spar 29
Hybrid <- dplyr::filter(Extended_expevol, strain=="2Spar")
Hybrid_30 <- dplyr::filter(Hybrid, rep==29)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")
Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                          ifelse(day==4,20,
                                                                                 ifelse(day==5,25,
                                                                                        ifelse(day==6,30,
                                                                                               ifelse(day==7,35,
                                                                                                      ifelse(day==8,40,
                                                                                                             ifelse(day==9,45,
                                                                                                                    ifelse(day==10,50,
                                                                                                                           ifelse(day==12,55,
                                                                                                                                  ifelse(day==13,65,
                                                                                                                                         ifelse(day==14,70,
                                                                                                                                                ifelse(day==15,75,
                                                                                                                                                       ifelse(day==16,80,
                                                                                                                                                              ifelse(day==17,85,
                                                                                                                                                                     ifelse(day==18,90,
                                                                                                                                                                            ifelse(day==19,95,
                                                                                                                                                                                   ifelse(day==20,100,
                                                                                                                                                                                          ifelse(day==21,105,NA)))))))))))))))))))
spar29_1 <- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
            geom_line()+
            geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
            ylab("Growth rate (OD/hour)")
          
spar29_1<-spar29_1 + scale_color_grey() + theme_classic()

loh_table<-Extended10_LOH_table_parents
loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- dplyr::filter(loh_table, Genotype=="Spar")
loh_tableB <- dplyr::filter(loh_tableB, Line==29)
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                                  ifelse(day==3,15,
                                                                         ifelse(day==4,20,
                                                                                ifelse(day==5,25,
                                                                                       ifelse(day==6,30,
                                                                                              ifelse(day==7,35,
                                                                                                     ifelse(day==8,40,
                                                                                                            ifelse(day==9,45,
                                                                                                                   ifelse(day==10,50,
                                                                                                                          ifelse(day==12,55,
                                                                                                                                 ifelse(day==13,65,
                                                                                                                                        ifelse(day==14,70,
                                                                                                                                               ifelse(day==15,75,
                                                                                                                                                      ifelse(day==16,80,
                                                                                                                                                             ifelse(day==17,85,
                                                                                                                                                                    ifelse(day==18,90,
                                                                                                                                                                           ifelse(day==19,95,
                                                                                                                                                                                  ifelse(day==20,100,
                                                                                                                                                                                         ifelse(day==21,105,NA))))))))))))))))))))
spar29_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
            geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
            theme_classic() +
            ylab("Relative Frequency")+xlab("Time (Cycle)")

spar29_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
            geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
            theme_classic() + 
            ylab("Relative Frequency")+xlab("Time (Cycle)")+
            scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

data<-dplyr::filter(Extended10_LOH_parents_table_proportion, Replicate==29)
          
spar29_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
            geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
            scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
            ylab("Relative Frequency")+xlab("Time (Cycles)") +
            scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 
          
Fig1A2<- plot_grid(spar29_1,labels="A",label_size=10)
Fig1B2<- plot_grid(spar29_2,labels="B",label_size=10)
Fig1C2<- plot_grid(spar29_3,labels="C",label_size=10)
Fig1D2<- plot_grid(spar29_4,labels="D",label_size=10)
          
Figure9<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("S. paradoxus Line 29", fontface='bold')
Figure9<-plot_grid(title,Figure9, ncol=1, rel_heights=c(0.1, 1))
#################
############Assemble and save raw Extended Figure 10############
###Save panels 
Fig_Extended_10<- plot_grid(Figure1,Figure2,Figure3,Figure6,Figure4,Figure5,Figure7,Figure8,Figure9,nrow=3)

ggsave (plot = Fig_Extended_10, filename = "Extended_Fig10_raw_low_quality.jpg", units = "cm", device = "jpg",width = 70, height =61, dpi = 300, bg = "white")
ggsave (plot = Fig_Extended_10, filename = "Extended_Fig10_raw.png", units = "cm", device = "png",width = 70, height =61, dpi = 1000, bg = "white")
ggsave (plot = Fig_Extended_10, filename = "Extended_Fig10_raw.jpg", units = "cm", device = "jpg",width = 70, height =61, dpi = 1000, bg = "white")
ggsave (plot = Fig_Extended_10, filename = "Extended_Fig10_raw.svg", units = "cm", device = "svg",width = 70, height =61, dpi = 1000, bg = "white")
ggsave (plot = Fig_Extended_10, filename = "Extended_Fig10_raw.pdf", units = "cm", device = "pdf",width = 70, height =61, dpi = 1000, bg = "white")
#################
############Assemble and save Extended Figure 10############
#We added some aesthetic modifications so we import it
gpp <- rasterGrob(Extended10_Sanger)

Extended_Fig10<-plot_grid(gpp)

ggsave (plot = Extended_Fig10, filename = "Extended_Fig10_low_quality.jpg", units = "cm", device = "jpg",width =72, height =52, dpi = 300, bg = "white")
ggsave (plot = Extended_Fig10, filename = "Extended_Fig10.png", units = "cm", device = "png",width =72, height =52, dpi = 1000, bg = "white")
ggsave (plot = Extended_Fig10, filename = "Extended_Fig10.jpg", units = "cm", device = "jpg",width =72, height =52, dpi = 1000, bg = "white")
ggsave (plot = Extended_Fig10, filename = "Extended_Fig10.svg", units = "cm", device = "svg",width =72, height =52, dpi = 1000, bg = "white")
ggsave (plot = Extended_Fig10, filename = "Extended_Fig10.pdf", units = "cm", device = "pdf",width =72, height =52, dpi = 1000, bg = "white")
#################

############Figure Extended 11############
loh_table_day2<-Extended_loh_table_day2
loh_table_day2<-dplyr::filter(loh_table_day2, Mutation!="White")
loh_table_day2<-dplyr::filter(loh_table_day2, Mutation!="White NQO")

order_manual <- c("No mutation", "Heterozygous", "Homozygous")
    
fdata_LOH<-Extended_fdata
fdata_LOH$day <- fdata_LOH$Cycles
fdata_LOH<-fdata_LOH%<>% mutate(Generations = ifelse(day==1,5,
                                                         ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA))))))))))))))))))))
fdata_LOH <- dplyr::filter(fdata_LOH,Day==2)
fdata_LOH <- dplyr::filter(fdata_LOH,Mutation!="White")
fdata_LOH <- na.omit(fdata_LOH)
fdata_LOH_NQO <- dplyr::filter(fdata_LOH,Condition.y=="NQO_4μM")

#Data everyone until 23h
Figure_Extended_A<- fdata_LOH_NQO %>% dplyr::filter(hour<23.5)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Replicate,Mutation,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Mutation, alpha=0.4)) +
  facet_grid(Line~Generations)+
  ylab("OD \n (595 nm)")+xlab("Time (hours)")+
  theme_bw() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=7)) +
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "No mutation" = "gray30"),
                         labels = c("Heterozygous" = "Heterozygous PDR1 mutation","No mutation" = "No PDR1 mutation", "Homozygous" = "Homozygous PDR1 mutation")) +
  theme(axis.title = element_text(size=20, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=8,face = "bold"),
        axis.text.y = element_text(size=16,face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12), 
        strip.background = element_rect(fill = "white"))
    
Figure_Extended_B<- fdata_LOH_NQO %>% dplyr::filter(hour<23.5)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Replicate,Mutation,Generations,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Mutation, alpha=0.4)) +
  ylab("OD \n (595 nm)")+xlab("Time (hours)")+
  theme_bw() +theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=7)) +
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "No mutation" = "gray30"),
                         labels = c("Heterozygous" = "Heterozygous PDR1 mutation","No mutation" = "No PDR1 mutation", "Homozygous" = "Homozygous PDR1 mutation")) +
  theme(axis.title = element_text(size=20, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=16,face = "bold"),
        axis.text.y = element_text(size=16,face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  theme(strip.text.x = element_text(size = 12),    
        strip.text.y = element_text(size = 12),      
        strip.background = element_rect(fill = "white"))
    
loh_table_day2b<-dplyr::filter(loh_table_day2,Codon!="White")
Figure_Extended_C <- loh_table_day2b %>% dplyr::filter(Condition.y=="NQO_4μM")%>% 
  ggplot(aes(x=as.numeric(Generations), y=rval)) +
  geom_boxplot(outlier.shape = NA,aes(col=Mutation,fill=Mutation,alpha=0.7, group=interaction(Generations,Mutation)))+
  scale_x_continuous(breaks = as.numeric(Extended_loh_table_day2$Generations))+xlab("Time (Generations)")+
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "No mutation" = "gray30"),
                         labels = c("Heterozygous" = "Heterozygous PDR1 mutation","No mutation" = "No PDR1 mutation", "Homozygous" = "Homozygous PDR1 mutation")) +
  scale_fill_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "No mutation" = "gray30"),
                        labels = c("Heterozygous" = "Heterozygous PDR1 mutation","No mutation" = "No PDR1 mutation", "Homozygous" = "Homozygous PDR1 mutation")) +
  ylab("Growth rate \n (OD/hour)") +theme_bw() +
  facet_grid(.~Line)+
  theme(axis.title = element_text(size=20, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=16,face = "bold"),
        axis.text.y = element_text(size=16,face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  theme(legend.position="none")
    
#Disclaimer data and graph are invented to get legend purposes
FigureE_B<-Extended_loh_table_day2%>% dplyr::filter(Mutation!="White")%>%
  ggplot(aes(x=as.numeric(Cycles), y=rval, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  ylab("Relative Frequency")+xlab("Time (Generations)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30","gray30")) +
  theme(legend.position = "top")+
  theme(axis.title = element_text(size=18, face = "bold"),
        strip.text = element_text(color = "black", face = "bold",size = 19))+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.x = element_text(size=16,face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y =  element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black", face = "bold"))+
  theme(strip.text.x = element_text(size = 19))+
  theme(strip.text = element_text(size = 14))+
  theme(legend.text = element_text(size = 14),  
        legend.key.size = unit(3, "lines"),
        legend.title = element_text(size = 14, face = "bold"))

leg <- get_legend(FigureE_B)
    
#################
############Assemble and save Extended Figure 11############
FigA_label<- plot_grid(Figure_Extended_A,labels="a",label_size=32)
FigB_label<- plot_grid(Figure_Extended_B,labels="b",label_size=32)
FigC_label<- plot_grid(Figure_Extended_C,labels="c",label_size=32)
Fig_Extended_Hyb_up <-plot_grid(FigA_label,FigC_label,nrow=2)
Fig_Extended_Hyb <-plot_grid(Fig_Extended_Hyb_up,FigB_label,nrow=1,rel_widths = c(1,0.5))  
Figure11leg<-plot_grid(leg, Fig_Extended_Hyb,nrow=2,rel_heights = c(0.6,7))
    
ggsave (plot = Figure11leg, filename = "Extended_Fig11_low_quality.jpg", units = "cm", device = "jpg",width = 43, height =25, dpi = 300, bg = "white")
ggsave (plot = Figure11leg, filename = "Extended_Fig11.png", units = "cm", device = "png",width = 43, height =25, dpi = 1000, bg = "white")
ggsave (plot = Figure11leg, filename = "Extended_Fig11.jpg", units = "cm", device = "jpg",width = 43, height =25, dpi = 1000, bg = "white")
ggsave (plot = Figure11leg, filename = "Extended_Fig11.svg", units = "cm", device = "svg",width = 43, height =25, dpi = 1000, bg = "white")
ggsave (plot = Figure11leg, filename = "Extended_Fig11.pdf", units = "cm", device = "pdf",width = 43, height =25, dpi = 1000, bg = "white")
#################
