#Main Figures of the Manuscript Bautista_2024

##Figure 4###

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
#Define the directory where you save all the data from Bautista_2024#################
setwd("")
#################

####1.Get data####
my_data <- read_excel("4_Table_SNPs_PDR1.xlsx")
genomics <- read_csv("4_growth.csv")
B_CRISPR_test<-read_csv("4B_CRISPR_test.csv")
C_LOH_track_parental_species <- read_csv("4C_LOH_track_parental_species2.csv")
C_LOH_track_hybrid<-read_csv("4C_LOH_track_hybrid2.csv")
E_Homozygous_Heterozygous_proportions<-read_csv("4E_Homozygous_Heterozygous_proportions.csv")
E_Homozygous_Heterozygous_proportions_medians<-read_csv("4E_Homozygous_Heterozygous_proportions_medians.csv")
F_Growth_Hybrids<-read_csv("4F_Growth_Hybrids.csv")
################

############Figure 4A############
Fig4A<- readPNG("Figure4A.png")
Fig4A <- rasterGrob(Fig4A,width = unit(1, "npc"), height = unit(1, "npc"))
Fig4A<-plot_grid(Fig4A)
##########

############Figure 4B############
#Statistics
B_CRISPR_test$Mutation <- as.factor(B_CRISPR_test$Mutation)
amod <- aov(rval~Mutation, data=B_CRISPR_test)
inter.test1 <- glht(amod,  mcp(Mutation = "Tukey"))
summary(inter.test1)
cld(inter.test1)

Fig4B<- B_CRISPR_test%>% 
  ggplot(aes(y=rval,x=Mutation)) +
  geom_boxplot(aes(alpha=0.4,fill=Mutation))+
  geom_point(colour="black",pch=21,size=2, aes(fill=Mutation))+
  facet_grid(Day~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  scale_fill_manual(values = c("gray30","mediumpurple3","mediumpurple3","gray","#21908CFF"))+
  scale_x_discrete(breaks=c("1_WTa_WTalpha","7c_M308Ia_WTalpha","3b_WTa_M308Ialpha","8c_M308Ia_STOPalpha","9d_M308Ia_M308Ialpha"),
                   limits=c("1_WTa_WTalpha","7c_M308Ia_WTalpha","3b_WTa_M308Ialpha","8c_M308Ia_STOPalpha","9d_M308Ia_M308Ialpha"))+
  ylab("Growth rate (OD/hour)")+
  theme_bw()+
  annotate("text", x = c(1,2,3,4,5), y=c(0.34,0.49, 0.5, 0.53,0.55), label = c("a","b","b","c","c"),size=4)+
  annotate("text",
           y = c(0.55),
           x = c(1),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=5) +
  theme(strip.placement = "outside",
        strip.text = element_blank(),
        axis.title = element_text(size=16, face = "bold"),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.title.y = element_text(size=15,face = "bold"),
        panel.background = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))

#We are importing the graph because it has undergone some modifications in BioRender:
Fig4B<- readPNG("Figure4B.png")
Fig4B <- rasterGrob(Fig4B,width = unit(1, "npc"), height = unit(1, "npc"))
Fig4B<-plot_grid(Fig4B)
##########

############Figure 4C############
Fig4C_A<- ggplot(C_LOH_track_parental_species, aes(x=as.numeric(Generations), y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(C_LOH_track_parental_species$Generations))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Generations)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30","gray30")) +
  facet_grid(.~ Genotype) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size=21, face = "bold"),
        axis.text.x = element_text(size=18,face = "bold"),
        axis.text.y = element_text(size=18,face = "bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.text.x = element_text(size = 20, face="bold"),
        strip.text = element_text(size = 14))

Fig4C_B<-C_LOH_track_hybrid%>% 
  ggplot(aes(x=as.numeric(Generations), y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = unique(C_LOH_track_hybrid$Generations))+
  theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Generations)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30","gray30")) +
  facet_grid(.~ Genotype) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size=21, face = "bold"),
        axis.text.x = element_text(size=20,face = "bold"),
        axis.text.y = element_text(size=20,face = "bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.text.x = element_text(size = 20, face="bold"),
        strip.text = element_text(size = 14))

Fig4C<- plot_grid(Fig4C_A, Fig4C_B,nrow=1,rel_widths = c(1,1))
##########

############Figure 4D############
Fig4D<- readPNG("Figure4D.png")
Fig4D <- rasterGrob(Fig4D,width = unit(1, "npc"), height = unit(1, "npc"))
Fig4D<-plot_grid(Fig4D)
##########

############Figure 4E############
E_Homozygous_Heterozygous_proportions_het <- dplyr::filter(E_Homozygous_Heterozygous_proportions,Mutation=="Heterozygous")
E_Homozygous_Heterozygous_proportions_hom <- dplyr::filter(E_Homozygous_Heterozygous_proportions,Mutation=="Homozygous")

#Statistics
#Heterozygous
E_Homozygous_Heterozygous_proportions_het$Genotype <- as.factor(E_Homozygous_Heterozygous_proportions_het$Genotype)
amod <- aov(median~Genotype, data=E_Homozygous_Heterozygous_proportions_het)
inter.test1 <- glht(amod,  mcp(Genotype = "Tukey"))
summary(amod)
summary(inter.test1)
cld(inter.test1)
#Homozygous
E_Homozygous_Heterozygous_proportions_hom$Genotype <- as.factor(E_Homozygous_Heterozygous_proportions_hom$Genotype)
amod <- aov(median~Genotype, data=E_Homozygous_Heterozygous_proportions_hom)
summary(amod)
inter.test1 <- glht(amod,  mcp(Genotype = "Tukey"))
summary(inter.test1)
cld(inter.test1)

Fig4E<- E_Homozygous_Heterozygous_proportions%>% 
  ggplot(aes(y=as.numeric(median),x=Genotype)) +
  scale_shape_identity() + 
  geom_point(colour="black",pch=21,size=2, aes(fill=Mutation))+
  geom_boxplot(aes(alpha=0.4,fill=Mutation))+
  geom_label(data = E_Homozygous_Heterozygous_proportions_medians, aes(x = Genotype, y =median, label =  median),fontface = "bold",
             size = 5)+
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF")) +
  scale_color_manual(values = c("mediumpurple3", "#21908CFF")) +
  facet_grid(.~Mutation)+
  ylab("Relative Frequency")+xlab("Genotype")+
  scale_x_discrete(breaks=c("Scer","Spar","Hybrid"),
                   limits=c("Scer","Spar","Hybrid"),
                   labels = expression(bold("Hybrid"), bolditalic("S. paradoxus"), bolditalic("S. cerevisiae")))+
  annotate("text",
          y = c(1.05),
          x = c(3),
          label = c("p < 0.05"),
          family = "", fontface = 3, size=5) +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(size=20, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 16,face = "bold"),
        axis.title.x =  element_blank(),
        legend.position = "none")


#We are importing the graph because it has undergone some modifications in BioRender:
Fig4E<- readPNG("Figure4E.png")
Fig4E <- rasterGrob(Fig4E,width = unit(1, "npc"), height = unit(1, "npc"))
Fig4E<-plot_grid(Fig4E)
##########

############Figure 4F############
#Statistics
hybrids_nqo<-dplyr::filter(F_Growth_Hybrids,Condition.y=="NQO_4μM")
model_anova <- aov(rval ~ Mutation, data = hybrids_nqo)
anova_resultados <- anova(model_anova)
print(anova_resultados)

tukey_resultados <- TukeyHSD(model_anova)
print(tukey_resultados)

hybrids_nqo$Mutation <- as.factor(hybrids_nqo$Mutation)
amod <- aov(rval~Mutation, data=hybrids_nqo)
summary(amod)
inter.test1 <- glht(amod,  mcp(Mutation= "Tukey"))
summary(amod)
summary(inter.test1)
cld(inter.test1)


hybrids_nqo<-dplyr::filter(F_Growth_Hybrids,Condition.y=="NQO_4μM")
hybrids_nqo<-filter(hybrids_nqo, Line==13)

hybrids_nqo$Mutation <- as.factor(hybrids_nqo$Mutation)
amod <- aov(rval~Mutation, data=hybrids_nqo)
summary(amod)
inter.test1 <- glht(amod,  mcp(Mutation= "Tukey"))
summary(amod)
summary(inter.test1)
cld(inter.test1)

hybrids_nqo<-dplyr::filter(F_Growth_Hybrids,Condition.y=="NQO_4μM")
hybrids_nqo<-filter(hybrids_nqo, Line==28)

hybrids_nqo$Mutation <- as.factor(hybrids_nqo$Mutation)
amod <- aov(rval~Mutation, data=hybrids_nqo)
summary(amod)
inter.test1 <- glht(amod,  mcp(Mutation= "Tukey"))
summary(amod)
summary(inter.test1)
cld(inter.test1)

hybrids_nqo<-dplyr::filter(F_Growth_Hybrids,Condition.y=="NQO_4μM")
hybrids_nqo<-filter(hybrids_nqo, Line==30)

hybrids_nqo$Mutation <- as.factor(hybrids_nqo$Mutation)
amod <- aov(rval~Mutation, data=hybrids_nqo)
summary(amod)
inter.test1 <- glht(amod,  mcp(Mutation= "Tukey"))
summary(amod)
summary(inter.test1)
cld(inter.test1)

Fig4F <- F_Growth_Hybrids %>% dplyr:::filter(Condition.y=="NQO_4μM")%>% 
  ggplot(aes(x=Mutation, y=rval)) +
  facet_grid(Line ~ .) +
  geom_point(colour="black",pch=21,size=2, aes(fill=Mutation))+
  geom_boxplot(aes(alpha=0.4,fill=Mutation))+
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3",
                                "Homozygous" = "#21908CFF",
                                "No mutation" = "gray30"),
                     labels = c("Heterozygous" = "Heterozygous PDR1 mutation",
                                "No mutation" = "No PDR1 mutation",
                                "Homozygous" = "Homozygous PDR1 mutation")) +
  scale_fill_manual(values = c("Heterozygous" = "mediumpurple3",
                               "Homozygous" = "#21908CFF",
                               "No mutation" = "gray30"),
                    labels = c("Heterozygous" = "Heterozygous PDR1 mutation",
                               "No mutation" = "No PDR1 mutation",
                               "Homozygous" = "Homozygous PDR1 mutation")) +
  ylab("Growth rate (OD/hour)") +theme_bw() +
  theme(legend.position = "none")+
  geom_text(data = data.frame(Mutation = 0.8, rval = 0.7, Line = c(13,28,30),
                              label = c("p < 0.0001", "p < 0.0001","p < 0.001")),
            aes(label = label,family = "", fontface = 3, size=5))+
  ylim(0,0.7)+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_text(size=16,face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        strip.text = element_text(size = 14))

#We are importing the graph because it has undergone some modifications in BioRender:
Fig4F<- readPNG("Figure4F.png")
Fig4F <- rasterGrob(Fig4F,width = unit(1, "npc"), height = unit(1, "npc"))
Fig4F<-plot_grid(Fig4F)
##########

############Create legend Figure 4############
#Disclaimer data and graph are invented to get legend purposes
Figure_leg<-C_LOH_track_hybrid%>% 
  ggplot(aes(x=as.numeric(Generations), y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(C_LOH_track_hybrid$Generations))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Generations)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30","gray30")) +
  facet_grid(.~ Genotype) +
  theme(legend.position = "top",
       legend.text = element_text(size = 30), 
       legend.key.size = unit(4, "lines"),
       legend.title = element_text(size = 40, face = "bold"))

leg <- get_legend(Figure_leg)
##########

############Assemble and save Figure 4############
Fig4A_label<- plot_grid(Fig4A,labels="a",label_size=38)
Fig4B_label<- plot_grid(Fig4B,labels="b",label_size=38)
Fig4C_label<- plot_grid(Fig4C,labels="c",label_size=38)
Fig4D_label<- plot_grid(Fig4D,labels="d",label_size=38)
Fig4E_label<- plot_grid(Fig4E,labels="e",label_size=38)
Fig4F_label<- plot_grid(Fig4F,labels="f",label_size=38)

Figure4top<-plot_grid(Fig4A_label,Fig4B_label,nrow=2,rel_heights = c(0.62,0.57))
Figure4med2<-plot_grid(Fig4E_label,Fig4F_label,nrow=1,rel_widths = c(1.2,1))
Figure4top2 <- plot_grid(Fig4D_label,Figure4med2,nrow=2,rel_heights = c(0.55,0.45))

Figure4<-plot_grid(Figure4top,Figure4top2,nrow=1,rel_widths = c(1,1.3))
Figure4c<-plot_grid(Figure4,Fig4C_label,nrow=2,rel_heights = c(1,0.35))
Figure4leg<-plot_grid(leg,Figure4c,nrow=2,rel_heights = c(0.4,7))

#Save the image in the previously set working directory
ggsave (plot = Figure4leg, filename = "Bautista2024_Figure4_low_quality.jpg", units = "cm", device = "jpg",width =65, height =55, dpi = 300)
ggsave (plot = Figure4leg, filename = "Bautista2024_Figure4.jpg", units = "cm", device = "jpg",width =65, height =55, dpi = 1000)
ggsave (plot = Figure4leg, filename = "Bautista2024_Figure4.png", units = "cm", device = "png",width =65, height =55, dpi = 1000, bg = "white")
ggsave (plot = Figure4leg, filename = "Bautista2024_Figure4.svg", units = "cm", device = "svg",width =65, height =55, dpi = 800)
ggsave (plot = Figure4leg, filename = "Bautista2024_Figure4.pdf", units = "cm", device = "pdf",width =65, height =55, dpi = 1000, bg = "white")
##########
