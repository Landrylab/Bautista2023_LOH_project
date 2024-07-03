#Main Figures of the Manuscript Bautista_2024

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
#install.packages("dunn.test")
library(dunn.test)
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
#Define the directory where you save all the data from Bautista_2024
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
img <-image_read_pdf("Fig1A.pdf")
genomics <- read_csv("1_growth.csv")
expevol <- read_csv("1_expevol.csv")
expevol_all <- read_csv("1_expevol_all_gen_rval.csv")
#################

############Figure 1A############
gpp <- rasterGrob(img, interpolate=TRUE)
Fig1A<-plot_grid(gpp)
#################

##########Figure 1B############

#Function to have italic font in legend
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

#Change column time by generations
#fdataAREA <- expevol %>% mutate(gen= day*5)

pd <- position_dodge(0.1)

#Fig 2B in cm
df <- expevol_all

#Asymtotic
# define drm function to use with map
drm.func <- function(x) {
  drm(rval ~ gen, 
      fct = LL2.2(),#names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), 
      data = x)
}

predict.fun <- function(x) {
  add_predictions(data.frame(gen = seq(15,105)), x)
}

coefs.fun <- function(x) {coef(x) %>% tidy}

df2 <- df %>% group_by(strain) %>% nest() %>%
  mutate(drmod = map(data, drm.func), 
         pred = map(drmod, predict.fun),
         coefs = map(drmod, coefs.fun))

Fig2C2<- df2 %>% unnest(data) %>% 
  ggplot() + 
  geom_line(aes(gen, pred, color = strain), data =
              df2 %>% unnest(pred), linetype = "solid", size=0.8) +
  geom_vline(aes(xintercept = log(x), color = strain),
             data = df2 %>% unnest(coefs) %>% filter(names == "ED50:(Intercept)")) +
  theme_bw() +ylim(0,0.8) +xlim(15,130) +
  theme(legend.position = "NONE", axis.title.x = element_text(size=10), axis.text.x = element_text(size = 9),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Time (generations)") + ylab("Growth rate (OD/hour)") +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + 
  scale_x_continuous(breaks = c(25,50,75,100)) +
  guides(color = guide_legend(title = "Genotype",
                              override.aes=list(fill=NA)), 
         linetype = guide_legend(title = "UV mimetic", 
                                 override.aes=list(fill=NA)),
         shape = guide_legend(title = "UV mimetic", override.aes=list(colour="black"))) +
  theme(legend.text=element_text(size=4), legend.title=element_text(size=5), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4)) +
  theme_prism()+
  theme_bw(base_size=24) +
  theme(legend.position = "none",
        axis.title = element_text(size=11, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust = 0.5, size = 12),
        panel.background = element_blank()) +theme(strip.text.x = element_blank())

#Filter by the growth in 4-NQO (UV mimetic)
genomics1 <- dplyr:::filter(expevol,NQO=="Yes")
genomics1 <- dplyr:::filter(genomics1,day==3)
genomics1$rep<-as.numeric(genomics1$rep)
genomics2<-genomics1 

genomics1 <- dplyr:::filter(expevol,NQO=="Yes")
genomics1 <- dplyr:::filter(genomics1,day==21)
genomics1$rep<-as.numeric(genomics1$rep)
genomics3<-genomics1

genomics4<- left_join(genomics2,genomics3, by=c("strain"="strain","rep"="rep"))

genomics4$fitness_gain<- (genomics4$rval.y - genomics4$rval.x) / genomics4$rval.x
genomics4$percentage <- genomics4$fitness_gain * 100

kruskal.test(percentage ~ strain, data = genomics4)
dunn.test(genomics4$percentage, genomics4$strain)

Fig1B <- genomics4 %>% 
  ggplot(aes(x=strain,y=percentage, col=strain, group=strain))+
  geom_boxplot(outlier.shape=NA,aes(fill=strain),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=strain), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotypic background", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("% increase in growth rate") + 
  scale_color_manual(values=c("green4", "dodgerblue1", "black"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  geom_segment(aes(x=1, xend=2, y=365, yend=365)) + 
  geom_segment(aes(x=2, xend=3, y=390, yend=390)) + 
  geom_segment(aes(x=1, xend=3, y=340, yend=340)) +
  annotate("text", x = 1, y = 400, size = 3,
           label = c("p < 0.0001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(375, 400, 350),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.0001", "p < 0.0001", "p < 0.05"),
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
        panel.background = element_blank()) +theme(strip.text.x = element_blank())+
  ylim(-30,400)

Fig1B<- Fig2C2
##########

############Figure 1C############
#try correlation
genomics1 <- dplyr:::filter(expevol,NQO=="Yes")
genomics1 <- dplyr:::filter(genomics1,day==3)
genomics1$rep<-as.numeric(genomics1$rep)
genomics2<-genomics1 

genomics1 <- dplyr:::filter(expevol,NQO=="Yes")
genomics1 <- dplyr:::filter(genomics1,day==21)
genomics1$rep<-as.numeric(genomics1$rep)
genomics3<-genomics1

genomics4<- left_join(genomics2,genomics3, by=c("strain"="strain","rep"="rep"))

genomics4$fitness_gain<- (genomics4$rval.y - genomics4$rval.x) / genomics4$rval.x
genomics4$percentage <- genomics4$fitness_gain * 100

genomics5<-genomics4

#Filter by the growth in 4-NQO (UV mimetic)
AV<- dplyr:::filter(genomics, condition=="NQO")
AV<- dplyr:::filter(AV, Evolved=="Evolved_NQO")
AV<- dplyr:::filter(AV, day==2)

#genomics <- AV
genomics1 <- dplyr:::filter(genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics2<-genomics1 %>% dplyr:::filter(Evolved=="Evolved_NQO")
genomics1 <- dplyr:::filter(genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics3<-genomics1 %>% dplyr:::filter(Evolved=="Ancestor")

genomics4<- left_join(genomics2,genomics3, by=c("Specie"="Specie","condition"="condition","Replicate"="Replicate"))

genomics4$fitness_gain<- (genomics4$rval.x - genomics4$rval.y) / genomics4$rval.y
genomics4$percentage <- genomics4$fitness_gain * 100

joined <- full_join(genomics5, genomics4, by=c("strain"="Specie","rep"="Replicate"))

legend_title <- " "

#install.packages("ggpmisc")
library(ggpmisc)

Fig1C <- joined %>% 
  ggplot(aes(x=percentage.x,y = percentage.y))+
  border()  +
  geom_point(pch=21, aes(fill=strain), col="black", size=2, alpha=0.5)+
  scale_fill_manual(legend_title,values=c("green4","dodgerblue1","#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus","Hybrid"))) +
  stat_smooth(method="lm", alpha=0.2,aes(fill = strain)) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               label.x = "right",
               formula = y ~ x,
               parse = TRUE,
               size = 3) +
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
  xlab("Growth rate (OD/hour) of entire populations") + 
  ylab("Growth rate (OD/hour) of isolated clones") +
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
  xlab("Growth rate (OD/hour) of populations") + 
  ylab("Growth rate (OD/hour) of isolated clones") +
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
ggsave (plot = Fig1_img, filename = "Bautista2024_Figure1_low_quality.jpg", units = "cm", device = "jpg",width =28, height =25, dpi = 300,bg = "white")
ggsave (plot = Fig1_img, filename = "Bautista2024_Figure1.png", units = "cm", device = "png",width =28, height =25, dpi = 1000,bg = "white")
ggsave (plot = Fig1_img, filename = "Bautista2024_Figure1.jpg", units = "cm", device = "jpg",width =28, height =25, dpi = 1000,bg = "white")
ggsave (plot = Fig1_img, filename = "Bautista2024_Figure1.svg", units = "cm", device = "svg",width =28, height =25, dpi = 1000,bg = "white")
ggsave (plot = Fig1_img, filename = "Bautista2024_Figure1.pdf", units = "cm", device = "pdf",width =28, height =25, dpi = 1000,bg = "white")
#################
