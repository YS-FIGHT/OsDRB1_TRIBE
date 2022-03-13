# TRIBEseq
#
###################################

library(ggplot2)
library(xlsx)
library(reshape2)
library(scales)
library(ggsignif)
library(RColorBrewer)
library(car)
library(agricolae)

setwd("C:/Users/Ys/Desktop/paper_loading/data and code")
dt <- read.xlsx('Supplemental table 1.xlsx',sheetIndex = 1,startRow = 2)

melt_dt <- melt(dt,c("site","distance_to_MS2_stemloop_nt"),variable.name = "repeats",value.name = "editing_efficiency")
melt_dt$site <- as.factor(melt_dt$site)
leveneTest(melt_dt$editing_efficiency~melt_dt$site,data = melt_dt)
oneway<-aov(melt_dt$editing_efficiency~melt_dt$site,data = melt_dt)
anova(oneway)
out <- LsD.test(oneway,"melt_dt$site",p.adj="none")
dt1 <- out$groups
dt1$site <- row.names(dt1)
dt1 <- dt1[order(dt1$site),]
dt1$yma <- out$means[,'Max']
ggplot() + 
geom_boxplot(data=melt_dt,aes(x= factor(site, level = c("A 61", "A 82", "A 146", "A 185", "A 213","A 236")),y=editing_efficiency),width = 0.5,size=0.6,outlier.alpha=0) + 
geom_jitter(data=melt_dt,aes(x= factor(site, level = c("A 61", "A 82", "A 146", "A 185", "A 213","A 236")),y=editing_efficiency,color=site),width=0.2,shape=19,alpha=0.6,size=4) + 
xlab('editing site in Ms2 transcript') + 
theme(legend.position = 'none') + 
scale_color_manual(values=c("#6495ED","#FFA500","#228B22","#FF4500","#9400D3","#828282")) + 
scale_y_continuous(labels=percent,limits = c(0,1)) + 
geom_text(data=dt1,aes(x=site,y=yma+0.1,label=groups),size=5,position= position_dodge(0.6))
ggsave("figure_1E.pdf",width = 8,height = 10)
####################################################################################################################################
ggplot(melt_dt,aes(x=distance_to_MS2_stemloop_nt,y=editing_efficiency)) + 
geom_point(aes(color=site),width=0.2,shape=19,alpha=0.6,size=4) + 
scale_color_manual(values=c("#6495ED","#FFA500","#228B22","#FF4500","#9400D3","#828282")) + 
geom_smooth(color='black') + 
scale_y_continuous(labels=percent,limits = c(0,1))
ggsave("figure_1G.pdf",width = 8,height = 10)
