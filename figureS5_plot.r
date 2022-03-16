# TRIBEseq
#
###########################################################################################
library(ggplot2)
library(scales)
library(ggpubr)
library(reshape2)
library(ggthemes)
library(patchwork)
setwd("C:/Users/Ys/Desktop/paper_loading/data and code")

###########################################################################################

#copy "uniq_map_summary.txt" file into getwd()

f1 <- read.table("uniq_map_summary.txt",sep="\t",header=F)
f1 <- f1[order(f1$V1),]
f1$tissues <- c(rep("leaf",6),rep("root",6))
f1$replicates <- c(rep(c("r1","r2"),6))
ggplot(f1,aes(x=V1,y=V3,label=replicates)) + 
		geom_bar(width = 0.6, position =position_dodge2(0.8),stat='identity',aes(fill=tissues)) + 
		ylim(0,100)+coord_flip()+scale_fill_manual(values=c(leaf='#84bd00',root='#c68143')) + 
		geom_text(position = position_dodge2(0.8),vjust = 0.5,hjust=-0.5) + 
		labs(x=NULL,y="mapping percentage (%)") + 
		theme_pander() + 
		theme(axis.line.x = element_line(colour = "black"),legend.title =element_blank(),
		axis.title = element_text(size = 12,color ="black"),
		axis.text = element_text(size= 12,color = "black"))
ggsave("figure_S4C.pdf",width=5,height=5)

###########################################################################################
# reference:https://www.jianshu.com/p/eece90bdddf9

countdata<-read.table("leaf_all_feature.txt",skip = 1,sep="\t",header = T,row.names = 1)
metadata <- countdata[,1:5]
countdata <- countdata[,6:ncol(countdata)]
prefix<-"couts"
#-----TPM Calculation------
options(scipen = 200)
kb <- metadata$Length / 1000
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
tpm <- as.data.frame(tpm)
colnames(tpm) <- c("OsDRB1-OE_r1",	"OsDRB1-OE_r2",	"ADARdd-OE_r1",	"ADARdd-OE_r2",	"OsDRB1-ADARdd-OE_r1",	"OsDRB1-ADARdd-OE_r2")
leaf_adar <- subset(tpm,select=c("ADARdd-OE_r1",	"ADARdd-OE_r2",	"OsDRB1-ADARdd-OE_r1",	"OsDRB1-ADARdd-OE_r2"))["ADAR",]

countdata<-read.table("root_all_feature.txt",skip = 1,sep="\t",header = T,row.names = 1)
metadata <- countdata[,1:5]
countdata <- countdata[,6:ncol(countdata)]
prefix<-"couts"
options(scipen = 200)
kb <- metadata$Length / 1000
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
tpm <- as.data.frame(tpm)
colnames(tpm) <- c("OsDRB1-OE_r1",	"OsDRB1-OE_r2",	"ADARdd-OE_r1",	"ADARdd-OE_r2",	"OsDRB1-ADARdd-OE_r1",	"OsDRB1-ADARdd-OE_r2")
root_adar <- subset(tpm,select=c("ADARdd-OE_r1",	"ADARdd-OE_r2",	"OsDRB1-ADARdd-OE_r1",	"OsDRB1-ADARdd-OE_r2"))["ADAR",]

rownames(leaf_adar) <- "leaf_ADAR"
rownames(root_adar) <- "root_ADAR"
tpm_adar <- rbind(leaf_adar,root_adar)
tpm_adar <- t(tpm_adar)
tpm_adar <- as.data.frame(tpm_adar)
tpm_adar$genotypes <- row.names(tpm_adar)
tpm_adar <- melt(tpm_adar,"genotypes",variable.name = "tissue",value.name = "tpm_value")
tpm_adar$replicates <- c("L1","L1","L2","L2","L3","L3","L4","L4")
ggplot(tpm_adar,aes(x=genotypes,y=tpm_value)) + 
		geom_bar(stat="identity",pos="dodge",aes(fill=replicates)) + 
		facet_wrap(~tissue)+theme_pander() + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1 ),legend.position = 'none', 
		axis.title = element_text(size = 12,color ="black"),axis.text = element_text(size= 12,color = "black")) + 
		labs(y="tpm_ADAR",x=NULL) + 
		scale_fill_manual(values = (c(L1="#67B83C",L2="#9BC82F",L3="#7F4F21",L4="#C9A063")))
ggsave("figure_S4D.pdf",width=5,height=5)

###########################################################################################
#copy "leaf_drb1adar_against_REF_DRB1_ADAR_mutations_retained_vs_removed.txt" file into getwd()
f1 <- read.table('leaf_drb1adar_against_REF_DRB1_ADAR_mutations_retained_vs_removed.txt',header=T,sep='\t')
f1_m <- melt(f1,"type_of_mutation",variable.name = "replicates",value.name = "counts")
f1_m$states <- c(rep("retained",12),rep("removed",12),rep("retained",12),rep("removed",12))
f1_m$replicates <- c(rep("r1",24),rep('r2',24))

f2 <- read.table('root_drb1adar_against_REF_DRB1_ADAR_mutations_retained_vs_removed.txt',header=T,sep='\t')
f2_m <- melt(f2,"type_of_mutation",variable.name = "replicates",value.name = "counts")
f2_m$states <- c(rep("retained",12),rep("removed",12),rep("retained",12),rep("removed",12))
f2_m$replicates <- c(rep("r1",24),rep('r2',24))

p1 <- ggplot(f1_m,aes(x=replicates,y=counts)) + 
		geom_bar(stat="identity",position='fill',aes(fill=states)) + 
		scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels = c("0","25","50","75","100")) + 
		facet_wrap(~type_of_mutation) + 
		theme_bw() + 
		scale_fill_manual(values = c(removed='#99aab5',retained='#0f723a')) + 
		labs(x="OsDRB1-ADARdd-OE vs. OsDRB1-OE & ADARdd-OE",y="percentage of mutations (%)") + 
		theme(legend.title = element_blank(),
		axis.title = element_text(size = 12,color ="black"),
		axis.text = element_text(size= 12,color = "black"),
		panel.grid=element_blank(),
		legend.position = "top",
		legend.text = element_text(size= 12))

p2 <- ggplot(f2_m,aes(x=replicates,y=counts)) + 
		geom_bar(stat="identity",position='fill',aes(fill=states)) + 
		scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels = c("0","25","50","75","100")) + 
		facet_wrap(~type_of_mutation) + 
		theme_bw() + 
		scale_fill_manual(values = c(removed='#b2a895',retained='#c98a30')) + 
		labs(x="OsDRB1-ADARdd-OE vs. OsDRB1-OE & ADARdd-OE",y="percentage of mutations (%)") + 
		theme(legend.title = element_blank(),
		axis.title = element_text(size = 12,color ="black"),
		axis.text = element_text(size= 12,color = "black"),
		panel.grid=element_blank(),
		legend.position = "top",
		legend.text = element_text(size= 12))
pout <- (p1+p2)+plot_annotation(tag_levels = list(c('SNPs in leaf', 'SNPs in root')))
ggsave("figure_S4E.pdf",width=8,height=5)
