setwd("C:/Users/Ys/Desktop/paper_loading/data and code")

library(circlize)
library(ComplexHeatmap)
library(ggprism)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(ggvenn)
options(scipen = 200)
db1ar_potential_r1 <- read.table('DB1AR_r1_raw.fq.gz_fpout.gz_cut3_trimout._18-34_genome_sorted_pileup_pile2base_no_strand_editing_MIRNA_filterDepthEfficiency.txt',header = F)
db1ar_potential_r2 <- read.table('DB1AR_r2_raw.fq.gz_fpout.gz_cut3_trimout._18-34_genome_sorted_pileup_pile2base_no_strand_editing_MIRNA_filterDepthEfficiency.txt',header = F)
adar_potential_r1 <- read.table('ADAR_r1_raw.fq.gz_fpout.gz_cut3_trimout._18-34_genome.sorted_pileup_pile2base_no_strand_editing_MIRNA_filterDepthEfficiency.txt',header = F)
adar_potential_r2 <- read.table('ADAR_r2_raw.fq.gz_fpout.gz_cut3_trimout._18-34_genome.sorted_pileup_pile2base_no_strand_editing_MIRNA_filterDepthEfficiency.txt',header = F)

db1ar <- c(nrow(db1ar_potential_r1),nrow(db1ar_potential_r2))
ar <- c(nrow(adar_potential_r1),nrow(adar_potential_r2)) 
mydata <- data.frame( 
    group = rep(c("OsDRB1-ADARdd", "ADARdd"), each = 2),
    number = c(db1ar,  ar)
)
ggbarplot(mydata,
          'group',
          'number',
          color = 'group',
          fill = 'group',
          palette = c('#93a1b2','#b2b08a'),
          add = "mean_sd",,xlab = F,ylab = 'Amount of editing sites',legend='none',
          ggtheme = theme_classic())+ 
    stat_compare_means(aes(label = ..p.signif..),
                       comparisons = list(c('OsDRB1-ADARdd','ADARdd')),
                       method = 't.test', tip.length=0,vjust=0,bracket.size=1)+
    scale_y_continuous(breaks = c(0,10,20,30,40,50)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
          legend.position = 'none',
          panel.grid=element_blank(),
          axis.text = element_text(size= 12,color = "black"),
          axis.title = element_text(size = 12,color ="black"))
ggsave("figure_4A.pdf",width=5,height=5)	
com_srna <- merge(srna_r1,srna_r2,by='editing_pos')
com_srna_eff <- com_srna[,c(16,31)]
colnames(com_srna_eff) <- c('r1','r2')
ggplot(data=com_srna_eff, aes(x=r1, y=r2)) + 
    geom_point(size=10,colour='#a38ebf',alpha=0.7) + 
    stat_smooth(method="lm",se=T,color='#b0dae0',fill='#daf47d') + 
    stat_cor(data=com_srna_eff, method = "pearson") + 
    scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1),labels = c("0","20","40","60","80","100")) + 
    scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1),labels = c("0","20","40","60","80","100")) + 
    labs(x='Common editing efficiency_r1 (%)',y='Common editing efficiency_r2 (%)') + 
    theme_classic() + 
    theme(axis.title = element_text(size = 12,color ="black"), 
          axis.text = element_text(size= 12,color = "black"), 
          panel.grid = element_blank())
ggsave("figure_4C.pdf",width=5,height=5)

com_editing <- read.table('common_editing_MIRNA.txt',header=F)
com_editing_mir <- com_editing[,10]
com_editing_mir_freq <- as.data.frame(table(com_editing_mir))
com_editing_mir_freq <-com_editing_mir_freq$Freq
com_editing_mir_freq <- as.data.frame(table(com_editing_mir))
com_editing_mir_count <-com_editing_mir_freq$Freq
com_editing_mir_count <- as.data.frame(table(com_editing_mir_count))
colnames(com_editing_mir_count) <- c('editing_counts','number_of_miRNAs')
ggplot(com_editing_mir_count,aes(x=editing_counts,y=number_of_miRNAs)) + 
    geom_bar(stat='identity',position = position_dodge2(0.5),width=0.8,fill='#c6dbbf') + 
    theme_classic() +
    theme(axis.text = element_text(size= 12,color = "black"),
          panel.grid=element_blank(),
          axis.title = element_text(size = 12,color ="black")) + 
    labs(x='Editing counts',y='Number of miRNAs')
ggsave("figure_4D.pdf",width=5,height=5)

srna_r1 <- read.table('DB1AR_r1_raw.fq.gz_fpout.gz_cut3_trimout._18-34_genome_sorted_pileup_pile2base_no_strand_editing_MIRNA_filterDepthEfficiency_against_adar.txt',header=F)
srna_r2 <- read.table('DB1AR_r2_raw.fq.gz_fpout.gz_cut3_trimout._18-34_genome_sorted_pileup_pile2base_no_strand_editing_MIRNA_filterDepthEfficiency_against_adar.txt',header=F)
srna_r1 <- unite(srna_r1,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
srna_r2 <- unite(srna_r2,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
srna_editing_sites <- list('OsDRB1-ADARdd_r1'=srna_r1[,1],'OsDRB1-ADARdd_r2'=srna_r2[,1])
ggvenn(srna_editing_sites,c('OsDRB1-ADARdd_r1','OsDRB1-ADARdd_r2'),
       show_percentage = T,
       fill_color =c('#d1a669','#99bdb8'),
       stroke_size = 2,
       set_name_color =c('black','black'),
       stroke_color = "white",fill_alpha = 0.7)
ggsave("figure_4B.pdf",width=5,height=5)
srna_editing_mir <- list('OsDRB1-ADARdd_r1'=srna_r1[,15],'OsDRB1-ADARdd_r2'=srna_r2[,15])
ggvenn(srna_editing,c('OsDRB1-ADARdd_r1','OsDRB1-ADARdd_r2'),
       show_percentage = T,
       fill_color =c('#88a5a7','#cfad9d'),
       stroke_size = 2,
       set_name_color =c('black','black'),
       stroke_color = "white",fill_alpha = 0.7)
ggsave("figure_4D_2.pdf",width=5,height=5)

pdf("figure_4E.pdf.pdf")
db1ar_ex1 <- read.table('DB1AR_r1_against_hairpin_count.txt',header=F)
total_reads_db1ar_ex1 <- sum(db1ar_ex1$V1)
db1ar_ex1$db1ar_ex1_cpm <- db1ar_ex1$V1/total_reads_db1ar_ex1*1000000
db1ar_ex2 <- read.table('DB1AR_r2_against_hairpin_count.txt',header=F)
total_reads_db1ar_ex2 <- sum(db1ar_ex2$V1)
db1ar_ex2$db1ar_ex2_cpm <- db1ar_ex2$V1/total_reads_db1ar_ex2*1000000

adar_ex1 <- read.table('ADAR_r1_against_hairpin_count.txt',header=F)
total_reads_adar_ex1 <- sum(adar_ex1$V1)
adar_ex1$adar_ex1_cpm <- adar_ex1$V1/total_reads_adar_ex1*1000000
adar_ex2 <- read.table('ADAR_r2_against_hairpin_count.txt',header=F)
total_reads_adar_ex2 <- sum(adar_ex2$V1)
adar_ex2$adar_ex2_cpm <- adar_ex2$V1/total_reads_adar_ex2*1000000

final_mir <- read.table('final_mirna_unique.txt',header=F)

db1ar_final_mir_count_r1 <- read.table('DB1AR_r1_common_editing_unique_mir_count',header = F)
db1ar_final_mir_count_r2 <- read.table('DB1AR_r2_common_editing_unique_mir_count',header = F)
mir_editing_counts <- merge(db1ar_final_mir_count_r1,db1ar_final_mir_count_r2,by='V1')
mir_editing_counts$adar_r1 <- rep(0,nrow(mir_editing_counts))
mir_editing_counts$adar_r2 <- rep(0,nrow(mir_editing_counts))
colnames(mir_editing_counts) <- c('mir_id','db1ar_r1','db1ar_r2','adar_r1','adar_r2')
rownames(mir_editing_counts) <- mir_editing_counts$mir_id
mir_editing_counts <- mir_editing_counts[order(mir_editing_counts$mir_id),]

final_mir_id <- final_mir[,1]
final_mir_db1ar_ex1 <- db1ar_ex1[db1ar_ex1$V2 %in% final_mir_id,]
final_mir_db1ar_ex2 <- db1ar_ex2[db1ar_ex2$V2 %in% final_mir_id,]
final_mir_adar_ex1 <- adar_ex1[adar_ex1$V2 %in% final_mir_id,]
final_mir_adar_ex2 <- adar_ex2[adar_ex2$V2 %in% final_mir_id,]
final_mir_exp <- cbind(final_mir_db1ar_ex1$V2,final_mir_db1ar_ex1$db1ar_ex1_cpm,final_mir_db1ar_ex2$db1ar_ex2_cpm,final_mir_adar_ex1$adar_ex1_cpm,final_mir_adar_ex2$adar_ex2_cpm)
final_mir_exp <- as.data.frame(final_mir_exp)
colnames(final_mir_exp) <- c('mir_id','mir_db1ar_ex1','mir_db1ar_ex2','mir_adar_ex1','mir_adar_ex2')
rownames_final_mir_exp <- final_mir_exp$mir_id


final_mir_exp <- final_mir_exp[,-1]

final_mir_exp <- as.data.frame(lapply(final_mir_exp,as.numeric))
rownames(final_mir_exp) <- rownames_final_mir_exp
exp <- apply(final_mir_exp, 1, scale)
rownames(exp) <- colnames(final_mir_exp)
exp <- t(exp)
exp <- as.data.frame(exp)
exp$mir_id <- row.names(exp)
exp <- exp[order(exp$mir_id),]
exp <- merge(exp,mir_editing_counts,by='mir_id')
rowname_exp <- exp$mir_id
exp <- exp[,-1]
exp <- as.data.frame(lapply(exp,as.numeric))
rownames(exp) <- rowname_exp

col_fun <- colorRamp2(c(-2, 0, 2),c("#ff948d", "white", "#4382a6"))
ht_srna <- Heatmap(exp[,seq(1,4)], name = "expression", col = col_fun,rect_gp = gpar(col='white',lwd=0.5),
			cell_fun = function(j, i, x, y, width, height, fill) { if(exp[i, j+4] > 0)
				grid.text(exp[i,j+4], x, y, gp = gpar(fontsize = 10))},
				heatmap_legend_param = list(
                                      legend_direction = "horizontal"
                                     ))
draw(ht_srna, heatmap_legend_side = "bottom")
dev.off()
