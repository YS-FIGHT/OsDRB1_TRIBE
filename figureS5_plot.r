setwd("C:/Users/Ys/Desktop/paper_loading/data and code")

library(ggplot2)
library(ggpubr)
library(ggvenn)
library(scales)
library(ggcor)

leaf_tpm<-read.table("leaf_all_feature.txt",header =T,row.names=1)
leaf_tpm <- leaf_tpm[,6:ncol(leaf_tpm)]
colnames(leaf_tpm) <- c("leaf_drb1_r1","leaf_drb1_r2","leaf_adar_r1","leaf_adar_r2","leaf_drb1adar_r1","leaf_drb1adar_r2")
p1 <- quickcor(log2(leaf_tpm+0.1), 
               type = "lower", 
               show.diag = F,
               cor.test=T,
               method='pearson',
               grid.size=0.1,
               grid.colour = "black") + 
    geom_ellipse2() + 
    scale_fill_gradientn(colours = c("#E9F2C0","#D4EBA8","#A2DD6F","#16AF00","#109F00","#0B8F00","#067F00","#036E00","#005E00","#003300")) + 
    geom_number(data = get_data(type = "lower"),aes(num = r),hjust=-0.2,vjust=3,size=2)
ggsave("figure_S5A_l.pdf",p1,width=5,height = 5)
root_tpm<-read.table("root_all_feature.txt",header =T,row.names=1)
root_tpm <- root_tpm[,6:ncol(root_tpm)]
colnames(root_tpm) <- c("root_drb1_r1","root_drb1_r2","root_adar_r1","root_adar_r2","root_drb1adar_r1","root_drb1adar_r2")
p2 <- quickcor(log2(root_tpm+0.1), 
			type = "upper", 
			show.diag = F,
			cor.test=T,
			method='pearson',
			grid.size=0.1,
			grid.colour = "black") + 
			geom_ellipse2() + 
			scale_fill_gradientn(colours = c("#B3AEA0","#B9B296","#C0B68C","#D5C96D","#C6B610","#B49C07","#8E6900","#7B5100","#552800","#421800")) + 
			geom_number(data = get_data(type = "upper"),aes(num = r),hjust=-0.2,vjust=3,size=2)
ggsave("figure_S5A_r.pdf",p2,width=5,height = 5)

#####################################################################################################

root_dp0_r1 <- read.table('root_drb1adar_against_REF_DRB1_ADAR_r1_potential_editing_dp0.vcf_efficiency',header=F)
root_dp0_r2 <- read.table('root_drb1adar_against_REF_DRB1_ADAR_r2_potential_editing_dp0.vcf_efficiency',header=F)
root_dp0_r1_pos <- root_dp0_r1[,1:3]
root_dp0_r2_pos <- root_dp0_r2[,1:3]
root_dp0_r1_pos <- unite(root_dp0_r1_pos,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
root_dp0_r2_pos <- unite(root_dp0_r2_pos,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
lt_root_dp0_pos <- list('root_dp0_r1'=root_dp0_r1_pos[,1],'root_dp0_r2'=root_dp0_r2_pos[,1])
p_root_dp0_venn <- ggvenn(lt_root_dp0_pos,c('root_dp0_r1','root_dp0_r2'),show_percentage = T,stroke_color = "black",fill_color =c('white','white'),set_name_color =c('black','black'))
ggsave('figure_S5B_r_dp0_venn.pdf',p_root_dp0_venn,width=5,height = 5)

root_dp10_r1 <- read.table('root_drb1adar_against_REF_DRB1_ADAR_r1_potential_editing_dp10.vcf_efficiency',header=F)
root_dp10_r2 <- read.table('root_drb1adar_against_REF_DRB1_ADAR_r2_potential_editing_dp10.vcf_efficiency',header=F)
root_dp10_r1_pos <- root_dp10_r1[,1:3]
root_dp10_r2_pos <- root_dp10_r2[,1:3]
root_dp10_r1_pos <- unite(root_dp10_r1_pos,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
root_dp10_r2_pos <- unite(root_dp10_r2_pos,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
lt_root_dp10_pos <- list('root_dp10_r1'=root_dp10_r1_pos[,1],'root_dp10_r2'=root_dp10_r2_pos[,1])
p_root_dp10_venn <- ggvenn(lt_root_dp10_pos,c('root_dp10_r1','root_dp10_r2'),show_percentage = T,stroke_color = "black",fill_color =c('white','white'),set_name_color =c('black','black'))
ggsave('figure_S5B_r_dp10_venn.pdf',p_root_dp10_venn,width=5,height = 5)

root_dp20_r1 <- read.table('root_drb1adar_against_REF_DRB1_ADAR_r1_potential_editing_dp20.vcf_efficiency',header=F)
root_dp20_r2 <- read.table('root_drb1adar_against_REF_DRB1_ADAR_r2_potential_editing_dp20.vcf_efficiency',header=F)
root_dp20_r1_pos <- root_dp20_r1[,1:3]
root_dp20_r2_pos <- root_dp20_r2[,1:3]
root_dp20_r1_pos <- unite(root_dp20_r1_pos,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
root_dp20_r2_pos <- unite(root_dp20_r2_pos,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
lt_root_dp20_pos <- list('root_dp20_r1'=root_dp20_r1_pos[,1],'root_dp20_r2'=root_dp20_r2_pos[,1])
p_root_dp20_venn <- ggvenn(lt_root_dp20_pos,c('root_dp20_r1','root_dp20_r2'),show_percentage = T,stroke_color = "black",fill_color =c('white','white'),set_name_color =c('black','black'))
ggsave('figure_S5B_r_dp20_venn.pdf',p_root_dp20_venn,width=5,height = 5)

leaf_dp0_r1 <- read.table('leaf_drb1adar_against_ref_drb1_adar_r1_potential_editing_dp0.vcf_efficiency',header=F)
leaf_dp0_r2 <- read.table('leaf_drb1adar_against_ref_drb1_adar_r2_potential_editing_dp0.vcf_efficiency',header=F)
leaf_dp0_r1_pos <- leaf_dp0_r1[,1:3]
leaf_dp0_r2_pos <- leaf_dp0_r2[,1:3]
leaf_dp0_r1_pos <- unite(leaf_dp0_r1_pos,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
leaf_dp0_r2_pos <- unite(leaf_dp0_r2_pos,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
lt_leaf_dp0_pos <- list('leaf_dp0_r1'=leaf_dp0_r1_pos[,1],'leaf_dp0_r2'=leaf_dp0_r2_pos[,1])
p_leaf_dp0_venn <- ggvenn(lt_leaf_dp0_pos,c('leaf_dp0_r1','leaf_dp0_r2'),show_percentage = T,stroke_color = "black",fill_color =c('white','white'),set_name_color =c('black','black'))
ggsave('figure_S5B_l_dp0_venn.pdf',p_leaf_dp0_venn,width=5,height = 5)

leaf_dp10_r1 <- read.table('leaf_drb1adar_against_ref_drb1_adar_r1_potential_editing_dp10.vcf_efficiency',header=F)
leaf_dp10_r2 <- read.table('leaf_drb1adar_against_ref_drb1_adar_r2_potential_editing_dp10.vcf_efficiency',header=F)
leaf_dp10_r1_pos <- leaf_dp10_r1[,1:3]
leaf_dp10_r2_pos <- leaf_dp10_r2[,1:3]
leaf_dp10_r1_pos <- unite(leaf_dp10_r1_pos,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
leaf_dp10_r2_pos <- unite(leaf_dp10_r2_pos,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
lt_leaf_dp10_pos <- list('leaf_dp10_r1'=leaf_dp10_r1_pos[,1],'leaf_dp10_r2'=leaf_dp10_r2_pos[,1])
p_leaf_dp10_venn <- ggvenn(lt_leaf_dp10_pos,c('leaf_dp10_r1','leaf_dp10_r2'),show_percentage = T,stroke_color = "black",fill_color =c('white','white'),set_name_color =c('black','black'))
ggsave('figure_S5B_l_dp10_venn.pdf',p_leaf_dp10_venn,width=5,height = 5)

leaf_dp20_r1 <- read.table('leaf_drb1adar_against_ref_drb1_adar_r1_potential_editing_dp20.vcf_efficiency',header=F)
leaf_dp20_r2 <- read.table('leaf_drb1adar_against_ref_drb1_adar_r2_potential_editing_dp20.vcf_efficiency',header=F)
leaf_dp20_r1_pos <- leaf_dp20_r1[,1:3]
leaf_dp20_r2_pos <- leaf_dp20_r2[,1:3]
leaf_dp20_r1_pos <- unite(leaf_dp20_r1_pos,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
leaf_dp20_r2_pos <- unite(leaf_dp20_r2_pos,"editing_pos",c("V1","V2","V3"), sep="_", remove = T)
lt_leaf_dp20_pos <- list('leaf_dp20_r1'=leaf_dp20_r1_pos[,1],'leaf_dp20_r2'=leaf_dp20_r2_pos[,1])
p_leaf_dp20_venn <- ggvenn(lt_leaf_dp20_pos,c('leaf_dp20_r1','leaf_dp20_r2'),show_percentage = T,stroke_color = "black",fill_color =c('white','white'),set_name_color =c('black','black'))
ggsave('figure_S5B_l_dp20_venn.pdf',p_leaf_dp20_venn,width=5,height = 5)


leaf_dp0_eff <- read.table('leaf_common_editing_information_dp0',header=F) 
leaf_dp0_eff <- leaf_dp0_eff[c('V8','V16')]
colnames(leaf_dp0_eff) <- c('leaf_dp0_r1','leaf_dp0_r2')
p_leaf_dp0_common_efficiency <- ggplot(data=leaf_dp0_eff, aes(x=leaf_dp0_eff$leaf_dp0_r1, y=leaf_dp0_eff$leaf_dp0_r2)) + 
										geom_point(color="#8FBC8F") + 
										stat_smooth(method="lm",se=T,color='#8FBC8F',fill='#8FBC8F') + 
										stat_cor(data=leaf_dp0_eff, method = "pearson") + 
										scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
										scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
										labs(x='leaf dp0 common editing efficiency_r1 (%)',y='leaf dp0 common editing efficiency_r2 (%)') + 
										theme_bw() + 
										theme(axis.title = element_text(size = 12,color ="black"), 
												axis.text = element_text(size= 12,color = "black"), 
												panel.grid = element_blank())
ggsave('figure_S5B_l_dp0_efficiency.pdf',p_leaf_dp0_common_efficiency,width=5,height = 5)

leaf_dp10_eff <- read.table('leaf_common_editing_information_dp10',header=F) 
leaf_dp10_eff <- leaf_dp10_eff[c('V8','V16')]
colnames(leaf_dp10_eff) <- c('leaf_dp10_r1','leaf_dp10_r2')
p_leaf_dp10_common_efficiency <- ggplot(data=leaf_dp10_eff, aes(x=leaf_dp10_eff$leaf_dp10_r1, y=leaf_dp10_eff$leaf_dp10_r2)) + 
										geom_point(color="#3CB371") + 
										stat_smooth(method="lm",se=T,color='#3CB371',fill='#3CB371') + 
										stat_cor(data=leaf_dp10_eff, method = "pearson") + 
										scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
										scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
										labs(x='leaf dp10 common editing efficiency_r1 (%)',y='leaf dp10 common editing efficiency_r2 (%)') + 
										theme_bw() + 
										theme(axis.title = element_text(size = 12,color ="black"), 
												axis.text = element_text(size= 12,color = "black"), 
												panel.grid = element_blank())
ggsave('figure_S5B_l_dp10_efficiency.pdf',p_leaf_dp10_common_efficiency,width=5,height = 5)

leaf_dp20_eff <- read.table('leaf_common_editing_information_dp20',header=F) 
leaf_dp20_eff <- leaf_dp20_eff[c('V8','V16')]
colnames(leaf_dp20_eff) <- c('leaf_dp20_r1','leaf_dp20_r2')
p_leaf_dp20_common_efficiency <- ggplot(data=leaf_dp20_eff, aes(x=leaf_dp20_eff$leaf_dp20_r1, y=leaf_dp20_eff$leaf_dp20_r2)) + 
										geom_point(color="#008000") + 
										stat_smooth(method="lm",se=T,color='#008000',fill='#008000') + 
										stat_cor(data=leaf_dp20_eff, method = "pearson") + 
										scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
										scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
										labs(x='leaf dp20 common editing efficiency_r1 (%)',y='leaf dp20 common editing efficiency_r2 (%)') + 
										theme_bw() + 
										theme(axis.title = element_text(size = 12,color ="black"), 
												axis.text = element_text(size= 12,color = "black"), 
												panel.grid = element_blank())
ggsave('figure_S5B_l_dp20_efficiency.pdf',p_leaf_dp20_common_efficiency,width=5,height = 5)


root_dp0_eff <- read.table('root_common_editing_information_dp0',header=F) 
root_dp0_eff <- root_dp0_eff[c('V8','V16')]
colnames(root_dp0_eff) <- c('root_dp0_r1','root_dp0_r2')
p_root_dp0_common_efficiency <- ggplot(data=root_dp0_eff, aes(x=root_dp0_eff$root_dp0_r1, y=root_dp0_eff$root_dp0_r2)) + 
										geom_point(color="#BDB76B") + 
										stat_smooth(method="lm",se=T,color='#BDB76B',fill='#BDB76B') + 
										stat_cor(data=root_dp0_eff, method = "pearson") + 
										scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
										scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
										labs(x='root dp0 common editing efficiency_r1 (%)',y='root dp0 common editing efficiency_r2 (%)') + 
										theme_bw() + 
										theme(axis.title = element_text(size = 12,color ="black"), 
												axis.text = element_text(size= 12,color = "black"), 
												panel.grid = element_blank())
ggsave('figure_S5B_r_dp0_efficiency.pdf',p_root_dp0_common_efficiency,width=5,height = 5)

root_dp10_eff <- read.table('root_common_editing_information_dp10',header=F) 
root_dp10_eff <- root_dp10_eff[c('V8','V16')]
colnames(root_dp10_eff) <- c('root_dp10_r1','root_dp10_r2')
p_root_dp10_common_efficiency <- ggplot(data=root_dp10_eff, aes(x=root_dp10_eff$root_dp10_r1, y=root_dp10_eff$root_dp10_r2)) + 
										geom_point(color="#808000") + 
										stat_smooth(method="lm",se=T,color='#808000',fill='#808000') + 
										stat_cor(data=root_dp10_eff, method = "pearson") + 
										scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
										scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
										labs(x='root dp10 common editing efficiency_r1 (%)',y='root dp10 common editing efficiency_r2 (%)') + 
										theme_bw() + 
										theme(axis.title = element_text(size = 12,color ="black"), 
												axis.text = element_text(size= 12,color = "black"), 
												panel.grid = element_blank())
ggsave('figure_S5B_r_dp10_efficiency.pdf',p_root_dp10_common_efficiency,width=5,height = 5)

root_dp20_eff <- read.table('root_common_editing_information_dp20',header=F) 
root_dp20_eff <- root_dp20_eff[c('V8','V16')]
colnames(root_dp20_eff) <- c('root_dp20_r1','root_dp20_r2')
	p_root_dp20_common_efficiency <- ggplot(data=root_dp20_eff, aes(x=root_dp20_eff$root_dp20_r1, y=root_dp20_eff$root_dp20_r2)) + 
											geom_point(color="#593600") + 
											stat_smooth(method="lm",se=T,color='#593600',fill='#593600') + 
											stat_cor(data=root_dp20_eff, method = "pearson") + 
											scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
											scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
											labs(x='root dp20 common editing efficiency_r1 (%)',y='root dp20 common editing efficiency_r2 (%)') + 
											theme_bw() + 
											theme(axis.title = element_text(size = 12,color ="black"), 
													axis.text = element_text(size= 12,color = "black"), 
													panel.grid = element_blank())
	ggsave('figure_S5B_r_dp20_efficiency.pdf',p_root_dp20_common_efficiency,width=5,height = 5)




p_root_dp0_venn/p_root_dp0_common_efficiency































