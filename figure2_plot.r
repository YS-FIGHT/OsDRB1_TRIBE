# TRIBEseq
#
###################################
library(ggplot2)
library(scales)
library(ggpubr)
library(reshape2)
library(ggsignif)
library(ggthemes)
library(agricolae)
library(car)
library(MASS)
library(viridis)
library(ggseqlogo) 
library(patchwork)
setwd("C:/Users/Ys/Desktop/paper_loading/data and code")

#################################################################################################

f1 <- read.table('leaf_replicates_nucleotide_mutation_summary.txt',header=T,sep='\t')
f2 <- read.table('root_replicates_nucleotide_mutation_summary.txt',header=T,sep='\t')
f1 <- f1[,c(1,14,16,2,4,10,12,6,8)]
f2 <- f2[,c(1,14,16,2,4,10,12,6,8)]
#data_a <- cbind(f1,f2[,-1])
data <- cbind(f1[,c(seq(1,7,1))],f2[,seq(2,7,1)])
# f1_m <- melt(f1,"type_of_mutation")
# f2_m <- melt(f2,"type_of_mutation")
colnames(data) <- c('mutation_type','drb1_l_1','drb1_l_2','adar_l_1','adar_l_2','drb1adar_l_1','drb1adar_l_2','drb1_r_1','drb1_r_2','adar_r_1','adar_r_2','drb1adar_r_1','drb1adar_r_2')
rownames(data) <- data[,1]
data <- data[,-1]
data_rownames <- rownames(data)
data_colnames <- colnames(data)
data$mutation <- data_rownames
data_m <- melt(data, id.vars=c("mutation"))
data_t <- t(data)
data_t <- t(data[,-ncol(data)])
data_t <- as.data.frame(data_t)
data_t$total_mutation <- apply(data_t,1,sum)
data_t$replicates <- row.names(data_t)
p <- ggplot() + 
    geom_bar(data=data_m, aes(x=factor(variable,level=c('adar_l_2','adar_l_1','adar_r_2','adar_r_1','drb1_l_2','drb1_l_1','drb1_r_2','drb1_r_1','drb1adar_l_2','drb1adar_l_1','drb1adar_r_2','drb1adar_r_1')),
				y=value,
				fill=factor(mutation,level=c("A>G","T>C","A>C","A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A","T>G")))
							,stat="identity",width = 0.7 ) + 
    scale_fill_manual(values=c('A>G'='#73c8dc', 'T>C'='#fcd72e', 'A>C'='#5A7FC1', 'A>T'='#C3DF0A', 'C>A'='#DEB982','C>G'='#f2cce1', 'C>T'='#f6ac5d','G>A'='#A3DCEC', 'G>C'='#6C9C44', 'G>T'='#CEB5BD', 'T>A'='#6C6C74', 'T>G'='#ec6d5c')) + 
    scale_y_continuous(limits=c(0, 11500))+
    coord_flip() + 
    labs(x="Number of mutations",y=NULL)+theme_pander() + 
    geom_text(data=data_t,aes(x=replicates,y=total_mutation,label=total_mutation),position = position_dodge2(0.8),vjust = 0.5,hjust=-0.05) + 
    theme(legend.title =element_blank(),
			axis.line.x = element_line(colour = 'black'), 
			axis.title = element_text(size = 12,color ="black"),
			axis.text = element_text(size= 12,color = "black"))
ggsave("figure_2B.pdf",width=6,height=5)

#################################################################################################

library(agricolae)
library(car)
f2c_d1 <- data[,c(5,6)]
f2c_d2 <- data[,c(11,12)]
f2c_d1$mean <- apply(f2c_d1,2,mean)
f2c_d1$mean <- apply(f2c_d1,1,mean)
f2c_d1 <- data[,c(5,6)]
f2c_d1$average <- apply(f2c_d1,1,mean)
f2c_d2 <- data[,c(11,12)]
f2c_d2$average <- apply(f2c_d2,1,mean)
f2c_d1$sd <- apply(f2c_d1[,c(1,2)],1,sd)
f2c_d2$sd <- apply(f2c_d2[,c(1,2)],1,sd)
f2c_d1$mutation <- rownames(f2c_d1)
f2c_d2$mutation <- rownames(f2c_d2)
f2c_d1_m <- melt(f2c_d1,c("mutation","average","sd"))
f2c_d2_m <- melt(f2c_d2,c("mutation","average","sd"))
one.way1 <- aov(value~mutation, data = f2c_d1_m)
out1 <- LSD.test(one.way1,"mutation",p.adj="none")
one.way2 <- aov(value~mutation, data = f2c_d2_m)
out2 <- LSD.test(one.way2,"mutation",p.adj="none")
d1 <- out1$groups
d1$mutation <- rownames(d1)
d1 <- d1[order(d1$mutation),]
f2c_d1 <- f2c_d1[order(f2c_d1$mutation),]
d1$sd <- f2c_d1$sd
d1$yma <- out1$means[,"Max"]
d2 <- out2$groups
d2$mutation <- rownames(d2)
d2 <- d2[order(d2$mutation),]
f2c_d2 <- f2c_d2[order(f2c_d2$mutation),]
d2$sd <- f2c_d2$sd
d2$yma <- out2$means[,"Max"]
p1 <- ggplot(d1,aes(x=mutation,y=value,fill=mutation)) + 
		geom_bar(stat="identity",position = position_dodge2(0.8)) + 
		geom_errorbar(aes(ymax=value+sd,ymin=value-sd),position=position_dodge(0.9),width=0.2,size=1.2,colour="#4C5454",alpha=0.7) + 
		scale_fill_manual(values=c('A>G'='#73c8dc', 'T>C'='#fcd72e', 'A>C'='#5A7FC1', 'A>T'='#C3DF0A', 'C>A'='#DEB982','C>G'='#f2cce1', 'C>T'='#f6ac5d','G>A'='#A3DCEC', 'G>C'='#6C9C44', 'G>T'='#CEB5BD', 'T>A'='#6C6C74', 'T>G'='#ec6d5c')) + 
		theme_bw() + 
		geom_text(aes(y=yma+400,label=groups),size=5) + 
		theme(axis.text.x=element_blank(),
				legend.position = 'none',
				panel.grid=element_blank(),
				axis.text.y = element_text(size= 12,color = "black")) + 
		labs(x=NULL,y=NULL)

p2 <- ggplot(d2,aes(x=mutation,y=value,fill=mutation)) + 
		geom_bar(stat="identity",position = position_dodge2(0.8)) + 
		geom_errorbar(aes(ymax=value+sd,ymin=value-sd),position=position_dodge(0.9),width=0.2,size=1.2,colour="#4C5454",alpha=0.7) + 
		scale_fill_manual(values=c('A>G'='#73c8dc', 'T>C'='#fcd72e', 'A>C'='#5A7FC1', 'A>T'='#C3DF0A', 'C>A'='#DEB982','C>G'='#f2cce1', 'C>T'='#f6ac5d','G>A'='#A3DCEC', 'G>C'='#6C9C44', 'G>T'='#CEB5BD', 'T>A'='#6C6C74', 'T>G'='#ec6d5c')) + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1 ),
				legend.position = 'none',
				panel.grid=element_blank(),
				axis.text = element_text(size= 12,color = "black")) + 
		geom_text(aes(y=yma+200,label=groups),size=5) + 
		labs(x=NULL,y=NULL)
pout <- (p1/p2)+plot_annotation(tag_levels = list(c('leaf', 'root')))	
ggsave("figure_2C.pdf",width=6,height=5)

#####################################################################################################

f2d1 <- f1[,c(1,8,9)]
f2d2 <- f2[,c(1,8,9)]
rownames(f2d1) <- f2d1[,1]
f2d1 <- f2d1[,-1]
rownames(f2d2) <- f2d2[,1]
f2d2 <- f2d2[,-1]
f2d1 <- f2d1[c('A>G','T>C'),]
f2d_l <- as.data.frame(t(apply(f2d1,2,sum)))
f2d2 <- f2d2[c('A>G','T>C'),]
f2d_r <- as.data.frame(t(apply(f2d2,2,sum)))
f2d_l$mean <- apply(f2d_l,1,mean)
f2d_l$sd <- apply(f2d_l,1,sd)
f2d_r$mean <- apply(f2d_r,1,mean)
f2d_r$sd <- apply(f2d_r,1,sd)
rownames(f2d_l) <- 'leaf_potential_editing'
rownames(f2d_r) <- 'root_potential_editing'
f2d <- rbind(f2d_l[,c(3,4)],f2d_r[,c(3,4)])
f2d$tissue <- rownames(f2d)
ggplot(f2d,aes(x=tissue,y=mean)) + 
		geom_bar(aes(fill=tissue),stat="identity",position=position_dodge(0.9)) + 
		geom_errorbar(aes(x=tissue,ymin=mean-sd,ymax=mean+sd),width=0.2,size=1,alpha=0.7,colour="#4C5454") + 
		scale_fill_manual(values = c(leaf_potential_editing="#6fba30",root_potential_editing="#c8a063")) + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1 ),
				legend.position = 'none',panel.grid=element_blank(),
				axis.text = element_text(size= 12,color = "black"),
				axis.title = element_text(size = 12,color ="black")) + 
		labs(x=NULL,y='Amount of editing sites')
ggsave("figure_2D.pdf",width=5,height=5)

#####################################################################################################

f2e_l_com <- read.table('leaf_common_editing_information_dp20',skip = 2)
f2e_l_c1 <- f2e_l_com[,'V8',drop = F]
f2e_l_c2 <- f2e_l_com[,'V16',drop = F]
f2e_l_c1$tissue <- rep('leaf',nrow(f2e_l_c1))
f2e_l_c1$replicates <- rep('leaf_HiCEs_r1',nrow(f2e_l_c1))
f2e_l_c2$tissue <- rep('leaf',nrow(f2e_l_c2))
f2e_l_c2$replicates <- rep('leaf_HiCEs_r2',nrow(f2e_l_c2))
f2e_l_sp1 <- read.table('leaf_r1_specific_editing_information_dp20',skip=1)
f2e_l_sp1 <- f2e_l_sp1[,'V8',drop=F]
f2e_l_sp1$tissue <- rep('leaf',nrow(f2e_l_sp1))
f2e_l_sp1$replicates <- rep('leaf_specific_r1',nrow(f2e_l_sp1))
f2e_l_sp2 <- read.table('leaf_r2_specific_editing_information_dp20',skip=1)
f2e_l_sp2 <- f2e_l_sp2[,'V8',drop=F]
f2e_l_sp2$tissue <- rep('leaf',nrow(f2e_l_sp2))
f2e_l_sp2$replicates <- rep('leaf_specific_r2',nrow(f2e_l_sp2))
colnames(f2e_l_c1) <- c('editing_efficiency','tissue','replicates')
colnames(f2e_l_c2) <- c('editing_efficiency','tissue','replicates')
colnames(f2e_l_sp1) <- c('editing_efficiency','tissue','replicates')
colnames(f2e_l_sp2) <- c('editing_efficiency','tissue','replicates')
f2e_l <- rbind(f2e_l_c1,f2e_l_c2,f2e_l_sp1,f2e_l_sp2)

f2e_r_com <- read.table('root_common_editing_information_dp10',skip = 2)
f2e_r_c1 <- f2e_r_com[,'V8',drop = F]
f2e_r_c2 <- f2e_r_com[,'V16',drop = F]
f2e_r_c1$tissue <- rep('root',nrow(f2e_r_c1))
f2e_r_c1$replicates <- rep('root_HiCEs_r1',nrow(f2e_r_c1))
f2e_r_c2$tissue <- rep('root',nrow(f2e_r_c2))
f2e_r_c2$replicates <- rep('root_HiCEs_r2',nrow(f2e_r_c2))
f2e_r_sp1 <- read.table('root_r1_specific_editing_information_dp10',skip=1)
f2e_r_sp1 <- f2e_r_sp1[,'V8',drop=F]
f2e_r_sp1$tissue <- rep('root',nrow(f2e_r_sp1))
f2e_r_sp1$replicates <- rep('root_specific_r1',nrow(f2e_r_sp1))
f2e_r_sp2 <- read.table('root_r2_specific_editing_information_dp10',skip=1)
f2e_r_sp2 <- f2e_r_sp2[,'V8',drop=F]
f2e_r_sp2$tissue <- rep('root',nrow(f2e_r_sp2))
f2e_r_sp2$replicates <- rep('root_specific_r2',nrow(f2e_r_sp2))
colnames(f2e_r_c1) <- c('editing_efficiency','tissue','replicates')
colnames(f2e_r_c2) <- c('editing_efficiency','tissue','replicates')
colnames(f2e_r_sp1) <- c('editing_efficiency','tissue','replicates')
colnames(f2e_r_sp2) <- c('editing_efficiency','tissue','replicates')
f2e_r <- rbind(f2e_r_c1,f2e_r_c2,f2e_r_sp1,f2e_r_sp2)

f2e_data <- rbind(f2e_l,f2e_r)
ggplot(f2e_data, aes(x =tissue,y=editing_efficiency)) + 
		geom_split_violin(aes(fill=replicates),alpha=0.6) + 
		scale_fill_manual(values = c(leaf_HiCEs_r1="#046C1C",
										leaf_HiCEs_r2="#7CFC04",
										leaf_specific_r1='#B06CFF',
										leaf_specific_r2='#1CACF4',
										root_HiCEs_r1='#7C3404',
										root_HiCEs_r2='#444444',
										root_specific_r1='#FCA706',
										root_specific_r2='#F7DB12')) + 
		theme_bw() + 
		scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
		theme(axis.title = element_text(size = 12,color ="black"),
				axis.text = element_text(size= 12,color = "black"),
				panel.grid = element_blank(),
				axis.text.x = element_text(angle = 45, hjust = 1 ),
				legend.position = "top",
				legend.text = element_text(size= 12),
				legend.title= element_blank()
				) + 
		labs(x=NULL,y='Editing efficiency (%)')
ggsave("figure_2E.pdf",width=8,height=5)

#####################################################################################################
library(ggseqlogo) 

f2g_dl <- read.table('leaf_hices_5nt_flank.fq',header=F)
f2g_dr <- read.table('root_hices_5nt_flank.fq',header=F)
f2g_dl <- as.vector(f2g_dl$V1)
f2g_dr <- as.vector(f2g_dr$V1) 

col_l = make_col_scheme(chars=c('A', 'U', 'C', 'G'), cols=c('#a5cc29', '#dba272', '#efb8d3', '#5bbdc1'))
col_r = make_col_scheme(chars=c('A', 'U', 'C', 'G'), cols=c('#787cb4', '#b65620', '#ea9320', '#578d57'))

p1 <- ggplot() + 
		annotate('rect', xmin = 5.5, xmax = 6.5, ymin = -0.02, ymax = 1, alpha = .1, col='black',fill='transparent') + 
		geom_logo(f2g_dl, stack_width = 0.90,method='prob',col_scheme = col_l) +
		theme_logo() + 
		theme_bw() + 
		scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels = c("0","25","50","75","100")) +
		labs(x=NULL,y="Probablity (%)") + 
		theme(axis.text = element_text(size= 12,color = "black"),
				panel.grid = element_blank(),
				axis.text.x=element_blank()
				)
p2 <- ggplot() + 
		annotate('rect', xmin = 5.5, xmax = 6.5, ymin = -0.02, ymax = 1, alpha = .1, col='black',fill='transparent') + 
		geom_logo(f2g_dr, stack_width = 0.90,method='prob',col_scheme = col_r) +
		theme_logo() + 
		theme_bw() + 
		scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels = c("0","25","50","75","100")) +
		labs(x=NULL,y="Probablity (%)") + 
		theme(axis.text = element_text(size= 12,color = "black"),
				panel.grid = element_blank()) + 
				scale_x_continuous(breaks=c(seq(1,11,1)),labels=c("-5","-4","-3","-2","-1","0","1","2","3","4","5")
				)
f2g_plot <- (p1/p2)+plot_annotation(tag_levels = list(c('leaf', 'root')))
ggsave("figure_2G.pdf",width=8,height=5)

#####################################################################################################

# reference:https://slowkow.com/notes/ggplot2-color-by-density/

library(MASS)
library(viridis)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

f2f_dl <- cbind(f2e_l_c1$editing_efficiency,f2e_l_c2$editing_efficiency)
f2f_dr <- cbind(f2e_r_c1$editing_efficiency,f2e_r_c2$editing_efficiency)
f2f_dl <- as.data.frame(f2f_dl)
f2f_dr <- as.data.frame(f2f_dr)
colnames(f2f_dl) <- c('leaf_hices_r1','leaf_hices_r2')
colnames(f2f_dr) <- c('root_hices_r1','root_hices_r2')
f2f_dl$density <- get_density(f2f_dl$leaf_hices_r1,f2f_dl$leaf_hices_r2, n = 100)
f2f_dr$density <- get_density(f2f_dr$root_hices_r1,f2f_dr$root_hices_r2 ,n = 100)
pl <- ggplot(f2f_dl) + 
			geom_jitter(aes(leaf_hices_r1,leaf_hices_r2, color = density),size=1) + 
			scale_color_viridis(direction = -1,option = 'G') + 
			stat_cor(method = 'pearson', aes(x = leaf_hices_r1, y = leaf_hices_r2)) + 
			theme_bw() + 
			scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
			scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
			theme(axis.line = element_line(colour = "black"),
					axis.text = element_text(size= 12,color = "black"),
					panel.border = element_blank(),
					legend.position = 'none',
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank()
					) + 
			labs(x='hices_r1_editing efficiency (%)',y='hices_r2_editing efficiency (%)')
p1 <- ggExtra::ggMarginal(pl,type = 'density',margins = 'both',size = 1,colour = '#000000',fill = '#E6E6E6')
pr <- ggplot(f2f_dr) + 
			geom_jitter(aes(root_hices_r1,root_hices_r2, color = density),size=1) + 
			scale_color_viridis(direction = -1,option = 'E') + 
			stat_cor(method = 'pearson', aes(x = root_hices_r1, y = root_hices_r2)) + 
			theme_bw() + 
			scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
			scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels = c("20","40","60","80","100")) + 
			theme(axis.line = element_line(colour = "black"),
			axis.text = element_text(size= 12,color = "black"),
			panel.border = element_blank(),
			legend.position = 'none',
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
			) + 
			labs(x='hices_r1_editing efficiency (%)',y='hices_r2_editing efficiency (%)')
p2 <- ggExtra::ggMarginal(pr,type = 'density',margins = 'both',size = 1,colour = '#000000',fill = '#E6E6E6')
f2f_plot <- wrap_elements(p1)+wrap_elements(p2)+plot_annotation(tag_levels = list(c('leaf', 'root')))
ggsave("figure_2F.pdf",width=10,height=5)




