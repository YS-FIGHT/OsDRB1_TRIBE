# TRIBEseq
#
###########################################################################################
setwd("C:/Users/Ys/Desktop/paper_loading/data and code")

library(ggplot2)
library(scales)
library(ggpubr)
library(ggthemes)
library(patchwork)
library("webshot")
library(networkD3)
library(dplyr)
library(circlize)
library(grid)
library(ComplexHeatmap)

###########################################################################################

f3bdl1 <- read.table('data_summary_of_leaf_usual_hices_dp20_in_rice_anno_genomic',header=F)
f3bdl1 <- cbind(f3bdl1[c(1,3,5,7),,drop=F],f3bdl1[c(2,4,6,8),,drop=F])
colnames(f3bdl1) <- c('group','counts')
f3bdl1$counts <- as.numeric(f3bdl1$counts)
f3bdl1$labs <- paste0(round(100 * f3bdl1$counts / sum(f3bdl1$counts),2),"%")
pl1 <- ggplot(f3bdl1,aes(x="",y=counts,fill=group))+
				geom_bar(stat="identity",position = "fill",color='white') + 
				scale_fill_manual(values = c('3utr'="#095e02",'5utr'="#5cb76d",'intron'="#a4cb3b",'cds'="#468bca")) + 
				geom_text(aes(x=1.15,label=labs),position=position_fill(vjust = 0.5),colour = "white",size=4) + 
				labs(x=NULL,y=NULL,title="leaf HiCEs distribution") + 
				theme(axis.ticks = element_blank(),
						legend.title = element_blank(),
						axis.text= element_blank(),
						axis.line= element_blank(),
						panel.border = element_blank(),
						panel.background = element_blank(),
						panel.grid.major = element_blank(),
						panel.grid.minor = element_blank(),
						legend.position = 'top'
						) + 
				coord_polar(theta = "y")
ggsave("figure_3A_l1.pdf",width=3,height=3)

f3bdr1 <- read.table('data_summary_of_root_usual_hices_dp10_in_rice_anno_genomic',header=F)
f3bdr1 <- cbind(f3bdr1[c(1,3,5,7),,drop=F],f3bdr1[c(2,4,6,8),,drop=F])
colnames(f3bdr1) <- c('group','counts')
f3bdr1$counts <- as.numeric(f3bdr1$counts)
f3bdr1$labs <- paste0(round(100 * f3bdr1$counts / sum(f3bdr1$counts),2),"%")
pr1 <- ggplot(f3bdr1,aes(x="",y=counts,fill=group))+
				geom_bar(stat="identity",position = "fill",color='white') + 
				scale_fill_manual(values = c('3utr'="#999637",'5utr'="#dc9322",'intron'="#946134",'cds'="#aca085")) + 
				geom_text(aes(x=1.15,label=labs),position=position_fill(vjust = 0.5),colour = "white",size=4) + 
				labs(x=NULL,y=NULL,title="root HiCEs distribution") + 
				theme(axis.ticks = element_blank(),
						legend.title = element_blank(),
						axis.text= element_blank(),
						axis.line= element_blank(),
						panel.border = element_blank(),
						panel.background = element_blank(),
						panel.grid.major = element_blank(),
						panel.grid.minor = element_blank(),
						legend.position = 'top'
						) + 
				coord_polar(theta = "y")
ggsave("figure_3A_r1.pdf",width=3,height=3)

f3bdl2 <- read.csv('leaf_usual_hices_dp20_in_rice_anno_gene',skip = 1,sep="\t",header = T,check.names = F)
f3bdr2 <- read.csv('root_usual_hices_dp10_in_rice_anno_gene',skip = 1,sep="\t",header = T,check.names = F)
f3bdl2 <- f3bdl2[,ncol(f3bdl2),drop=F]
f3bdr2 <- f3bdr2[,ncol(f3bdr2),drop=F]
l_freq <- as.data.frame(table(f3bdl2))
l_freq <- l_freq$Freq
l_freq <- as.data.frame(table(l_freq))
r_freq <- as.data.frame(table(f3bdr2))
r_freq <- r_freq$Freq
r_freq <- as.data.frame(table(r_freq))
colnames(l_freq) <- c('editing_counts','number_of_genes')
colnames(r_freq) <- c('editing_counts','number_of_genes')
l_freq$editing_counts <- as.numeric(as.character(l_freq$editing_counts))
r_freq$editing_counts <- as.numeric(as.character(r_freq$editing_counts))
pl2 <- ggplot(l_freq,aes(x=editing_counts,y=number_of_genes)) + 
				geom_bar(stat='identity',position = position_dodge2(0.5),width=0.8,fill='#60ab34') + 
				theme_bw()+xlim(0,41) +
				scale_y_continuous(limits=c(0,100)) + 
				theme(axis.text = element_text(size= 12,color = "black")) + 
				labs(x='Editing counts',y='Number of genes')
ggsave("figure_3A_l2.pdf",width=5,height=5)

pr2 <- ggplot(r_freq,aes(x=editing_counts,y=number_of_genes)) + 
				geom_bar(stat='identity',position = position_dodge2(0.5),width=0.8,fill='#755a41') + 
				theme_bw() + xlim(0,50) +
				scale_y_continuous(limits=c(0,120)) + 
				theme(axis.text = element_text(size= 12,color = "black")) + 
				labs(x='Editing counts',y='Number of genes')
ggsave("figure_3A_r2.pdf",width=5,height=5)

pl1/pl2|pr1/pr2
ggsave("figure_3A.pdf",width=10,height=5)

#################################################################################################

dlg <- read.table('data_summary_of_leaf_usual_hices_dp20_in_rice_anno_gene_repeat',header=F,sep="\t",fill=T)
dlg <- dlg[,1:2]
colnames(dlg) <- c('target','value')
dlg$source <- c(rep('genomic',nrow(dlg)))
dlg <- dlg[,c('source','target','value')]
dlug <- read.table('data_summary_of_leaf_usual_hices_dp20_in_rice_anno_uncommented_repeat',header=F,sep="\t",fill=T)
dlug <- dlug[,1:2]
colnames(dlug) <- c('target','value')
dlug$source <- c(rep('uncommented',nrow(dlug)))
dlug <- dlug[,c('source','target','value')]

drg <- read.table('data_summary_of_root_usual_hices_dp10_in_rice_anno_gene_repeat',header=F,sep="\t",fill=T)
drg <- drg[,1:2]
colnames(drg) <- c('target','value')
drg$source <- c(rep('genomic',nrow(drg)))
drg <- drg[,c('source','target','value')]
drug <- read.table('data_summary_of_root_usual_hices_dp10_in_rice_anno_uncommented_repeat',header=F,sep="\t",fill=T)
drug <- drug[,1:2]
colnames(drug) <- c('target','value')
drug$source <- c(rep('uncommented',nrow(drug)))
drug <- drug[,c('source','target','value')]

dt <- rbind(dlg,dlug,drg,drug)
dt_temp <- dt[grep(pattern = "total",dt[,2]),]
dt <- dt[-grep(pattern = "total",dt[,2]),]
dt_temp$source <- c('leaf','leaf','root','root')
dt_temp$target <- rep(c('genomic','uncommented'),2)

links <- rbind(dt,dt_temp)
links$group <- c(rep('l',nrow(dlg)*2-2),rep('r',nrow(drg)*2-2),rep('l',2),rep('r',2))
nodes <- data.frame(
    name=c(as.character(links$source), 
           as.character(links$target)) %>% unique()
)
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

nodes$group <- as.factor(c(nodes$name))
my_color <- 'd3.scaleOrdinal() .domain(["l",
										"r",
										"leaf",
										"root",
										"genomic",
										"uncommented",
										"other_repeat",
										"mite",
										"ltr",
										"mudr",
										"line",
										"sine",
										"En_Spm",
										"TYPEU",
										"hAT",
										"no_repeat"]) 
								.range(["#8FBC8F",
										"#EEE8AA",
										"#006400",
										"#DAA520",
										"#1E90FF",
										"#B0E0E6",
										"#C0C0C0",
										"#FFD700",
										"#7B68EE",
										"#FF7F50",
										"#dbc7ac",
										"#b89985",
										"#a54d5e",
										"#004545",
										"#253493",
										"#000000"])'
p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                   Value = "value", NodeID = "name", 
                   colourScale=my_color, LinkGroup="group",fontSize = 15,sinksRight = F,fontFamily = 'Arial', nodePadding = 20, iterations = 0)
saveNetwork(p,"sankey.html")
webshot("sankey.html" , "fig_3C.pdf")


#################################################################################################################


df <- read.table('karyotype',header=T,sep='\t')
dt_l <- read.table('leaf_common_editing_information_dp20',skip=2,sep='\t')
dt_r <- read.table('root_common_editing_information_dp10',skip=2,sep='\t')
dt_l$mean_efficiency <- apply(dt_l[,c('V8','V16')],1,mean)
dt_l$mean_depth <- apply(dt_l[,c('V7','V15')],1,mean)
dt_r$mean_efficiency <- apply(dt_r[,c('V8','V16')],1,mean)
dt_r$mean_depth <- apply(dt_r[,c('V7','V15')],1,mean)
dt_l <- dt_l[c('V1','V2','V3','mean_efficiency','mean_depth')]
dt_r <- dt_r[c('V1','V2','V3','mean_efficiency','mean_depth')]
colnames(dt_l) <- c('Chr','Start','End','mean_efficiency','mean_depth')
colnames(dt_r) <- c('Chr','Start','End','mean_efficiency','mean_depth')

dt_r[dt_r=='8546372'] <- NA  # This site has extremely high depth (r1=7702; r2=7585), which will affect the display of the depth of other sites, here we omit this site.
dt_r <- na.omit(dt_r)

pdf("figure_3E.pdf.pdf")
circos.clear()
circos.par(gap.degree=c(rep(5,11),15),track.height=0.1,start.degree = 90)
circos.genomicInitialize(df,labels.cex = 0.7 ,axis.labels.cex = 0.4,sector.names = c(paste0('Chr ',1:12)))
circos.trackPlotRegion(sectors = df$Chr,ylim = c(0.2,1),track.height=0.15,track.margin=c(0.008, 0.001),bg.border=NA)
circos.yaxis(side = "left",
             at = c(0.25,0.5,0.75,1),
             labels = c('25%','50%','75%','100%'),
             labels.cex = 0.5, 
             lwd = 0.5,
             sector.index = 1,
             track.index=2)
col_fun = colorRamp2(seq(0.2,1,0.1), c('#3faf27', '#3bb81e', '#38bc19', '#35c111', '#00bd6a', '#00b983', '#00a9a8', '#007cd1', '#2b5cbe'))
circos.trackPoints(dt_l$Chr, dt_l$End, dt_l$mean_efficiency,pch = 16,cex = 0.15,col =col_fun(dt_l$mean_efficiency) )
col_fun = colorRamp2(seq(0.2,1,0.1), c('#d9b903', '#e2aa00', '#e6a300', '#e99b00', '#e98500', '#e97900', '#e56600','#dc3a00', '#d61604'))
circos.trackPoints(dt_r$Chr, dt_r$End, dt_r$mean_efficiency,pch = 16,cex = 0.15,col =col_fun(dt_r$mean_efficiency) )
circos.genomicDensity(dt_l[,1:3], col = '#008A00',bg.border=NA,track.margin=c(0.001, 0.001))
circos.genomicDensity(dt_r[,1:3], col = '#F0A30A',bg.border=NA)
circos.genomicTrack(
    dt_l[c(1:3,5)], track.height = 0.1, ylim = c(0, (max(dt_l$mean_depth) + 1)), bg.col = '#EEEEEE6E', bg.border = NA,
    panel.fun = function(region,value, ...) {
        circos.genomicLines(region, value, ytop.column = 1, ybottom = 0,col='#008A00',...)
        circos.yaxis(at = c(0,500,1000),
                     labels = c('0','500','1000'),
                     labels.cex = 0.3, 
                     lwd = 0.1, 
                     tick.length = convert_x(0.15, 'mm'),
                     sector.index = 1)
    }
)
circos.genomicTrack(
    dt_r[c(1:3,5)], track.height = 0.1, ylim = c(0, (max(dt_r$mean_depth) + 1)), bg.col = '#EEEEEE6E', bg.border = NA,
    panel.fun = function(region,value, ...) {
        circos.genomicLines(region, value, ytop.column = 1, ybottom = 0,col='#F0A30A',...)
        circos.yaxis(at = c(0,500),
                     labels = c('0','500'),
                     labels.cex = 0.3, 
                     lwd = 0.1, 
                     tick.length = convert_x(0.15, 'mm'),
                     sector.index = 1)
    }
)

title("Genome-wide distribution of high-confident editing sites")
pt <- Legend(labels =c("leaf-editing","root-editing"),
             type = "points",background = "white", 
             title = "Tissue",title_position = 'topcenter', 
             pch=c(16,16),
             title_gp = gpar(col = "black", fontsize = 8,font=2),
			 labels_gp=gpar(cex=0.8),
             legend_gp = gpar(col= c('#00b983', '#e97900')),ncol = 2
)
draw(pt,x = unit(x = 1, "npc"), y = unit(0.9, "npc"), just = c("right", "center"))
col_fun_l = colorRamp2(seq(0.2,1,0.1), c('#3faf27', '#3bb81e', '#38bc19', '#35c111', '#00bd6a', '#00b983', '#00a9a8', '#007cd1', '#2b5cbe'))
col_fun_r = colorRamp2(seq(0.2,1,0.1), c('#d9b903', '#e2aa00', '#e6a300', '#e99b00', '#e98500', '#e97900', '#e56600','#dc3a00', '#d61604'))
lgd_l <- Legend(col_fun = col_fun_l,
                at=c(0.2,0.4,0.6,0.8,1),
                labels = c("20","40","60","80",'100'),
                title = "leaf_editing_(%)",
                title_position = "leftcenter-rot",
                title_gp = gpar(cex=0.8),
                grid_width = unit(0.2, "cm")
)
lgd_r <- Legend(col_fun = col_fun_r,
                at=c(0.2,0.4,0.6,0.8,1),
                labels = c("20","40","60","80",'100'),
                title="root_editing_(%)",
                title_position = "leftcenter-rot",
                title_gp = gpar(cex=0.8),
                grid_width = unit(0.2, "cm")
)
draw(lgd_l,x = unit(x = 0.86, "npc"), y = unit(0.75, "npc"), just = c("right", "center"))
draw(lgd_r,x = unit(x = 0.98, "npc"), y = unit(0.75, "npc"), just = c("right", "center"))
dp <- Legend(labels =c("leaf","root"),
				title='Editing site depth',
				labels_gp=gpar(cex=0.8), 
				type = "lines",
				legend_gp = gpar(col =c('#008A00','#F0A30A')),
				grid_width = unit(0.3, "cm"),
				title_position = 'topcenter',
				nrow = 1,
				background = "white",
				title_gp = gpar(col = "black", fontsize =8,font=2),gap = unit(1, "cm"))

draw(dp,x = unit(x =0.96, "npc"), y = unit(0.6, "npc"), just = c("right", "center"))
dev.off()















