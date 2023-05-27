rm(list = ls())
library(tidyverse)
library(ConsensusClusterPlus)
load('Database/GTEx_expr.rda')
load('Database/Meta_RNAseq.rda')
load('Database/GRCh38_v39.gtf.rda')
GBM <- merge(gtf,Meta_RNAseq_expr,by.x=1,by.y=0)%>%
  dplyr::select(-gene_id)%>%
  group_by(gene_name)%>%
  summarise_all(max)%>%
  tibble::column_to_rownames('gene_name')

load("Database/net.rda")#Reading the gene interaction file
all_PPI_net <- net
all_PPI_net <- distinct(all_PPI_net,SYMBOL1,SYMBOL2,.keep_all = T)
library(tidyverse)
library(ggplot2)

#Construction of edge perturbation matrix-------------------------------------
source('Resource/lzq_genenetwork.R')#Loading the R package used to construct the edge perturbation matrix
network <- lzq_genenetwork(net = all_PPI_net,#Constructing the edge perturbation matrix
                           tumor_expr_data =GBM,
                           normal_expr_data = GTEx_expr)

save(network,file = "RData/03.network-clustering/network.rdata")

#Filtering edges for clustering and generate edge perturbation matrix for clustering-------------------------
sf <- lzq_selectfeature(network,count = 4000)
save(sf,file = "RData/03.network-clustering/sf.rdata")

#Consensus Clustering-------------------------
dir.create('Figure/03.network-clustering/ConsensusCluster')
results = ConsensusClusterPlus(d = as.matrix(sf$tumor_feature_matrix),
                               maxK=10,
                               reps=1000,
                               pItem=0.8,
                               pFeature=1,
                               title='Figure/03.network-clustering/ConsensusCluster',
                               innerLinkage="complete",
                               finalLinkage="complete",
                               clusterAlg="pam",
                               distance="euclidean", 
                               seed=123456,
                               plot="png")
icl <- calcICL(results,title = 'Figure/03.network-clustering/ConsensusCluster',plot = 'png')
save(results,file = "RData/03.network-clustering/consensus-results.rdata")

#Generating the files needed to draw the clustering heatmap----------------------------
rm(list = ls())
library(survminer)
library(survival)
load("RData/03.network-clustering/consensus-results.rdata")
clusterNum=4
cluster=results[[clusterNum]][["consensusClass"]]

sub <- data.frame(Sample=names(cluster),Cluster=cluster)
sub$Cluster <- paste0('C',sub$Cluster)
table(sub$Cluster)

head(sub)
cc <- sub$Cluster
names(cc) <- sub$Sample
cc2 <- sort(cc)

my <- results[[4]][["ml"]]
rownames(my) <- sub$Sample
colnames(my) <- sub$Sample
my2 <- my[names(cc2),names(cc2)]
save(cc2,my2,file = 'RData/03.network-clustering/consensus-matrix.rda')

#Drawing Clustering heatmap-------------------
rm(list = ls())
library(pheatmap)
library(tidyverse)
load('RData/03.network-clustering/Consensus-matrix.rda')
pdf(file='Figure/Figure2/A_pheat-cluster.pdf',width=5,height=4.5)
col <- c('#438F5D','#E19412','#2C6589',"#BC3C29FF")
pheatmap(my2,show_colnames = F,show_rownames = F,
         cluster_rows = F,cluster_cols = F,
         clustering_method = 'complete',border_color = 'black',
         color = colorRampPalette(c("white","#C75D30"))(50),
         annotation_names_row = F,annotation_names_col = F,
         annotation_row = data.frame(Subtype=cc2),
         annotation_col = data.frame(Subtype=cc2),
         annotation_colors = list(Subtype=c('C1'=col[1],'C2'=col[2],
                                            'C3'=col[3],'C4'=col[4])))
dev.off()

#Survival Analysis----------------------
rm(list = ls())
library(survminer)
library(survival)
load("RData/03.network-clustering/consensus-results.rdata")
clusterNum=4
cluster=results[[clusterNum]][["consensusClass"]]
sub <- data.frame(Sample=names(cluster),Cluster=cluster)
sub$Cluster <- paste0('C',sub$Cluster)
table(sub$Cluster)
load("Database/Meta_RNAseq_clin.rda")
clin <- merge(sub,Meta_RNAseq_clin,by=1)
save(clin,file = "RData/phenotype_RNAseq.rda")
fit <- survfit(Surv(OS.time,OS)~Cluster,clin)
surv_pvalue(fit)$pval
summary(coxph(Surv(OS.time,OS)~Cluster,data=clin))

#Drawing survival curves---------------------------
ggsurvplot(fit,clin,pval = T)
sdf <- survdiff(Surv(OS.time,OS)~Cluster,data = clin)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}

library("ggquickeda")
ggplot(clin, aes(time = OS.time, status = OS,color=Cluster,linetype=Cluster)) +
  geom_km(size=1.2) + geom_kmticks()+
  theme_bw(base_rect_size = 1.5)+
  scale_color_manual(values = c('#438F5D','#E19412','#2C6589',"#BC3C29FF"))+  
  scale_linetype_manual(values = c("solid","dashed","twodash","longdash"))+
  labs(y='Overall survival',x='Time (years)')+
  theme_bw(base_rect_size = 2)+
  theme(legend.text = element_text(size = 10,colour = 'black'),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "#F5F8FF", color = NA),
    axis.text = element_text(size = 12,colour = 'black'),
    axis.title = element_text(colour = "black",size = 18),
    panel.background = element_rect(fill = "#F5F8FF", color = NA),
    panel.grid = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
    panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
    axis.line = element_line(color = "#606F7B"),
    plot.title = element_text(colour = "black",size = 18,hjust = 0.5))+
  annotate('text',x = 7.25,y=0.75,label=paste0('Log-rank\n',p),hjust=0,size=6,fontface = 'italic')+
  ggtitle('Meta-RNAseq')
ggsave(filename = "Figure/Figure2/B_survival_RNAseq.pdf",width = 7,height = 5)


#Validation in microarray cohort----------
rm(list = ls())
library(tidyverse)
library(CMScaller)
load('Database/Meta_RNAseq.rda')
load('RData/phenotype_RNAseq.rda')
source('Resource/lzq_genenetwork.R')
load("Database/Meta_Array.rda")
test <- Meta_Array_expr
train <- Meta_RNAseq_expr
train <- t(scale(t(train)))%>%as.data.frame()
test <- t(scale(t(test)))%>%as.data.frame()%>%filter(!if_all(.fns = is.na))

#Filtering of signature genes -------------------------------------------------------------------------
results <- lzq_get_ntp_template(phenotype_data = clin,Sample_column_index = 1,Subtype_column_index = 2,
                                expression_data = train,count = 300)

ENSG <- results$templates
load("Database/GRCh38_v39.gtf.rda")
gene <- results$templates
gene <- merge(gtf,gene,by.x=1,by.y=2)%>%
  dplyr::select(-gene_id)

save(ENSG,file = "RData/04.validation/signature_ENSG.rdata")
save(gene,file = "RData/04.validation/signature_gene.rdata")
writexl::write_xlsx(ENSG,path = "Table/signature_ENSG.xlsx")
#Validation with NTP-------------------------------
load("RData/04.validation/signature_ENSG.rdata")
library(clusterProfiler)
###COC
ll <- lzq_COC_normalize(data1 = train,
                        data2 = test,
                        ABS_Cor_cutoff = 0.5)
emat <- t(scale(t(ll$data2)))
save(emat,file="RData/04.validation/emat_microarray.rda")
ntp_res <- ntp(emat      = emat,
               templates = ENSG,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 1,
               seed      = 2020104,
               verbose   = T)

ntp_res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')
save(ntp_res,file = "RData/04.validation/ntp_heatmapinput_microarray.rda")
#Survival Analysis----------------
rm(list = ls())
load("RData/04.validation/ntp_heatmapinput_microarray.rda")
load("Database/Meta_Array_clin.rda")

table(ntp_res$prediction)
dd <- merge(ntp_res[,1:2],Meta_Array_clin,by=1)
colnames(dd)[2] <- 'Cluster'
phenotype_array <- dd
save(phenotype_array,file = "RData/04.validation/phenotype_microarray.rdata")

library(survminer)
library(survival)
#绘制生存图
fit <- survfit(Surv(OS.time,OS)~Cluster,data = dd)
ggsurvplot(fit,dd,pval = T)
sdf <- survdiff(Surv(OS.time,OS)~Cluster,data = dd)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}

library("ggquickeda")
ggplot(dd, aes(time = OS.time, status = OS,color=Cluster,linetype=Cluster)) +
  geom_km(size=1.2) + geom_kmticks()+
  theme_bw(base_rect_size = 1.5)+
  scale_color_manual(values = c('#438F5D','#E19412','#2C6589',"#BC3C29FF"))+    #c(pal_npg("nrc", alpha =0.9)(9)[c(1)],pal_npg("nrc", alpha =0.9)(9)[c(2)]))+
  scale_linetype_manual(values = c("solid","dashed","twodash","longdash"))+
  labs(y='Overall survival',x='Time (years)')+
  theme_bw(base_rect_size = 2)+
  theme(#legend.position = c(0.82,0.84),#定义位置可以把注释放在图内
    legend.text = element_text(size = 10,colour = 'black'),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "#F5F8FF", color = NA),##图注的背景色
    axis.text = element_text(size = 12,colour = 'black'),
    axis.title = element_text(colour = "black",size = 18),
    panel.background = element_rect(fill = "#F5F8FF", color = NA),
    panel.grid = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
    panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
    axis.line = element_line(color = "#606F7B"),
    plot.title = element_text(colour = "black",size = 18,hjust = 0.5))+
  annotate('text',x = 7.25,y=0.75,label=paste0('Log-rank\n',p),hjust=0,size=6,fontface = 'italic')+
  ggtitle('Meta-microarray')

ggsave(filename = "Figure/Figure2/G_survival_microarray.pdf",width = 7,height = 4.8)

#Drawing Clustering heatmap----------------------------------------
rm(list = ls())
load("RData/04.validation/ntp_heatmapinput_microarray.rda")
load("RData/04.validation/emat_microarray.rda")
load("RData/04.validation/signature_ENSG.rdata")
expr <- as.data.frame(emat)
expr <- expr[rownames(expr)%in%ENSG$probe,]
expr <- t(expr)
expr <- merge(expr,ntp_res[,1:2],by.x=0,by.y=1,all=F)
expr <- expr[,c(1,582,2:581)]

gene <- rownames(as.data.frame(emat))
gene <- gene[gene%in%ENSG$probe]
GENE <- ENSG[ENSG$probe%in%gene,]
library(tibble)
a <- expr[,c(1,3:582)] 
rownames(a) <- a[,1]
a <- a[,-1]
a <- a[,GENE$probe]
expr <- merge(expr[,1:2],a,by.x = 1,by.y = 0)
expr <- expr[order(expr$prediction),]

Cluster <- c('#438F5D','#E19412','#2C6589',"#BC3C29FF")
names(Cluster) <- c('C1','C2',"C3","C4")
library(ComplexHeatmap)
library(circlize)
Top = HeatmapAnnotation(
  border = T,
  show_legend = F,
  foo = anno_block(gp = gpar(fill = Cluster),
                   labels = c("C1", "C2", "C3","C4"), 
                   labels_gp = gpar(col = "black", fontsize = 10),
                   height = unit(0.5,"cm")),
  show_annotation_name = F )
Left = HeatmapAnnotation(border = T,
                         which = "row",
                         foo = anno_block(gp = gpar(fill = Cluster),
                                          labels = c("C1", "C2", "C3","C4"), 
                                          labels_gp = gpar(col = "black", fontsize = 10),
                                          width = unit(0.5,"cm")),#注释条的宽度 
                         show_annotation_name = F )
pdf(file = 'Figure/Figure2/C_NTPheatmap_microarray.pdf',width=5,height = 4.5)
Heatmap(t(scale(expr[,c(3:582)])),
        name='Meta-microarray',
        top_annotation = Top,
        left_annotation = Left,
        cluster_rows = FALSE,
        col=colorRamp2(c(-2,0,2),c('#82B0C7','white','#FF5768')),
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        column_order=NULL,
        show_column_names = FALSE,
        show_row_names = FALSE,  
        row_names_gp = gpar(fontsize = 7),
        column_split = c(rep(1,225),rep(2,140),rep(3,208),rep(4,170)),
        row_title = "Gene features",row_title_rot = 90,
        column_title_rot = 0,column_title_side = "top",column_title="Class predictions",
        row_split = factor(rep(c('C1','C2','C3','C4'),
                               times=c(172,96,150,162)),levels = c('C1','C2','C3','C4')),
        
        gap = unit(1.5, "mm"),
        column_title_gp = gpar(fontsize = 12),
        row_title_gp = gpar(fontsize=12),
        show_heatmap_legend = F,
) 
dev.off()

#Validation in ZZU cohort----------------------------------------
rm(list = ls())
library(tidyverse)
library(CMScaller)
library(survminer)
library(survival)
load('Database/Meta_RNAseq.rda')
source('Resource/lzq_genenetwork.R')
load("Database/ZGBM_expr.rda")
test <- ZGBM_expr
train <- Meta_RNAseq_expr
train <- t(scale(t(train)))%>%as.data.frame()
test <- t(scale(t(test)))%>%as.data.frame()%>%filter(!if_all(.fns = is.na))
load("RData/04.validation/signature_ENSG.rdata")
library(clusterProfiler)
#COC
ll <- lzq_COC_normalize(data1 = train,
                        data2 = test,
                        ABS_Cor_cutoff = 0.5)
emat <- t(scale(t(ll$data2)))
save(emat,file = "RData/04.validation/emat_ZGBM.rda")
#Validation with NTP----------------
ntp_res <- ntp(emat      = emat,
               templates = ENSG,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 1,
               seed      = 2020104,
               verbose   = T)
ntp_res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')
table(ntp_res$prediction)
save(ntp_res,file = "RData/04.validation/ntp_heatmapinput_ZGBM.rda")

rm(list = ls())
load("RData/04.validation/ntp_heatmapinput_ZGBM.rda")
load("Database/ZGBM_clin.rda")
ZGBM_clin <- ZGBM_clin[,1:5]
dd <- merge(ntp_res[,1:2],ZGBM_clin,by.x=1,by.y = 0)
colnames(dd)[2] <- 'Cluster'
save(dd,file = "RData/04.validation/phenotype_ZGBM.rdata")
library(survminer)
library(survival)
#Survival Analysis------------------------
fit <- survfit(Surv(OS.time,OS)~Cluster,data = dd)
ggsurvplot(fit,dd,pval = T)
sdf <- survdiff(Surv(OS.time,OS)~Cluster,data = dd)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}

library("ggquickeda")
ggplot(dd, aes(time = OS.time, status = OS,color=Cluster,linetype=Cluster)) +
  geom_km(size=1.2) + geom_kmticks()+
  theme_bw(base_rect_size = 1.5)+
  scale_color_manual(values = c('#438F5D','#E19412','#2C6589',"#BC3C29FF"))+    
  scale_linetype_manual(values = c("solid","dashed","twodash","longdash"))+
  labs(y='Overall survival',x='Time (years)')+
  theme_bw(base_rect_size = 2)+
  theme(legend.text = element_text(size = 10,colour = 'black'),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "#F5F8FF", color = NA),
    axis.text = element_text(size = 12,colour = 'black'),
    axis.title = element_text(colour = "black",size = 18),
    panel.background = element_rect(fill = "#F5F8FF", color = NA),
    panel.grid = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
    panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
    axis.line = element_line(color = "#606F7B"),
    plot.title = element_text(colour = "black",size = 18,hjust = 0.5))+
  annotate('text',x = 2.25,y=0.75,label=paste0('Log-rank\n',p),hjust=0,size=6,fontface = 'italic')+
  ggtitle('ZZU in-house')

ggsave(filename = "Figure/Figure2/H_survival_ZGBM.pdf.pdf",width = 7,height = 4.8)

#Drawing Clustering heatmap---------------------------------------
rm(list = ls())
load("RData/04.validation/ntp_heatmapinput_ZGBM.rda")
load("RData/04.validation/emat_ZGBM.rda")
load("RData/04.validation/signature_ENSG.rdata")
expr <- as.data.frame(emat)
expr <- expr[rownames(expr)%in%ENSG$probe,]
expr <- t(expr)
expr <- merge(expr,ntp_res[,1:2],by.x=0,by.y=1,all=F)
expr <- expr[,c(1,1011,2:1010)]
gene <- rownames(as.data.frame(emat))
gene <- gene[gene%in%ENSG$probe]
GENE <- ENSG[ENSG$probe%in%gene,]
library(tibble)
a <- expr[,c(1,3:1011)] 
rownames(a) <- a[,1]
a <- a[,-1]
a<- a[,GENE$probe]
expr <- merge(expr[,1:2],a,by.x = 1,by.y = 0)
expr <- expr[order(expr$prediction),]
Cluster <- c('#438F5D','#E19412','#2C6589',"#BC3C29FF")
names(Cluster) <- c('C1','C2',"C3","C4")
library(ComplexHeatmap)
library(circlize)
Top = HeatmapAnnotation(border = T,
                        show_legend = F,
                        foo = anno_block(gp = gpar(fill = Cluster),
                                         labels = c("C1", "C2", "C3","C4"), 
                                         labels_gp = gpar(col = "black", fontsize = 10),
                                         height = unit(0.5,"cm")),
                        show_annotation_name = F )
Left = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = Cluster),
                                          labels = c("C1", "C2", "C3","C4"), 
                                          labels_gp = gpar(col = "black", fontsize = 10),
                                          width = unit(0.5,"cm")),
                         border = T,
                         which = "row",
                         show_annotation_name = T )

pdf(file = 'Figure/Figure2/D_NTPheatmap_ZGBM.pdf',width=5,height = 4.5)
Heatmap(t(scale(expr[,c(3:1011)])),
        top_annotation = Top, 
        left_annotation = Left,
        cluster_rows = FALSE,
        col=colorRamp2(c(-2,0,2),c('#82B0C7','white','#FF5768')),
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        column_order=NULL,
        show_column_names = FALSE,
        show_row_names = FALSE,  
        row_names_gp = gpar(fontsize = 7),
        column_split = c(rep(1,63),rep(2,40),rep(3,58),rep(4,49)),
        row_title = "Gene features",row_title_rot = 90,
        column_title_rot = 0,column_title_side = "top",column_title="Class predictions",
        row_split = factor(rep(c('C1','C2','C3','C4'),
                               times=c(241,270,266,232)),levels = c('C1','C2','C3','C4')),
        
        gap = unit(1.5, "mm"),
        column_title_gp = gpar(fontsize = 12),
        row_title_gp = gpar(fontsize=12),
        show_heatmap_legend = F,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 9), 
                                  title_gp = gpar(fontsize = 9, fontface = "bold"))
) 

dev.off()

#Performing SubMap analysis--------------------------------
rm(list = ls())
library(tibble)
library(dplyr)

generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

load("RData/04.validation/signature_ENSG.rdata")
#discovery cohort
load('Database/Meta_RNAseq_expr.rda')
load('RData/phenotype_RNAseq.rda')

discovery.logNC <- Meta_RNAseq_expr
rm(Meta_RNAseq_expr)
discovery.info <- clin
rm(clin)
discovery.info <- discovery.info[order(discovery.info$Cluster),1:2]
discovery.info$rank <- rep(c(1,2,3,4),times=as.character(table(discovery.info$Cluster))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R
rownames(discovery.info) <- discovery.info$Sample
discovery.info <- discovery.info[,2:3]


#validation cohort
##Meta_array
load("Database/Meta_Array.rda")
rm(Meta_Array_clin)
array.logNC <- Meta_Array_expr
rm(Meta_Array_expr)
load("RData/04.validation/phenotype_microarray.rdata")
array.info <- phenotype_array
rm(phenotype_array)
array.info <- array.info[order(array.info$Cluster),1:2]
array.info$rank <- rep(c(1,2,3,4),times=as.character(table(array.info$Cluster))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R
rownames(array.info) <- array.info$ID
array.info <- array.info[,2:3]

##ZGBM
load("Database/ZGBM_expr.rda")
ZGBM.logNC <- ZGBM_expr
rm(ZGBM_expr)
load("RData/04.validation/phenotype_ZGBM.rdata")
ZGBM.info <- dd
rm(dd)
ZGBM.info <- ZGBM.info[order(ZGBM.info$Cluster),1:2]
ZGBM.info$rank <- rep(c(1,2,3,4),times=as.character(table(ZGBM.info$Cluster))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R
rownames(ZGBM.info) <- ZGBM.info$ID
ZGBM.info <- ZGBM.info[,2:3]

#discovery VS array
GENELIST <- ENSG$probe
GENELIST <- intersect(GENELIST,rownames(array.logNC))
sam_info <- discovery.info
in_gct <- discovery.logNC[GENELIST,rownames(discovery.info)]
gct_file <- "Reference/04.validation/discovery/discovery_array.gct"
cls_file <- "Reference/04.validation/discovery/discovery_array.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

sam_info <- array.info
in_gct <- array.logNC[GENELIST,rownames(array.info)]
gct_file <- "Reference/04.validation/validation/array.gct"
cls_file <- "Reference/04.validation/validation/array.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

#discovery VS ZGBM
GENELIST <- intersect(ENSG$probe,rownames(ZGBM.logNC))
sam_info <- discovery.info
in_gct <- discovery.logNC[GENELIST,rownames(discovery.info)]
gct_file <- "Reference/04.validation/discovery/discovery_ZGBM.gct"
cls_file <- "Reference/04.validation/discovery/discovery_ZGBM.cls"

generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")
sam_info <- ZGBM.info
in_gct <- ZGBM.logNC[GENELIST,rownames(ZGBM.info)]
in_gct[is.na(in_gct)] <- 0

gct_file <- "Reference/04.validation/validation/ZGBM.gct"
cls_file <- "Reference/04.validation/validation/ZGBM.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

###Drawing submap heatmap
rm(list=ls())
library(pheatmap)
library(export)
library(ggsci)
tmp <- matrix(c(0.000999001, 0.261738262, 0.991008991, 0.944055944,
                0.844155844, 0.000999001, 0.737262737, 0.987012987,
                1.000000000, 0.702297702, 0.000999001, 0.903096903,
                0.893106893, 0.997002997, 0.894105894, 0.000999001), 
              nrow =4,byrow = T,dimnames = list(c("C1A","C2A","C3A","C4A"),
                                                c("C1R","C2R","C3R","C4R")))
library(tidyverse)
library(data.table)
library(pROC)
library(caret)
library(survival)
library(survminer)
library(ggkm)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(export)
dd1 <- as.data.frame(t(tmp))
dd2 <- dd1
rownames(dd2) <- substr(rownames(dd2),1,2)
colnames(dd2) <- substr(colnames(dd2),1,2)
cell_fun = function(j, i, x, y, width, height, fill) {
  if(dd2[i, j]==1){
    grid.text(sprintf("%.0f", dd2[i, j]), x, y, gp = gpar(fontsize = 10))}
  if(dd2[i, j]<0.05&dd2[i, j]>0.001){
    grid.text(sprintf("%.3f", dd2[i, j]), x, y, gp = gpar(fontsize = 10))}
  if(dd2[i, j]<0.001){
    grid.text(sprintf("<0.001", dd2[i, j]), x, y, gp = gpar(fontsize = 10))}
  if(dd2[i, j]>0.05&dd2[i, j]!=1){
    grid.text(sprintf("%.2f", dd2[i, j]), x, y, gp = gpar(fontsize = 10))}
}

Top = HeatmapAnnotation(Subtype=paste0('C',1:4),border = F,
                        show_legend = F,show_annotation_name = F,
                        simple_anno_size =unit(3,'mm'),
                        col = list(Subtype = c('C1' = '#438F5D',
                                               'C2' = '#E19412',
                                               'C3' = '#2C6589',
                                               'C4' = "#BC3C29FF")))
left <- HeatmapAnnotation(Subtype=paste0('C',1:4),border = F,
                          show_legend = F,show_annotation_name = F,
                          simple_anno_size =unit(3,'mm'),
                          which = 'row',
                          col = list(Subtype = c('C1' = '#438F5D',
                                                 'C2' = '#E19412',
                                                 'C3' = '#2C6589',
                                                 'C4' = "#BC3C29FF")))
dd2 <- as.matrix(dd2)
pdf(file="Figure/Figure2/E_Submap_microarray.pdf",width=4.2,height=4.2)
Heatmap(dd2,name = 'Sig',
        top_annotation = Top,left_annotation = left,
        col = colorRamp2(c(0,0.33,0.66,1),c("#58739C","#C9D5E8","#D0DCF0","#DBE7FC")),
        cell_fun = cell_fun,
        show_row_names = F,show_column_names = F,
        cluster_rows = F,border = T,
        cluster_columns = F,show_heatmap_legend = F,
        column_title = "Subtypes in Meta-microarray Cohort",
        row_title = 'Subtypes in Discovery Cohort',
        gap = unit(1.2,'mm'),
        column_gap = unit(1.5,'mm'))
dev.off()
rm(list=ls())
library(pheatmap)
library(export)
library(ggsci)

tmp <- matrix(c(0.000999001, 0.359640360, 0.971028971, 0.929070929,
                0.533466533, 0.005994006, 0.390609391, 0.982017982,
                0.999000999, 0.730269730, 0.000999001, 0.768231768,
                0.724275724, 0.897102897, 0.921078921, 0.000999001), # 校正p值
              nrow = 4,byrow = T,dimnames = list(c("C1Z","C2Z","C3Z","C4Z"),
                                                 c("C1R","C2R","C3R","C4R")))

library(tidyverse)
library(data.table)
library(pROC)
library(caret)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(export)
dd1 <- as.data.frame(t(tmp))
dd2 <- dd1
rownames(dd2) <- substr(rownames(dd2),1,2)
colnames(dd2) <- substr(colnames(dd2),1,2)
cell_fun = function(j, i, x, y, width, height, fill) {
  if(dd2[i, j]==1){
    grid.text(sprintf("%.0f", dd2[i, j]), x, y, gp = gpar(fontsize = 10))}
  if(dd2[i, j]<0.05&dd2[i, j]>0.001){
    grid.text(sprintf("%.3f", dd2[i, j]), x, y, gp = gpar(fontsize = 10))}
  if(dd2[i, j]<0.001){
    grid.text(sprintf("<0.001", dd2[i, j]), x, y, gp = gpar(fontsize = 10))}
  if(dd2[i, j]>0.05&dd2[i, j]!=1){
    grid.text(sprintf("%.2f", dd2[i, j]), x, y, gp = gpar(fontsize = 10))}
}

Top = HeatmapAnnotation(Subtype=paste0('C',1:4),border = F,
                        show_legend = F,show_annotation_name = F,
                        simple_anno_size =unit(3,'mm'),
                        col = list(Subtype = c('C1' = '#438F5D',
                                               'C2' = '#E19412',
                                               'C3' = '#2C6589',
                                               'C4' = "#BC3C29FF")))
left <- HeatmapAnnotation(Subtype=paste0('C',1:4),border = F,
                          show_legend = F,show_annotation_name = F,
                          simple_anno_size =unit(3,'mm'),
                          which = 'row',
                          col = list(Subtype = c('C1' = '#438F5D',
                                                 'C2' = '#E19412',
                                                 'C3' = '#2C6589',
                                                 'C4' = "#BC3C29FF")))
dd2 <- as.matrix(dd2)
pdf(file="Figure/Figure2/F_Submap_ZGBM.pdf",width=4.2,height=4.2)
Heatmap(dd2,name = 'Sig',
        top_annotation = Top,left_annotation = left,
        col = colorRamp2(c(0,0.33,0.66,1),c("#58739C","#C9D5E8","#D0DCF0","#DBE7FC")),#c("#CC6600","#F0D8C0","#F5DDC4","#FFEEE3")
        cell_fun = cell_fun,
        show_row_names = F,show_column_names = F,
        cluster_rows = F,border = T,
        cluster_columns = F,show_heatmap_legend = F,
        column_title = "Subtypes in ZZU in-house Cohort",
        row_title = 'Subtypes in Discovery Cohort',
        gap = unit(1.2,'mm'),
        column_gap = unit(1.5,'mm'))
dev.off()

#Comparing the distribution of the four subtypes in the three cohorts-----------------
rm(list = ls())
library(ggalluvial)
library(readxl)
library(dplyr)
library(reshape2)
library(ggsci)
load("RData/phenotype_RNAseq.rda")
table(clin$Cluster)
load("RData/04.validation/phenotype_microarray.rdata")
table(phenotype_array$Cluster)
load("RData/04.validation/phenotype_ZGBM.rdata")
table(dd$Cluster)
df <- data.frame(subtype <- c("C1","C2","C3","C4"),
                 discovery <- c(131,150,139,112),
                 array <- c(225,140,208,170),
                 zzu <- c(63,40,58,49))
colnames(df) <- c("subtype","Discovery","Meta-array","ZZU in-house")
df$Meta_RNAseq <- df$Discovery/sum(df$Discovery)*100
df$Meta_microarray <- df$`Meta-array`/sum(df$`Meta-array`)*100
df$ZZU <- df$`ZZU in-house`/sum(df$`ZZU in-house`)*100

sr <- df[,c(1,5:7)]
df_final <- melt(sr)
head(df_final)
df_final <- rename(df_final, Time=variable, Ratio=value)
df_final$Subject <- c(1:4, 1:4, 1:4) 
df_final <- transform(df_final, Time = factor(Time, levels = c("Meta_RNAseq","Meta_microarray","ZZU")))
head(df_final)
save(df_final,file = "RData/04.validation/proportion_input.rda")

#Drawing subtype scale maps
rm(list = ls())
load("RData/04.validation/proportion_input.rda")
library(ggplot2)
ggplot(df_final, aes(x = Time, stratum = subtype, alluvium = Subject, y = Ratio , 
                     fill = subtype,label = subtype, position = "fill")) +
  scale_x_discrete(expand = c(.1, .1)) +
  theme_bw()+
  scale_fill_manual(values = c('#438F5D','#E19412','#2C6589',"#BC3C29FF"),guide = guide_legend(reverse = TRUE))+
  geom_flow(alpha = 0.5, width = 1/3,) +
  theme_classic() + 
  geom_stratum(alpha = 1, width = 1/3) +
  geom_text(stat = "stratum", size = 3) +
  scale_color_discrete(limits =c("C1","C2","C3","C4"))+
  theme(legend.position = "right") + ylab("Ratio") + xlab("") 
ggsave(filename = "Figure/Figure2/I_proportion.pdf",width = 5,height = 2.5)


#GO analysis of C1 subtype----------------------------------------------------------
rm(list = ls())
library( org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
load("RData/04.validation/signature_gene.rdata")
cordata1 <- gene[gene$class=="C1",]
ID <- bitr(cordata1$gene_name,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
GO_pos <- enrichGO(gene = ID$ENTREZID,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = 'all',
                   pvalueCutoff = 0.05,qvalueCutoff = 0.05,minGSSize = 1,maxGSSize = 50000) 
d1 <- GO_pos@result%>%filter(ONTOLOGY=='BP')
dd1 <- d1[(1:10),]
tmp1 <- dd1%>%dplyr::select(ID,Description,GeneRatio,ONTOLOGY)%>%
  tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
  dplyr::rename(Ontology=ONTOLOGY)%>%
  dplyr::mutate(Enriched = as.numeric(Enriched),
                Detected = as.numeric(Detected),
                Ontology = ifelse(Ontology=='BP','Biological Process',
                                  ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
  dplyr::arrange(Ontology) %>%
  dplyr::mutate(ID = factor(ID, rev(unique(ID))))
d2 <- GO_pos@result%>%filter(ONTOLOGY=='CC')
dd2 <- d2[c(1:10),]
tmp2 <- dd2%>%dplyr::select(ID,Description,GeneRatio,ONTOLOGY)%>%
  tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
  dplyr::rename(Ontology=ONTOLOGY)%>%
  dplyr::mutate(Enriched = as.numeric(Enriched),
                Detected = as.numeric(Detected),
                Ontology = ifelse(Ontology=='BP','Biological Process',
                                  ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
  dplyr::arrange(Ontology) %>%
  dplyr::mutate(ID = factor(ID, rev(unique(ID))))
d3 <- GO_pos@result%>%filter(ONTOLOGY=='MF')
dd3 <- d3[c(1:10),]
tmp3 <- dd3%>%dplyr::select(ID,Description,GeneRatio,ONTOLOGY)%>%
  tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
  dplyr::rename(Ontology=ONTOLOGY)%>%
  dplyr::mutate(Enriched = as.numeric(Enriched),
                Detected = as.numeric(Detected),
                Ontology = ifelse(Ontology=='BP','Biological Process',
                                  ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
  dplyr::arrange(Ontology) %>%
  dplyr::distinct(Description,.keep_all = T)%>%
  dplyr::mutate(ID = factor(ID, rev(unique(ID))))
tmp <- rbind(tmp1,tmp2,tmp3)%>%dplyr::distinct(ID,.keep_all = T)
tmp$ID <- factor(tmp$ID,rev(tmp$ID))
rm(d1,d2,d3,dd1,dd2,dd3,tmp1,tmp2,tmp3)
gc()
tmp_l1 <- data.frame(nrow(tmp) - (cumsum(table(tmp$Ontology)) - table(tmp$Ontology)/2))
tmp_l1$start <- tmp_l1$Freq - table(tmp$Ontology)/2
tmp_l1$end <- tmp_l1$Freq + table(tmp$Ontology)/2
m1 <- ifelse(log(max(tmp$Detected),base = 20)*1000<1000,1000,log(max(tmp$Detected),base = 20)*1000)
m2 <- ifelse(log(max(tmp$Enriched),base = 20)*1000<1000,1000,log(max(tmp$Enriched),base = 20)*1000)
p1 <- ggplot(tmp) +
  geom_col(mapping = aes(ID, Detected, fill= Ontology),
           color = "black", fill = "#e6b8a2",
           width = 0.75, 
           show.legend = F) +
  geom_text(mapping = aes(ID, Detected, label = Detected),hjust=-0.3, size = 2.5) +
  scale_y_log10(limits = c(1, m1),expand = c(0,0),breaks=c(1,100),labels=c(1,100)) +
  coord_flip() + 
  theme_classic() +
  labs(x = NULL, y = NULL, title = "Detected Genes") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,2,1,0), "mm"),
        axis.text.x = element_text(size=7),
        plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',color='black',vjust=-1))

p2 <- ggplot(tmp) +
  geom_col(mapping = aes(ID, Enriched, fill= Ontology),color = "black", width = 0.75, show.legend = F) +
  geom_text(mapping = aes(ID, Enriched, label = Enriched),hjust=-0.3, size = 2.5) +
  scale_y_log10(limits = c(1, m2),expand = expansion(),breaks=c(1,100,1000),labels=c(1,100,1000)) +
  scale_fill_manual(values = alpha(c('#00b4d8','#ff477e','#FFA500'),0.9)) +
  coord_flip() + 
  theme_classic() +
  labs(x = NULL, y = NULL, title = "Enriched Genes",fill=NULL) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size=7),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',vjust=-1))

p0 <- ggplot(tmp) + 
  geom_text(mapping = aes(ID, 0, label = ID ,color = Ontology),size = 3, show.legend = F, hjust = 0) +
  scale_color_manual(values = alpha(c('#00b4d8','#ff477e','#FFA500'),0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.02)) +
  coord_flip() + theme_void() +
  labs(x = NULL, y = NULL, title = "Identifiers") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,0,1,3), "mm"),
        plot.title = element_text(hjust = 0.01, size = 10,face = 'bold',vjust=1.5))

p3 <- ggplot(tmp) +
  geom_text(mapping = aes(ID, 0, label = Description,color = Ontology), 
            size = 3, show.legend = F, hjust = 0) +
  scale_color_manual(values = alpha(c('#00b4d8','#ff477e','#FFA500'),0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.1)) +
  coord_flip() + 
  theme_void() +
  labs(x = NULL, y = NULL, title = "Description") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,0.2,1,1), "mm"),
        plot.title = element_text(hjust = 0, size = 10,face = 'bold',vjust=1.5))

p4 <- ggplot(tmp_l1) +
  geom_segment(mapping = aes(x = start+0.1, xend = end-0.1, y = -0.1, yend = -0.1), size = 2)+
  geom_text(mapping = aes(Freq, 0, label = Var1), size = 3, show.legend = F, hjust = 0) +
  scale_y_continuous(expand = expansion(), limits = c(-0.1,1)) +
  scale_x_continuous(expand = expansion()) +
  coord_flip() + 
  theme_void() +
  labs(x = NULL, y = NULL, title = "Ontoloty") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,1,1,0.2), "mm"),
        plot.title = element_text(hjust = 0.1,size = 10,face = 'bold',vjust=1.5))
cowplot::plot_grid(p0,p1,p2,p3,p4,align = "h", nrow = 1, 
                   rel_widths = c(0.1,0.15,0.15,0.35,0.2))

#GO analysis of C2 subtype--------------------------------------------------------
rm(list = ls())
library( org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
load("RData/04.validation/signature_gene.rdata")
cordata2 <- gene[gene$class=="C2",]
ID <- bitr(cordata2$gene_name,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
GO_pos <- enrichGO(gene = ID$ENTREZID,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = 'all',
                   pvalueCutoff = 0.05,qvalueCutoff = 0.05,minGSSize = 1,maxGSSize = 50000) 
d1 <- GO_pos@result%>%filter(ONTOLOGY=='BP')
a <- c("GO:0072359","GO:0035295","GO:0001944","GO:0035239","GO:0030334","GO:0001568","GO:0001525","GO:0007155","GO:0022610","GO:0033627")
dd1 <- d1[a,]
tmp1 <- dd1%>%dplyr::select(ID,Description,GeneRatio,ONTOLOGY)%>%
  tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
  dplyr::rename(Ontology=ONTOLOGY)%>%
  dplyr::mutate(Enriched = as.numeric(Enriched),
                Detected = as.numeric(Detected),
                Ontology = ifelse(Ontology=='BP','Biological Process',
                                  ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
  dplyr::arrange(Ontology) %>%
  dplyr::mutate(ID = factor(ID, rev(unique(ID))))
d2 <- GO_pos@result%>%filter(ONTOLOGY=='CC')
dd2 <- d2[c(1:10),]
tmp2 <- dd2%>%dplyr::select(ID,Description,GeneRatio,ONTOLOGY)%>%
  tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
  dplyr::rename(Ontology=ONTOLOGY)%>%
  dplyr::mutate(Enriched = as.numeric(Enriched),
                Detected = as.numeric(Detected),
                Ontology = ifelse(Ontology=='BP','Biological Process',
                                  ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
  dplyr::arrange(Ontology) %>%
  dplyr::mutate(ID = factor(ID, rev(unique(ID))))

d3 <- GO_pos@result%>%filter(ONTOLOGY=='MF')
b <- c("GO:0140272","GO:0050839","GO:0005509","GO:0019900","GO:0019901","GO:0001618","GO:0046872","GO:0043167","GO:0043169","GO:0045296")
dd3 <- d3[b,]
tmp3 <- dd3%>%dplyr::select(ID,Description,GeneRatio,ONTOLOGY)%>%
  tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
  dplyr::rename(Ontology=ONTOLOGY)%>%
  dplyr::mutate(Enriched = as.numeric(Enriched),
                Detected = as.numeric(Detected),
                Ontology = ifelse(Ontology=='BP','Biological Process',
                                  ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
  dplyr::arrange(Ontology) %>%
  dplyr::distinct(Description,.keep_all = T)%>%
  dplyr::mutate(ID = factor(ID, rev(unique(ID))))

tmp <- rbind(tmp1,tmp2,tmp3)%>%dplyr::distinct(ID,.keep_all = T)
tmp$ID <- factor(tmp$ID,rev(tmp$ID))

rm(d1,d2,d3,dd1,dd2,dd3,tmp1,tmp2,tmp3)
gc()

tmp_l1 <- data.frame(nrow(tmp) - (cumsum(table(tmp$Ontology)) - table(tmp$Ontology)/2))
tmp_l1$start <- tmp_l1$Freq - table(tmp$Ontology)/2
tmp_l1$end <- tmp_l1$Freq + table(tmp$Ontology)/2

m1 <- ifelse(log(max(tmp$Detected),base = 20)*1000<1000,1000,log(max(tmp$Detected),base = 20)*1000)
m2 <- ifelse(log(max(tmp$Enriched),base = 20)*1000<1000,1000,log(max(tmp$Enriched),base = 20)*1000)

p1 <- ggplot(tmp) +
  geom_col(mapping = aes(ID, Detected, fill= Ontology),
           color = "black", fill = "#e6b8a2",
           width = 0.75, 
           show.legend = F) +
  geom_text(mapping = aes(ID, Detected, label = Detected),hjust=-0.3, size = 2.5) +
  scale_y_log10(limits = c(1, m1),expand = c(0,0),breaks=c(1,100),labels=c(1,100)) +
  coord_flip() + 
  theme_classic() +
  labs(x = NULL, y = NULL, title = "Detected Genes") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,2,1,0), "mm"),
        axis.text.x = element_text(size=7),
        plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',color='black',vjust=-1))

p2 <- ggplot(tmp) +
  geom_col(mapping = aes(ID, Enriched, fill= Ontology),color = "black", width = 0.75, show.legend = F) +
  geom_text(mapping = aes(ID, Enriched, label = Enriched),hjust=-0.3, size = 2.5) +
  scale_y_log10(limits = c(1, m2),expand = expansion(),breaks=c(1,100,1000),labels=c(1,100,1000)) +
  scale_fill_manual(values = alpha(c('#00b4d8','#ff477e','#FFA500'),0.9)) +
  coord_flip() + 
  theme_classic() +
  labs(x = NULL, y = NULL, title = "Enriched Genes",fill=NULL) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size=7),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',vjust=-1))

p0 <- ggplot(tmp) + 
  geom_text(mapping = aes(ID, 0, label = ID ,color = Ontology),size = 3, show.legend = F, hjust = 0) +
  scale_color_manual(values = alpha(c('#00b4d8','#ff477e','#FFA500'),0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.02)) +
  coord_flip() + theme_void() +
  labs(x = NULL, y = NULL, title = "Identifiers") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,0,1,3), "mm"),
        plot.title = element_text(hjust = 0.01, size = 10,face = 'bold',vjust=1.5))

p3 <- ggplot(tmp) +
  geom_text(mapping = aes(ID, 0, label = Description,color = Ontology), 
            size = 3, show.legend = F, hjust = 0) +
  scale_color_manual(values = alpha(c('#00b4d8','#ff477e','#FFA500'),0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.1)) +
  coord_flip() + 
  theme_void() +
  labs(x = NULL, y = NULL, title = "Description") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,0.2,1,1), "mm"),
        plot.title = element_text(hjust = 0, size = 10,face = 'bold',vjust=1.5))

p4 <- ggplot(tmp_l1) +
  geom_segment(mapping = aes(x = start+0.1, xend = end-0.1, y = -0.1, yend = -0.1), size = 2)+
  geom_text(mapping = aes(Freq, 0, label = Var1), size = 3, show.legend = F, hjust = 0) +
  scale_y_continuous(expand = expansion(), limits = c(-0.1,1)) +
  scale_x_continuous(expand = expansion()) +
  coord_flip() + 
  theme_void() +
  labs(x = NULL, y = NULL, title = "Ontoloty") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,1,1,0.2), "mm"),
        plot.title = element_text(hjust = 0.1,size = 10,face = 'bold',vjust=1.5))
cowplot::plot_grid(p0,p1,p2,p3,p4,align = "h", nrow = 1, 
                   rel_widths = c(0.1,0.15,0.15,0.35,0.2))

#GO analysis of C3 subtype---------------------------------------------------------
rm(list = ls())
library( org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
load("RData/04.validation/signature_gene.rdata")
cordata3 <- gene[gene$class=="C3",]
ID <- bitr(cordata3$gene_name,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
GO_pos <- enrichGO(gene = ID$ENTREZID,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = 'all',
                   pvalueCutoff = 0.05,qvalueCutoff = 0.05,minGSSize = 1,maxGSSize = 50000) 

d1 <- GO_pos@result%>%filter(ONTOLOGY=='BP')
dd1 <- d1[(1:10),]
tmp1 <- dd1%>%dplyr::select(ID,Description,GeneRatio,ONTOLOGY)%>%
  tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
  dplyr::rename(Ontology=ONTOLOGY)%>%
  dplyr::mutate(Enriched = as.numeric(Enriched),
                Detected = as.numeric(Detected),
                Ontology = ifelse(Ontology=='BP','Biological Process',
                                  ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
  dplyr::arrange(Ontology) %>%
  dplyr::mutate(ID = factor(ID, rev(unique(ID))))

d2 <- GO_pos@result%>%filter(ONTOLOGY=='CC')
dd2 <- d2[c(1:10),]
tmp2 <- dd2%>%dplyr::select(ID,Description,GeneRatio,ONTOLOGY)%>%
  tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
  dplyr::rename(Ontology=ONTOLOGY)%>%
  dplyr::mutate(Enriched = as.numeric(Enriched),
                Detected = as.numeric(Detected),
                Ontology = ifelse(Ontology=='BP','Biological Process',
                                  ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
  dplyr::arrange(Ontology) %>%
  dplyr::mutate(ID = factor(ID, rev(unique(ID))))

d3 <- GO_pos@result%>%filter(ONTOLOGY=='MF')
dd3 <- d3[c(1:10),]
tmp3 <- dd3%>%dplyr::select(ID,Description,GeneRatio,ONTOLOGY)%>%
  tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
  dplyr::rename(Ontology=ONTOLOGY)%>%
  dplyr::mutate(Enriched = as.numeric(Enriched),
                Detected = as.numeric(Detected),
                Ontology = ifelse(Ontology=='BP','Biological Process',
                                  ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
  dplyr::arrange(Ontology) %>%
  dplyr::distinct(Description,.keep_all = T)%>%
  dplyr::mutate(ID = factor(ID, rev(unique(ID))))

tmp <- rbind(tmp1,tmp2,tmp3)%>%dplyr::distinct(ID,.keep_all = T)
tmp$ID <- factor(tmp$ID,rev(tmp$ID))

rm(d1,d2,d3,dd1,dd2,dd3,tmp1,tmp2,tmp3)
gc()

tmp_l1 <- data.frame(nrow(tmp) - (cumsum(table(tmp$Ontology)) - table(tmp$Ontology)/2))
tmp_l1$start <- tmp_l1$Freq - table(tmp$Ontology)/2
tmp_l1$end <- tmp_l1$Freq + table(tmp$Ontology)/2

m1 <- ifelse(log(max(tmp$Detected),base = 20)*1000<1000,1000,log(max(tmp$Detected),base = 20)*1000)
m2 <- ifelse(log(max(tmp$Enriched),base = 20)*1000<1000,1000,log(max(tmp$Enriched),base = 20)*1000)

p1 <- ggplot(tmp) +
  geom_col(mapping = aes(ID, Detected, fill= Ontology),
           color = "black", fill = "#e6b8a2",
           width = 0.75, 
           show.legend = F) +
  geom_text(mapping = aes(ID, Detected, label = Detected),hjust=-0.3, size = 2.5) +
  scale_y_log10(limits = c(1, m1),expand = c(0,0),breaks=c(1,100),labels=c(1,100)) +
  coord_flip() + 
  theme_classic() +
  labs(x = NULL, y = NULL, title = "Detected Genes") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,2,1,0), "mm"),
        axis.text.x = element_text(size=7),
        plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',color='black',vjust=-1))

p2 <- ggplot(tmp) +
  geom_col(mapping = aes(ID, Enriched, fill= Ontology),color = "black", width = 0.75, show.legend = F) +
  geom_text(mapping = aes(ID, Enriched, label = Enriched),hjust=-0.3, size = 2.5) +
  scale_y_log10(limits = c(1, m2),expand = expansion(),breaks=c(1,100,1000),labels=c(1,100,1000)) +
  scale_fill_manual(values = alpha(c('#00b4d8','#ff477e','#FFA500'),0.9)) +
  coord_flip() + 
  theme_classic() +
  labs(x = NULL, y = NULL, title = "Enriched Genes",fill=NULL) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size=7),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',vjust=-1))

p0 <- ggplot(tmp) + 
  geom_text(mapping = aes(ID, 0, label = ID ,color = Ontology),size = 3, show.legend = F, hjust = 0) +
  scale_color_manual(values = alpha(c('#00b4d8','#ff477e','#FFA500'),0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.02)) +
  coord_flip() + theme_void() +
  labs(x = NULL, y = NULL, title = "Identifiers") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,0,1,3), "mm"),
        plot.title = element_text(hjust = 0.01, size = 10,face = 'bold',vjust=1.5))

p3 <- ggplot(tmp) +
  geom_text(mapping = aes(ID, 0, label = Description,color = Ontology), 
            size = 3, show.legend = F, hjust = 0) +
  scale_color_manual(values = alpha(c('#00b4d8','#ff477e','#FFA500'),0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.1)) +
  coord_flip() + 
  theme_void() +
  labs(x = NULL, y = NULL, title = "Description") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,0.2,1,1), "mm"),
        plot.title = element_text(hjust = 0, size = 10,face = 'bold',vjust=1.5))

p4 <- ggplot(tmp_l1) +
  geom_segment(mapping = aes(x = start+0.1, xend = end-0.1, y = -0.1, yend = -0.1), size = 2)+
  geom_text(mapping = aes(Freq, 0, label = Var1), size = 3, show.legend = F, hjust = 0) +
  scale_y_continuous(expand = expansion(), limits = c(-0.1,1)) +
  scale_x_continuous(expand = expansion()) +
  coord_flip() + 
  theme_void() +
  labs(x = NULL, y = NULL, title = "Ontoloty") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,1,1,0.2), "mm"),
        plot.title = element_text(hjust = 0.1,size = 10,face = 'bold',vjust=1.5))
cowplot::plot_grid(p0,p1,p2,p3,p4,align = "h", nrow = 1, 
                   rel_widths = c(0.1,0.15,0.15,0.35,0.2))

#GO analysis of C1 subtype--------------------------------------------------------------
rm(list = ls())
library( org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
load("RData/04.validation/signature_gene.rdata")
cordata4 <- gene[gene$class=="C4",]


ID <- bitr(cordata4$gene_name,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
GO_pos <- enrichGO(gene = ID$ENTREZID,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = 'all',
                   pvalueCutoff = 0.05,qvalueCutoff = 0.05,minGSSize = 1,maxGSSize = 50000) 

d1 <- GO_pos@result%>%filter(ONTOLOGY=='BP')
dd1 <- d1[(1:10),]
tmp1 <- dd1%>%dplyr::select(ID,Description,GeneRatio,ONTOLOGY)%>%
  tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
  dplyr::rename(Ontology=ONTOLOGY)%>%
  dplyr::mutate(Enriched = as.numeric(Enriched),
                Detected = as.numeric(Detected),
                Ontology = ifelse(Ontology=='BP','Biological Process',
                                  ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
  dplyr::arrange(Ontology) %>%
  dplyr::mutate(ID = factor(ID, rev(unique(ID))))

d2 <- GO_pos@result%>%filter(ONTOLOGY=='CC')
dd2 <- d2[c(1:10),]
tmp2 <- dd2%>%dplyr::select(ID,Description,GeneRatio,ONTOLOGY)%>%
  tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
  dplyr::rename(Ontology=ONTOLOGY)%>%
  dplyr::mutate(Enriched = as.numeric(Enriched),
                Detected = as.numeric(Detected),
                Ontology = ifelse(Ontology=='BP','Biological Process',
                                  ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
  dplyr::arrange(Ontology) %>%
  dplyr::mutate(ID = factor(ID, rev(unique(ID))))

d3 <- GO_pos@result%>%filter(ONTOLOGY=='MF')
dd3 <- d3[c(1:10),]
tmp3 <- dd3%>%dplyr::select(ID,Description,GeneRatio,ONTOLOGY)%>%
  tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
  dplyr::rename(Ontology=ONTOLOGY)%>%
  dplyr::mutate(Enriched = as.numeric(Enriched),
                Detected = as.numeric(Detected),
                Ontology = ifelse(Ontology=='BP','Biological Process',
                                  ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
  dplyr::arrange(Ontology) %>%
  dplyr::distinct(Description,.keep_all = T)%>%
  dplyr::mutate(ID = factor(ID, rev(unique(ID))))

tmp <- rbind(tmp1,tmp2,tmp3)%>%dplyr::distinct(ID,.keep_all = T)
tmp$ID <- factor(tmp$ID,rev(tmp$ID))

rm(d1,d2,d3,dd1,dd2,dd3,tmp1,tmp2,tmp3)
gc()

tmp_l1 <- data.frame(nrow(tmp) - (cumsum(table(tmp$Ontology)) - table(tmp$Ontology)/2))
tmp_l1$start <- tmp_l1$Freq - table(tmp$Ontology)/2
tmp_l1$end <- tmp_l1$Freq + table(tmp$Ontology)/2

m1 <- ifelse(log(max(tmp$Detected),base = 20)*1000<1000,1000,log(max(tmp$Detected),base = 20)*1000)
m2 <- ifelse(log(max(tmp$Enriched),base = 20)*1000<1000,1000,log(max(tmp$Enriched),base = 20)*1000)

p1 <- ggplot(tmp) +
  geom_col(mapping = aes(ID, Detected, fill= Ontology),
           color = "black", fill = "#e6b8a2",
           width = 0.75, 
           show.legend = F) +
  geom_text(mapping = aes(ID, Detected, label = Detected),hjust=-0.3, size = 2.5) +
  scale_y_log10(limits = c(1, m1),expand = c(0,0),breaks=c(1,100),labels=c(1,100)) +
  coord_flip() + 
  theme_classic() +
  labs(x = NULL, y = NULL, title = "Detected Genes") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,2,1,0), "mm"),
        axis.text.x = element_text(size=7),
        plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',color='black',vjust=-1))

p2 <- ggplot(tmp) +
  geom_col(mapping = aes(ID, Enriched, fill= Ontology),color = "black", width = 0.75, show.legend = F) +
  geom_text(mapping = aes(ID, Enriched, label = Enriched),hjust=-0.3, size = 2.5) +
  scale_y_log10(limits = c(1, m2),expand = expansion(),breaks=c(1,100,1000),labels=c(1,100,1000)) +
  scale_fill_manual(values = alpha(c('#00b4d8','#ff477e','#FFA500'),0.9)) +
  coord_flip() + 
  theme_classic() +
  labs(x = NULL, y = NULL, title = "Enriched Genes",fill=NULL) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size=7),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',vjust=-1))

p0 <- ggplot(tmp) + 
  geom_text(mapping = aes(ID, 0, label = ID ,color = Ontology),size = 3, show.legend = F, hjust = 0) +
  scale_color_manual(values = alpha(c('#00b4d8','#ff477e','#FFA500'),0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.02)) +
  coord_flip() + theme_void() +
  labs(x = NULL, y = NULL, title = "Identifiers") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,0,1,3), "mm"),
        plot.title = element_text(hjust = 0.01, size = 10,face = 'bold',vjust=1.5))

p3 <- ggplot(tmp) +
  geom_text(mapping = aes(ID, 0, label = Description,color = Ontology), 
            size = 3, show.legend = F, hjust = 0) +
  scale_color_manual(values = alpha(c('#00b4d8','#ff477e','#FFA500'),0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.1)) +
  coord_flip() + 
  theme_void() +
  labs(x = NULL, y = NULL, title = "Description") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,0.2,1,1), "mm"),
        plot.title = element_text(hjust = 0, size = 10,face = 'bold',vjust=1.5))

p4 <- ggplot(tmp_l1) +
  geom_segment(mapping = aes(x = start+0.1, xend = end-0.1, y = -0.1, yend = -0.1), size = 2)+
  geom_text(mapping = aes(Freq, 0, label = Var1), size = 3, show.legend = F, hjust = 0) +
  scale_y_continuous(expand = expansion(), limits = c(-0.1,1)) +
  scale_x_continuous(expand = expansion()) +
  coord_flip() + 
  theme_void() +
  labs(x = NULL, y = NULL, title = "Ontoloty") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,1,1,0.2), "mm"),
        plot.title = element_text(hjust = 0.1,size = 10,face = 'bold',vjust=1.5))
cowplot::plot_grid(p0,p1,p2,p3,p4,align = "h", nrow = 1, 
                   rel_widths = c(0.1,0.15,0.15,0.35,0.2))

#Performing GSVA analysis---------------------------
rm(list = ls())
library(tidyverse)
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)

load("Database/Meta_RNAseq_expr.rda")
load('RData/phenotype_RNAseq.rda')
load("Database/GRCh38_v39.gtf.rda")
exp <- merge(gtf,Meta_RNAseq_expr,by.x=1,by.y=0)%>%
  dplyr::select(-gene_id)%>%
  group_by(gene_name)%>%
  summarise_all(max)%>%
  tibble::column_to_rownames('gene_name')
# Hallmarker list
library(GSEABase)
library(GSVA)
gs <- getGmt('Reference/05.function/h.all.v7.4.symbols.gmt')
#gsva
gsva_es <- gsva(as.matrix(exp), gs)
write.csv(gsva_es, "Table/05.function/gsva_output.csv", quote = F)
save(gsva_es, file = "RData/05.function/gsva_output.rda")

load("RData/05.function/gsva_output.rda")

DD1 <- clin[order(clin$Cluster),1:2]
rownames(DD1) <- DD1$Sample
res_es1 <- scale(gsva_es)   
res_es1 <- res_es1[,rownames(DD1)]
rownames(res_es1) <- tolower(rownames(res_es1))

##Drawing
library(ggsci)
library(ComplexHeatmap)

ccol <- c('#438F5D','#E19412','#2C6589',"#BC3C29FF")
names(ccol) <- c('C1','C2',"C3","C4")
Top = HeatmapAnnotation(border = T,
                        show_legend = F,
                        foo = anno_block(gp = gpar(fill = ccol),
                                         labels = c("C1", "C2", "C3","C4"), 
                                         labels_gp = gpar(col = "black", fontsize = 10),
                                         height = unit(0.5,"cm")),
                        show_annotation_name = F )
library('circlize')

library(stringr)
x <- str_split(rownames(res_es1), "_",2)
x1 <- as.data.frame(unlist(x)) 
x2 <- x1[x1$`unlist(x)`!='hallmark',]
rownames(res_es1) <- x2
rownames(res_es1) <- gsub("_"," ",rownames(res_es1))
table(DD1$Cluster)

pdf(file = 'Figure/Figure3/E_GSVA.pdf',width=11,height = 10)
Heatmap(res_es1,
        name='GSVA',
        top_annotation = Top,  
        cluster_rows = T,
        col=colorRamp2(c(-2,0,2),c('seagreen','white',"orange")),
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,row_names_side = 'right',
        column_order=NULL,
        show_column_names = FALSE,
        show_row_names = T, 
        show_row_dend = F,
        row_names_gp = gpar(fontsize = 12),
        column_split = c(rep(1,131),rep(2,150),rep(3,139),rep(4,112)),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 16),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 16), 
                                  title_gp = gpar(fontsize = 16, fontface = "bold"))) 

dev.off()


#Immunoscore and tumor purity score---------------------------------------------
rm(list = ls())
library(IOBR)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggsci)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
rm(list = ls())
load('Database/Meta_RNAseq_expr.rda')
load('RData/phenotype_RNAseq.rda')
logtpm <- Meta_RNAseq_expr
zz <- rownames(logtpm)
zz1 <- bitr(zz, fromType = "ENSEMBL",toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
a <- !duplicated(zz1$SYMBOL)
zz1 <- zz1[a,]
logtpm <- logtpm[zz1$ENSEMBL,]
rownames(logtpm) <- zz1$SYMBOL
TCGA_mRNA=logtpm
estimate <- deconvo_tme(eset = TCGA_mRNA, method = "estimate")[,-5]
save(estimate,file = "RData/06.estimate/estimatescore.rda")

rm(list = ls())
library(IOBR)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggsci)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
load('RData/phenotype_RNAseq.rda')
load("RData/06.estimate/estimatescore.rda")
#immune score
dat <-  estimate[,c(1,3)] %>% 
  merge(.,clin[,c("Sample","Cluster"),drop=F], by.x = 1,by.y=1) %>% 
  column_to_rownames("ID")
dat$Cluster <- factor(dat$Cluster)
dat$ImmuneScore_estimate <- scale(dat$ImmuneScore_estimate)

library(ggsci)
cols <- c('#438F5D','#E19412','#2C6589',"#BC3C29FF")
ggplot(dat,aes_string("Cluster","ImmuneScore_estimate"))+
  geom_jitter(shape = 21,size=2,width = 0.2,aes_string(fill="Cluster",color="Cluster"))+
  geom_boxplot(outlier.colour = NA,aes_string(fill="Cluster"),color='black',size=0.6,alpha=0.65)+
  geom_violin(alpha=0.5,aes_string(fill="Cluster"),color=NA,trim = T)+
  stat_compare_means(label.y = max(dat[,1])*1.05,color='black',size=5,hjust = 0.1)+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  expand_limits(y=max(dat[,1])*1.1)+
  theme_bw(base_rect_size = 2)+
  labs(title ='Immune Score',x=NULL,y=NULL)+
  theme(axis.text.y = element_text(size = 11, colour = 'black'),
        axis.text.x = element_text(size = 14,colour = 'black',angle = 0),
        axis.title.x = element_text(size = 15,colour = 'darkred',face='bold'),
        axis.title.y = element_text(size = 15,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 15,colour = 'black',face='bold'),
        legend.position = 'none',
        axis.ticks.x = element_blank())
ggsave(filename = "Figure/Figure4/A_immunescore.pdf",width = 4.5,height = 4.5)


#tumorpurity
TumorPurity = cos(0.6049872018+0.0001467884 * estimate[,3])
lj <- as.data.frame(TumorPurity)
out <- merge(estimate,lj,by=0)
dat <-  out[,c(2,6)] %>% 
  merge(.,clin[,c("Sample","Cluster"),drop=F], by.x = 1,by.y=1) %>% 
  column_to_rownames("ID")
dat$Cluster <- factor(dat$Cluster)
colnames(dat)[1] <- "TumorPurity"
dat$TumorPurity <- scale(dat$TumorPurity)
colnames(dat)[1] <- "TumorPurity"

library(ggsci)
cols <- c('#438F5D','#E19412','#2C6589',"#BC3C29FF")
ggplot(dat,aes_string("Cluster","TumorPurity"))+
  geom_jitter(shape = 21,size=2,width = 0.2,aes_string(fill="Cluster",color="Cluster"))+
  geom_boxplot(outlier.colour = NA,aes_string(fill="Cluster"),color='black',size=0.6,alpha=0.65)+
  geom_violin(alpha=0.5,aes_string(fill="Cluster"),color=NA,trim = T)+
  stat_compare_means(label.y = max(dat[,1])*1.05,color='black',size=5,hjust = 0.1)+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  expand_limits(y=max(dat[,1])*1.1)+
  theme_bw(base_rect_size = 2)+
  labs(title ='Tumor Purity',x=NULL,y=NULL)+
  theme(axis.text.y = element_text(size = 11, colour = 'black'),
        axis.text.x = element_text(size = 14,colour = 'black',angle = 0),
        axis.title.x = element_text(size = 15,colour = 'darkred',face='bold'),
        axis.title.y = element_text(size = 15,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 15,colour = 'black',face='bold'),
        legend.position = 'none',
        axis.ticks.x = element_blank())
ggsave(filename = "Figure/Figure3/F_TumorPurity.pdf",width = 4.5,height = 4.5)


#Proliferation score----------------------------------------------
rm(list = ls())
library(ggplot2)
library(data.table)
library(tidyverse)
options(digits = 3)
load("Reference/06.score/immune_landscape.RData")
ABSOLUTE_scores$ID <- substr(rownames(ABSOLUTE_scores),1,12)
ABSOLUTE_scores <- ABSOLUTE_scores[!duplicated(ABSOLUTE_scores$ID),-5]
rownames(ABSOLUTE_scores) <- substr(rownames(ABSOLUTE_scores),1,12)

GBM_im <- immune_landscape[immune_landscape$TCGA.Study == "GBM",]
GBM_ab <- ABSOLUTE_scores[rownames(GBM_im),]
GBM_il<- merge(GBM_im,GBM_ab,by=0) %>%  column_to_rownames("Row.names")
rm(immune_landscape,ABSOLUTE_scores,GBM_im,GBM_ab)

load("RData/phenotype_RNAseq.rda")
Cluster <- clin[,1:2,drop=F]
Cluster <- Cluster[389:532,]
Cluster$Sample <- substr(Cluster$Sample,1,12)
metadata <- merge(Cluster,GBM_il,by.x = 1,by.y = 0)%>% 
  column_to_rownames("Sample") 
ggdata <- metadata[,-c(2:4,8)] 
ggdata$Cluster <- substr(ggdata$Cluster,2,2)
ggdata <- apply(ggdata[,-2],2,as.numeric) %>% 
  as.data.frame()
ggdata$Cluster <- factor(ggdata$Cluster)

library(compareGroups)
library(ggpubr)

library(ggsci)
mm <- ggdata[,c(1,4)]
mm$Cluster <- paste("C",mm$Cluster)
mm$Cluster <- gsub(" ","",mm$Cluster)
cols <- c('#438F5D','#E19412','#2C6589',"#BC3C29FF")
ggplot(mm,aes_string("Cluster","Proliferation"))+
  geom_jitter(shape = 21,size=2,width = 0.2,aes_string(fill="Cluster",color="Cluster"))+
  geom_boxplot(outlier.colour = NA,aes_string(fill="Cluster"),color='black',size=0.6,alpha=0.65)+
  geom_violin(alpha=0.5,aes_string(fill="Cluster"),color=NA,trim = T)+
  stat_compare_means(label.y = max(mm[,2])*1.15,color='black',size=5,hjust = 0.1)+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  expand_limits(y=max(mm[,2])*1.2)+
  theme_bw(base_rect_size = 2)+
  labs(title ='Proliferation',x=NULL,y=NULL)+
  theme(axis.text.y = element_text(size = 11, colour = 'black'),
        axis.text.x = element_text(size = 14,colour = 'black',angle = 0),
        axis.title.x = element_text(size = 15,colour = 'darkred',face='bold'),
        axis.title.y = element_text(size = 15,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 15,colour = 'black',face='bold'),
        legend.position = 'none',
        axis.ticks.x = element_blank())
ggsave(filename = "Figure/Figure3/G_proliferation.pdf",width = 4.5,height = 4.5)

#
rm(list=ls())
library(GSVA)
library(limma)
library(GSEABase)
library(magrittr)
library(dplyr)

load("Database/Meta_RNAseq_expr.rda")
load("Database/GRCh38_v39.gtf.rda")
Meta_RNAseq_expr <- merge(gtf,Meta_RNAseq_expr,by.x=1,by.y=0)%>%
  dplyr::select(-gene_id)%>%
  group_by(gene_name)%>%
  summarise_all(max)%>%
  tibble::column_to_rownames('gene_name')
load("Reference/07.immune/cellMarker_ssGSEA.Rdata")
ssGSEA_Score=gsva(as.matrix(Meta_RNAseq_expr), as.list(cellMarker), method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)#ssGSEA计算

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
norm_ssGSEA_Score=normalize(ssGSEA_Score)
norm_ssGSEA_Score=rbind(id=colnames(norm_ssGSEA_Score),norm_ssGSEA_Score)
ssGSEA=norm_ssGSEA_Score
class(ssGSEA)
write.table(norm_ssGSEA_Score,file="Table/07.immune/ICP_ssGSEA_output.txt",sep="\t",quote=F,col.names=F)#此处的输出文件即为ssGSEA富集得分文件

#Comparison of immune cells and immune checkpoint expression in four subtypes-----------------
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(data.table)
library(tibble)
library(readxl)
library(glmnet)
library(timeROC)
library(survival)
library(survminer)
library(ggsci)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyr)

load('RData/phenotype_RNAseq.rda')
Cluster <- clin[,1:2,drop=F]

#immune cells
ssGSEA_output <- fread("Table/07.immune/ICP_ssGSEA_output.txt",sep = "\t")
ssGSEA_output <- as.data.frame(ssGSEA_output)
rownames(ssGSEA_output) <- ssGSEA_output[,1]
ssGSEA_output <- ssGSEA_output[,-1]
ssGSEA_output <- t(ssGSEA_output)
colnames(ssGSEA_output)
ic <- ssGSEA_output[,c("Macrophage","Neutrophil",
                       "Natural killer cell","Activated CD4 T cell",
                       "Activated CD8 T cell","Gamma delta T cell")]
dd <- merge(Cluster,ic,by.x = 1,by.y =0) %>%column_to_rownames("Sample")

#immune checkpoint
load("Database/Meta_RNAseq_expr.rda")
exp <- Meta_RNAseq_expr
load("Database/GRCh38_v39.gtf.rda")
exp <- merge(gtf,exp,by.x=1,by.y=0)%>%
  dplyr::select(-gene_id)%>%
  group_by(gene_name)%>%
  summarise_all(max)%>%
  tibble::column_to_rownames('gene_name')
rm(Meta_RNAseq_expr,gtf)

rr2 <- clin[,1:2,drop=F] %>% 
  .[order(.$Cluster),,drop=F]
TCGA_ee <- exp
IC <- rbind(data.frame(Type='B7-CD28',
                       Gene=c('CD274','CD276','CTLA4','HHLA2','ICOS','ICOSLG','PDCD1','PDCD1LG2','TMIGD2','VTCN1')),
            data.frame(Type='TNF superfamily',
                       Gene=c('BTLA','CD27','CD40','CD40LG','CD70','TNFRSF18','TNFRSF4','TNFRSF9','TNFSF14','TNSF4','TNSF9')),
            data.frame(Type='Others',
                       Gene=c('C10ORF54','ENTPD1','FGL1','HAVCR2','IDO1','LAG3','NCR3','NT5E','SIGLEC15')))
IC$Type <- factor(IC$Type,levels = unique(IC$Type))
x <- intersect(IC$Gene,rownames(TCGA_ee))
IC <- IC[IC$Gene%in%x,]
IC <- IC[order(IC$Type,IC$Gene),]
ee <- TCGA_ee[IC$Gene,]
dd2 <- merge(dd,t(ee),by = 0)
save(dd2,file = 'Rdata/07.immune/ICP_input.Rda')
Cluster <- c('#438F5D','#E19412','#2C6589',"#BC3C29FF")
names(Cluster) <- c('C1','C2',"C3","C4")
dd2 <- dd2[order(dd2$Cluster),]
Top = HeatmapAnnotation(Cluster=dd2$Cluster,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 9),
                                                     title_gp = gpar(fontsize = 9, fontface = "bold"),
                                                     ncol=1),
                        border = T,
                        col=list(Cluster=Cluster),
                        show_annotation_name = F)
table(dd2$Cluster)
pdf(file = 'Figure/Figure4/D_ICP_heatmap.pdf',width=7.7,height = 5.8)
Heatmap(t(scale(dd2[,c(3:35)])),
        name='Z-score',
        top_annotation = Top,
        cluster_rows = FALSE,
        col=colorRamp2(c(-2,0,2),c('#82B0C7','white','#DF5525')),
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 7),
        column_split = c(rep(1,131),rep(2,150),rep(3,139),rep(4,112)),
        row_split = factor(rep(c('Immune cell','B7-CD28','TNF superfamily','Others'),
                               times=c(6,10,9,8)),levels = c('Immune cell','B7-CD28','TNF superfamily','Others')),
        
        gap = unit(1.5, "mm"),
        column_title = "Immune Checkpoint",
        column_title_gp = gpar(fontsize = 15),
        row_title_gp = gpar(fontsize=9),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 9), 
                                  title_gp = gpar(fontsize = 9, fontface = "bold"))
) 
dev.off()

#Comparison of APS in four subtypes--------------------
rm(list = ls())
library(xlsx)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggsci)

imsig <- read.xlsx("Reference/07.immune/APS.xlsx", sheetIndex = 1,header = T)
imsig <- imsig[1:78,1:8]
ap <- imsig[imsig$Super.Category == "Antigen presentation",]
ap_gene <- ap$Gene

load("Database/Meta_RNAseq_expr.rda")
load("Database/GRCh38_v39.gtf.rda")
apExpr <- Meta_RNAseq_expr
apExpr <- merge(gtf,apExpr,by.x=1,by.y=0)%>%
  dplyr::select(-gene_id)%>%
  group_by(gene_name)%>%
  summarise_all(max)%>%
  tibble::column_to_rownames('gene_name')
rm(Meta_RNAseq_expr,gtf)
expr <- apExpr
apExpr <- apExpr[rownames(apExpr) %in% ap_gene,]

load("RData/phenotype_RNAseq.rda")
metadata <-  as.data.frame(t(apExpr)) %>% 
  merge(.,clin[,c("Sample","Cluster"),drop=F], by.x = 0,by.y=1) %>% 
  column_to_rownames("Row.names")
metadata$Cluster <- factor(metadata$Cluster)
plotinput <- pivot_longer(metadata,cols = 1:13, names_to = "signature", values_to = "value")
#Calculating APS scores via gsva
library(GSVA)
expr <-expr 
APS <- list('APS'= ap_gene)
GSVA <- gsva(expr = as.matrix(expr),
             gset.idx.list = APS,
             method="gsva")
APS <- data.frame(APS=t(GSVA))%>%rownames_to_column('ID')
save(APS,file = 'RData/07.immune/APS.Rdata')
rm(list = ls())
load("RData/phenotype_RNAseq.rda")
load('RData/07.immune/APS.Rdata')
dat <-  APS %>% 
  merge(.,clin[,c("Sample","Cluster"),drop=F], by.x = 1,by.y=1) %>% 
  column_to_rownames("ID")
dat$Cluster <- factor(dat$Cluster)

library(ggsci)
cols <- c('#438F5D','#E19412','#2C6589',"#BC3C29FF")
ggplot(dat,aes_string("Cluster","APS"))+
  geom_jitter(shape = 21,size=2,width = 0.2,aes_string(fill="Cluster",color="Cluster"))+
  geom_boxplot(outlier.colour = NA,aes_string(fill="Cluster"),color='black',size=0.6,alpha=0.65)+
  geom_violin(alpha=0.5,aes_string(fill="Cluster"),color=NA,trim = T)+
  stat_compare_means(label.y = max(dat[,1])*1.05,color='black',size=5,hjust = 0.1)+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  expand_limits(y=max(dat[,1])*1.1)+
  theme_bw(base_rect_size = 2)+
  labs(title ='APS',x=NULL,y=NULL)+
  theme(axis.text.y = element_text(size = 11, colour = 'black'),
        axis.text.x = element_text(size = 14,colour = 'black',angle = 0),
        axis.title.x = element_text(size = 15,colour = 'darkred',face='bold'),
        axis.title.y = element_text(size = 15,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 15,colour = 'black',face='bold'),
        legend.position = 'none',
        axis.ticks.x = element_blank())
ggsave(filename = "Figure/Figure4/B_APS.pdf",width = 4.5,height = 4.5)

#Comparison of TIS in four subtypes-----------------------------
rm(list = ls())
options(stringsAsFactors = F)
library(MCPcounter)
library(GSVA)
library(tibble)
library(tidyverse)
library(ggpubr)
library(ggsci)
load("Database/Meta_RNAseq_expr.rda")
expr <- Meta_RNAseq_expr
rm(Meta_RNAseq_expr)
load("Database/GRCh38_v39.gtf.rda")
expr <- merge(gtf,expr,by.x=1,by.y=0)%>%
  dplyr::select(-gene_id)%>%
  group_by(gene_name)%>%
  summarise_all(max)%>%
  tibble::column_to_rownames('gene_name')
TISG <- c("CCL5", "CD27", "CD274", "CD276", 
          "CD8A", "CMKLR1", "CXCL9", "HLA-DQA1",
          "HLA-DRB1", "HLA-E", "IDO1",
          "LAG3", "NKG7", "PDCD1LG2", 
          "PSMB10", "STAT1", "TIGIT")
GSVA <- gsva(expr=as.matrix(expr), 
             gset.idx.list=list(TISG), method="gsva")
TIS <- data.frame(TIS=t(GSVA))%>%rownames_to_column('ID')
save(TIS,file = 'RData/07.immune/TIS.Rdata')

rm(list = ls())
load("RData/phenotype_RNAseq.rda")
load('RData/07.immune/TIS.Rdata')
dat <-  TIS %>% 
  merge(.,clin, by.x=1,by.y=1) %>% 
  column_to_rownames("ID")
dat <- dat[,1:2]
dat$Cluster <- factor(dat$Cluster)
dat %>%
  ggplot(aes(Cluster,TIS)) +
  geom_boxplot(aes(color=Cluster), outlier.colour = NA, size=1, fill=NA) +
  geom_jitter(aes(fill=Cluster,color=Cluster), width = 0.2, shape=21, size=3, alpha=0.7) +
  stat_compare_means( label = 'p.format', size=5,label.x = 2.4,label.y = 0.9) +
  theme_bw(base_rect_size = 1.5) +
  labs(x=NULL,y='Relative TIS', title = 'TIS score') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5,size = 15,colour = 'black',face = "bold"),
        axis.title.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'black'),
        axis.text.y = element_text(size = 12,colour = 'black'),
        axis.ticks = element_line(size = 1),
        panel.grid = element_blank(),#
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6')) +
  scale_color_manual(values = c('#438F5D','#E19412','#2C6589',"#BC3C29FF")) +
  scale_fill_manual(values = c('#438F5D','#E19412','#2C6589',"#BC3C29FF"))
ggsave(filename = 'Figure/Figure4/F_TIS.pdf',height = 4.5,width = 4.5)


#Comparison of MHC expression in four subtypes-----------------------
rm(list = ls())
library(data.table)
library(xlsx)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(compareGroups)
library(reshape2)
mm <- read.xlsx2("Reference/07.immune/immune_related_molecular.xlsx",sheetIndex = 1)
load("Database/Meta_RNAseq_expr.rda")
load("RData/phenotype_RNAseq.rda")
load("Database/GRCh38_v39.gtf.rda")
df <- Meta_RNAseq_expr
df <- merge(gtf,df,by.x=1,by.y=0)%>%
  dplyr::select(-gene_id)%>%
  group_by(gene_name)%>%
  summarise_all(max)%>%
  tibble::column_to_rownames('gene_name')
Cluster <- clin[,1:2,drop=F]
rm(Meta_RNAseq_expr,gtf,clin)

gene <- intersect(rownames(df),mm$GENE.SYMBLE)
mm <- mm[mm$GENE.SYMBLE %in% gene,]
sp <- intersect(colnames(df),Cluster$Sample)
mm_df <- df[gene,sp] %>% 
  t()%>%
  as.data.frame() %>% 
  merge(Cluster,.,by.x = 1,by.y = 0)%>% 
  column_to_rownames("Sample") 


mm_df$Cluster <- factor(mm_df$Cluster)
dd <- descrTable(Cluster~.,data = mm_df)
print(dd)
pvals <- getResults(dd, "p.overall")
p.adjust(pvals, method = "BH")


gene <-  colnames(mm_df)[which(pvals<0.01)+1]
mm1 <- mm[mm$GENE.SYMBLE %in% gene,]
mm_df1 <- cbind(Cluster=mm_df$Cluster,mm_df[,gene])
ggdata <- melt(mm_df1,id.vars = "Cluster",
               variable.name = "Gene",value.name = "expr") %>% 
  merge(.,mm1,by.x = 2,by.y = 1,all = T) 
ggdata <- ggdata[order(ggdata$Type,ggdata$Gene),]
ggdata$Gene <- factor(ggdata$Gene)
ggdata_mhc <- ggdata[ggdata$Type %in% c("MHC Class II" ,"MHC Class I"),]

c = ggdata_mhc %>% 
  group_by(Gene) %>% 
  summarise(m = median(expr)) %>% 
  arrange(desc(m)) %>% 
  pull(Gene)
ggdata_mhc$Gene = factor(ggdata_mhc$Gene,levels = c)

gg_mhc <- ggplot(ggdata_mhc,aes(Gene,expr,fill=Cluster))+
  geom_boxplot(outlier.size = 0.8)+
  scale_fill_manual(values = c('#438F5D','#E19412','#2C6589',"#BC3C29FF"))+
  theme_bw()+
  labs(title="MHC molecules",x=NULL,y="Expression")+
  theme(axis.text = element_text(size = 8),
    axis.title.y = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust = 0.9),
    plot.title = element_text(size = 12,face = "bold",hjust = 0.5,vjust=0.8),
    panel.grid = element_blank(),#
    panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),#
    panel.background = element_rect(fill='#f3f6f6'),
    legend.title = element_text(size = 12,face = "bold"))+
  stat_compare_means(aes(label=..p.signif..))
gg_mhc
ggsave(gg_mhc,filename = "Figure/Figure4/C_mhc ICTs.pdf",width = 7.5,height = 3.5)

#Comparison of the correlation of four subtypes with immune cells and CIC--------------
rm(list = ls())
library(data.table)
library(GSVA)
library(ggplot2)
library(ggcor)
library(dplyr)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)

gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}


load("Database/Meta_RNAseq_expr.rda")
load("Database/GRCh38_v39.gtf.rda")
expr <- merge(gtf,Meta_RNAseq_expr,by.x=1,by.y=0)%>%
  dplyr::select(-gene_id)%>%
  group_by(gene_name)%>%
  summarise_all(max)%>%
  tibble::column_to_rownames('gene_name')
rm(Meta_RNAseq_expr)
# Loading immunotherapy-predicted pathways
load('Reference/07.immune/CIC_Marker.rda')
names(CIC_Marker) <- c(paste("step1:","Release of cancer cell antigens"),
                       paste("step2:","Cancer antigen presentation"),
                       paste("step3:","Priming and activation"),
                       paste("step4:","Basophil recruiting"),
                       paste("step4:","B cell recruiting"),
                       paste("step4:","CD4 T cell recruiting"),
                       paste("step4:","CD8 T cell recruiting"),
                       paste("step4:","Dendritic cell recruiting"),
                       paste("step4:","Eosinophil recruiting"),
                       paste("step4:","Macrophage recruiting"),
                       paste("step4:","MDSC recruiting"),
                       paste("step4:","Monocyte recruiting"),
                       paste("step4:","Neutrophil recruiting"),
                       paste("step4:","NK cell recruiting"),
                       paste("step4:","T cell recruiting"),
                       paste("step4:","Th1 cell recruiting"),
                       paste("step4:","Th22 cell recruiting"),
                       paste("step4:","Th2 cell recruiting"),
                       paste("step4:","Treg cell recruiting"),
                       paste("step5:","Infiltration of immune cells into tumors"),
                       paste("step6:","Recognition of cancer cells by T cells"),
                       paste("step7:","Killing of cancer cells"))

immPath.score <- gsva(expr = as.matrix(expr),
                      CIC_Marker, 
                      method = "ssgsea")

CIC.score <- immPath.score

write.table(CIC.score, "Table/07.immune/butter_easy_input_CIC.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

load("RData/phenotype_RNAseq.rda")
a <- as.data.frame(t(CIC.score))

clin <- clin[,c(1,2)]
clin <- merge(clin,a,by.x = 1,by.y = 0)

rownames(clin) <- clin$Sample
clin$Subtype1 <- ifelse(clin$Cluster=='C1',1,0)
clin$Subtype2 <- ifelse(clin$Cluster=='C2',1,0)
clin$Subtype3 <- ifelse(clin$Cluster=='C3',1,0)
clin$Subtype4 <- ifelse(clin$Cluster=='C4',1,0)
clin <- clin[,-c(1,2)]
clin <- as.data.frame(t(clin))


outTab=data.frame()
corFilter=0             
pvalueFilter=1    
for(i in row.names(clin)){
  for(j in row.names(clin)){
    x=as.numeric(clin[i,])
    y=as.numeric(clin[j,])
    corT=cor.test(x,y)
    cor=corT$estimate
    pvalue=corT$p.value
    if((cor>corFilter) & (pvalue<pvalueFilter)){
      outTab=rbind(outTab,cbind(d1=j,d2=i,cor,pvalue,Regulation="postive"))
    }
    if((cor< -corFilter) & (pvalue<pvalueFilter)){
      outTab=rbind(outTab,cbind(d1=j,d2=i,cor,pvalue,Regulation="negative"))
    }
  }
}

table(outTab$d1)

sub <- outTab[outTab$d1%in%c('Subtype1','Subtype2','Subtype3','Subtype4'),]
sub <- sub[order(sub$d1,decreasing = F),]
sub <- sub[!sub$d2%in%c('Subtype1','Subtype2','Subtype3','Subtype4'),]
colnames(sub)[3:4] <- c('r','p.value')
sub$p.value <- as.numeric(sub$p.value)
sub$pd <- ifelse(sub$p.value<0.05,'< 0.05','>= 0.05')
sub$r <- as.numeric(sub$r)
sub$rd <- cut(sub$r, breaks = c(-Inf,-0.4,-0.2,0.2,0.4, Inf),
              labels = c("<= -0.4",'-0.4 - -0.2','-0.2 - 0.2', "0.2 - 0.4", ">= 0.4"))
table(sub$rd)
table(sub$d1)
sub$d1[sub$d1=='Subtype1'] <- 'C1'
sub$d1[sub$d1=='Subtype2'] <- 'C2'
sub$d1[sub$d1=='Subtype3'] <- 'C3'
sub$d1[sub$d1=='Subtype4'] <- 'C4'
str(sub)
varechem <- as.data.frame(a)
save(varechem,file = "RData/07.immune/butter_varechem.rdata")

quickcor(varechem, type = "upper") + 
  geom_colour() +
  anno_link(aes(colour = rd, size = pd), data = sub) +
  scale_size_manual(values = c(0.4, 0.2)) +
  scale_color_manual(values = c("#26B1D0","#DA6003","#A2A2A288","#E8288E","#65A818")) +
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#E11953",midpoint=0) +
  remove_axis("x")+
  guides(colour = guide_legend(title = "Subtype's r", 
                               override.aes = list(size = 3), 
                               order = 1),
         size = guide_legend(title = "Subtype's p",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         fill = guide_colorbar(title = "Pearson's r", order =3))

ggsave(filename = "Figure/Figure4/E_ggcor plot in top right.pdf", width = 10,height = 8)

rm(list = ls())
library(data.table)
library(GSVA)
library(ggplot2)
library(ggcor)
library(dplyr)
library(readxl)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}

load("Database/Meta_RNAseq_expr.rda")
load("Database/GRCh38_v39.gtf.rda")
expr <- merge(gtf,Meta_RNAseq_expr,by.x=1,by.y=0)%>%
  dplyr::select(-gene_id)%>%
  group_by(gene_name)%>%
  summarise_all(max)%>%
  tibble::column_to_rownames('gene_name')
rm(Meta_RNAseq_expr)

immPath <- read.csv("Reference/07.immune/ICP.txt",sep = "\t")
a <- list()
for (i in unique(immPath$Cell.type)) {
  tmp <- immPath[immPath$Cell.type==i,]
  a[[i]] <- tmp
}

b <- list()
for (i in names(a)) {
  b[[i]] <- as.list(a[[i]])
  
}

immPath.list= list()
for (i in names(b)) {
  immPath.list[[i]] <- b[[i]]$Metagene
  
}

ICP.score <- gsva(expr = as.matrix(expr),
                  immPath.list, 
                  method = "ssgsea")


write.table(ICP.score, "Table/07.immune/butter_easy_input_ICP.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

load("RData/phenotype_RNAseq.rda")
a <- as.data.frame(t(ICP.score))

clin <- clin[,c(1,2)]
clin <- merge(clin,a,by.x = 1,by.y = 0)

rownames(clin) <- clin$Sample
clin$Subtype1 <- ifelse(clin$Cluster=='C1',1,0)
clin$Subtype2 <- ifelse(clin$Cluster=='C2',1,0)
clin$Subtype3 <- ifelse(clin$Cluster=='C3',1,0)
clin$Subtype4 <- ifelse(clin$Cluster=='C4',1,0)
clin <- clin[,-c(1,2)]
clin <- as.data.frame(t(clin))


outTab=data.frame()
corFilter=0             
pvalueFilter=1    
for(i in row.names(clin)){
  for(j in row.names(clin)){
    x=as.numeric(clin[i,])
    y=as.numeric(clin[j,])
    corT=cor.test(x,y)
    cor=corT$estimate
    pvalue=corT$p.value
    if((cor>corFilter) & (pvalue<pvalueFilter)){
      outTab=rbind(outTab,cbind(d1=j,d2=i,cor,pvalue,Regulation="postive"))
    }
    if((cor< -corFilter) & (pvalue<pvalueFilter)){
      outTab=rbind(outTab,cbind(d1=j,d2=i,cor,pvalue,Regulation="negative"))
    }
  }
}

table(outTab$d1)

sub <- outTab[outTab$d1%in%c('Subtype1','Subtype2','Subtype3','Subtype4'),]
sub <- sub[order(sub$d1,decreasing = F),]
sub <- sub[!sub$d2%in%c('Subtype1','Subtype2','Subtype3','Subtype4'),]
colnames(sub)[3:4] <- c('r','p.value')
sub$p.value <- as.numeric(sub$p.value)
sub$pd <- ifelse(sub$p.value<0.05,'< 0.05','>= 0.05')
sub$r <- as.numeric(sub$r)
sub$rd <- cut(sub$r, breaks = c(-Inf,-0.4,-0.2,0.2,0.4, Inf),
              labels = c("<=-0.4",'-0.4 - -0.2','-0.2 - 0.2', "0.2 - 0.4", ">= 0.4"))
table(sub$rd)
table(sub$d1)
sub$d1[sub$d1=='Subtype1'] <- 'C1'
sub$d1[sub$d1=='Subtype2'] <- 'C2'
sub$d1[sub$d1=='Subtype3'] <- 'C3'
sub$d1[sub$d1=='Subtype4'] <- 'C4'
varechem <- as.data.frame(a)
save(varechem,file = "RData/07.immune/butter_varechem_bottom.rdata")
quickcor(varechem, type = "lower") + 
  geom_colour() +
  anno_link(aes(colour = rd, size = pd), data = sub) +
  scale_size_manual(values = c(0.4, 0.2)) +
  scale_color_manual(values = c("#26B1D0","#DA6003","#A2A2A288","#E8288E","#65A818")) +
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#E11953",midpoint=0) +
  remove_axis("x")+
  guides(colour = guide_legend(title = "Subtype's r", 
                               override.aes = list(size = 3), 
                               order = 1),
         size = guide_legend(title = "Subtype's p",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         fill = guide_colorbar(title = "Pearson's r", order =3))

ggsave(filename = "Figure/Figure4/E_ggcor plot in bottom left.pdf", width = 10,height = 8)

#Performing submap immunotherapy prediction------------------------------
rm(list = ls())
library(tibble)
library(dplyr)

generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

### Creating the data format required to perform submap
load("Reference/07.immune/ICI_data.rda")
load("Database/GRCh38_v39.gtf.rda")
gene <- read.table("Reference/07.immune/skcm.immunotherapy.47samples.log2CountsNorm.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
gene <- rownames(gene)
skcm.immunotherapy.logNC <- list()
for (i in names(ICI_data)) {
  skcm.immunotherapy.logNC[[i]] <- ICI_data[[i]]$expr 
  skcm.immunotherapy.logNC[[i]] <- merge(gtf,skcm.immunotherapy.logNC[[i]],by.x=1,by.y=0)%>%
    dplyr::select(-gene_id)%>%
    group_by(gene_name)%>%
    summarise_all(max)%>%
    tibble::column_to_rownames('gene_name')
  skcm.immunotherapy.logNC[[i]] <- skcm.immunotherapy.logNC[[i]][gene,]
  
  rownames(skcm.immunotherapy.logNC[[i]]) <- toupper(rownames(skcm.immunotherapy.logNC[[i]])) # 基因大写，因为我使用的数据是把基因名都大写的
  skcm.immunotherapy.logNC[[i]] <- na.omit(skcm.immunotherapy.logNC[[i]])
  
}
save(skcm.immunotherapy.logNC,file = "RData/07.immune/submap/skcm.immunotherapy.logNC.rda")

for (iii in names(ICI_data)) {
  a <- colnames(ICI_data[[iii]]$clin)
  print(a)
}

skcm.immunotherapy.info <- list()
for (ii in names(ICI_data)) {
  skcm.immunotherapy.info[[ii]] <- ICI_data[[ii]]$clin[,1:2]
  skcm.immunotherapy.info[[ii]] <- skcm.immunotherapy.info[[ii]][order(skcm.immunotherapy.info[[ii]]$Response),]
  skcm.immunotherapy.info[[ii]] <- na.omit(skcm.immunotherapy.info[[ii]])
  skcm.immunotherapy.info[[ii]]$rank <- rep(c(1,2),times=as.character(table(skcm.immunotherapy.info[[ii]]$Response))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R
  skcm.immunotherapy.info[[ii]]$ID <- as.character(skcm.immunotherapy.info[[ii]]$ID)
}

save(skcm.immunotherapy.info,file = "RData/07.immune/submap/skcm.immunotherapy.info.rda")


library(tibble)
library(dplyr)
load('Database/Meta_RNAseq_expr.rda')
load("Database/GRCh38_v39.gtf.rda")
tmp <- merge(gtf,Meta_RNAseq_expr,by.x=1,by.y=0)%>%
  dplyr::select(-gene_id)%>%
  group_by(gene_name)%>%
  summarise_all(max)%>%
  tibble::column_to_rownames('gene_name')
load('RData/phenotype_RNAseq.rda')
rm(Meta_RNAseq_expr,gtf) 

GENELIST <- list()
in_gct <- list()
for (i in names(ICI_data)) {
  GENELIST[[i]] <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC[[i]]))
  in_gct[[i]] <- skcm.immunotherapy.logNC[[i]][GENELIST[[i]],skcm.immunotherapy.info[[i]]$ID]
  
}

save(GENELIST,file = "RData/07.immune/submap/GENELIST.rda")
sam_info <- skcm.immunotherapy.info

for (i in names(ICI_data)) {
  gct_file <- file.path("Table/07.immune/submap/",paste0(i,"for.SubMap.gct"))
  cls_file <- file.path("Table/07.immune/submap/",paste0(i,"for.SubMap.cls"))
  generateInputFileForSubMap(in_gct = in_gct[[i]], gct_file = gct_file, cls_file = cls_file, sam_info = sam_info[[i]], type_name = "rank")
  
}

 
samples.C1 <- clin$Sample[which(clin$Cluster == "C1")]
samples.C2 <- clin$Sample[which(clin$Cluster == "C2")]
samples.C3 <- clin$Sample[which(clin$Cluster == "C3")]
samples.C4 <- clin$Sample[which(clin$Cluster == "C4")]

sam_info <- data.frame("ImmClust"=c(samples.C1,samples.C2,samples.C3,samples.C4),row.names = c(samples.C1,samples.C2,samples.C3,samples.C4))
sam_info$rank <- rep(c(1,2,3,4),times=c(length(samples.C1),length(samples.C2),length(samples.C3),length(samples.C4))) 
in_gct <- tmp[gene,rownames(sam_info)] 

gct_file <- "Table/07.immune/submap/Immune2.for.SubMap.gct"
cls_file <- "Table/07.immune/submap/Immune2.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")
###GSE35640
rm(list = ls())
library(pheatmap)
library(export)
library(ggsci)
heatmap.YlGnPe <- c("#7AA6DCFF","#7AA6DC99","#7AA6DC66","#7AA6DC33")
cherry    <- pal_jco("default",alpha = 0.8)(8)[6]
lightgrey <- pal_jco("default",alpha = 0.4)(8)[3]
tmp <- matrix(c(1,0.008,1,1,0.032,1,1,1,1,0.001,1,1,0.008,1,1,1), 
              nrow = 8,byrow = T,dimnames = list(c("C1-P","C2-P","C3-P","C4-P",'C1 adj_P',"C2 adj_P","C3 adj_P","C4 adj_P"),
                                                 c("NoR","R")))
pheatmap(tmp, cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[1:4],
         gaps_row = 4,
         display_numbers = matrix(ifelse(tmp < 0.1, paste0("P=",tmp), ""),nrow(tmp)),number_color = "black",
         annotation_row = data.frame(pvalue=c("P value","P value","P value","P value",
                                              "FDR","FDR","FDR","FDR"),
                                     row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("P value"=lightgrey,
                                           "FDR"=cherry)),
         legend_breaks = c(0.1,0.5,1),
         main = "GSE35640")
graph2pdf(file = "Figure/Figure4/G_GSE35640.pdf",width=4,height=8)
###GSE91061
rm(list = ls())
library(pheatmap)
library(export)
library(ggsci)
heatmap.YlGnPe <- c("#7AA6DCFF","#7AA6DC99","#7AA6DC66","#7AA6DC33")
cherry    <- pal_jco("default",alpha = 0.8)(8)[6]
lightgrey <- pal_jco("default",alpha = 0.4)(8)[3]

tmp <- matrix(c(1,0.016,1,1,0.024,1,1,1,1,0.002,1,1,0.006,1,1,1),
              nrow = 8,byrow = T,dimnames = list(c("C1-P","C2-P","C3-P","C4-P",'C1 adj_P',"C2adj_P","C3adj_P","C4adj_P"),
                                                 c("PD1-NoR","PD1-R")))
pheatmap(tmp, cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[1:4],
         gaps_row = 4,
         display_numbers = matrix(ifelse(tmp < 0.1, paste0("P=",tmp), ""),nrow(tmp)),number_color = "black",
         annotation_row = data.frame(pvalue=c("P value","P value","P value","P value",
                                              "FDR","FDR","FDR","FDR"),
                                     row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("P value"=lightgrey,
                                           "FDR"=cherry)),
         legend_breaks = c(0.1,0.5,1),
         main = "GSE91061")
graph2pdf(file = "Figure/Figure4/G_GSE91061.pdf",width=4,height=8)
###GSE100797
rm(list = ls())
library(pheatmap)
library(export)
library(ggsci)
heatmap.YlGnPe <- c("#7AA6DCFF","#7AA6DC99","#7AA6DC66","#7AA6DC33")
cherry    <- pal_jco("default",alpha = 0.8)(8)[6]
lightgrey <- pal_jco("default",alpha = 0.4)(8)[3]

tmp <- matrix(c(1,0.008,1,1,0.040,1,1,1,1,0.001,1,1,0.005,1,1,1),
              nrow = 8,byrow = T,dimnames = list(c("C1-P","C2-P","C3-P","C4-P",'C1 adj_P',"C2 adj_P","C3 adj_P","C4 adj_P"),
                                                 c("CAR-T-NoR","CAR-T-R")))
pheatmap(tmp, cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[1:4],
         gaps_row = 4,
         display_numbers = matrix(ifelse(tmp < 0.1, paste0("P=",tmp), ""),nrow(tmp)),number_color = "black",
         annotation_row = data.frame(pvalue=c("P value","P value","P value","P value",
                                              "FDR","FDR","FDR","FDR"),
                                     row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("P value"=lightgrey,
                                           "FDR"=cherry)),
         legend_breaks = c(0.1,0.5,1),
         main = "GSE100797")
graph2pdf(file = "Figure/Figure4/G_GSE100797.pdf",width=4,height=8)
###GSE126044
rm(list = ls())
library(pheatmap)
library(export)
library(ggsci)
heatmap.YlGnPe <- c("#7AA6DCFF","#7AA6DC99","#7AA6DC66","#7AA6DC33")
cherry    <- pal_jco("default",alpha = 0.8)(8)[6]
lightgrey <- pal_jco("default",alpha = 0.4)(8)[3]

tmp <- matrix(c(1,0.008,1,1,0.104,1,1,1,1,0.001,1,1,0.026,1,1,1), 
              nrow = 8,byrow = T,dimnames = list(c("C1-P","C2-P","C3-P","C4-P",'C1 adj_P',"C2 adj_P","C3 adj_P","C4 adj_P"),
                                                 c("PD1-NoR","PD1-R")))
pheatmap(tmp, cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[1:4],
         gaps_row = 4,
         display_numbers = matrix(ifelse(tmp < 0.1, paste0("P=",tmp), ""),nrow(tmp)),number_color = "black",
         annotation_row = data.frame(pvalue=c("P value","P value","P value","P value",
                                              "FDR","FDR","FDR","FDR"),
                                     row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("P value"=lightgrey,
                                           "FDR"=cherry)),
         legend_breaks = c(0.1,0.5,1),
         main = "GSE126044")
graph2pdf(file = "Figure/Figure4/G_GSE126044.pdf",width=4,height=8)
###GSE135222
rm(list = ls())
library(pheatmap)
library(export)
library(ggsci)
heatmap.YlGnPe <- c("#7AA6DCFF","#7AA6DC99","#7AA6DC66","#7AA6DC33")
cherry    <- pal_jco("default",alpha = 0.8)(8)[6]
lightgrey <- pal_jco("default",alpha = 0.4)(8)[3]

tmp <- matrix(c(1,0.016,1,1,0.028,1,1,1,1,0.002,1,1,0.007,1,1,1),
              nrow = 8,byrow = T,dimnames = list(c("C1-P","C2-P","C3-P","C4-P",'C1 adj_P',"C2 adj_P","C3 adj_P","C4 adj_P"),
                                                 c("PD1-NoR","PD1-R")))
pheatmap(tmp, cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[1:4],
         gaps_row = 4,
         display_numbers = matrix(ifelse(tmp < 0.1, paste0("P=",tmp), ""),nrow(tmp)),number_color = "black",
         annotation_row = data.frame(pvalue=c("P value","P value","P value","P value",
                                              "FDR","FDR","FDR","FDR"),
                                     row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("P value"=lightgrey,
                                           "FDR"=cherry)),
         legend_breaks = c(0.1,0.5,1),
         main = "GSE135222")
graph2pdf(file = "Figure/Figure4/G_GSE135222.pdf",width=4,height=8)
###NATH
rm(list = ls())
library(pheatmap)
library(export)
library(ggsci)
heatmap.YlGnPe <- c("#7AA6DCFF","#7AA6DC99","#7AA6DC66","#7AA6DC33")
cherry    <- pal_jco("default",alpha = 0.8)(8)[6]
lightgrey <- pal_jco("default",alpha = 0.4)(8)[3]
tmp <- matrix(c(1,0.012,1,1,0.044,1,1,1,1,0.001,1,1,0.011,1,1,1),
              nrow = 8,byrow = T,dimnames = list(c("C1-P","C2-P","C3-P","C4-P",'C1 adj_P',"C2 adj_P","C3 adj_P","C4 adj_P"),
                                                 c("CTLA4-NoR","CTLA4-R")))
pheatmap(tmp, cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[1:4],
         gaps_row = 4,
         display_numbers = matrix(ifelse(tmp < 0.1, paste0("P=",tmp), ""),nrow(tmp)),number_color = "black",
         annotation_row = data.frame(pvalue=c("P value","P value","P value","P value",
                                              "FDR","FDR","FDR","FDR"),
                                     row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("P value"=lightgrey,
                                           "FDR"=cherry)),
         legend_breaks = c(0.1,0.5,1),
         main = "Nathanson")
graph2pdf(file = "Figure/Figure4/G_Nathanson.pdf",width=4,height=8)

#Comparison of gene mutations in four groups-----------------------
rm(list = ls())
load('RData/phenotype_RNAseq.rda')
options(stringsAsFactors = F)
library(data.table)
library(tibble)
library(dplyr)
library(maftools)
library(ggsci)
maf <- read.maf('Reference/08.mutation/All_maf.txt',isTCGA = T)
Cluster <- clin[389:532,1:2]
Cluster$sample <- substr(Cluster$Sample,1,12)
Cluster <- Cluster[,-1]
gg <- maf@variant.type.summary %>% as.data.frame()
gg$INDEL <- gg$INS+gg$DEL
colnames(gg)
gg <- select(gg,Tumor_Sample_Barcode,SNP,INDEL,total)
gg$TMB <- log2(gg$total+1)
gg$SNP <- log2(gg$SNP+1)
gg$INDEL <- log2(gg$INDEL+1)
gg <- gg[,-4]
gg <- merge(gg,Cluster,by.x=1,by.y = 2)
maf2 <- read.maf('Reference/08.mutation/All_maf.txt',isTCGA = T,clinicalData = gg)
c1 <- as.character(gg[gg$Cluster=='C1',][,1])###41
c2 <- as.character(gg[gg$Cluster=='C2',][,1])###35
c3 <- as.character(gg[gg$Cluster=='C3',][,1])###41
c4 <- as.character(gg[gg$Cluster=='C4',][,1])###22
maf3 <- subsetMaf(maf2,tsb = gg$Tumor_Sample_Barcode)
gene <- c("TP53","RB1","PIK3R1","PTEN","EGFR","ERBB2","PIK3CA","DNAH17","ATRX","NF1","PDGFRA","COL6A3","TTN")#,"TERT","IDH1","MUC16","SPTA1","RYR2","FLG","AHNAK2","COL6A3","DNAH17","DOCK5","HYDIN","PCLO","PKHD1"
oncoplot(maf = maf3,genes = gene,writeMatrix = T,removeNonMutated = F)

pp <- t(read.table('onco_matrix.txt',h=T,row.names = 1,sep = '\t',check.names = F))%>%as.data.frame()
pp[pp=='0'] <- 'x'
pp[pp==''] <- 'x'
save(pp,file = "RData/08.mutation/mutate_gene.rda")

rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(maftools)
library(stringr)
library(tidyverse)
load('RData/08.mutation/mutate_gene.rda')
load('RData/phenotype_RNAseq.rda')
cluster <- clin[389:532,]
cluster$Sample <- substr(cluster$Sample,1,12)
rownames(cluster) <- cluster$Sample
cluster <- cluster[rownames(pp),]
cluster <- cluster[order(cluster$Cluster),1:2]
mut <- as.data.frame(t(pp))
mut <- mut[,cluster$Sample]
col <- c('#438F5D','#E19412','#2C6589',"#BC3C29FF")
str_length(rownames(mut))
tmp <- mut
tmp[!tmp=='x'] <- 'Mutated'
cumsum(table(cluster$Cluster))
l <- apply(tmp[,1:41],1,function(x){sum(x=='Mutated')/length(x)})%>%as.numeric()%>%round(2)
l <- paste0(l*100,'%')
right <- rowAnnotation(foo = anno_text(l,gp = gpar(fontsize = 8)))
Top = HeatmapAnnotation(Cluster=cluster$Cluster[1:41],border = F,
                        show_legend = F,show_annotation_name = F,
                        simple_anno_size =unit(4,'mm'),
                        col = list(Cluster = c('C1' = col[1],
                                               'C2' = col[2],
                                               'C3' = col[3],
                                               'C4' = col[4])))
h1 <- Heatmap(as.matrix(tmp[,1:41]),
              col = c('Mutated'='#67AB9F','x'='WhiteSmoke'),
              right_annotation = right,
              row_names_gp = gpar(fontsize = 9,fontface='italic'),
              top_annotation = Top,
              column_split = cluster$Cluster[1:41],
              cluster_rows = F,cluster_columns = F,row_names_side = 'left',
              show_column_names = F,column_title = NULL,
              show_heatmap_legend = F)
h1
dev.off()

cumsum(table(cluster$Cluster))
l <- apply(tmp[,42:76],1,function(x){sum(x=='Mutated')/length(x)})%>%as.numeric()%>%round(2)
l <- paste0(l*100,'%')
right <- rowAnnotation(foo = anno_text(l,gp = gpar(fontsize = 8)))
Top = HeatmapAnnotation(Cluster=cluster$Cluster[42:76],border = F,
                        show_legend = F,show_annotation_name = F,
                        simple_anno_size =unit(4,'mm'),
                        col = list(Cluster = c('C1' = col[1],
                                               'C2' = col[2],
                                               'C3' = col[3],
                                               'C4' = col[4])))
h2 <- Heatmap(as.matrix(tmp[,42:76]),
              col = c('Mutated'='#67AB9F','x'='WhiteSmoke'),
              right_annotation = right,
              row_names_gp = gpar(fontsize = 9),
              top_annotation = Top,
              column_split = cluster$Cluster[42:76],
              cluster_rows = F,cluster_columns = F,row_names_side = 'left',
              show_column_names = F,column_title = NULL,
              show_row_names = F,
              show_heatmap_legend = F)
h2
dev.off()

cumsum(table(cluster$Cluster))
l <- apply(tmp[,77:117],1,function(x){sum(x=='Mutated')/length(x)})%>%as.numeric()%>%round(2)
l <- paste0(l*100,'%')
right <- rowAnnotation(foo = anno_text(l,gp = gpar(fontsize = 8)))
Top = HeatmapAnnotation(Cluster=cluster$Cluster[77:117],border = F,
                        show_legend = F,show_annotation_name = F,
                        simple_anno_size =unit(4,'mm'),
                        col = list(Cluster = c('C1' = col[1],
                                               'C2' = col[2],
                                               'C3' = col[3],
                                               'C4' = col[4])))
h3 <- Heatmap(as.matrix(tmp[,77:117]),
              col = c('Mutated'='#67AB9F','x'='WhiteSmoke'),
              right_annotation = right,
              row_names_gp = gpar(fontsize = 9),
              top_annotation = Top,
              column_split = cluster$Cluster[77:117],
              cluster_rows = F,cluster_columns = F,row_names_side = 'left',
              show_column_names = F,column_title = NULL,
              show_row_names = F,
              show_heatmap_legend = F)
h3
dev.off()

cumsum(table(cluster$Cluster))
l <- apply(tmp[,118:139],1,function(x){sum(x=='Mutated')/length(x)})%>%as.numeric()%>%round(2)
l <- paste0(l*100,'%')
right <- rowAnnotation(foo = anno_text(l,gp = gpar(fontsize = 8)))
Top = HeatmapAnnotation(Cluster=cluster$Cluster[118:139],border = F,
                        show_legend = F,show_annotation_name = F,
                        simple_anno_size =unit(4,'mm'),
                        col = list(Cluster = c('C1' = col[1],
                                               'C2' = col[2],
                                               'C3' = col[3],
                                               'C4' = col[4])))
h4 <- Heatmap(as.matrix(tmp[,118:139]),
              col = c('Mutated'='#67AB9F','x'='WhiteSmoke'),
              right_annotation = right,
              row_names_gp = gpar(fontsize = 9),
              top_annotation = Top,
              column_split = cluster$Cluster[118:139],
              cluster_rows = F,cluster_columns = F,row_names_side = 'left',
              show_column_names = F,column_title = NULL,
              show_row_names = F,
              show_heatmap_legend = F)
h4
dev.off()

mut.ht <- h1+h2+h3+h4
pdf("Figure/Figure5/A_mut_gene.pdf",7,2.5)
draw(mut.ht)
dev.off()


rm(list = ls())
library(tidyr)
library(ggplot2)
library(export)
library(ggsci)
load("Database/RNAseq_clin_new.rda")
#Comparison of IDH mutation-------------------------------------------------
clin <- RNAseq_clin[,c(2,10)]
table(clin$IDH)
table(clin$Cluster,clin$IDH)  

data <- data.frame(C1=c(11,115),C2=c(12,138),C3=c(35,100),C4=c(42,64),row.names = c("Mut","WT"))
fisher.test(data)#p-value = 2.351e-11
data2 <- data.frame(C1=c(11/126,115/126),
                    C2=c(12/150,138/150),
                    C3=c(35/135,100/135),
                    C4=c(42/106,64/106),
                    row.names = c("Mut","WT"))
data3 <- gather(data2,'Cluster','Value')
data3$IDH <- rep(c("Mut","WT"),4)


library(ggsci)
ggplot(data3,aes(Cluster,Value,fill=IDH))+
  geom_bar(stat = 'identity',color='#303030',width = 0.65)+
  theme_classic()+
  scale_fill_brewer(palette = 'Set2')+
  ggtitle('IDH Mutation****')+theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title = NULL))+
  scale_y_continuous(breaks=seq(0, 1, 0.5))+
  xlab(NULL)+ylab(NULL)+
  scale_fill_manual(values = c("#305883",'#ffc858'))+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size = 16,hjust = 0.5))

graph2pdf(file='Figure/Figure5/B_IDH_RNAseq.pdf',height=2.9,width=2.5)

#Comparison of G-CIMP------------------------------------------------
clin <- RNAseq_clin[389:532,c(2,14)]
table(clin$`G-CIMP
      methylation`)
table(clin$Cluster,clin$`G-CIMP
      methylation`)
clin$`G-CIMP
methylation`[clin$`G-CIMP
methylation`=="non-G-CIMP"] <- "non"
data <- data.frame(C1=c(0,41),C2=c(0,35),C3=c(4,38),C4=c(3,19),row.names = c('G-CIMP','non-G-CIMP'))
fisher.test(data)
data2 <- data.frame(C1=c(0/41,41/41),C2=c(0/35,35/35),C3=c(4/42,38/42),C4=c(3/22,19/22),row.names = c('G-CIMP','non-G-CIMP'))
data3 <- gather(data2,'Cluster','Value')
data3$G_CIMP_methylation <- rep(c('G-CIMP','non'),4)


ggplot(data3,aes(Cluster,Value,fill=G_CIMP_methylation))+
  geom_bar(stat = 'identity',color='#303030', width = 0.65)+
  theme_classic()+
  scale_fill_brewer(palette = 'Set2')+
  scale_y_continuous(breaks=seq(0, 1, 0.5))+
  ggtitle('G-CIMP Status**')+theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title = NULL))+
  xlab(NULL)+ylab(NULL)+
  scale_fill_manual(values = c("#05007f","#87ceeb"))+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size = 16,hjust = 0.5))

graph2pdf(file='Figure/Figure5/E_G_CIMP_RNAseq.pdf',height=2.9,width=2.5)

#Comparison of Chr_1p19q-----------------------------------------------------------
clin <- RNAseq_clin[,c(2,12)]
table(clin$Chr_1p19q)
clin$Chr_1p19q[clin$Chr_1p19q=="Non_codel"] <- "non-codel"
clin$Chr_1p19q[clin$Chr_1p19q=="non_codel"] <- "non"
table(clin$Cluster,clin$Chr_1p19q)  

data <- data.frame(C1=c(1,127),C2=c(1,116),C3=c(9,127),C4=c(9,101),row.names = c("Codel","non-codel"))
fisher.test(data)
data2 <- data.frame(C1=c(1/128,127/128),C2=c(1/117,116/117),C3=c(9/136,127/136),C4=c(9/110,101/110),row.names = c("Codel","non-codel"))
data3 <- gather(data2,'Cluster','Value')
data3$Chr_1p19q <- rep(c("Codel","Non-codel"),4)

library(ggsci)
ggplot(data3,aes(Cluster,Value,fill=Chr_1p19q))+
  geom_bar(stat = 'identity',color='#303030',width = 0.65)+
  theme_classic()+
  scale_fill_brewer(palette = 'Set2')+
  ggtitle('Chr_1p19q**')+theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title = NULL))+
  xlab(NULL)+ylab(NULL)+
  scale_y_continuous(breaks=seq(0, 1, 0.5))+
  scale_fill_manual(values = c("#85110f","#ff8b00"))+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size = 16,hjust = 0.5))

graph2pdf(file='Figure/Figure5/F_Chr_1p19q_RNAseq.pdf',height=2.9,width=2.5)

#ZZU cohort
rm(list = ls())
load("RData/04.validation/phenotype_ZGBM.rdata")
load("Database/ZGBM_clin.rda")
ZGBM_clin <- ZGBM_clin[dd$ID,]

clin <- ZGBM_clin[,c(5,6,8,10,11,17,18,20,21)]
clin <- merge(dd[,1:2],clin,by.x = 1,by.y = 0)


#Comparison of IDHmutation in ZZU cohort---------------------------------------------
table(clin$IDH)
table(clin$Cluster,clin$IDH)
data <- data.frame(C1=c(1,62),C2=c(0,40),C3=c(7,51),C4=c(7,42),row.names = c('Mut','WT'))
fisher.test(data)
data2 <- data.frame(C1=c(1/63,62/63),C2=c(0/40,40/40),C3=c(7/58,51/58),C4=c(7/49,42/49),row.names = c('Mut','WT'))
data3 <- gather(data2,'Cluster','Value')
data3$IDH <- rep(c('Mut','WT'),4)

library(ggsci)
library(export)
ggplot(data3,aes(Cluster,Value,fill=IDH))+
  geom_bar(stat = 'identity',color='#303030',width = 0.65)+
  theme_classic()+
  scale_fill_brewer(palette = 'Set2')+
  ggtitle('IDH Mutation**')+theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title = NULL))+
  xlab(NULL)+ylab(NULL)+
  scale_y_continuous(breaks=seq(0, 1, 0.5))+
  scale_fill_manual(values = c("#305883",'#ffc858'))+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size = 16,hjust = 0.5))

graph2pdf(file='Figure/Figure5/C_IDH_ZGBM.pdf',height=2.9,width=2.5)

#Comparison of TERTp--------------------------------------
table(clin$TERTp)
table(clin$Cluster,clin$TERTp)
data <- data.frame(C1=c(32,13),C2=c(26,7),C3=c(20,22),C4=c(21,18),row.names = c('Mut','WT'))
fisher.test(data)
data2 <- data.frame(C1=c(32/45,13/45),
                    C2=c(26/33,7/33),
                    C3=c(20/42,22/42),
                    C4=c(21/39,18/39),
                    row.names = c('Mut','WT'))
data3 <- gather(data2,'Cluster','Value')
data3$TERTp <- rep(c('Mut','WT'),4)

library(ggsci)
library(export)
ggplot(data3,aes(Cluster,Value,fill=TERTp))+
  geom_bar(stat = 'identity',color='#303030',width = 0.65)+
  theme_classic()+
  scale_fill_brewer(palette = 'Set2')+
  ggtitle('TERTp Mutation**')+theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title = NULL))+
  xlab(NULL)+ylab(NULL)+
  scale_y_continuous(breaks=seq(0, 1, 0.5))+
  scale_fill_manual(values = c("#305883",'#ffc858'))+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size = 16,hjust = 0.5))

graph2pdf(file='Figure/Figure5/D_TERTp_ZGBM.pdf',height=2.9,width=2.5)

