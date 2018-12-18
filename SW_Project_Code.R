#Stephanie Wankowicz#
setwd('/Users/stephaniewankowicz/Dropbox/Biostats/CM+iPSC_Hi-C/')

#packages
library(stringr)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(bedr)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

##########################LOAD IN AND PREPROCESS DATA###################################
#1)gene lists
#2)RNA
#3)Full Gene Lists
#4)Signifianct Promoters
#5)Histone Marks

###Gene List###
setwd('/Users/stephaniewankowicz/Dropbox/Biostats/Final_Project/Individual_Project')
gene_lists<-read.csv('CHD_gene_list.csv',header = T)
CHD_only<-(gene_lists$OMIM_CHDGenes_Only)
CHD_only<-CHD_only[CHD_only != ""]

ASD_only<-gene_lists$ASD_Genes.1
ASD_only<-ASD_only[ASD_only != ""]
nrow(high_expression_cm_RNA[high_expression_cm_RNA$Name %in% ASD_only, ])

###RNA Data###
#function for reading in & processing RNA
load_preprocess_RNA_data<-function(file_list){
  df=data.frame()
  for (i in 1:length(file_list)){
    df_tmp<-read.table(file_list[i], sep = '\t', header = TRUE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))
    df <- rbind(df, df_tmp)
    rm(df_tmp)
  }
  gene_names<-unique(df$Name)
  collapsed_gene_expression<-data.frame()
  n=1
  for(i in gene_names){
    df_tmp<-df[df$Name==i, ]
    collapsed_gene_expression[n,'Name']<-i
    collapsed_gene_expression[n, 'Mean_TPM']<-mean(df_tmp$TPM)
    collapsed_gene_expression[n,'SD_TPM']<-sd(df_tmp$TPM)
    n=n+1
    rm(df_tmp)
  }
  for(i in 1:nrow(collapsed_gene_expression)){
    expression=collapsed_gene_expression[i,'Mean_TPM']
    if(expression==0){
      collapsed_gene_expression[i,'TPM_Cateogry']<-'TPM=0'}
    else if(expression>0&expression<=3){
      collapsed_gene_expression[i,'TPM_Cateogry']<-'TPM=[0-3]'}
    else if(expression>4&expression<=25){
      collapsed_gene_expression[i,'TPM_Cateogry']<-'TPM=[4-25]'}
    else if(expression>25&expression<=150){
      collapsed_gene_expression[i,'TPM_Cateogry']<-'TPM=[26-150]'}
    else{
      collapsed_gene_expression[i,'TPM_Cateogry']<-'TPM>150'}
    rm(expression)
  }
  print(head(collapsed_gene_expression))
  return(collapsed_gene_expression)
}

setwd('/Users/stephaniewankowicz/Dropbox/Biostats/Final_Project/RNA/')
file_list_ipsc <- list.files(pattern = 'iPSC*')
ipsc_processed_rna<-load_preprocess_RNA_data(file_list_ipsc)
file_list_cm <-list.files(pattern='cardiomyocyte*')
cm_processed_rna<-load_preprocess_RNA_data(file_list_cm)

###All Possible Genes###
#rearrange sup. table 9.1 (all promoters captured in experiment)
setwd('/Users/stephaniewankowicz/Dropbox/Biostats/Final_Project/')
table_9.1<-read.csv('elife-35788-supp9-v2.txt', sep = '\t', header = TRUE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))
table_9.1$chr<-str_split_fixed(table_9.1$Hg19.coordinate, fixed(":"), 2)[, 1]
table_9.1$Gene_Name2<-str_split_fixed(table_9.1$Gene.name, fixed("*"), 2)[, 1]

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene_list<-getBM(attributes= "hgnc_symbol", mart=ensembl)
#subsetting genes down to only those with promoters tested in this experiment.
gene_list_w_promoters<-intersect(gene_list$hgnc_symbol,table_9.1$Gene_Name2)

####Significant Promoters####
setwd('/Users/stephaniewankowicz/Dropbox/Biostats/Final_Project/Sign_Promoters/')
CM_sign_promoters<-read.csv('CM_elife-35788-supp2-v2.txt', sep = '\t', header = FALSE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))
iPSC_sign_promoters<-read.csv('iPSC_elife-35788-supp1-v2.txt', sep = '\t', header = FALSE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))
sign_promoter_col_names<-c('Chr_Prom', 'Start_Pos_Prom','End_Pos_Prom', 'Chr_Inter', 'Start_Pos_Inter', 'End_Pos_Inter', 'Score', 'Gene_ID_Prom')
colnames(CM_sign_promoters)<-sign_promoter_col_names
colnames(iPSC_sign_promoters)<-sign_promoter_col_names
CM_sign_promoters$Gene_Name<-str_split_fixed(CM_sign_promoters$Gene_ID_Prom, fixed("*"), 2)[, 1]
CM_sign_promoters$Gene_Name_Distal<-str_split_fixed(CM_sign_promoters$Gene_ID_Prom, fixed("*"), 6)[, 5]
iPSC_sign_promoters$Gene_Name<-str_split_fixed(iPSC_sign_promoters$Gene_ID_Prom, fixed("*"), 2)[, 1]
iPSC_sign_promoters$Gene_Name_Distal<-str_split_fixed(iPSC_sign_promoters$Gene_ID_Prom, fixed("*"), 6)[, 5]
CM_sign_promoters$Gene_Name_Distal<-str_split_fixed(CM_sign_promoters$Gene_ID_Prom, fixed("*"), 6)[, 5]
#output bedfile for BEDTOOLS
write.table(cm_bed_file, '181125_CM_BedFile.txt', sep = '\t', row.names = FALSE, quote = FALSE, col.names =FALSE)
write.table(ipsc_bed_file, '181125_ipsc_BedFile.txt', sep = '\t', row.names = FALSE, quote = FALSE, col.names =FALSE)

####Histone Marks####
setwd('/Users/stephaniewankowicz/Dropbox/Biostats/Final_Project/Histone_Marks')
#CM_H3K27ac<-read.csv('E095-H3K27ac.narrowPeak', sep = '\t', header = FALSE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))
CM_H3K27me3<-read.csv('E095-H3K27me3.narrowPeak', sep = '\t', header = FALSE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))
CM_H3K4me1<-read.csv('E095-H3K4me1.narrowPeak', sep = '\t', header = FALSE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))

#female fetal brain
FFB_H3K27me3<-read.csv('E082-H3K27me3.narrowPeak',sep = '\t', header = FALSE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))
FFB_H3K4me1<-read.csv('E082-H3K4me1.narrowPeak',sep = '\t', header = FALSE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))
#FFB_H3K27ac<-read.csv('',sep = '\t', header = FALSE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))

#male fetal brain
FMB_H3K27me3<-read.csv('E081-H3K27me3.narrowPeak',sep = '\t', header = FALSE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))
FMB_H3K4me1<-read.csv('E081-H3K4me1.narrowPeak',sep = '\t', header = FALSE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))

#creating bedfiles for promoter regions
ipsc_bed_file<-iPSC_sign_promoters[,1:3]
cm_bed_file<-CM_sign_promoters[,1:3]
CM_H3K27me3_bed<-CM_H3K27me3[,1:3]
CM_H3K4me1_bed<-CM_H3K4me1[,1:3]
write.table(CM_H3K27me3_bed, 'CM_H3K27me3_bed',sep = '\t', row.names = FALSE, quote = FALSE, col.names =FALSE)
write.table(CM_H3K4me1_bed, 'CM_H3K4me1_bed',sep = '\t', row.names = FALSE, quote = FALSE, col.names =FALSE)

##############ANALYSIS & PLOTTING######################

###Expression Plot & Analysis###
CHD_ASD<-c(CHD_only,ASD_only)
expression_CHD_ASD_CM<-cm_processed_rna[cm_processed_rna$Name %in% CHD_ASD,]
expression_CHD_ASD_CM$Disease<-ifelse(expression_CHD_ASD$Name %in% CHD_only, 'CHD', 'ASD')
expression_CHD_ASD_i<-ipsc_processed_rna[ipsc_processed_rna$Name %in% CHD_ASD,]
expression_CHD_ASD_i$Disease<-ifelse(expression_CHD_ASD$Name %in% CHD_only, 'CHD', 'ASD')

ggplot(data=expression_CHD_ASD_i, aes(x=Disease, y=log(Mean_TPM), fill=Disease)) +geom_boxplot()+scale_fill_brewer()
ggplot(data=expression_CHD_ASD_CM, aes(x=Disease, y=log(Mean_TPM), fill=Disease)) +geom_boxplot()+scale_fill_brewer(palette = "Reds")

wilcox.test(log(expression_CHD_ASD[expression_CHD_ASD$Disease=='CHD',]$Mean_TPM),log(expression_CHD_ASD[expression_CHD_ASD$Disease=='ASD',]$Mean_TPM))

##Obtaining Bivalent Promoters##
#bedtools intersect -a CM_H3K4me1_bed -b CM_H3K27me3_bed -u -f 0.2 -r #20% overlap
#bedtools intersect -a CM_bivalent_bed -b CM_sign_promoter.bed -u #any overlap

##random genes##
#random genes#
#getting background/expected proportions
random_gene = DataFrame()
for (i in 1:1000){
  random_genes<-sample(gene_list_w_promoters, 200, replace=F)
  CM_bi<-CM_bivalent_Promoter2[CM_bivalent_Promoter2$Gene_Name %in% random_genes, ]
  iPSC_bi<-iPSC_bivalent_Promoter2[iPSC_bivalent_Promoter2$Gene_Name %in% random_genes, ]
  CM<-CM_sign_promoters[CM_sign_promoters$Gene_Name %in% random_genes, ]
  iPSC<-iPSC_sign_promoters[iPSC_sign_promoters$Gene_Name %in% random_genes, ]
  CM_bi_table<-CM_bi %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
  iPSC_bi_table<-iPSC_bi %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
  CM_table<-CM %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
  iPSC_table<-iPSC %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
  random_gene[i,'CM_bivalent_Promoter']<-mean(CM_bi_table$`n_distinct(Start_Pos_Inter)`)
  random_gene[i,'iPSC_bivalent_Promoter']<-mean(iPSC_bi_table$'n_distinct(Start_Pos_Inter)')
  random_gene[i,'CM_all_sign_Promoters']<-mean(CM_table$'n_distinct(Start_Pos_Inter)')
  random_gene[i,'iPSC_all_sign_Promoters']<-mean(iPSC_table$'n_distinct(Start_Pos_Inter)')
  random_gene[i, 'iPSC_bi_number']<-length(unique(iPSC_bi$Gene_Name))
  random_gene[i, 'CM_bi_number']<-length(unique(CM_bi$Gene_Name))
  random_gene[i, 'iPSC_number']<-length(unique(iPSC$Gene_Name))
  random_gene[i, 'CM_number']<-length(unique(CM$Gene_Name))
}
random_gene<-as.data.frame(random_gene)
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
random_gene[is.nan(random_gene)]<-0

###Number of CM/iPSC significant interactions in CHD & ASD gene lists###
CM_table<-CM_sign_promoters %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
mean_CM_inter=mean(CM_table$`n_distinct(Start_Pos_Inter)`)

iPSC_table<-iPSC_sign_promoters %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
mean_iPSC_inter=mean(iPSC_table$`n_distinct(Start_Pos_Inter)`)

#CHD#
CM_interactions_CHD<-CM_sign_promoters[CM_sign_promoters$Gene_Name %in% gene_lists$OMIM_CHDGenes_Only, ]
iPSC_interactions_CHD<-iPSC_sign_promoters[iPSC_sign_promoters$Gene_Name %in% gene_lists$OMIM_CHDGenes_Only, ]
CM_table_CHD<-CM_interactions_CHD %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
c=mean(CM_table_CHD$`n_distinct(Start_Pos_Inter)`)
iPSC_table_CHD<-iPSC_interactions_CHD %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
mean_iPSC_CHD_inter=mean(iPSC_table_CHD$`n_distinct(Start_Pos_Inter)`)

#ASD#
CM_interactions_ASD<-CM_sign_promoters[CM_sign_promoters$Gene_Name %in% gene_lists$ASD_Genes, ]
iPSC_interactions_ASD<-iPSC_sign_promoters[iPSC_sign_promoters$Gene_Name %in% gene_lists$ASD_Genes, ]
CM_table_ASD<-CM_interactions_ASD %>%group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
mean_CM_ASD_inter=mean(CM_table_ASD$`n_distinct(Start_Pos_Inter)`)
iPSC_table_ASD<-iPSC_interactions_ASD %>%group_by(Gene_Name) %>%summarise(n_distinct(Start_Pos_Inter))
mean_iPSC_ASD_inter=mean(iPSC_table_ASD$`n_distinct(Start_Pos_Inter)`)

##stats (random versus CHD or ASD)
(sum(random_gene$CM_all_sign_Promoters > mean_CM_CHD_inter) + 1) / (1000 + 1) #CHD CM
(sum(random_gene$CM_all_sign_Promoters > mean_CM_ASD_inter) + 1) / (1000 + 1) #ASD CM
(sum(random_gene$iPSC_all_sign_Promoters > mean_iPSC_ASD_inter) + 1) / (1000 + 1) #ASD iPSC
(sum(random_gene$iPSC_all_sign_Promoters > mean_iPSC_CHD_inter) + 1) / (1000 + 1) #CHD iPSC


#bivalent
iPSC_bivalent_Promoter<-read.csv('iPSC_bivalent_interactions',sep = '\t', header = FALSE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))
CM_bivalent_Promoter<-read.csv('CM_bivalent_interactions',sep = '\t', header = FALSE, comment.char = '#', skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings=c("","NA"))
iPSC_bivalent_Promoter2<-iPSC_sign_promoters[iPSC_sign_promoters$Promoter %in% iPSC_bivalent_Promoter$Promoter, ]
CM_bivalent_Promoter2<-CM_sign_promoters[CM_sign_promoters$Promoter %in% CM_bivalent_Promoter$Promoter, ]
#bivalent all
CM_table_bivalent<-CM_bivalent_Promoter2 %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
mean_CM_inter_bivalent=mean(CM_table_bivalent$`n_distinct(Start_Pos_Inter)`)

iPSC_table_bivalent<-iPSC_bivalent_Promoter2 %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
mean_iPSC_inter_bivalent=mean(iPSC_table_bivalent$`n_distinct(Start_Pos_Inter)`)

#bivalent CHD
CM_inter_bivalent_CHD<-CM_bivalent_Promoter2[CM_bivalent_Promoter2$Gene_Name %in% gene_lists$OMIM_CHDGenes_Only, ]
iPSC_inter_bivalent_CHD<-iPSC_bivalent_Promoter2[iPSC_bivalent_Promoter2$Gene_Name %in% gene_lists$OMIM_CHDGenes_Only, ]
CM_table_CHD_bivalent<-CM_inter_bivalent_CHD %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
mean_CM_CHD_inter_bivalent=mean(CM_table_CHD_bivalent$`n_distinct(Start_Pos_Inter)`)
iPSC_table_CHD_bivalent<-iPSC_inter_bivalent_CHD %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
c=mean(iPSC_table_CHD_bivalent$`n_distinct(Start_Pos_Inter)`)

#bivalent ASD
CM_inter_bivalent_ASD<-CM_bivalent_Promoter2[CM_bivalent_Promoter2$Gene_Name %in% ASD_only, ]
iPSC_inter_bivalent_ASD<-iPSC_bivalent_Promoter2[iPSC_bivalent_Promoter2$Gene_Name %in% ASD_only, ]
CM_table_ASD_bivalent<-CM_inter_bivalent_ASD %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
mean_CM_ASD_inter_bivalent=mean(CM_table_ASD_bivalent$`n_distinct(Start_Pos_Inter)`)
iPSC_table_ASD_bivalent<-iPSC_inter_bivalent_ASD %>% group_by(Gene_Name) %>% summarise(n_distinct(Start_Pos_Inter))
mean_iPSC_ASD_inter_bivalent=mean(iPSC_table_ASD_bivalent$`n_distinct(Start_Pos_Inter)`)

#random genes 
mean_random_CM=mean(random_gene$CM_all_sign_Promoters)
mean_random_iPSC=mean(random_gene$iPSC_all_sign_Promoters)
#random genes bivalent
mean_random_CM_bi=mean(random_gene$CM_bivalent_Promoter)
mean_random_iPSC_bi=mean(random_gene$iPSC_bivalent_Promoter)

##stats
(sum(random_gene$CM_bivalent_Promoter > mean_CM_CHD_inter_bivalent) + 1) / (1000 + 1) #CHD CM
(sum(random_gene$CM_bivalent_Promoter < mean_CM_ASD_inter_bivalent) + 1) / (1000 + 1) #ASD CM
(sum(random_gene$iPSC_bivalent_Promoter > mean_iPSC_ASD_inter_bivalent) + 1) / (1000 + 1) #ASD iPSC
(sum(random_gene$iPSC_bivalent_Promoter < mean_iPSC_CHD_inter_bivalent) + 1) / (1000 + 1) #CHD iPSC


###bargraphs### (plot B&C)

dist_inter<-read.table('Distance_Interactions.txt', sep='\t', header = TRUE)
ggplot(data=dist_inter, aes(x=Cell_Type, y=Freq_Interactions)) + geom_bar(stat='Identity', position='Dodge',aes(fill=Disease)) +scale_fill_brewer(palette = "Reds") +theme(axis.text=element_text(size=16),axis.title=element_text(size=14,face="bold"))
ggplot(data=dist_inter, aes(x=Cell_Type, y=Bivalent_Promoter_Freq)) + geom_bar(stat='Identity', position='Dodge',aes(fill=Disease)) +scale_fill_brewer(palette = "Reds")+theme(axis.text=element_text(size=16),axis.title=element_text(size=14,face="bold"))


#second expression boxplot (Bivalent Promoters)
expression_bivalent<-cm_processed_rna[cm_processed_rna$Name %in% CM_bivalent_Promoter2$Gene_Name,]
expression_bivalent$Disease<-ifelse(expression_bivalent$Name %in% CHD_only, 'CHD', 'No_Association')
ggplot(data=expression_bivalent, aes(x=Disease, y=log(Mean_TPM), fill=Disease)) +geom_boxplot()+scale_fill_brewer(palette = 'Reds')

expression_bivalent_i<-ipsc_processed_rna[ipsc_processed_rna$Name %in% iPSC_bivalent_Promoter2$Gene_Name,]
expression_bivalent_i$Disease<-ifelse(expression_bivalent_i$Name %in% CHD_only, 'CHD', 'No_Association')
ggplot(data=expression_bivalent_i, aes(x=Disease, y=log(Mean_TPM), fill=Disease)) +geom_boxplot()+scale_fill_brewer()

