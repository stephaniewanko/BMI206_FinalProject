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
