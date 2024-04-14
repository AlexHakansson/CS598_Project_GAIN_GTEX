library(readxl)
library(tidyr)
library(tidyverse)
library(dplyr)




#### Read in lung and colon data to combine

lung1 <- read.table("data/gene_reads_2017-06-05_v8_lung.gct/gene_reads_lung.gct",
                    sep="\t",skip = 2,header = 1)
colon <- read.table("data/gene_reads_2017-06-05_v8_colon_transverse.gct/gene_reads_colon_transverse.gct",
                    sep="\t",skip = 2,header = 1)

### make function to process data frame so that it is more readable to users
df_process <- function(gtex_df,gene_list = NULL){
  
  if (!is.null(gene_list)){
    gtex_df <- gtex_df %>% filter(Description%in% gene_list)
  }
  
  fix_df <- gtex_df %>% column_to_rownames("Name") %>% select(-Description,-id)
  fix_df <- fix_df %>% t()
  # row.names(fix_df) <- substr(rownames(fix_df),1,10) %>% gsub(".","-",.,fixed=T)
  row.names(fix_df) <- rownames(fix_df) %>% gsub(".","-",.,fixed=T)
  return(fix_df)
}

## load in p53 genes
p53_genes <-  c('AIFM2', 'APAF1', 'ATM', 'ATR', 'BAX', 'BBC3', 'BCL2', 'BCL2L1',
                'BID', 'CASP3', 'CASP8', 'CASP9', 'CCNB1', 'CCNB2', 'CCND1',
                'CCND2', 'CCND3', 'CCNE1', 'CCNE2', 'CCNG1', 'CCNG2', 'CD82',
                'CDK1', 'CDK2', 'CDK4', 'CDK6', 'CDKN1A', 'CDKN2A', 'CHEK1',
                'CHEK2', 'CYCS', 'DDB2', 'EI24', 'FAS', 'GADD45A', 'GADD45B',
                'GADD45G', 'GORAB', 'GTSE1', 'IGFBP3', 'MDM2', 'MDM4', 'PERP',
                'PMAIP1', 'PPM1D', 'PTEN', 'RCHY1', 'RRM2B', 'SERPINE1', 'SESN1',
                'SESN2', 'SESN3', 'SFN', 'SHISA5', 'SIAH1', 'SIVA1', 'STEAP3',
                'THBS1', 'TNFRSF10A', 'TNFRSF10B', 'TP53', 'TP53I3', 'TP73',
                'TSC2', 'ZMAT3')
### process lung and colon dat and add tissue column
lungp <- df_process(lung1,p53_genes) %>% data.frame(check.rows = F)
colp <- df_process(colon,p53_genes)%>% data.frame(check.rows = F)
lungp$tissue <- "lung"
colp$tissue <- "colon"

com_cols <- intersect(colnames(lungp),colnames(colp))

### save combined data
comb_dat <- rbind(lungp[,com_cols],colp[,com_cols])
comb_dat %>% write.csv("data/gtex_lung_colon.csv")


### read in meta dataframe so that we can fix sample ID to match gene expression sample ID
meta_dat <- read.table("data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                       sep="\t",header=1)


meta_dat %>% head()

rname_id <- data.frame(FULL_ID = rownames(comb_dat),
                       SUBJID= substr(rownames(comb_dat),1,10) )

meta_dat_comb <- full_join(meta_dat,rname_id,by="SUBJID")

dim(meta_dat_comb)

meta_dat_lc <- meta_dat_comb %>% filter(!is.na(FULL_ID))
meta_dat_lc$SUBJID <- meta_dat_lc$FULL_ID

### fix age so that it is a number not a character
meta_dat_lc$AGE_chr <- meta_dat_lc$AGE
meta_dat_lc$AGE <- meta_dat_lc$AGE_chr %>% substr(1,2) %>% gsub("0","5",.) %>% as.numeric()
meta_dat_lc %>% write.table("data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS_Lung_Colon.txt",row.names = F,
                            sep = "\t")

### check that all samples are mathing
setdiff(rownames(comb_dat),meta_dat_lc$SUBJID) %>% length()


