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
p53_genes <-  c('AIFM2', 'APAF1', 'ATM', 'ATR', 'BAX',  'BBC3', 'BCL2', 'BCL2L1',
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
                       SUBJID= substr(rownames(comb_dat),1,10)%>% gsub("-$","",.) )

substr(rownames(comb_dat),1,10) %>% gsub("-$","",.) %>% tail()

meta_dat_comb <- full_join(meta_dat,rname_id,by="SUBJID")

dim(meta_dat_comb)

meta_dat_lc <- meta_dat_comb %>% filter(!is.na(FULL_ID))
meta_dat_lc$SUBJID <- meta_dat_lc$FULL_ID

### fix age so that it is a number not a character
meta_dat_lc$AGE_chr <- meta_dat_lc$AGE
meta_dat_lc$AGE <- meta_dat_lc$AGE_chr %>% substr(1,2) %>% gsub("0","5",.) %>% as.numeric()

meta_dat_lc_nn <- meta_dat_lc %>% filter(!is.na(AGE))
# meta_dat_lc <- meta_dat_lc %>% mutate(
#   AGE=meta_dat_lc_nn[match(SUBJID,meta_dat_lc_nn$SUBJID),"AGE"],
#   SEX=meta_dat_lc_nn[match(SUBJID,meta_dat_lc_nn$SUBJID),"SEX"],
#   DTHHRDY=meta_dat_lc_nn[match(SUBJID,meta_dat_lc_nn$SUBJID),"DTHHRDY"]
# )

meta_dat_lc %>% tail()

meta_dat_lc %>% write.table("data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS_Lung_Colon.txt",row.names = F,
                            sep = "\t")

### check that all samples are mathing
setdiff(rownames(comb_dat),meta_dat_lc$SUBJID) %>% length()


### now we clean up tcga dat

### read in lung tcga mrna data file
lung_tcga <- read.table("data/luad_tcga.tar/luad_tcga/data_mrna_seq_v2_rsem.txt", 
                        header = 1)
### check that all genes are in the table
setdiff(p53_genes,lung_tcga$Hugo_Symbol)

### they are not so we will have to correct some gene names
  
rename_genes <-  c("MDM2","PTEN","GADD45B","RCHY1")
names(rename_genes) <- c("MGC5370","TEP1","DKFZP566B133", "PRO1996")
lung_tcga <- lung_tcga %>% select(-Entrez_Gene_Id)
p53_genes_tcga <- c(p53_genes,names(rename_genes))
lung_tcga <- lung_tcga %>% filter(Hugo_Symbol%in%p53_genes_tcga)

### rename the hugo genes so that they match GTEX
lung_tcga <- lung_tcga %>% mutate(
  Hugo_Symbol_corrected=ifelse(Hugo_Symbol%in%names(rename_genes),
                               rename_genes[Hugo_Symbol],Hugo_Symbol))
lung_tcga_cleaned <- lung_tcga %>% select(-Hugo_Symbol) %>% column_to_rownames("Hugo_Symbol_corrected")

### get only one sample per patient to make things cleaner
lung_samps <- lung_tcga_cleaned %>% colnames() %>%substr(1,12) %>% gsub(".","-",.,
                                                                        fixed=T)


# lung_tcga_cleaned %>%t()%>% write.csv("data/lung_tcga_p53_cleaned.csv")
lung_tcga_cleaned %>%# t()%>%
  write.csv("data/LUAD_tcga_p53_cleaned.csv")

### now lets do the same for colon cancer

colon_tcga <- read.table("data/coadread_tcga.tar/coadread_tcga/data_mrna_seq_v2_rsem.txt",
                         header = 1)

### double check we are missing the same genes
setdiff(p53_genes,colon_tcga$Hugo_Symbol)


colon_tcga <- colon_tcga %>% filter(Hugo_Symbol%in%p53_genes_tcga) %>% 
  select(-Entrez_Gene_Id)

colon_tcga <- colon_tcga %>% mutate(
  Hugo_Symbol_corrected=ifelse(Hugo_Symbol%in%names(rename_genes),
                               rename_genes[Hugo_Symbol],Hugo_Symbol))
colon_tcga_cleaned <- colon_tcga %>% select(-Hugo_Symbol) %>% column_to_rownames("Hugo_Symbol_corrected")
colon_tcga_cleaned 

# colon_tcga_cleaned %>% t() %>% write.csv("data/colon_tcga_p53_cleaned.csv")
colon_tcga_cleaned %>% #t() %>% 
  write.csv("data/COAD_tcga_p53_cleaned.csv")

### read in clinical data for tcga and clean it

lung_clinical_samp <- read.table("data/luad_tcga.tar/luad_tcga/data_clinical_sample.txt",skip=4,
                            sep="\t",header=1)
lung_clinical_pat <- read.table("data/luad_tcga.tar/luad_tcga/data_clinical_patient.txt",skip=4,
                                 sep="\t",header=1)

lung_clinical <- left_join(lung_clinical_samp,lung_clinical_pat,by="PATIENT_ID")

### impute age as average to make things easier
lung_average <- mean(as.numeric(lung_clinical$AGE),na.rm = T) %>% as.integer()
lung_clinical$AGE <- lung_clinical$AGE %>% as.numeric()
lung_clinical[is.na(lung_clinical$AGE),"AGE"] <- lung_average

lung_clinical_clean <- lung_clinical %>% mutate(
  SUBJID=SAMPLE_ID %>% gsub("-",".",.,fixed = T)) %>%
  select(SUBJID,SEX,AGE)



# lung_clinical_clean %>% write.csv("data/lung_clinical_clean.csv")
lung_clinical_clean %>% write.csv("data/LUAD_clinical_clean.csv")


### clean colon clinical

colon_clinical_samp <- read.table("data/coadread_tcga.tar/coadread_tcga/data_clinical_sample.txt",skip=4,
                                 sep="\t",header=1)
colon_clinical_pat <- read.table("data/coadread_tcga.tar/coadread_tcga/data_clinical_patient.txt",skip=4,
                                sep="\t",header=1)

colon_clinical <- left_join(colon_clinical_samp,colon_clinical_pat,by="PATIENT_ID")
colon_average <- mean(as.numeric(colon_clinical$AGE),na.rm = T) %>% as.integer()
colon_clinical$AGE <- colon_clinical$AGE %>% as.numeric()
colon_clinical[is.na(colon_clinical$AGE),"AGE"] <- colon_average

colon_clinical_clean <- colon_clinical %>% mutate(SUBJID=SAMPLE_ID %>% gsub("-",".",.,fixed = T)) %>%
  select(SUBJID,SEX,AGE)

# colon_clinical_clean %>% write.csv("data/colon_clinical_clean.csv")
colon_clinical_clean %>% write.csv("data/COAD_clinical_clean.csv")


  
