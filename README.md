# Data download:

All necessary data can be found in the data folder and is already pre-processed.

The script used to preprocess data can be found here:
https://github.com/AlexHakansson/CS598_Project_GAIN_GTEX/blob/main/Data_Cleaning/project_preprocessing_clean.R

### Original publicly availible data files:

GTEx Data: URL: https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression

The data files used are: gene_reads_2017-06-05_v8_lung.gct gene_reads_2017-06-05_v8_colon_transverse.gct GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt

TCGA: The cancer genome atlas, a database of genomic and clinical data for various cancers. The data is readily available at cbioportal where it has been made easier to download. URL: https://www.cbioportal.org/

The data sets used are Colorectal Adenocarcinoma (TCGA, Firehose Legacy) Lung Adenocarcinoma (TCGA, Firehose Legacy)

# Installing dependencies
```
pip install pyyaml
pip install keras
pip install wandb
pip install tensorflow
```

To run the notebook ensure that you make this your current working directory.
