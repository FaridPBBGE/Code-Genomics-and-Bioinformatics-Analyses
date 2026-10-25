## Installing karyoploteR package and loading the packages
BiocManager::install("karyoploteR",force = TRUE)
BiocManager::install("GenomicFeatures", force = TRUE)

library(karyoploteR)
library(GenomicFeatures)

## Creating TxDb object from gff3 file
txdb_rice<- makeTxDbFromGFF("osa1_r7.all_models.gff3", format = "gff3")

genes_data_rice<- makeGenesDataFromTxDb(txdb = txdb_rice, )