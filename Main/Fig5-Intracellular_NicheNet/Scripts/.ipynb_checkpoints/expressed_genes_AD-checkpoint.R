## Import Libraries
library(nichenetr)
library(Seurat)
library(tidyverse)
options(stringsAsFactors = FALSE)
library(patchwork)
library(dplyr)

library(SeuratDisk)

## Load Seurat Object
seurat_object <- LoadH5Seurat("/extra/zhanglab0/CommonData/AMP-AD/AMP-AD_restricted/AMP-AD_ROSMAP_annotated.h5seurat",assays = "RNA" , meta.data = FALSE, misc = FALSE)

## Add meta data to seurat_object
meta <- read.csv(file = '/extra/zhanglab0/CommonData/AMP-AD/AMP-AD_restricted/meta_withDisorderInfo.csv',sep = ',',row.names='barcodekey')
meta <- meta %>% mutate(subclass = recode(subclass, "L5 ET" = 'L5', "L5 IT" = 'L5')) #Update specific column categories in meta data
seurat_object <- AddMetaData(object = seurat_object, metadata = meta)

## Debug 
# unique(seurat_object@meta.data["subclass"]$subclass)

## Set Idents(seurat_object)
Idents(seurat_object) <- seurat_object@meta.data["subclass"]$subclass

## Sender Neurons
sender_celltypes = c("Astro", "Micro/PVM")

## All Neurons
receiver_celltypes = c("Chandelier","Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Sncg","Sst", "Sst Chodl","Vip",
                       "L2/3 IT","L4 IT","L5","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")

list_expressed_genes_receiver = receiver_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, 0.10)
expressed_genes_receiver= list_expressed_genes_receiver %>% unlist() %>% unique()
saveRDS(expressed_genes_receiver, file = "../Notebooks/Alz_Disease/AD_expressed_genes_receiver-all.rds")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, 0.10)
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
saveRDS(expressed_genes_sender, file = "../Notebooks/Alz_Disease/AD_expressed_genes_sender-all.rds")

## EXC Neurons
receiver_celltypes = c("L2/3 IT","L4 IT","L5","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")

list_expressed_genes_receiver = receiver_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, 0.10)
expressed_genes_receiver= list_expressed_genes_receiver %>% unlist() %>% unique()
saveRDS(expressed_genes_receiver, file = "../Notebooks/Alz_Disease/AD_expressed_genes_receiver-EXC.rds")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, 0.10)
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
saveRDS(expressed_genes_sender, file = "../Notebooks/Alz_Disease/AD_expressed_genes_sender-EXC.rds")

## INH Neurons
receiver_celltypes = c("Chandelier","Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Sncg","Sst", "Sst Chodl","Vip")

list_expressed_genes_receiver = receiver_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, 0.10)
expressed_genes_receiver= list_expressed_genes_receiver %>% unlist() %>% unique()
saveRDS(expressed_genes_receiver, file = "../Notebooks/Alz_Disease/AD_expressed_genes_receiver-INH.rds")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, 0.10)
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
saveRDS(expressed_genes_sender, file = "../Notebooks/Alz_Disease/AD_expressed_genes_sender-INH.rds")

## Session Info
# sessionInfo()