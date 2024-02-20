library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
options(stringsAsFactors = FALSE)
library(patchwork)

library(SeuratDisk)


seurat_object <- LoadH5Seurat("/gpfs/gibbs/pi/girgenti/JZhang/CL/cellChat/brain_call/new/SZBD-Kellis_annotated-BD_CON.h5seurat", meta.data = FALSE, misc = FALSE)
meta <- read.csv(file = '/gpfs/gibbs/pi/girgenti/JZhang/CL/cellChat/brain_call/new/SZBD-Kellis_annotated-BD_CON.csv',sep = ',',row.names='barcodekey')
seurat_object <- AddMetaData(object = seurat_object, metadata = meta)
print(seurat_object)


# Idents(seurat_object) <- seurat_object@meta.data["Disorder"]
# seurat_object <- subset(seurat_object, idents = c("BD","CON"))
Idents(seurat_object) <- seurat_object@meta.data["subclass"]


ligand_target_matrix = readRDS("/gpfs/gibbs/pi/girgenti/JZhang/CL/cellChat/brain_call/NicheNet_db/ligand_target_matrix.rds")
lr_network = readRDS("/gpfs/gibbs/pi/girgenti/JZhang/CL/cellChat/brain_call/NicheNet_db/lr_network.rds")
weighted_networks = readRDS("/gpfs/gibbs/pi/girgenti/JZhang/CL/cellChat/brain_call/NicheNet_db/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))





sender_celltypes = c("Astro","Endo","Immune","Micro","Oligo","OPC","PC","SMC","VLMC")




## All Neurons
receiver_celltypes = c("Chandelier","Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Sncg","Sst", "Sst Chodl","Vip",
                       "L2/3 IT","L4 IT","L5 ET","L5 IT","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")


list_expressed_genes_receiver = receiver_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, 0.10)
expressed_genes_receiver= list_expressed_genes_receiver %>% unlist() %>% unique()
saveRDS(expressed_genes_receiver, file = "/gpfs/gibbs/pi/girgenti/JZhang/CL/cellChat/brain_call/SZBD_bipolar_control_expressed_genes_receiver-all.rds")


list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, 0.10)
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
saveRDS(expressed_genes_sender, file = "/gpfs/gibbs/pi/girgenti/JZhang/CL/cellChat/brain_call/SZBD_bipolar_control_expressed_genes_sender-all.rds")




## EXC Neurons
receiver_celltypes = c("L2/3 IT","L4 IT","L5 ET","L5 IT","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")


list_expressed_genes_receiver = receiver_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, 0.10)
expressed_genes_receiver= list_expressed_genes_receiver %>% unlist() %>% unique()
saveRDS(expressed_genes_receiver, file = "/gpfs/gibbs/pi/girgenti/JZhang/CL/cellChat/brain_call/SZBD_bipolar_control_expressed_genes_receiver-EXC.rds")


list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, 0.10)
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
saveRDS(expressed_genes_sender, file = "/gpfs/gibbs/pi/girgenti/JZhang/CL/cellChat/brain_call/SZBD_bipolar_control_expressed_genes_sender-EXC.rds")



## INH Neurons
receiver_celltypes = c("Chandelier","Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Sncg","Sst", "Sst Chodl","Vip")


list_expressed_genes_receiver = receiver_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, 0.10)
expressed_genes_receiver= list_expressed_genes_receiver %>% unlist() %>% unique()
saveRDS(expressed_genes_receiver, file = "/gpfs/gibbs/pi/girgenti/JZhang/CL/cellChat/brain_call/SZBD_bipolar_control_expressed_genes_receiver-INH.rds")


list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, 0.10)
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
saveRDS(expressed_genes_sender, file = "/gpfs/gibbs/pi/girgenti/JZhang/CL/cellChat/brain_call/SZBD_bipolar_control_expressed_genes_sender-INH.rds")
