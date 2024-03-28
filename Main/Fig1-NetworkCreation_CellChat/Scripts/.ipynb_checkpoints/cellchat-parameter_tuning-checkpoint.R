library(Seurat)
library(SeuratDisk)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)
library(CellChat)
library(patchwork)
library(graphics)
library(Matrix)
options(stringsAsFactors = FALSE)
library(glue)
args = commandArgs(trailingOnly=TRUE)

# Parameter Tuning
condition <- args[1]               # c("CON", "AD")
type <- args[2]                    # c("triMean", "truncatedMean", "thresholdedMean", "median")
trim <- as.numeric(args[3])        # c(0.05, 0.10, 0.15, 0.20, 0.25)

#loading count matrix and meta data
seurat_object <- LoadH5Seurat(glue("/extra/zhanglab0/CommonData/AMP-AD/AMP-AD_restricted/AMP-AD_ROSMAP_annotated.h5seurat"), meta.data = FALSE, misc = FALSE)
meta <- read.csv(file = glue('/extra/zhanglab0/CommonData/AMP-AD/AMP-AD_restricted/meta_withDisorderInfo.csv'),sep = ',',row.names='barcodekey')
seurat_object <- AddMetaData(object = seurat_object, metadata = meta)
print(seurat_object)

#get count matrix
data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data")

#create cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "subclass_renamed")

#choose database
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
#set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

#do parallel
library(future)
options(future.globals.maxSize = +Inf)
future::plan("multicore", workers = 16)

#analysis
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

######################################
## KEY LINE ##########################
cellchat <- computeCommunProb(cellchat, type = type, trim = trim)
######################################

#Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

#Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
nPatterns = 3  #default, set to 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

saveRDS(cellchat, file = "subclass_cellchat_alz.rds")

print("\n\n\n")
print(paste0("Condition: ", condition))
print(paste0("Type: ", type))
print(paste0("Trim: ", trim))
print("DONE!")