library(NeuronChat)
library(CellChat)
library(ggalluvial)
library(glue)
library(Seurat)
library(SeuratDisk)
library(ComplexHeatmap)
library(circlize)

NMF_and_river_OUT <- function(neuronchat_object, study, disease, inhibitory_celltypes, excitatory_celltypes, support_celltypes, chosen_version, quant_threshold, k, parallel='p16') {
    for (chosen_seed in seeds){
        for (chosen_nrun in nruns){
            for (chosen_method in methods){
                try({
                slot.name = "net"
                pattern = "outgoing"
                object = neuronchat_object

                prob <- methods::slot(object, slot.name)$prob
                # prob <- simplify2array(methods::slot(object, slot.name))                                       # new

                data_sender <- apply(prob, c(1,3), sum)

                ########################################
                ### Normalization
                ########################################
                data_sender = sweep(data_sender, 2L, apply(data_sender, 2, function(x) max(x, na.rm = TRUE)), '/', check.margin = FALSE)
                # if (normalization_by > 0) {
                #     data_sender = sweep(data_sender, 2L, apply(data_sender, 2, function(x) quantile(x[x != 0], probs = normalization_by)), '/', check.margin = FALSE)
                # }
                ########################################
                data_sender[is.na(data_sender)] <- 0                                                           # new
                data0 = as.matrix(data_sender)

                options(warn = -1)
                data <- data0
                data <- data[rowSums(data)!=0,colSums(data)!=0]                                                # new

                ########################################
                ### Diff Quant Threshold to filter out Signaling
                ########################################
                data <- data[,colSums(data) >= quantile(colSums(data), quant_threshold)]                       # new
                ########################################
                                                           
                #########################################
                ### Removing ODC Pairs
                #########################################
                # data <- data[,!colnames(data) %in% to_remove]
                #########################################
                                                           
                data <- data[rowSums(data)!=0,]
                                                           
                print("beginning NMF")
                outs_NMF <- NMF::nmf(data, rank = k, method = chosen_method, seed = chosen_seed, nrun=chosen_nrun, .opt=parallel)  
                W <- scaleMat(outs_NMF@fit@W, 'r1')
                H <- scaleMat(outs_NMF@fit@H, 'c1')
                print("Finished NMF")
                colnames(W) <- paste0("Pattern ", seq(1,ncol(W))); rownames(H) <- paste0("Pattern ", seq(1,nrow(H)));

                data_W <- as.data.frame(as.table(W)); colnames(data_W) <- c("CellGroup","Pattern","Contribution")
                data_H <- as.data.frame(as.table(H)); colnames(data_H) <- c("Pattern","Signaling","Contribution")

                res.pattern = list("cell" = data_W, "signaling" = data_H)
                # methods::slot(object, 'net_analysis')$pattern[[pattern]] <- list(data = data0, pattern = res.pattern)
                methods::slot(object, 'net')$pattern[[pattern]] <- list(data = data0, pattern = res.pattern)

                })
            }
        }
    }
    return(object)
}