cellType_pathways_heatmap <- function() {
    data_whole.CT  <- cellchat.CT@netP$prob
    data_whole.DIS <- cellchat.DIS@netP$prob * normalization
    
    if (comm_type == "receiver") {
        comm_slice.CT  <- data_whole.CT[, cell_type, ]
        comm_slice.DIS <- data_whole.DIS[, cell_type, ]
        row_title = paste0("Sender Cell Types to ", cell_type)
    } else if (comm_type == "sender") {
        comm_slice.CT  <- data_whole.CT[cell_type, , ]
        comm_slice.DIS <- data_whole.DIS[cell_type, , ]
        row_title = paste0(cell_type, " to Receiver Cell Types")
    }
    
    if (length(setdiff(colnames(comm_slice.CT), colnames(comm_slice.DIS))) > 0) {
        original_names <- colnames(comm_slice.DIS)
        for (i in 1:length(setdiff(colnames(comm_slice.CT), colnames(comm_slice.DIS)))) {
            comm_slice.DIS <- cbind(comm_slice.DIS, 0)
        }
        colnames(comm_slice.DIS) <- c(original_names, setdiff(colnames(comm_slice.CT), colnames(comm_slice.DIS)))
    }
    
    if (length(setdiff(colnames(comm_slice.DIS), colnames(comm_slice.CT))) > 0) {
        original_names <- colnames(comm_slice.CT)
        for (i in 1:length(setdiff(colnames(comm_slice.DIS), colnames(comm_slice.CT)))) {
            comm_slice.CT <- cbind(comm_slice.CT, 0)
        }
        colnames(comm_slice.CT) <- c(original_names, setdiff(colnames(comm_slice.DIS), colnames(comm_slice.CT)))
    }
    
    comm_slice.DIS <- comm_slice.DIS[, colnames(comm_slice.CT)]
    
    net.diff_comm_slice <- comm_slice.DIS - comm_slice.CT

    net.diff_comm_slice <- net.diff_comm_slice[celltype_name, colSums(net.diff_comm_slice) != 0]
    return(net.diff_comm_slice)
}