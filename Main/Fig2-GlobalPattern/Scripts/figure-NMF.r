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

library(circlize)
library(colorspace)
options(repr.plot.width = 12, repr.plot.height = 9, repr.plot.res = 300)

library(pracma)
library(glue)
args = commandArgs(trailingOnly=TRUE)

# Parameter Tuning
type <- args[1]                    # c("triMean", "truncatedMean", "thresholdedMean", "median")
trim <- as.numeric(args[2])        # c(0.05, 0.10, 0.15, 0.20, 0.25)
numVersions <- as.integer(args[3])


################
## Parameters ##
################
plot_dir <- "/gpfs/gibbs/project/girgenti/cl2553/C2C/PTSD_Call/figures/NMF_last2"

# NMF
#versions that will be made for each combination (resulting in that many figures being made)
versions = seq(from = 1, to = numVersions, by = 1) 
# NMF methods to be used, can be "brunet" or "lee"...
methods = list("brunet")
# NMF seeding methods to be used, can be "random" or "nndsvd"...
seeds = list("random")
# NMF number of runs to be used, can be 200, 500, 1000, 2000, 5000, 10000...
nrums = list(200)

inhibitory_celltypes = c('Inh KCNG1', 'Inh LAMP5', 'Inh PVALB', 'Inh SST', 'Inh VIP')
excitatory_celltypes = c('Exc CUX2', 'Exc FEZF2', 'Exc OPRK1', 'Exc RORB')
support_celltypes = c('Astro', 'Endo', 'MG', 'Oligo', 'OPC')

# these are manually ligand ordering, if left bank it will be done automatically using patterns from control
incoming_ligand_order <- c()
outgoing_ligand_order <- c()
#################





###########################
## Dylan's NMF Functions ##
###########################

# ----------------------------
# this is the function to make
# We include this for comparision purposes
# ----------------------------
my_netAnalysis_river <- function(object, slot.name = "netP", pattern = c("outgoing","incoming"), cutoff = 0.5,
                              sources.use = NULL, targets.use = NULL, signaling = NULL,
                              color.use = NULL, color.use.pattern = NULL, color.use.signaling = "#c5c7c9",
                              do.order = FALSE, main.title = NULL,
                              font.size = 5, font.size.title = 12){
  message("Please make sure you have load `library(ggalluvial)` when running this function")
  requireNamespace("ggalluvial")
  #  suppressMessages(require(ggalluvial))
  res.pattern <- methods::slot(object, slot.name)$pattern[[pattern]]
  data1 = res.pattern$pattern$cell
  data2 = res.pattern$pattern$signaling
  if (is.null(color.use.pattern)) {
    nPatterns <- length(unique(data1$Pattern))
    if (pattern == "outgoing") {
      color.use.pattern = ggPalette(nPatterns*2)[seq(1,nPatterns*2, by = 2)]
    } else if (pattern == "incoming") {
      color.use.pattern = ggPalette(nPatterns*2)[seq(2,nPatterns*2, by = 2)]
    }
  }
  if (is.null(main.title)) {
    if (pattern == "outgoing") {
      main.title = "Outgoing communication patterns of secreting cells"
    } else if (pattern == "incoming") {
      main.title = "Incoming communication patterns of target cells"
    }
  }

  if (is.null(data2)) {
    data1$Contribution[data1$Contribution < cutoff] <- 0
    plot.data <- data1
    nPatterns<-length(unique(plot.data$Pattern))
    nCellGroup<-length(unique(plot.data$CellGroup))
    if (is.null(color.use)) {
      color.use <- scPalette(nCellGroup)
    }
    if (is.null(color.use.pattern)){
      color.use.pattern <- ggPalette(nPatterns)
    }

    plot.data.long <- to_lodes_form(plot.data, axes = 1:2, id = "connection")
    if (do.order) {
      mat = tapply(plot.data[["Contribution"]], list(plot.data[["CellGroup"]], plot.data[["Pattern"]]), sum)
      d <- dist(as.matrix(mat))
      hc <- hclust(d, "ave")
      k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
      cluster <- hc %>% cutree(k)
      order.name <- order(cluster)
      plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(names(cluster)[order.name], colnames(mat)))
      color.use <- color.use[order.name]
    }
    color.use.all <- c(color.use, color.use.pattern)
    gg <- ggplot(plot.data.long,aes(x = factor(x, levels = c("CellGroup", "Pattern")),y=Contribution,
                                    stratum = stratum, alluvium = connection,
                                    fill = stratum, label = stratum)) +
      geom_flow(width = 1/3,aes.flow = "backward") +
      geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
      geom_text(stat = "stratum", size = font.size) +
      scale_x_discrete(limits = c(),  labels=c("Cell groups", "Patterns")) +
      scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size=10))+
      ggtitle(main.title)

  } else {
    data1$Contribution[data1$Contribution < cutoff] <- 0
    plot.data <- data1
    nPatterns<-length(unique(plot.data$Pattern))
    nCellGroup<-length(unique(plot.data$CellGroup))
    cells.level = levels(object@idents)
    if (is.null(color.use)) {
      color.use <- scPalette(length(cells.level))[cells.level %in% unique(plot.data$CellGroup)]
    }
    if (is.null(color.use.pattern)){
      color.use.pattern <- ggPalette(nPatterns)
    }
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      plot.data <- subset(plot.data, CellGroup %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      plot.data <- subset(plot.data, CellGroup %in% targets.use)
    }
    ## connect cell groups with patterns
    plot.data.long <- to_lodes_form(plot.data, axes = 1:2, id = "connection")
    if (do.order) {
      mat = tapply(plot.data[["Contribution"]], list(plot.data[["CellGroup"]], plot.data[["Pattern"]]), sum)
      d <- dist(as.matrix(mat))
      hc <- hclust(d, "ave")
      k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
      cluster <- hc %>% cutree(k)
      order.name <- order(cluster)
      plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(names(cluster)[order.name], colnames(mat)))
      color.use <- color.use[order.name]
    }


    data_subset <- subset(data1, data1$Contribution > 0.4)

    data_subset

    data_subset[["Pattern"]]
    # length(patterns)

    colors = list()


    if (pattern == "outgoing") {
        for (cell_type in data_subset[["Pattern"]]){
        if(cell_type == "Pattern 1"){
        colors <- append(colors, "#f9918a")
        }
        if(cell_type == "Pattern 2"){
        colors <- append(colors,"#33c860")
        }
        if(cell_type == "Pattern 3"){
        colors <- append(colors,"#81b0ff")
        }
    }
    colors <- append(colors, "#f9918a")
    colors <- append(colors,"#33c860")
    colors <- append(colors,"#81b0ff")
    } 
    else if (pattern == "incoming") {
        for (cell_type in data_subset[["Pattern"]]){
            if(cell_type == "Pattern 1"){
                colors <- append(colors, "#c5b233")
            }
            if(cell_type == "Pattern 2"){
                colors <- append(colors,"#33ccd0")
            }
            if(cell_type == "Pattern 3"){
                colors <- append(colors,"#f783e9")
            }
        }
        colors <- append(colors, "#c5b233")
        colors <- append(colors,"#33ccd0")
        colors <- append(colors,"#f783e9")
    }


    names(colors) <- data_subset[["CellGroup"]]
    names(colors)[length(colors)-2] <- "Pattern 1"
    names(colors)[length(colors)-1] <- "Pattern 2"
    names(colors)[length(colors)] <- "Pattern 3"
    # color.use <- colors

    color.use.all <- c(color.use,color.use.pattern)

    color.use.all <- colors

    message(color.use.all)

    # color.use.all <- colors


    StatStratum <- ggalluvial::StatStratum
    gg1 <- ggplot(plot.data.long,aes(x = factor(x, levels = c("CellGroup", "Pattern")),y=Contribution,
                                     stratum = stratum, alluvium = connection,
                                     fill = stratum, label = stratum)) +
      geom_flow(width = 1/3,aes.flow = "forward") +
      geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
      geom_text(stat = "stratum", size = font.size) +
      scale_x_discrete(limits = c(),  labels=c("Cell groups", "Patterns")) +

      scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size=10)) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    ## connect patterns with signaling
    data2$Contribution[data2$Contribution < cutoff] <- 0
    plot.data <- data2
    nPatterns<-length(unique(plot.data$Pattern))
    nSignaling<-length(unique(plot.data$Signaling))
    if (length(color.use.signaling) == 1) {
      color.use.all <- c(color.use.pattern, rep(color.use.signaling, nSignaling))
    } else {
      color.use.all <- c(color.use.pattern, color.use.signaling)
    }

    if (!is.null(signaling)) {
      plot.data <- plot.data[plot.data$Signaling %in% signaling, ]
    }

    plot.data.long <- ggalluvial::to_lodes_form(plot.data, axes = 1:2, id = "connection")

    gg2 <- ggplot(plot.data.long,aes(x = factor(x, levels = c("Pattern", "Signaling")),y= Contribution,
                                     stratum = stratum, alluvium = connection,
                                     fill = stratum, label = stratum)) +
      geom_flow(width = 1/3,aes.flow = "forward") +
      geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
      geom_text(stat = "stratum", size = font.size) + # 2.5
      scale_x_discrete(limits = c(),  labels=c("Patterns", "Signaling")) +
      scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size= 10))+
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

    gg <- cowplot::plot_grid(gg1, gg2,align = "h", nrow = 1)
    title <- cowplot::ggdraw() + cowplot::draw_label(main.title,size = font.size.title)
    gg <- cowplot::plot_grid(title, gg, ncol=1, rel_heights=c(0.1, 1))
  }
  return(gg)
}





#--------------------------------------------
# Custom merged figured.
# This is what we will use on the file paper
#--------------------------------------------
my_merged_river <- function(object,pattern, disease, inhibitory_celltypes, excitatory_celltypes, support_celltypes){

    if (pattern == "outgoing") {
      main.title = "Outgoing communication patterns of secreting cells"
      # The color palette is defined in the function ggPalette
      # color.use just needs three colors you can define them in a list with hex codes if you want 
      color.use <- ggPalette(3*2)[seq(1,3*2, by = 2)]
      
    } else if (pattern == "incoming") {
      main.title = "Incoming communication patterns of target cells"
      # The color palette is defined in the function ggPalette
      # color.use just needs three colors you can define them in a list with hex codes if you want 
      color.use <- ggPalette(3*2)[seq(2,3*2, by = 2)]
    }

    slot.name = "netP"
    res.pattern <- methods::slot(object, slot.name)$pattern[[pattern]]

    data1 = res.pattern$pattern$cell
    data2 = res.pattern$pattern$signaling

    sorted_group <- data2[order(data2$Signaling, -data2$Contribution),]
    data2 <- sorted_group[!duplicated(sorted_group$Signaling),]

    sorted_group <- data1[order(data1$CellGroup, -data1$Contribution),]
    data1 <- sorted_group[!duplicated(sorted_group$CellGroup),]

    plot.data <- merge(data1,data2,by="Pattern")

    plot.data$Signaling <- as.character(plot.data$Signaling)
    print(unique(plot.data$Signaling))
    plot.data <- plot.data[order(plot.data$Pattern, plot.data$Signaling),]

    # -- this is code for pulling out custom ligands
    # ligand_in_order = ligand_in_order[!ligand_in_order %in% c("FGF", "NRG", "PSAP")]
    # ligand_in_order = c(ligand_in_order, "FGF", "NRG", "PSAP")
    # plot.data$Signaling <- factor(plot.data$Signaling, levels = ligand_in_order)

    # order ligands by pattern
    if(pattern == "outgoing"){
      if(length(outgoing_ligand_order) > 0){
          plot.data$Signaling <- factor(plot.data$Signaling, levels = outgoing_ligand_order) 
      } 
      else {
        if(disease == "control"){
            outgoing_ligand_order <<- c()
            for (pattern in unique(plot.data$Pattern)){
              outgoing_ligand_order <<- c(outgoing_ligand_order, unique(plot.data$Signaling[plot.data$Pattern == pattern]))
            }
            plot.data$Signaling <- factor(plot.data$Signaling, levels = outgoing_ligand_order) 
        }
      }
    }
    else if(pattern == "incoming"){
      if(length(incoming_ligand_order) > 0){
          plot.data$Signaling <- factor(plot.data$Signaling, levels = incoming_ligand_order) 
      } 
      else {
        if(disease == "control"){
            incoming_ligand_order <<- c()
            for (pattern in unique(plot.data$Pattern)){
              incoming_ligand_order <<- c(incoming_ligand_order, unique(plot.data$Signaling[plot.data$Pattern == pattern]))
            }
            plot.data$Signaling <- factor(plot.data$Signaling, levels = incoming_ligand_order) 
        }
      }
    }

    #order cell types by cell type groups
    # inhibitory_celltypes <- c("Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sst", "Sst Chodl", "Vip", "Sncg")
    inhibitory_celltypes <- inhibitory_celltypes[order(inhibitory_celltypes)]
    # excitatory_celltypes <- c("L2/3 IT", "L4 IT", "L5", "L5 ET", "L5/6 ET", "L6 CT", "L6 IT", "L6 IT Car3","L5 IT","L5/6 NP","L6b")
    excitatory_celltypes <- excitatory_celltypes[order(excitatory_celltypes)]
    # support_celltypes <- c("Astro","Endo","Micro/PVM", "Oligo", "OPC", "VLMC")

    celltypesOrder <- c(support_celltypes, inhibitory_celltypes, excitatory_celltypes)
    plot.data$CellGroup <- factor(plot.data$CellGroup, levels = celltypesOrder)

    #scale by cellgroup
    goal_sum <- 1
    for(celltype in unique(plot.data$CellGroup)){
        goal_sum <- pracma::Lcm(
            goal_sum,
            sum(plot.data$CellGroup == celltype)
        )
    }

    print(goal_sum)

    for(celltype in unique(plot.data$CellGroup)){
        index_to_copy <- which(match(plot.data$CellGroup, celltype) == 1)
        while(sum(plot.data$CellGroup == celltype) != goal_sum){
            for (index in index_to_copy){
                    plot.data <- rbind(plot.data, plot.data[rep(index, 1), ])

            }
        }
    }

    # Some signaling is NA, this will cause angles to be different length.
    plot.data <- subset(plot.data, !is.na(plot.data$Signaling))

    #make Pattern text Vertical
    angles <- rep(0, length(unique(plot.data$Pattern)) + length(unique(plot.data$CellGroup)) + length(unique(plot.data$Signaling)))
    for (i in seq(length(unique(plot.data$CellGroup))+1,length(unique(plot.data$CellGroup))+3)){
        angles[i] <- 90
    }

    #make our plot
    gg <- ggplot(plot.data,
        aes(axis1 = CellGroup, axis2 = Pattern, axis3 = Signaling)) +
    geom_flow(width = 1/3, aes.flow = "forward", aes(fill = Pattern)) + 
    scale_x_discrete(limits = c("Cell groups", "Pattern", "Signaling" )) +
    geom_stratum(alpha = 0.8, aes(fill = Pattern),width = 1/3, size=0.1) + 
    scale_fill_manual(values = color.use) +
    geom_text(angle=angles, size = 4, stat = "stratum", aes(label = after_stat(stratum))) +
    theme_bw()+
        theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size= 10)) + 
            theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(main.title)
    return(gg)
}







# ----------------------------
# This function takes in cellchat object and is responsibule for doing all the runs
# and making all the figures
# ----------------------------
NMF_and_river <- function(cellchat_object,study,disease, inhibitory_celltypes, excitatory_celltypes, support_celltypes, chosen_version, parallel='p16') {
    for (chosen_seed in seeds){
        for (chosen_nrun in nrums){
            for (chosen_method in methods){
                try({
                slot.name = "netP"
                pattern = "outgoing"
                object = cellchat_object
                k = 3 #three patterns

                prob <- methods::slot(object, slot.name)$prob

                data_sender <- apply(prob, c(1,3), sum)
                data_sender = sweep(data_sender, 2L, apply(data_sender, 2, function(x) max(x, na.rm = TRUE)), '/', check.margin = FALSE)
                data0 = as.matrix(data_sender)

                options(warn = -1)
                data <- data0
                data <- data[rowSums(data)!=0,]
                print("beginning NMF")
                outs_NMF <- NMF::nmf(data, rank = 3, method = chosen_method, seed = chosen_seed, nrun=chosen_nrun, .opt=parallel) 
                W <- scaleMat(outs_NMF@fit@W, 'r1')
                H <- scaleMat(outs_NMF@fit@H, 'c1')
                print("Finished NMF")
                colnames(W) <- paste0("Pattern ", seq(1,ncol(W))); rownames(H) <- paste0("Pattern ", seq(1,nrow(H)));

                data_W <- as.data.frame(as.table(W)); colnames(data_W) <- c("CellGroup","Pattern","Contribution")
                data_H <- as.data.frame(as.table(H)); colnames(data_H) <- c("Pattern","Signaling","Contribution")

                res.pattern = list("cell" = data_W, "signaling" = data_H)
                methods::slot(object, slot.name)$pattern[[pattern]] <- list(data = data0, pattern = res.pattern)

                base_file_name = paste0("outgoing_", chosen_nrun, "runs-", glue("type_{type}-trim_{trim}-"), "v", chosen_version, ".pdf")
                rds_name = paste0("outgoing", chosen_nrun, "_runs_", chosen_version, ".rds")
                folder_name = paste(study,disease,sep="/")
                file_name = paste(chosen_method,chosen_seed,base_file_name,sep="_")
                print("creating directory")
                file_dir = paste(folder_name,chosen_method,chosen_seed,sep="/")
                dir.create(folder_name, recursive = TRUE)
                rds_path = paste(folder_name,chosen_method,chosen_seed,rds_name,sep="/")
                file_path = paste(folder_name, file_name, sep="/")
                # png(filename=file_path)
                colors = c("green","green","green","green","green","green","green","green","green","green","green","green","green")
                my_merged_river(object, pattern = "outgoing", disease, 
                                inhibitory_celltypes=inhibitory_celltypes, 
                                excitatory_celltypes=excitatory_celltypes, 
                                support_celltypes=support_celltypes)
                print("saving fig:")
                print(file_path)
                ggsave(file_path, device="pdf")
                # saveRDS(object, file =rds_path)
                print("done saving fig")


                # base_file_name = paste0("ORG_outgoing_", chosen_nrun, "runs-", glue("type_{type}-trim_{trim}-"), "v", chosen_version, ".pdf")
                # folder_name = paste(study,disease,sep="/")
                # file_name = paste(chosen_method,chosen_seed,base_file_name,sep="_")
                # file_path = paste(folder_name, file_name, sep="/")
                # # png(filename=file_path)
                # colors = c("green","green","green","green","green","green","green","green","green","green","green","green","green")
                # my_netAnalysis_river(object, pattern = "outgoing",cutoff=0.5,do.order=TRUE)
                # print("saving fig:")
                # print(file_path)
                # ggsave(file_path, device="pdf")
                # print("done saving fig")

                slot.name = "netP"
                pattern = "incoming"
                object = cellchat_object
                k = 3

                prob <- methods::slot(object, slot.name)$prob

                data_sender <- apply(prob, c(1,3), sum)
                data_sender = sweep(data_sender, 2L, apply(data_sender, 2, function(x) max(x, na.rm = TRUE)), '/', check.margin = FALSE)
                data0 = as.matrix(data_sender)

                options(warn = -1)
                data <- data0
                data <- data[rowSums(data)!=0,]

                outs_NMF <- NMF::nmf(data, rank = 3, method = chosen_method, seed = chosen_seed, nrun=chosen_nrun, .opt=parallel) 
                W <- scaleMat(outs_NMF@fit@W, 'r1')
                H <- scaleMat(outs_NMF@fit@H, 'c1')

                colnames(W) <- paste0("Pattern ", seq(1,ncol(W))); rownames(H) <- paste0("Pattern ", seq(1,nrow(H)));

                data_W <- as.data.frame(as.table(W)); colnames(data_W) <- c("CellGroup","Pattern","Contribution")
                data_H <- as.data.frame(as.table(H)); colnames(data_H) <- c("Pattern","Signaling","Contribution")

                res.pattern = list("cell" = data_W, "signaling" = data_H)
                methods::slot(object, slot.name)$pattern[[pattern]] <- list(data = data0, pattern = res.pattern)


                base_file_name = paste0("incoming_", chosen_nrun, "runs-", glue("type_{type}-trim_{trim}-"), "v", chosen_version, ".pdf")
                rds_name = paste0("incoming", chosen_nrun, "_runs_", chosen_version, ".rds")
                folder_name = paste(study,disease,sep="/")
                file_name = paste(chosen_method,chosen_seed,base_file_name,sep="_")

                file_dir = paste(folder_name,chosen_method,chosen_seed,sep="/")
                dir.create(folder_name, recursive = TRUE) 
                print("making directory")
                file_path = paste(folder_name, file_name, sep="/")
                rds_path = paste(folder_name,chosen_method,chosen_seed,rds_name,sep="/")
                # png(filename=file_path)
                colors = c("green","green","green","green","green","green","green","green","green","green","green")
                my_merged_river(object, pattern = "incoming", disease, 
                                inhibitory_celltypes=inhibitory_celltypes, 
                                excitatory_celltypes=excitatory_celltypes, 
                                support_celltypes=support_celltypes)
                print("saving fig:")
                print(file_path)
                ggsave(file_path, device="pdf")
                # saveRDS(object, file =rds_path)
                print("done saving fig")

                # base_file_name = paste0("ORG_incoming_", chosen_nrun, "runs-", glue("type_{type}-trim_{trim}-"), "v", chosen_version, ".pdf")
                # folder_name = paste(study,disease,sep="/")
                # file_name = paste(chosen_method,chosen_seed,base_file_name,sep="_")
                # file_path = paste(folder_name, file_name, sep="/")
                # # png(filename=file_path)
                # my_netAnalysis_river(object, pattern = "incoming",cutoff=0.5,do.order=TRUE)
                # print("saving fig:")
                # print(file_path)
                # ggsave(file_path, device="pdf")
                # print("done saving fig")

                })
            }
        }
    }

}

###########################





###############
## Execution ##
###############
                                        
# --------------------------------------------
# run the function
# --------------------------------------------
# !! IMPORTANT !!
# RIGHT NOW FOR LIGAND ORDERING TO WORK YOU MUST RUN 
# CONTROL FIRST (i'm working on fixing this)
# --------------------------------------------

print(glue("/gpfs/gibbs/pi/girgenti/JZhang/CL/C2C/PTSD_Call/Obj-CellChat/RNA_FINAL-CON-type_{type}-trim_{trim}_cellchat.rds"))
cellchat.CT <- readRDS(file = glue("/gpfs/gibbs/pi/girgenti/JZhang/CL/C2C/PTSD_Call/Obj-CellChat/RNA_FINAL-CON-type_{type}-trim_{trim}_cellchat.rds"))
print(cellchat.CT)
cat("\n")

print(glue("/gpfs/gibbs/pi/girgenti/JZhang/CL/C2C/PTSD_Call/Obj-CellChat/RNA_FINAL-PTSD-type_{type}-trim_{trim}_cellchat.rds"))
cellchat.PTSD <- readRDS(file = glue("/gpfs/gibbs/pi/girgenti/JZhang/CL/C2C/PTSD_Call/Obj-CellChat/RNA_FINAL-PTSD-type_{type}-trim_{trim}_cellchat.rds"))
print(cellchat.PTSD)
cat("\n")

print(glue("/gpfs/gibbs/pi/girgenti/JZhang/CL/C2C/PTSD_Call/Obj-CellChat/RNA_FINAL-MDD-type_{type}-trim_{trim}_cellchat.rds"))
cellchat.MDD <- readRDS(file = glue("/gpfs/gibbs/pi/girgenti/JZhang/CL/C2C/PTSD_Call/Obj-CellChat/RNA_FINAL-MDD-type_{type}-trim_{trim}_cellchat.rds"))
print(cellchat.MDD)
cat("\n")

                                        
for (chosen_version in versions){
    incoming_ligand_order <- c()
    outgoing_ligand_order <- c()
    NMF_and_river(cellchat.CT, plot_dir, "control", inhibitory_celltypes, excitatory_celltypes, support_celltypes, chosen_version, parallel='p16')
    NMF_and_river(cellchat.PTSD, plot_dir, "PTSD", inhibitory_celltypes, excitatory_celltypes, support_celltypes, chosen_version, parallel='p16')
    NMF_and_river(cellchat.MDD, plot_dir, "MDD", inhibitory_celltypes, excitatory_celltypes, support_celltypes, chosen_version, parallel='p16')
}
