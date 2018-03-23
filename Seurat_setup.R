########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
sample1 <- c("SO1A","SO1B","SO2A","SO2B")
sample2 <- c("S11","S12","S21","S22","S31","S32","S1A","S1B","S2A","S2B","S3A","S3B")
sample3 <- c("SYM1A","SYM1B","SYM2A","SYM2B")
group <- c("aged_male1","aged_male2","young_female1","young_female2",
           "young_male")
condition <- c("aged_male","young_female","young_male")
groups <- c(rep(group[1:2],each =4),rep(group[3:4],each =6),rep(group[5],each =4))
conditions <- c(rep(condition[1],each =8),rep(condition[2],each =12),rep(condition[3],each =4))
sample <- c(rep(sample1,time =2),sample2,sample3)
length(groups);length(sample);length(conditions)
pwd <- getwd();paths <- c(); samples <- c()
for(i in c(1:24)){
        paths[i] <- paste0(pwd,"/data/", groups[i],"/out_",sample[i],
                           "/out_gene_exon_tagged.dge_",sample[i],".txt.gz")
        samples[i] <- paste0(groups[i],"_",sample[i])
}

Thymus_raw <- lapply(paths, read.table, header=T,row.names = "GENE",stringsAsFactors=F)

for(i in 1:length(samples)) colnames(Thymus_raw[[i]]) <- paste0(samples[i],
                                                            ".",colnames(Thymus_raw[[i]])) #. is easier for sub
Thymus_list <- lapply(Thymus_raw, CreateSeuratObject,
                 min.cells = 3,
                 min.genes = 200,
                 project = "DropSeq")
Thymus_list <- lapply(Thymus_list, FilterCells, 
                 subset.names = "nGene", 
                 low.thresholds = 500, 
                 high.thresholds = Inf)
Thymus_list <- lapply(Thymus_list, NormalizeData)
Thymus_list <- lapply(Thymus_list, FindVariableGenes, do.plot = FALSE)
Thymus_list <- lapply(Thymus_list, ScaleData)
for(i in 1:length(samples)) Thymus_list[[i]]@meta.data$conditions <- conditions[i]

# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
genes.use <- lapply(Thymus_list, function(x) head(rownames(x@hvg.info), 500))
genes.use <- unique(unlist(genes.use))
for(i in 1:length(samples)){
        genes.use <- intersect(genes.use, rownames(Thymus_list[[i]]@scale.data))
}
length(genes.use) # 1/10 of total sample size

#======1.2 Perform a canonical correlation analysis (CCA) =========================
# run a canonical correlation analysis to identify common sampleurces
# of variation between the two datasets.
Thymus <- RunMultiCCA(object.list = Thymus_list, 
                      genes.use = genes.use,
                      niter = 25, num.ccs = 30,
                      standardize =TRUE)
save(Thymus, file = paste0(pwd,"/data/Thymus_alignment.Rda"))
remove(Thymus_raw)
remove(Thymus_list)

# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = Thymus, reduction.use = "cca", 
              group.by = "conditions", pt.size =1, 
              do.return = TRUE)
p2 <- VlnPlot(object = Thymus, features.plot = "CC1", 
              group.by = "conditions", do.return = TRUE)
png(paste0(pwd,'/output/Dim_violin_plots.png'))
plot_grid(p1, p2)
dev.off()
png(paste0(pwd,'/output/MetageneBicorPlot.png'))
p3 <- MetageneBicorPlot(Thymus, grouping.var = "conditions", dims.eval = 1:30, 
                        display.progress = FALSE)
dev.off()

#======1.3 QC ==================================
# Run rare non-overlapping filtering
Thymus <- CalcVarExpRatio(object = Thymus, reduction.type = "pca",
                      grouping.var = "conditions", dims.use = 1:15)
Thymus <- SubsetData(Thymus, subset.name = "var.ratio.pca",accept.low = 0.5)


#======1.4 align seurat objects =========================
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned

Thymus <- AlignSubspace(object = Thymus, reduction.type = "cca", grouping.var = "conditions", 
                        dims.align = 1:15)
#Now we can run a single integrated analysis on all cells!
Thymus <- FindClusters(object = Thymus, reduction.type = "cca.aligned", dims.use = 1:15, 
                       resolution = 0.8, force.recalc = TRUE, save.SNN = TRUE)
Thymus <- RunTSNE(object = Thymus, reduction.use = "cca.aligned", dims.use = 1:15, 
                  dim.embed = 2, do.fast = TRUE)
p1 <- TSNEPlot(Thymus, do.return = T, pt.size = 1, group.by = "conditions")
p2 <- TSNEPlot(Thymus, do.label = F, do.return = T, pt.size = 1)
png(paste0(pwd,'/output/TSNES_plots.png'))
plot_grid(p1, p2)
dev.off()
#Now, we annotate the clusters as before based on canonical markers.
png(paste0(pwd,'/output/TheClusters.png'))
TSNEPlot(object = Thymus,do.label = TRUE, group.by = "ident", 
               do.return = TRUE, no.legend = TRUE,
               pt.size = 1,label.size = 8 )+
        ggtitle("Tsne Plot of all clusters")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
dev.off()
save(Thymus, file = paste0(pwd,"/data/Thymus_alignment.Rda"))

# Compare clusters for each dataset
cell.all <- FetchData(Thymus,"conditions")
cell.subsets <- lapply(condition, function(x) 
        rownames(cell.all)[cell.all$conditions == x])

Thymus.subsets <- list()
for(i in 1:length(condition)){
        Thymus.subsets[[i]] <- SubsetData(Thymus, cells.use =cell.subsets[[i]])
}

table(Thymus.subsets[[1]]@ident)

p <- list()
for(i in 1:length(condition)){
        p[[i]] <- TSNEPlot(object = Thymus.subsets[[i]],do.label = F, group.by = "ident", 
                           do.return = TRUE, no.legend = TRUE,
                           pt.size = 1,label.size = 4 )+
                ggtitle(condition[i])+
                theme(text = element_text(size=20),     #larger text including legend title							
                      plot.title = element_text(hjust = 0.5)) #title in middle
}
do.call(plot_grid, p)
