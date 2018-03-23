library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/Thymus_alignment.Rda")
lnames

# Adipocytes
Adipocytes <- MouseGenes(Thymus,c("SLC36A2","P2RX5","TRIP4","ASCC1","MYF5","UCP1"))
FeaturePlot(object = Thymus, 
            features.plot = Adipocytes, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)

# Endothelial Cells
Endothelial <- MouseGenes(Thymus,c("Cdh5","Pecam1","Flt1","Vwf","Plvap","Kdr",
                          "EMCN","Car4","ptprb"))
FeaturePlot(object = Thymus, 
            features.plot = Endothelial, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)

# Epithelium
Epithelium <- MouseGenes(Thymus,c("KRT19","Epcam","KRT5","MUC1","SCGB3A2","SCGB1A1",
              "SCGB3A1","SFTPB","FOXJ1","Rpe65","Rlbp1","Msln","Upk3b","Lrrn4"))
TECs <- MouseGenes(Thymus,c("Ccl25","Cxcl12","Ctsl","Psmb11","Aire","HLA-DMA","Krt5",
                   "Gas1","Plet1","Ly6d","Spink5","Reg3g","Bpifa1"))
FeaturePlot(object = Thymus, 
            features.plot = TECs, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)

#fibroblasts
fibroblasts <- MouseGenes(Thymus,c("FGF1","FGF9","SFRP1"))
FeaturePlot(object = Thymus, 
            features.plot = fibroblasts, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)

# Hematopoietic cells
Hematopoietic <- MouseGenes(Thymus,c("PTPRC","LAPTM5","SRGN"))
myeloid <- MouseGenes(Thymus,c("PPBP","GNG11","HBA2","HBB","Cma1",
        "Mcpt4","Tpsb2","Cpa3","LYZ","S100A9","CD14","CCL2","Emr1","CD68",
        "MARCO","LYZ","FCGR3A","MS4A7","VMO1","Itgax","GPR183","CST3","HLA-DQA1",
        "FCER1A","TSPAN13","IL3RA","IGJ"))
Tcell <- MouseGenes(Thymus,c("CD3G","CD3D","CD2","CREM","Cd4","CD62L",
        "IL7R","IL2RG","GIMAP5","SELL","CREM","Cd8a","CCL5",
        "CACYBP","Foxp3"))
Bcell <- MouseGenes(Thymus,c("CD19","CD79A","MS4A1","MIR155HG","NME1",
                             "HLA-DQA1","CD27","SDC1","IL6R","SLAMF7"))
NK <- MouseGenes(Thymus,c("GNLY","Ncr1","CCL5","KLRD1","NKG7"))
FeaturePlot(object = Thymus, 
            features.plot = NK, min.cutoff = 1, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)

#Mesenchymal Cells
Mesenchymal <- MouseGenes(Thymus,c("Pdgfrb","Vim","Has2","Dcn"))
Neurons <- MouseGenes(Thymus,c("Rbfox3"))
Pericytes <- MouseGenes(Thymus,c("Pdgfrb","Cspg4","Anpep","Rgs5",
                                 "Myh11","Mylk","Sost","Des","Vtn","Ifitm1"))
Smooth_muscle_cells <- MouseGenes(Thymus,c("Acta2","Myh11"))
FeaturePlot(object = Thymus, 
            features.plot = Smooth_muscle_cells, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)

# Stem cell
Embryonic_Stem_Cell <- MouseGenes(Thymus,c("POU5F1","FUT4"))
Hematopoietic_Stem_Cell <- MouseGenes(Thymus,c("CD34","PROM1","ABCG2","ATXN1","RUNX1"))
Mesenchymal_Stem_Cells <- MouseGenes(Thymus,c("Eng","CD44","Nt5e","Alcam","Kit",
                        "Ly6a","Thy1","Itgb1","Vcam1","Icam2","CD72","CD2"))
Neural_Stem_Cell <- MouseGenes(Thymus,c("NES","NCAM","NGFR"))
FeaturePlot(object = Thymus, 
            features.plot = Neural_Stem_Cell, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)

# Rename ident
table(Thymus@ident)
idents <- as.data.frame(table(Thymus@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("0.Mesenchymal Cells",
                     "1.T cells",
                     "2.Medullary TECs (mTEChi)",
                     "3.Medullary TECs (mTEClo)",
                     "4.T cells",
                     "5.Mesenchymal Cells",
                     "6.Endothelial cells",
                     "7.T cells",
                     "8.T cells & Cortical TECs (cTECs)",
                     "9.Mesothelial Cells",
                     "10.Hematopoietic cells & TEC subset",
                     "11.Mesenchymal Cells",
                     "12.Epithelium",
                     "13.Mesenchymal Cells")

Thymus@ident <- plyr::mapvalues(x = Thymus@ident,
                            from = old.ident.ids,
                            to = new.cluster.ids)
TSNEPlot(object = Thymus, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
        ggtitle("Major cell types")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle

# The SplitDotPlotGG function can be useful for viewing conserved cell type markers
# across conditions, showing both the expression level and the percentage of cells
# in a cluster expressing any given gene. 
# Here we plot 1-3 strong marker genes for each of our 13 clusters.
markers.to.plot <- c("Cdh5","Pecam1","Flt1","Vwf","Plvap","Kdr","EMCN","ptprb","KRT19",
                     "Epcam","PTPRC","LAPTM5","SRGN","CD14","CD68","CD3D","CCL5","Pdgfrb",
                     "Dcn")
markers.to.plot <- MouseGenes(Thymus,markers.to.plot, unique =T)
sort(markers.to.plot)
sdp <- SplitDotPlotGG(Thymus, genes.plot = rev(markers.to.plot),
                      cols.use = c("grey","blue"), x.lab.rot = T, plot.legend = T,
                      dot.scale = 8, do.return = T, grouping.var = "conditions")


# Thymus <- RenameIdentBack(Thymus)
# How many cells are in each cluster
table(Thymus@ident)
TSNEPlot(object = Thymus, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
  ggtitle("Major cell types")+
  theme(text = element_text(size=20),     #larger text including legend title							
        plot.title = element_text(hjust = 0.5)) #title in middle

#=====2.2 - A table with the number of cells of each cluster and subcluster, for both B6 and 129_B6 strains.
# We can also compare proportional shifts in the data. As can be seen in the barplot, 
freq_table <- prop.table(x = table(Thymus@ident, Thymus@meta.data[, "conditions"]), 
                         margin = 2)
barplot(height = freq_table)

freq_table


#====== 2.3 Compare cell type changes across conditions  ==========================================
# the two patients profiled have very different composition
lnames = load(file = "./data/Thymus_alignment.Rda")
idents <- as.data.frame(table(Thymus@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Unknown",
                     "Unknown",
                     "Unknown",
                     "Unknown",
                     "Myeloid cells",
                     "Unknown",
                     "Unknown",
                     "Unknown",
                     "Unknown",
                     "Epithelium",
                     "Mesenchymal Cells",
                     "Unknown",
                     "Myeloid cells",
                     "Lymphoid cells",
                     "Endothelium",
                     "Epithelium",
                     "Endothelium")
Thymus@ident <- plyr::mapvalues(x = Thymus@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
TSNEPlot(object = Thymus, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 6,
         colors.use = c("#DB72FB","#619CFF","#D39200","#00B9E3",
                        "#D39200","#F8766D"))+
        ggtitle("Major cell types")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
# Compare clusters for each dataset
cell.all <- FetchData(Thymus,"conditions")
conditions <- c("HF-L","HF-R","NC-B","NC-R")
cell.subsets <- lapply(conditions, function(x) 
        rownames(cell.all)[cell.all$conditions == x])

Thymus.subsets <- list()
for(i in 1:length(conditions)){
        Thymus.subsets[[i]] <- SubsetData(Thymus, cells.use =cell.subsets[[i]])
}

table(Thymus.subsets[[1]]@ident)

p <- list()
for(i in 1:length(conditions)){
        p[[i]] <- TSNEPlot(object = Thymus.subsets[[i]],do.label = TRUE, group.by = "ident", 
                           do.return = TRUE, no.legend = TRUE,
                           pt.size = 1,label.size = 4)+
                ggtitle(samples[i])+
                theme(text = element_text(size=20),     #larger text including legend title							
                      plot.title = element_text(hjust = 0.5)) #title in middle
}
do.call(plot_grid, p)
