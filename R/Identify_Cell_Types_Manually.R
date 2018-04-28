library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/Thymus_alignment.Rda")
lnames
Featureplot <- function(x){
        p <- FeaturePlot(object = Thymus, 
                         reduction.use = "tsne",
                         features.plot = x, min.cutoff = NA, 
                         cols.use = c("lightgrey","blue"), pt.size = 0.5)
        return(p)
}

#------
Adipocytes <- MouseGenes(Thymus,c("SLC36A2","P2RX5","TRIP4","ASCC1","MYF5","UCP1"))
Endothelium <- MouseGenes(Thymus,c("Cdh5","Pecam1","Flt1","Vwf","Plvap","Kdr",
                          "EMCN","Car4","ptprb"))
Mesothelial <- MouseGenes(Thymus,c("Msln","Upk3b","Lrrn4"))
cTECs <- MouseGenes(Thymus,c("Ccl25","Cxcl12","Ctsl","Psmb11"))
mTEChi <- MouseGenes(Thymus,c("Aire","HLA-DMB"))
mTEClo <- MouseGenes(Thymus,c("Krt5","Gas1"))
TEC_subset <- MouseGenes(Thymus,c("Plet1","Ly6d","Spink5","Reg3g","Bpifa1"))
TECs <- c(cTECs,mTEChi,mTEClo,TEC_subset)
Epithelium <- MouseGenes(Thymus,c("Epcam","KRT19"))
#--Hematopoietic----
Hematopoietic <- MouseGenes(Thymus,c("PTPRC","LAPTM5","SRGN"))
#------Myeloid----
megakaryocytes <-  MouseGenes(Thymus,c("PPBP","GNG11"))
erythrocyte <-  MouseGenes(Thymus,c("HBA2","HBB"))
MastCells <- MouseGenes(Thymus,c("Cma1","Mcpt4","Tpsb2","Cpa3"))
Monocytes <-  MouseGenes(Thymus,c("LYZ","S100A9","CD14","CCL2","FCGR3A","MS4A7","VMO1"))
Macrophages <- MouseGenes(Thymus,c("LYZ","CD68","MARCO","Emr1"))
DendriticCells <- MouseGenes(Thymus,c("Itgax","GPR183","CST3","HLA-DQA1","FCER1A","TSPAN13",
                                     "IL3RA","IGJ"))
Myeloid <-  c(megakaryocytes,erythrocyte,MastCells,Monocytes,Macrophages,DendriticCells)
#------Lymphoid----
T_Cell <- MouseGenes(Thymus,c("CD3G","CD3D","CD2"))
CD4_Naive_T <- MouseGenes(Thymus,c("CD4","IL7R","GIMAP5","SELL","IL2RG"))
NK <- MouseGenes(Thymus,c("GNLY","Ncr1","CCL5","KLRD1","NKG7"))
B_Cell <- MouseGenes(Thymus,c("CD19","MS4A1","CD38"))
B_Stem_Cell <- MouseGenes(Thymus,c("SPN","CD20"))
Lymphoid <- c(T_Cell,CD4_Naive_T,NK,B_Cell)
#Mesenchymal Cells
Mesenchymal <- MouseGenes(Thymus,c("Pdgfrb","Vim","Has2","Dcn"))
Neurons <- MouseGenes(Thymus,c("Rbfox3"))
Pericytes <- MouseGenes(Thymus,c("Pdgfrb","Cspg4","Anpep","Rgs5",
                                 "Myh11","Mylk","Sost","Des","Vtn","Ifitm1"))
Smooth_muscle_cells <- MouseGenes(Thymus,c("Acta2","Myh11"))
# Stem cell
Embryonic_Stem_Cell <- MouseGenes(Thymus,c("POU5F1","FUT4"))
Hematopoietic_Stem_Cell <- MouseGenes(Thymus,c("CD34","PROM1","ABCG2","ATXN1","RUNX1"))
Mesenchymal_Stem_Cells <- MouseGenes(Thymus,c("Eng","CD44","Nt5e","Alcam","Kit",
                        "Ly6a","Thy1","Itgb1","Vcam1","Icam2","CD72","CD2"))
Neural_Stem_Cell <- MouseGenes(Thymus,c("NES","NCAM","NGFR"))
# Featureplot
Featureplot(Adipocytes) # Adipocytes
Featureplot(Endothelium) # Endothelial Cells
Featureplot(Epithelium) # Epithelium
Featureplot(TECs) # TECs
Featureplot(Fibroblast) # Fibroblasts
Featureplot(Hematopoietic) # Hematopoietic cells
Featureplot(Myeloid) # Myeloid cells
Featureplot(T_Cell) # RPE, Melanocytes, Myelinating Schwann cells
eatureplot(B_Cell) # RPE, Melanocytes, Myelinating Schwann cells
Featureplot(MastCells)
Featureplot(Monocytes)
Featureplot(Macrophages)
Featureplot(DendriticCells)
Featureplot(Lymphoid) # Lymphoid cells
Featureplot(Mesenchymal)
Featureplot(Pericytes)
Featureplot(Smooth_muscle_cells)

markers.to.plot <- c(Mesenchymal[c(2,4)],Smooth_muscle_cells,Mesothelial,
                     Epithelium[1],TECs[-3],
                     Hematopoietic,T_Cell,Endothelium[c(1:3,9)])
markers.to.plot <- MouseGenes(Thymus,markers.to.plot, unique =T)
DotPlot(Thymus, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)
# Rename ident
table(Thymus@ident)
idents <- as.data.frame(table(Thymus@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("mTEChi 0",
                     "Mesenchymal Cells 1",
                     "T cells 2",
                     "T cells 3",
                     "Endothelial cells 4",
                     "Endothelial cells 5",
                     "T cells 6",
                     "cTECs & mTEClo 7",
                     "mTEClo 8",
                     "T cells & cTECs 9",
                     "Mesenchymal Cells 10",
                     "T cells 11",
                     "Mesothelial cells 12",
                     "Smooth muscle cells 13",
                     "mTEClo & mTEChi 14",
                     "Smooth muscle cells 15",
                     "Epithelium 16")

Thymus@ident <- plyr::mapvalues(x = Thymus@ident,
                            from = old.ident.ids,
                            to = new.cluster.ids)
DotPlot(Thymus, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)
new.cluster.ids <- c("mTEChi",
                     "Mesenchymal\n
                     cells",
                     "T cells",
                     "T cells",
                     "Endothelial\n
                     cells",
                     "Endothelial\n
                     cells",
                     "T cells",
                     "cTECs &\n
                     mTEClo",
                     "mTEClo",
                     "T cells &\n
                     cTECs",
                     "Mesenchymal\n
                     Cells",
                     "T cells",
                     "Mesothelial\n
                     cells",
                     "Smooth muscle\n
                     cells",
                     "mTEClo &\n
                     mTEChi",
                     "Smooth muscle\n
                     cells",
                     "Epithelium")
Thymus@ident <- plyr::mapvalues(x = Thymus@ident,
                                from = old.ident.ids,
                                to = new.cluster.ids)
table(Thymus@ident)
DotPlot(Thymus, genes.plo = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)

TSNEPlot(object = Thymus, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
        ggtitle("Major cell types")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle

SplitTSNEPlot(Thymus, "conditions",do.label = F,no.legend = T)
SplitTSNEPlot(Thymus, "conditions",do.label = F,no.legend = F)
#2.1.1 - A table with the number of cells of each cluster and subcluster, for both B6 and 129_B6 strains.
# We can also compare proportional shifts in the data. As can be seen in the barplot, 
freq_table <- prop.table(x = table(Thymus@ident, Thymus@meta.data[, "conditions"]), 
                         margin = 2)
barplot(height = freq_table)
freq_table
table(Thymus@meta.data[, "conditions"])

#====== 2.2 T Cell depletion  ==========================================
lnames = load(file = "./data/Thymus_alignment.Rda")
lnames
Thymus_noT <- SubsetData(Thymus,ident.remove = 
                                 (old.ident.ids[(new.cluster.ids %in% "T cells")]))
save(Thymus_noT, file = "./output/Thymus_noT.Rda")
# How many cells are in each cluster
table(Thymus_noT@ident)
Thymus_noT@ident <- plyr::mapvalues(x = Thymus_noT@ident,
                                from = old.ident.ids[!(new.cluster.ids %in% "T cells")],
                                to = new.cluster.ids[!(new.cluster.ids %in% "T cells")])
#Now we can run a single integrated analysis on all cells!
TSNEPlot(object = Thymus_noT, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
  ggtitle("Major cell types after T cell depletion")+
  theme(text = element_text(size=20),     #larger text including legend title							
        plot.title = element_text(hjust = 0.5)) #title in middle
# SplitTSNEPlot
SplitTSNEPlot(Thymus_noT, split.by = "conditions",do.label = F,no.legend = T)
SplitTSNEPlot(Thymus_noT, split.by = "conditions",do.label = F,no.legend = F)
# Population of different cell types after T cell depletion
freq_table <- prop.table(x = table(Thymus_noT@ident, Thymus_noT@meta.data[, "conditions"]), 
                         margin = 2)
barplot(height = freq_table)
freq_table

table(Thymus_noT@meta.data[, "conditions"])

p1 <- VlnPlot(object = Thymus_noT, features.plot = "Cd38", group.by = "conditions")
p2 <- VlnPlot(object = Thymus_noT, features.plot = "Cd38", group.by = "ident")
plot_grid(p1, p2, nrow=2)
