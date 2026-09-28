library(dplyr)
library(Seurat)
library(patchwork)

#########################################################################################
############################## Control Samples ##########################################
#########################################################################################

control_matrix <- ReadMtx(
  mtx = "../control_aggr/outs/count/filtered_feature_bc_matrix/matrix.mtx.gz",
  features = "../control_aggr/outs/count/filtered_feature_bc_matrix/features.tsv.gz",
  cells = "../control_aggr/outs/count/filtered_feature_bc_matrix/barcodes.tsv.gz"
)
control <- CreateSeuratObject(counts = control_matrix)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(control, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(control, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
par(mar = c(10,10,10,10))
# Filter
control <- subset(control, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
control <- NormalizeData(control, normalization.method = "LogNormalize", scale.factor = 10000)

control <- NormalizeData(control)

control <- FindVariableFeatures(control, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(control), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(control)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling data
all.genes <- rownames(control)
control <- ScaleData(control, features = all.genes)

# Perform linear dimensional reduction
control <- RunPCA(control, features = VariableFeatures(object = control))

# Examine and visualize PCA results a few different ways
print(control[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(control, dims = 1:2, reduction = "pca")

DimPlot(control, reduction = "pca") + NoLegend()

DimHeatmap(control, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(control, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset (use to choose # of PCs <-- 17?)
ElbowPlot(control)

# Cluster the cells
control <- FindNeighbors(control, dims = 1:17)

control <- FindClusters(control, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(control), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
control <- RunUMAP(control, dims = 1:17)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(control, reduction = "umap")


#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 0
cluster0.markers <- FindMarkers(control, ident.1 = 0)
head(cluster0.markers, n=10)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(control, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
control.markers <- FindAllMarkers(control, only.pos = TRUE)
control.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# ROC test
cluster0.markers <- FindMarkers(control, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster0.markers)

#expression probability distributions across clusters (*choose your genes of interest*)
VlnPlot(control, features = c("MYRF", "SLC44A1"))

# you can plot raw counts as well (*choose your genes of interest*)
VlnPlot(control, features = c("NKG7", "FMN1"), slot = "counts", log = TRUE)

# Multiple UMAPs (*choose your genes of interest*)
FeaturePlot(control, features = c("MYRF", "GNLY", "CD3E", "CD14", "NKG7", "FCGR3A", "LYZ", "FMN1", "CD8A"))

# expression heatmap for given cells and features
control.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(control, features = top10$gene) + NoLegend()

# Automatic cell type annotation
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/kris-nader/sc-type/master/R/sctype_wrapper.R")
control <- run_sctype(control,known_tissue_type=c("Brain"),
                  custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",name="sctype_classification")
DimPlot(control, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')

# You can save the object at this point so that it can easily be loaded back in without having to rerun the computationally intensive steps performed above, or easily shared with collaborators.
saveRDS(control, file = "control_seurat.rds")

##############
## CellChat ##
##############

library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

# Use previously saved Seurat object
seurat_object <- readRDS(file = "control_seurat.rds")
Idents(seurat_object)
Idents(object = seurat_object) <- "sctype_classification"
data.input <- control[["RNA"]]$data # normalized data matrix
Idents(object = seurat_object) <- "sctype_classification"
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(seurat_object)
meta <- data.frame(labels = labels, row.names = names(labels))

cellchat <- createCellChat(object = seurat_object, group.by = "ident", assay = "RNA")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
 CellChatDB.use <- CellChatDB # simply use the default CellChatDB.
# We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat@DB <- CellChatDB.use

showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
interactions <- CellChatDB$interaction
#write.csv(interactions, file="interactions.csv")

# subset the expression data of signaling genes for saving computation cost
# This step is necessary even if using the whole database
cellchat <- subsetData(cellchat)

future::plan("multisession", workers = 4) # do parallel processing
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#cellchat <- projectData(cellchat, PPI.human)

# Compute the communication probability and infer cellular communication network (*takes a few minutes)
cellchat <- computeCommunProb(cellchat, type = "triMean")

#filter out the cell-cell communication if there are only few cells in certain cell groups (default is 10)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
# returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors.
df.net <- subsetCommunication(cellchat)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot
# Thicker edge line indicates a stronger signal
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# examine the signaling sent from each cell group
mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# #1="Unknown" #2="Radial glial cells" #3="Oligodendrocytes" #4="Microglial cells" #5="Immune system cells"

# All the signaling pathways showing significant communications can be accessed by cellchat@netP$pathways
cellchat@netP$pathways
#  "PTN"           "EGF"           "MIF"           "SPP1"          "CD99"          
#  "CADM"          "APP"           "CypA"          "NCAM"           "MAG"
#  "CLDN"          "MK"            "CLEC"          "ADGRE"         "LAMININ"
#  "ApoE"          "ADGRL"         "ANNEXIN"       "CDH"           "MHC-I"
#  "ICAM"          "JAM"           "PLAU"          "GRN"        "Prostaglandin"
#  "SIRP"          "CNTN"           "CD86"
# *choose which pathway for below; what do they do?

pathways.show <- c("EGF") # <-- choose pathway of interest

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to Radial glial cells and the right portion shows signaling to immune cells 
vertex.receiver = seq(2,5) # a numeric vector. 
#dev.off() # if below doesn't plot, uncomment here, run to reset plotting device, and run below line again
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
#dev.off() #if below doesn't plot, uncomment here, run to reset plotting device, and run below line again
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
#dev.off()
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Chord diagram
# we can define a named char vector group to create multiple-group chord diagram, e.g., grouping cell clusters into different cell types
group.cellType <- c(rep("Radial glial cells", 4), rep("Microglial cells", 4)) # grouping cell clusters into glial cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType,
                     title.name = paste0(pathways.show, " Signaling Network"),
                     )

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair.
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.PTN <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.PTN[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,5) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
# We can also show all the significant interactions (L-R pairs) from some cell groups to other cell groups using netVisual_bubble.
# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')

# #1="Unknown" #2="Radial glial cells" #3="Oligodendrocytes" #4="Microglial cells" #5="Immune system cells"

netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:4), remove.isolate = FALSE)

# (2) show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:3,5), signaling = c("PTN","SPP1","EGF"), remove.isolate = FALSE) #<-- choose pathways of interest

# (3) show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("PTN","SPP1","EGF")) #<-- choose pathways of interest
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(1,2,5), pairLR.use = pairLR.use, remove.isolate = TRUE)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Immune cells (#5)
netVisual_chord_gene(cellchat, sources.use = 5, targets.use = c(1:4), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions received by Immune system cells (#5)
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 5, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5), signaling = c("PTN","SPP1"),legend.pos.x = 8) #<-- choose pathways of interest

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5), slot.name = "netP")

# Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "PTN", enriched.only = TRUE, type = "violin")

#show the expression of all signaling genes related to one signaling pathway
plotGeneExpression(cellchat, signaling = "PTN", enriched.only = FALSE)


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
#dev.off() #if below doesn't plot, uncomment here, run to reset plotting device, and run below line again
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Visualize dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("PTN", "SPP1")) #<-- choose pathwasy of interest
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

#  Identify signals contributing the most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("PTN", "SPP1", "EGF")) # <-- genes of interest
ht

# Identify and visualize outgoing communication pattern of secreting cells
library(NMF)
library(ggalluvial)
# Here we run selectK to infer the number of patterns
selectK(cellchat, pattern = "outgoing") # takes some time...
# Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3
nPatterns = 3
dev.off()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# Please make sure you have load `library(ggalluvial)` when running this function

netAnalysis_dot(cellchat, pattern = "outgoing")

# Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming") # takes some time...
 # Cophenetic values begin to drop when the number of incoming patterns is 5
nPatterns = 5
dev.off()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

# Save the CellChat object
saveRDS(cellchat, file = "cellchat_control.rds")

#Explore the cell-cell communication through the Interactive CellChat Explorer
runCellChatApp(cellchat)

###########################################################################################
############################## Treatment Samples ##########################################
###########################################################################################

treatment_matrix <- ReadMtx(
  mtx = "../treatment_aggr/outs/count/filtered_feature_bc_matrix/matrix.mtx.gz",
  features = "../treatment_aggr/outs/count/filtered_feature_bc_matrix/features.tsv.gz",
  cells = "../treatment_aggr/outs/count/filtered_feature_bc_matrix/barcodes.tsv.gz"
)
treatment <- CreateSeuratObject(counts = treatment_matrix)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
treatment[["percent.mt"]] <- PercentageFeatureSet(treatment, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(treatment, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(treatment, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(treatment, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter
treatment <- subset(treatment, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
treatment <- NormalizeData(treatment, normalization.method = "LogNormalize", scale.factor = 10000)

treatment <- NormalizeData(treatment)

treatment <- FindVariableFeatures(treatment, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(treatment), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(treatment)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling data
all.genes <- rownames(treatment)
treatment <- ScaleData(treatment, features = all.genes)

# Perform linear dimensional reduction
treatment <- RunPCA(treatment, features = VariableFeatures(object = treatment))

# Examine and visualize PCA results a few different ways
print(treatment[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(treatment, dims = 1:2, reduction = "pca")

DimPlot(treatment, reduction = "pca") + NoLegend()

DimHeatmap(treatment, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(treatment, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset (use to choose # of PCs <-- 18?)
ElbowPlot(treatment)

# Cluster the cells
treatment <- FindNeighbors(treatment, dims = 1:18)

treatment <- FindClusters(treatment, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(treatment), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
treatment <- RunUMAP(treatment, dims = 1:18)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(treatment, reduction = "umap")

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 0
cluster0.markers <- FindMarkers(treatment, ident.1 = 0)
head(cluster0.markers, n = 10)

# find all markers distinguishing cluster 4 from clusters 0 and 3
cluster4.markers <- FindMarkers(treatment, ident.1 = 4, ident.2 = c(0, 3))
head(cluster4.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
treatment.markers <- FindAllMarkers(treatment, only.pos = TRUE)
treatment.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# ROC test
cluster0.markers <- FindMarkers(treatment, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster0.markers

# Expression probability distributions across clusters (*choose your genes of interest*)
VlnPlot(treatment, features = c("MYRF", "SLC44A1"))

# you can plot raw counts as well (*choose your genes of interest*)
VlnPlot(treatment, features = c("NKG7", "FMN1"), slot = "counts", log = TRUE)

# Multiple UMAPs (*choose your genes of interest*)
FeaturePlot(treatment, features = c("MYRF", "GNLY", "CD3E", "CD14", "NKG7", "FCGR3A", "LYZ", "FMN1", "CD8A"))

# expression heatmap for given cells and features
treatment.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(treatment, features = top10$gene) + NoLegend()

# Automatic cell type annotation
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
treatment <- run_sctype(treatment,known_tissue_type="Brain",
                  custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",name="sctype_classification")
DimPlot(treatment, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')

treatment <- run_sctype(treatment,known_tissue_type=c("Brain", "Immune system"),
                        custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",name="sctype_classification")
DimPlot(treatment, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')

# You can save the object at this point so that it can easily be loaded back in without having to rerun the computationally intensive steps performed above, or easily shared with collaborators.
saveRDS(treatment, file = "treatment_seurat.rds")

##############
## CellChat ##
##############

library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

seurat_object <- readRDS(file = "treatment_seurat.rds")
Idents(seurat_object)
Idents(object = seurat_object) <- "sctype_classification"

data.input <- control[["RNA"]]$data # normalized data matrix
Idents(object = seurat_object) <- "sctype_classification"
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(seurat_object)
meta <- data.frame(labels = labels, row.names = names(labels))

cellchat <- createCellChat(object = seurat_object, group.by = "ident", assay = "RNA")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
 CellChatDB.use <- CellChatDB # simply use the default CellChatDB.
# We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat@DB <- CellChatDB.use


showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
interactions <- CellChatDB$interaction
#write.csv(interactions, file="interactions_all.csv")

# subset the expression data of signaling genes for saving computation cost
# This step is necessary even if using the whole database
cellchat <- subsetData(cellchat)

future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#cellchat <- projectData(cellchat, PPI.human)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")

#filter out the cell-cell communication if there are only few cells in certain cell groups (default is 10)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
# returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors.
df.net <- subsetCommunication(cellchat)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

#  showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# examine the signaling sent from each cell group
mat <- cellchat@net$weight
par(mfrow = c(1,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# All the signaling pathways showing significant communications can be accessed by cellchat@netP$pathways
cellchat@netP$pathways
#  "PTN"  "EGF"  "SPP1" <-- choose which one for below; what do they do?
pathways.show <- c("PTN") # <-- choose pathway of interest

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
#dev.off()
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
#dev.off()
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Chord diagram
# we can define a named char vector group to create multiple-group chord diagram, e.g., grouping cell clusters into different cell types
group.cellType <- c(rep("Radial glial cells", 4), rep("Microglial cells", 4)) # grouping cell clusters into glial cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair.
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.PTN <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.PTN[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
# We can also show all the significant interactions (L-R pairs) from some cell groups to other cell groups using netVisual_bubble.
# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# #1="Glutamatergic neurons" #2="Microglial cells" #3="Oligodendrocytes" #4="Radial glial cells"
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:3), remove.isolate = FALSE)

# (2) show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:3), signaling = c("PTN","SPP1","EGF"), remove.isolate = FALSE) #<-- choose pathways of interest

# (3) show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("PTN","SPP1","EGF")) #<-- choose pathways of interest
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(1,2), pairLR.use = pairLR.use, remove.isolate = TRUE)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Radial glial cells (#4)
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(1:3), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions received by Radial glial cells (#4)
netVisual_chord_gene(cellchat, sources.use = c(1:3), targets.use = 4, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1:3), targets.use = c(4), signaling = c("PTN","SPP1"),legend.pos.x = 8) #<-- choose pathways of interest

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1:3), targets.use = c(4), slot.name = "netP")

# Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "PTN", enriched.only = TRUE, type = "violin")

#show the expression of all signaling genes related to one signaling pathway
plotGeneExpression(cellchat, signaling = "PTN", enriched.only = FALSE)


# Compute the network centrality scores
#dev.off()
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Visualize dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("PTN", "SPP1")) #<-- choose pathwasy of interest
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

#  Identify signals contributing the most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("PTN", "SPP1", "EGF")) # <-- genes of interest
ht

# Identify and visualize outgoing communication pattern of secreting cells
library(NMF)
library(ggalluvial)
# Here we run selectK to infer the number of patterns
selectK(cellchat, pattern = "outgoing") # takes some time...
# Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3
nPatterns = 3
dev.off()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

netAnalysis_dot(cellchat, pattern = "outgoing")

# Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming") # takes some time...
# Cophenetic values begin to drop when the number of incoming patterns is 3
nPatterns = 3
dev.off()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

# Save the CellChat object
saveRDS(cellchat, file = "cellchat_treatment.rds")

#Explore the cell-cell communication through the Interactive CellChat Explorer
runCellChatApp(cellchat)

