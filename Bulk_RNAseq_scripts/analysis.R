# Load the required packages
library(DESeq2)
library(tximport)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Get the sample names from the file names
sample_names <- list.files(path = "../final/", pattern = "*.clean", full.names = F)
sample_names <- gsub('.genes.results.clean', '', sample_names)

# Get the actual files names
file_names <- list.files(path = "../final", pattern = "*.clean")

# Make a list of the file names
files <- file.path("../final", file_names)

# Add sample names to the files
names(files) <- sample_names

# Use tximport to import the files and create a matrix
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

# Change genes with length zero to 0.01
txi$length[txi$length == 0] <- 0.01

# check matrix
head(txi$counts)


# Create a dataframe for sample names and treatments
treatment <- c("control", "control", "control", "control", "disease", "disease", "disease", "disease", "disease")
colData <-data.frame(treatment, row.names = sample_names)

# Check our dataframe
head(colData)

# Check that the matrix and dataframe have the same information
all(rownames(colData) %in% colnames(txi$counts))
all(rownames(colData) == colnames(txi$counts))
which(rownames(colData) != colnames(txi$counts))

# Create a DESeq dataset
dds <- DESeqDataSetFromTximport(txi, colData, ~treatment)

# Check how many rows (genes) should be 78,932
nrow(dds) 

# Filter out low counts (keep only rows that have a count of at least 10 for a minimal number of samples)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Check  number of genes remaining (should be 21,782)
nrow(dds) 

# Set control 
dds$treatment <- relevel(dds$treatment, ref = "control")

# Run DESeq
dds <- DESeq(dds)
res <- results(dds)

# Add gene symbols to our results
res$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# Get a summary of the results
summary(res)

# reorder the results by p-value (most statistically significant)
resOrdered <- res[order(res$pvalue),]

# Most significant gene is WT1-AS
head(resOrdered)

# Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes
resLFC <- lfcShrink(dds,coef = "treatment_disease_vs_control", type = "apeglm")

# Use the plotCounts function to plot the counts for WT1-AS for both treatments
plotCounts(dds, gene=which.min(res$padj), intgroup="treatment", main = "WT1-AS")

# Use the variance stabilizing transformation (vst) function (good for visualization or clustering)
vsd <- vst(dds, blind=FALSE)

# Install and load the heatmap package "pheatmap"
#BiocManager::install("pheatmap")
library(pheatmap)

# Create a heatmap of the top 20 highest normalized counts
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,"treatment"])
rownames(df) <- colnames(dds)
colnames(df) = "Treatment"
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, annotation_names_col=FALSE)

# Create a heatmap of the distance matrix to see similarities and dissimilarities between samples
library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#Create a principal component analysis (PCA) plot
# useful for visualizing the overall effect of experimental covariates and batch effects
plotPCA(vsd, intgroup="treatment")

#################################################################################################
########################### Recreate figures from paper #########################################
#################################################################################################

#### Heatmap for 9hr treatments # adjusted p-value < 0.05 and |log2fc|> 1 ####
# Filter the results for adjusted p-value less than 0.05 and absolute value of log2 fold change greater than 1
resLFC.05fc1 <- as.data.frame(subset(resLFC, padj <0.05 & abs(log2FoldChange)>1))
# Check # of genes left (should be 1015)
nrow(resLFC.05fc1)
# Convert the VST transformed data into a data frame
vsd_df <- as.data.frame(assay(vsd))
# Filter the matrix based on above criteria for adjusted p-value and log2fc
vsd_df <- vsd_df[row.names(vsd_df) %in% row.names(resLFC.05fc1),]
# Check if filtered genes is correct (should be )
nrow(vsd_df)
# Add symbols to the matrix
vsd_df$symbol <- mapIds(org.Hs.eg.db, keys = row.names(vsd_df),
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
# Look at symbols - some are NA - these are long non-coding RNAs (lncRNAs) that don't have official gene symbols
vsd_df$symbol
# create a vector of gene symbols to use for the heatmap row labels
row_labels <- vsd_df$symbol
# remove the 'symbol' column from the data frame (need to do to make a matrix of just numbers)
vsd_df2 <- subset(vsd_df, select=-(symbol))
# Convert the data frame to a matrix
vsd_matrix <- as.matrix(vsd_df2)
# Scale (Z-score) the gene counts for better visualization
vsd_scaled <- t(scale(t(vsd_matrix)))
# Change column names to be more informative (*need to know order*)
colnames(vsd_scaled) <- c("Control_1", "Control_2", "Control_3", "Control_4",
                          "Disease_1", "Disease_2", "Disease_3", "Disease_4", "Disease_5")
# Create a color palette
myCol <- colorRampPalette(c("dodgerblue3","white","red2"))(n = 100)
# Load packages
library(ComplexHeatmap)
library(circlize)
library(colorRamp2)
# Create heatmap
Heatmap(vsd_scaled, col = myCol, column_names_rot=45,show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 8), row_names_max_width = unit(3, "cm"),
        heatmap_legend_param = list(title = "Z-score"))


#### Volcano Plot ####
# Add symbols to the adjusted log2fc data
resLFC$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
# Enhanced Volcano plot # default settings 
library(EnhancedVolcano)
EnhancedVolcano(resLFC,
                lab = resLFC$symbol,
                x = 'log2FoldChange',
                y = 'padj')


# Enhanced Volcano plot similar to paper's figure
EnhancedVolcano(resLFC,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                gridlines.major = FALSE, gridlines.minor = FALSE,
                ylab = bquote(~-Log[10] ~ "(FDR)"),
                cutoffLineType = "blank",
                subtitle=NULL,
                colAlpha = 1,
                caption = NULL,
                col = c("black", "black", "red", "green"),
                pCutoff = 5e-02,
                xlim=c(-4,4),
                legendPosition = 'none') + theme(plot.title=element_text(hjust=0.5))


#### Exporting a table ####
# Sort by adjusted p-vlue
resLFC.05fc1 <- arrange(resLFC.05fc1, padj)
# Add gene symbols
resLFC.05fc1$symbol <- mapIds(org.Hs.eg.db, keys = row.names(resLFC.05fc1),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
# Save as CSV
write.csv(resLFC.05fc1, file = "DEG_results.csv")


#### Protein-protein interaction (PPI) Networks using STRING #### Fig. 2 ####
# load the STRING database package
library(STRINGdb)
# Fetch the Homo sapiens database (#9606)
string_db <- STRINGdb$new(species=9606)

# map the gene names to the STRING database identifiers using the "map" method
vsd_df <- as.data.frame(assay(vsd))

# Filter the matrix based on above criteria for adjusted p-value and log2fc
vsd_df <- vsd_df[row.names(vsd_df) %in% row.names(resLFC.05fc1),]

# add gene symbols
vsd_df$symbol <- mapIds(org.Hs.eg.db, keys = row.names(vsd_df),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

# map to STRING database (takes about 1 minute)
mapped <- string_db$map(vsd_df, "symbol", removeUnmappedRows = TRUE)

# remove NAs
mapped <- mapped %>% na.exclude()

# Extract mapped genes
hits <- mapped$STRING_id

# Plot network
string_db$plot_network(hits)

# Cluster interactions and plot first 4 clusters (takes about 1 minute)
clustersList<- string_db$get_clusters(mapped$STRING_id)                
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]])
}

# Get a list of proteins known to interact with a protein (or proteins) of interest
VASH2 = string_db$mp("VASH2")
string_db$get_neighbors(VASH2)


#### Functional Enrichment with clusterProfiler ####
# Load package
library(clusterProfiler)
# Create gene list from previous data frame (Ensembl IDs)
genes <- row.names(resLFC.05fc1)

## Gene Ontology (GO) over-representation analysis (ORA)
## Biological Process (BP) ## (also CC and MF available)
go_bp_over <- enrichGO(gene = genes,
                    keyType = 'ENSEMBL',
                    OrgDb = org.Hs.eg.db,
                    universe = row.names(res),
                    ont = "BP",
                    readable = TRUE)
# Check results
head(go_bp_over)

# Visualize enriched GO terms with a Dot Plot
dotplot(go_bp_over) + ggtitle("GO Biological Process")


#### Functional Enrichment with EnrichR #### https://maayanlab.cloud/Enrichr/#libraries <-- see available libraries here
# Load package
library(enrichR)
# Use genes with adjusted p-value < 0.05
resLFC.05 <- as.data.frame(subset(resLFC, padj <0.05 ))
# Fetch the RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO database
db_disease <- c("RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO") # <-- choose your preferred library here
# Run analysis
enriched <- enrichr(resLFC.05$symbol, db_disease, background = res$symbol)
# Plot results
plotEnrich(enriched[["RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO"]], showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "Adjusted.P.value", title = "RNA-Seq Disease Gene and Drug Signatures from GEO")

#### fgsea package ####
library(fgsea)
library(msigdbr)
# Create fgsea input
res_filtered <- as.data.frame(res) %>% dplyr::select(symbol,stat) %>% na.omit() %>% distinct %>% group_by(symbol) %>% summarise(stat=mean(stat)) %>% filter(!duplicated(stat))
fgsea_input <- res_filtered$stat
names(fgsea_input) <- res_filtered$symbol
# Run fgsea
hallmarks = msigdbr(species = "human", collection = "H")
hallmarks_list = split(hallmarks$gene_symbol, f = hallmarks$gs_name)
fgsea_res <- fgsea(pathways=hallmarks_list, fgsea_input, minSize = 5, maxSize = 500)
fgsea_res_ordered <- fgsea_res[order(padj), ]
# The below line will save the results as a .csv file (uncomment to run)
#fwrite(fgsea_res_ordered, file = "fgsea_table_GSE231409.csv", sep =",", sep2 = c("", " ",""))
topPathwaysUp <- fgsea_res[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgsea_res[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
# The below line will save the plot as a .tiff file (uncomment to run)
#tiff("fgsea_hallmarks.tiff", width = 800, height = 800)
plotGseaTable(hallmarks_list[topPathways], fgsea_input, fgsea_res)


#### GSEA (BROAD) ####
# Process results to use with GSEA
gsea_input <- as.data.frame(res) %>%
  rownames_to_column(var = "gene_id") %>%
  dplyr::select(gene_id,stat) %>%
  na.omit() %>%
  distinct %>%
  group_by(gene_id) %>%
  summarise(stat=mean(stat)) %>%
  filter(!duplicated(stat))
#save file
write.table(gsea_input, file = "gsea_input.rnk", quote =F, sep = "\t", row.names=F, col.names = F)

#### Request a Linux Desktop session to run Broad GSEA GUI Application ####
#### see PowerPoint for details ####
