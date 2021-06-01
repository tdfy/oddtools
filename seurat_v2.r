#1. Make directory for all samples filtered feature barcodes
#2. Copy recursive from orignal directories and rename respective directory after sample
#3. 

library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(future)
library(stringr)

load("~/bin/extras/cycle.rda")

plan("multicore", workers = 10)
options(future.globals.maxSize = 4000 * 1024^2)

setwd("/home/yoderto/workspace/SC_v2/")

data <- "/home/yoderto/workspace/SC_v2/"

# for (file in c("IP", "D10","M4")){
#   seurat_data <- Read10X(data.dir = paste0(data, file))
#   seurat_obj <- CreateSeuratObject(counts = seurat_data,
#                                    min.features = 100,
#                                    project = file)
#   assign(file, seurat_obj)
# }
# 
# 
# P39.big <- merge(IP, y = c(D10, M4), add.cell.ids = c("IP", "D10", "M4"), project = "seurat_obj")
# 
# 
# P39.big$log10GenesPerUMI <- log10(P39.big$nFeature_RNA) / log10(P39.big$nCount_RNA)
# 
# P39.big$mitoRatio <- PercentageFeatureSet(object = P39.big,pattern = "^MT-")
# P39.big$mitoRatio <- P39.big@meta.data$mitoRatio / 100
# 
# 
# ###-------QC-----------------##
# 
# metadata <- P39.big@meta.data
# 
# metadata$cells <- rownames(metadata)
# 
# metadata$sample <- NA
# metadata$sample[which(str_detect(metadata$cells, "^IP_"))] <- "IP"
# metadata$sample[which(str_detect(metadata$cells, "^D10_"))] <- "D10"
# metadata$sample[which(str_detect(metadata$cells, "^M4_"))] <- "M4"
# 
# metadata <- metadata %>%
#   dplyr::rename(seq_folder = orig.ident,
#                 nUMI = nCount_RNA,
#                 nGene = nFeature_RNA)
# 
# # Add metadata back to Seurat object
# P39.big@meta.data <- metadata
# 
# # Create .RData object to load at any time
# save(merged_seurat, file="data/merged_filtered_seurat.RData")
# 
# 
# metadata %>%
#   ggplot(aes(x=sample, fill=sample)) +
#   geom_bar() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   theme(plot.title = element_text(hjust=0.5, face="bold")) +
#   ggtitle("NCells")
# 
# hope2 = paste("cellnum", ".png", sep="")
# ggsave(hope2,units = c("in"), height = 10, width = 10)
# 
# 
# 
# metadata %>%
#   ggplot(aes(color=sample, x=nUMI, fill= sample)) +
#   geom_density(alpha = 0.2) +
#   scale_x_log10() +
#   theme_classic() +
#   ylab("Cell density") +
#   geom_vline(xintercept = 500)
# hope2 = paste("umipercell", ".png", sep="")
# ggsave(hope2,units = c("in"), height = 10, width = 10)
# 
# 
# # Visualize the distribution of genes detected per cell via histogram
# metadata %>%
#   ggplot(aes(color=sample, x=nGene, fill= sample)) +
#   geom_density(alpha = 0.2) +
#   theme_classic() +
#   scale_x_log10() +
#   geom_vline(xintercept = 300)
# 
# hope2 = paste("genesdis", ".png", sep="")
# ggsave(hope2,units = c("in"), height = 10, width = 10)
# 
# # Visualize the distribution of genes detected per cell via boxplot
# metadata %>%
#   ggplot(aes(x=sample, y=log10(nGene), fill=sample)) +
#   geom_boxplot() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   theme(plot.title = element_text(hjust=0.5, face="bold")) +
#   ggtitle("NCells vs NGenes")
# 
# 
# hope2 = paste("genespercelldis", ".png", sep="")
# ggsave(hope2,units = c("in"), height = 10, width = 10)
# 
# # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
# metadata %>%
#   ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) +
#   geom_point() +
#   scale_colour_gradient(low = "gray90", high = "black") +
#   stat_smooth(method=lm) +
#   scale_x_log10() +
#   scale_y_log10() +
#   theme_classic() +
#   geom_vline(xintercept = 500) +
#   geom_hline(yintercept = 250) +
#   facet_wrap(~sample)
# 
# hope2 = paste("geneumicorr2", ".png", sep="")
# ggsave(hope2,units = c("in"), height = 10, width = 10)
# 
# 
# # Visualize the distribution of mitochondrial gene expression detected per cell
# metadata %>%
#   ggplot(aes(color=sample, x=mitoRatio, fill=sample)) +
#   geom_density(alpha = 0.2) +
#   scale_x_log10() +
#   theme_classic() +
#   geom_vline(xintercept = 0.2)
# 
# hope2 = paste("mito", ".png", sep="")
# ggsave(hope2,units = c("in"), height = 10, width = 10)
# 
# # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
# metadata %>%
#   ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
#   geom_density(alpha = 0.2) +
#   theme_classic() +
#   geom_vline(xintercept = 0.8)
# 
# hope2 = paste("comp", ".png", sep="")
# ggsave(hope2,units = c("in"), height = 10, width = 10)
# 
# 
# 
# 
# 
# # Filter out low quality cells using selected thresholds - these will change with experiment
# filtered_seurat <- subset(x = P39.big,
#                           subset= (nUMI >= 500) &
#                             (nGene >= 250) &
#                             (log10GenesPerUMI > 0.80) &
#                             (mitoRatio < 0.20))
# 
# # Extract counts
# counts <- GetAssayData(object = filtered_seurat, slot = "counts")
# 
# # Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
# nonzero <- counts > 0
# 
# # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
# keep_genes <- Matrix::rowSums(nonzero) >= 10
# 
# # Only keeping those genes expressed in more than 10 cells
# filtered_counts <- counts[keep_genes, ]
# 
# # Reassign to filtered Seurat object
# filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
# 
# 
# # Save filtered subset to new metadata
# metadata_clean <- filtered_seurat@meta.data



load("seurat_filtered.RData")
     
# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)

seurat_phase <- CellCycleScoring(seurat_phase, 
                                g2m.features = g2m_genes, 
                                s.features = s_genes)

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                    selection.method = "vst",
                                    nfeatures = 2000, 
                                    verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
       reduction = "pca",
       group.by= "Phase",
       split.by = "Phase")

# Regress cycling cells by non-cycling (preferred for differentiating cells)
# https://satijalab.org/seurat/archive/v3.0/cell_cycle_vignette.html




seurat_phase@meta.data$celldif <- seurat_phase@meta.data$S.Score - seurat_phase@meta.data$G2M.Score


# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                    breaks=c(-Inf, 0.03749, 0.04833, 0.06425, Inf), 
                                    labels=c("Low","Medium","Medium high", "High"))


# Plot the PCA colored by mitoFr
DimPlot(seurat_phase,
       reduction = "pca",
       group.by= "mitoFr",
       split.by = "mitoFr")
hope2 = paste("mito", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

split_seurat <- split_seurat[c("IP", "D10","M4")]
 
 
 
for (i in 1:length(split_seurat)) {
   split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("celldif","mitoRatio"))
 }
 
 
# saveRDS(split_seurat, "split_seurat.rds")
#  
# IP <- DimPlot(split_seurat$IP)
# D10 <- DimPlot(split_seurat$D10)
# M4 <- DimPlot(split_seurat$M4)
# IP + D10 + M4
# hope2 = paste("dim_not_int", ".png", sep="")
# ggsave(hope2,units = c("in"), height = 10, width = 10)

readRDS("split_seurat.rds")

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")


# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
PCAPlot(seurat_integrated,
        split.by = "sample")  

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated, split.by = "orig.ident")
hope2 = paste("dim_int", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)




DimHeatmap(int_seurat, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

print(x = int_seurat[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)


# Plot the elbow plot
ElbowPlot(object = int_seurat, 
          ndims = 40)


# Determine the K-nearest neighbor graph
int_seurat <- FindNeighbors(object = int_seurat, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
int_seurat <- FindClusters(object = int_seurat,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))


# Assign identity of clusters
Idents(object = int_seurat) <- "integrated_snn_res.0.8"



# Plot the UMAP
DimPlot(int_seurat,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
hope2 = paste("res_int", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)

# Assign identity of clusters
Idents(object = int_seurat) <- "integrated_snn_res.0.6"

# Plot the UMAP
DimPlot(int_seurat,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
hope2 = paste("res_int6", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)






##-----------------------------QC------------------------##

DimPlot(int_seurat,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()
hope2 = paste("phase_intck", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)



# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(int_seurat, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
hope2 = paste("compchk_int", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)

int_seuat <- readRDS("post_in.rds")

FeaturePlot(int_seuat, 
            reduction = "umap", 
            features = c("CD14", "LYZ","CD19","CD20","BLK","TNFRSF17"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
hope2 = paste("bcells", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)

FeaturePlot(int_seuat, 
            reduction = "umap", 
            features = c("CD8A","CD8B","CD244","EOMES","CD223"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
hope2 = paste("CD8cells", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)


FeaturePlot(int_seuat, 
            reduction = "umap", 
            features = c("CD4","FOXP3","CXCR3","CCR5","CXCR4","CCR4","CXCR5"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
hope2 = paste("CD4cells", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)


FeaturePlot(int_seuat, 
            reduction = "umap", 
            features = c("NYCE"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
hope2 = paste("NYCEcells", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)





###-----------------Conserved Markers----------------------##

DefaultAssay(int_seuat) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(int_seuat,
                                                   ident.1 = 4,
                                                   grouping.var = "sample",
                                                   only.pos = FALSE,
                                                   logfc.threshold = 0.25)

GOI <- rownames(cluster0_conserved_markers)


test.subset <- subset(x = int_seuat, subset = integrated_snn_res.0.6 == 4)

ScaleData(test.subset)

DefaultAssay(int_seuat) <- "integrated"

# Single cell heatmap of feature expression
DoHeatmap(subset(int_seuat, subset = integrated_snn_res.0.6 == 4), features = GOI, size = 3,group.by = "sample")
hope2 = paste("clus4cells", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)
