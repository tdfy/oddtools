#1. Make directory for all samples filtered feature barcodes
#2. Copy recursive from orignal directories and rename respective directory after sample
#3. workflow: https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html

library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(future)
library(stringr)

load("~/bin/extras/cycle.rda")

plan("multicore", workers = 10)
options(future.globals.maxSize = 2000 * 1024^2)

setwd("/home/yoderto/workspace/SC_v2/")

data <- "/home/yoderto/workspace/SC_v2/"

for (file in c("IP", "D10","M4")){
        seurat_data <- Read10X(data.dir = paste0(data, file))
        seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                         min.features = 100, 
                                         project = file)
        assign(file, seurat_obj)
}


P39.big <- merge(IP, y = c(D10, M4), add.cell.ids = c("IP", "D10", "M4"), project = "seurat_obj")


P39.big$log10GenesPerUMI <- log10P39.big$nFeature_RNA) / log10P39.big$nCount_RNA)

P39.big$mitoRatio <- PercentageFeatureSet(object = P39.big,pattern = "^MT-")
P39.big$mitoRatio <- P39.big@meta.data$mitoRatio / 100


###-------QC-----------------##

metadata <- P39.big@meta.data

metadata$cells <- rownames(metadata)

metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^IP_"))] <- "IP"
metadata$sample[which(str_detect(metadata$cells, "^D10_"))] <- "D10"
metadata$sample[which(str_detect(metadata$cells, "^M4_"))] <- "M4"

metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)

# Add metadata back to Seurat object
P39.big@meta.data <- metadata
                           
# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")


metadata %>% 
  	ggplot(aes(x=sample, fill=sample)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")

hope2 = paste("cellnum", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)



metadata %>% 
  	ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500)
hope2 = paste("umipercell", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)


# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)

hope2 = paste("genesdis", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  	ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  	geom_boxplot() + 
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells vs NGenes")


hope2 = paste("genespercelldis", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~sample)

hope2 = paste("geneumicorr", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)


# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  	ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)

hope2 = paste("mito", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)

hope2 = paste("comp", ".png", sep="")
ggsave(hope2,units = c("in"), height = 10, width = 10)

