library(ggrepel)
library(ggplot2)
library(data.table)
library(supclust)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dplyr)


tdn <- read.table("C:/Export/Betsy/TDN/PEAK_AUC/G2xG3_peak.csv", header=TRUE, sep=",")

div <- read.table("C:/Export/Betsy/APH/diversity.csv", header=TRUE, sep=",")


tdn$gene_name <- row.names(tdn)

tdn <- subset(tdn, !(gene_name %in% div$Gene))


tdn$diffexpressed <- "NO"
tdn$diffexpressed[tdn$logFC > 0.5 & tdn$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
tdn$diffexpressed[tdn$logFC < -0.5 & tdn$pvalue < 0.05] <- "DOWN"

# tdn <- subset(tdn, Gene %in% rownames(query4)) #<------------------------ filter
ranger <- "NA"
tdn$ranger[tdn$logFC > -7 & tdn$logFC < 7] <- "select" #<------------------------ filter
tdn <- subset(tdn, ranger == "select")

tdn$delabel <- NA
tdn$delabel[tdn$diffexpressed != "NO"] <- tdn$Gene[tdn$diffexpressed != "NO"]


tdn$delabel[tdn$diffexpressed != "NO"] 

tdn$adj <- "q-value>=0.10"
tdn$adj[tdn$qvalue < 0.10] <- "q-value=<0.10"


ggplot(data=tdn, aes(x=logFC, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point(aes(shape = adj), size = 2) + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_text_repel(aes(label=ifelse(pvalue < 0.05 & abs(logFC) >= 0.5,
                                   tdn$gene_name, '')))+ggtitle("NanoStringDiff Differential Expression for UPCC35313  High and Non-Responders")




qual <- tdn[(tdn$adj == "q-value=<0.10") & (tdn$diffexpressed != "NO"),]

write.table(qual,"C:/Export/Betsy/G1xG3_candidate.tsv", sep="\t",row.names=TRUE,quote = FALSE,col.names=TRUE)


diff_gene <- qual$gene_name

q2 <- read.table("C:/Export/Betsy/TDN/PEAK_AUC/Peak_AUC_supclust_G1xG2.tsv", header=TRUE, sep="\t") #<---------- nSolver cols easy 

pheat_list <- c(unique(c(diff_gene,q2$Probe.Name)))


###----------------nSolver Normalized---------------------------##
###----------------______---------------------------##

path <- "C:/Export/Betsy/TDN/PEAK_AUC"

probe <- read.table("C:/Export/Betsy/TDN/20210106_TDN_Normalized.csv", header=TRUE, sep=",") #<---------- nSolver cols easy 

# probe <- subset(probe, !(Probe.Name %in% div$Gene))


# colnames(probe)[5:14] <-substr(colnames(probe)[5:14], 30, 32)
colnames(probe)[2:11] <-substr(colnames(probe)[2:11], 30, 32)



##Peak_G1xG2-------
samp_list <- c("202","207","211","201","204","213")

nq2 <- select(probe,c("Probe.Name",samp_list))

nq4 <- select(probe,samp_list)


# nq2$SD <- rowSds(as.matrix(nq2[,c(2:8)]))

# q3 <- q2[q2$SD != 0,] ##<----------------- for gene names
# 
# q4 <- q3[,c(2:8)]

# q5 <- select(probe2,"202","204","207","211","201","213","216")
design <- c(0,0,0,1,1,1)



iris2 <- transpose(nq4)

iris_mat2 <- as.matrix(iris2)

log_iris <- log10(iris_mat2)

fit_p <- pelora(log_iris, design, noc = 2, trace = 1)

fit_n  <- wilma(log_iris, design, noc = 2, trace = 1)


nq2[c(531,589,389,78),]

cat <- rbind(pfs_p,pfs_n)


write.table(cat, file=paste(path,"/","Peak_AUC_supclust_G1xG3_APH.tsv",sep=""), sep="\t",row.names=TRUE,quote = FALSE,col.names=TRUE)

##-------------new PHEAT ------------------------###

library(pheatmap)
library(RColorBrewer)
library(viridis)

probe <- read.table("C:/Export/Betsy/TDN/20210106_TDN_Normalized.csv", header=TRUE, sep=",") #<---------- nSolver cols easy 

tdn <- subset(tdn, !(gene_name %in% div$Gene))

pheat_list <- c("MAML3","NOTCH1","NOTCH2","IGF1R","AKT1","AKT2","PIK3C3")
##-------------###
# peak_list <- unique(cat$Probe.Name)

# qual <- tdn[tdn$adj == "q-value=<0.10",]
# 
# dif_peak <- unique(qual$gene_name)
# 
# GOI <- c(as.list(peak_list),as.list(dif_peak))


probe2 <- subset(probe, probe$Probe.Name %in% pheat_list)


# new_list5<- new_list4[!str_detect(new_list4,pattern="TR.")]

rownames(probe2) <- probe2$Probe.Name

colnames(probe2)[2:11] <-substr(colnames(probe2)[2:11], 30, 32)



pheat_probe <- select(probe2,"202","207","211","201","204","213","205","209","216","217")

pheat_probe <- select(probe2,"201","204","213","205","209","216","217")

# pheat_probe <- select(probe2,samp_list)


my_sample_col <- data.frame(Group = rep(c("G1","G2","G3"), c(3,3,4)))
row.names(my_sample_col) <- colnames(pheat_probe)


#------------------------------------------------------##

pheatmap(
  mat               = log10(pheat_probe),
  color             = inferno(20),
  # breaks            = mat_breaks,
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  annotation_col    = my_sample_col,
  annotation_colors = NA,
  drop_levels       = TRUE,
  fontsize          = 8,
  cluster_cols =  FALSE,
  main              = "UPCC35313 TDN log10 Transformed nSolver Normalized Expression Data\n Expression Profile Subset by MAML3-NOTCH-AKT Signaling"
  
)


write.table(probe2, file=paste(path,"/","pheat_DEG_G1xG3_APH.tsv",sep=""), sep="\t",row.names=TRUE,quote = FALSE,col.names=TRUE)



vio <- read.table("C:/Export/Betsy/TDN/PEAK_AUC/violin.csv", header=TRUE, sep=",")

vio_gene <- vio[vio$Gene == "TERT",]

p <- ggplot(vio_gene, aes(x=Group, y=Expr)) + geom_violin()

fin <- p + stat_summary(fun=median, geom="point", size=2, color="red") + geom_boxplot(width=0.1) + ggtitle("TERT")





##-----------------PATHway-------------------------------------##
library(org.Hs.eg.db)


gen_path <- read.table("C:/Export/Betsy/APH/G1xG3_peak_APH.csv", header=TRUE, sep=",")

entr <- read.table("C:/Export/Betsy/APH/sym_convert.csv", header=TRUE, sep=",")

entr <- subset(entr, (gene_name %in% pheat_list))

combo <- merge(entr,gen_path, by = "gene_name")

gen_path$gene_name <- row.names(gen_path)

gen_path <- subset(gen_path, (gene_name %in% pheat_list))

## assume that 1st column is ID
## 2nd column is fold change

## feature 1: numeric vector
geneList <- combo[,5]

## feature 2: named vector
names(geneList) <- as.character(combo[,2])

## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# eg = bitr("TNFRSF17", fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# 
# 
gene <- combo$ENTREZ
# 
# gene <-eg = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


edo <- enrichDGN(names(geneList))


kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

write.table(head(kk,30), file=paste(path,"/","ENRICH_DEG_G1xG3_APH.csv",sep=""), sep=",",row.names=TRUE,quote = FALSE,col.names=TRUE)



library("pathview")
try <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04060",
                     species    = "hsa",
                      limit      = 3)

library(pathview)
data(gse16873.d)
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",species = "hsa", out.suffix = "gse16873")
