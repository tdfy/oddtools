library(dplyr)
library(tidyr)
library(tibble)
library(NACHO)


excel <- load_rcc(data_directory = "C:/Export/Betsy/TDN/20201120_208802640621_RCC", ssheet_csv = "sample_sheet.csv", id_colname = "File",housekeeping_genes = "ABCF1,G6PD,NRDE2,OAZ1,POLR2A,SDHA,STK11IP,TBC1D10B,TBP,UBB", housekeeping_predict = FALSE,housekeeping_norm = TRUE, normalisation_method = "GEO",n_comp = 10)

expr_counts <- excel[["nacho"]] %>% 
  filter(grepl("Endogenous", CodeClass)) %>% 
  select(File, Name, Count_Norm) %>% 
  pivot_wider(names_from = "Name", values_from = "Count_Norm") %>% 
  column_to_rownames("File") %>% 
  t()

get_counts <- function(
  nacho, 
  codeclass = "Endogenous", 
  rownames = "IDFILE", 
  colnames = c("Name", "Accession")
) {
  nacho[["nacho"]] %>% 
    dplyr::filter(grepl(codeclass, .data[["CodeClass"]])) %>% 
    dplyr::select(c("IDFILE", "Name", "Count_Norm")) %>% 
    tidyr::pivot_wider(names_from = colnames[1], values_from = "Count_Norm") %>%
    tibble::column_to_rownames(rownames) %>% 
    t()
}


##---------------------------------------------------------------------------------------------###

library(vsn)
library(NanoStringNorm)

tdn <- read.table("C:/Export/Betsy/TDN2.csv", header=TRUE, sep=",")


tdn_norm <- NanoStringNorm(tdn,anno = NA, header = NA, Probe.Correction.Factor = 'adjust', CodeCount = 'geo.mean', Background = 'none', SampleContent = 'housekeeping.geo.mea',OtherNorm = 'none',CodeCount.summary.target = NA,SampleContent.summary.target = NA,round.values = FALSE,is.log = FALSE,take.log = FALSE,return.matrix.of.endogenous.probes = TRUE)
#                            traits = NA,
#                            predict.conc = FALSE,
#                            verbose = TRUE,
#                            genes.to.fit,
#                            genes.to.predict,
#                            guess.cartridge = TRUE,
#                            ...
# ))




write.table(tdn_norm, file=paste(path,"/","TDN_Norm.tsv",sep=""), sep="\t",row.names=TRUE,quote = FALSE,col.names=TRUE)

tdn_norm2 <- NanoStringNorm(tdn, anno = NA, header = NA, Probe.Correction.Factor = 'adjust', CodeCount = 'geo.mean', Background = 'none', SampleContent = 'housekeeping.geo.mean',OtherNorm = 'none',CodeCount.summary.target = NA,SampleContent.summary.target = NA,round.values = FALSE,is.log = FALSE,take.log = FALSE,return.matrix.of.endogenous.probes = FALSE)

Plot.NanoStringNorm(tdn_norm2,label.best.guess = TRUE,label.ids = list(),label.as.legend = TRUE,plot.type='all')


Plot.NanoStringNorm(x = tdn_norm2,label.best.guess = FALSE,label.ids = list(genes = rownames(tdn_norm2$gene.summary.stats.norm),
                                                                            samples = rownames(tdn_norm2$sample.summary.stats)),plot.type = c('norm.factors'))

                    
                    
                    
png('NanoStringNorm_4_Plots_%03d.png', units = 'in', height = 6,width = 6, res = 250, pointsize = 10)
Plot.NanoStringNorm(x = tdn_norm3,label.best.guess = TRUE,plot.type = c('cv', 'mean.sd', 'RNA.estimates', 'volcano', 'missing','norm.factors', 'positive.controls', 'batch.effects'))
dev.off()



png('test_%03d.png', units = 'in', height = 6,width = 6, res = 250, pointsize = 10)
Plot.NanoStringNorm(x = tdn_norm2,label.best.guess = TRUE,plot.type = c('volcano'))
dev.off()


tdn_norm3 <- NanoStringNorm(tdn, anno = NA, header = NA, Probe.Correction.Factor = 'adjust', CodeCount = 'geo.mean', Background = 'mean.2sd', SampleContent = 'housekeeping.geo.mean',OtherNorm = 'none',CodeCount.summary.target = NA,SampleContent.summary.target = NA,round.values = FALSE,is.log = FALSE,take.log = FALSE,return.matrix.of.endogenous.probes = FALSE)

















##-------------------------------------------------------##
library(ggrepel)
library(ggplot2)


tdn <- read.table("C:/Export/Betsy/TDN/G1xG3_indeed.csv", header=TRUE, sep=",")

ggplot(data=tdn, aes(x=logFC, y=pvalue)) + geom_point()

p <- ggplot(data=tdn, aes(x=logFC, y=-log10(pvalue))) + geom_point()+ theme_minimal()


p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

tdn$logFC <- tdn$logFC*log10(2)


tdn$diffexpressed <- "NO"
tdn$diffexpressed[tdn$logFC > 0.6 & tdn$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
tdn$diffexpressed[tdn$logFC < -0.6 & tdn$pvalue < 0.05] <- "DOWN"

# tdn <- subset(tdn, Gene %in% rownames(query4)) #<------------------------ filter
ranger <- "NA"
tdn$ranger[tdn$logFC > -5 & tdn$logFC < 5] <- "select" #<------------------------ filter
tdn <- subset(tdn, ranger == "select")

# 
# p <- ggplot(data=tdn, aes(x=logFC, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
# 
# p2 <- p + geom_vline(xintercept=c(-1.0, 1,0), col="red") +
#   geom_hline(yintercept=-log10(0.05), col="red")
# 
# p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

tdn$delabel <- NA
tdn$delabel[tdn$diffexpressed != "NO"] <- tdn$Gene[tdn$diffexpressed != "NO"]


tdn$delabel[tdn$diffexpressed != "NO"] 

# ggplot(data=tdn, aes(x=logFC, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
#   geom_point() + 
#   theme_minimal() +
#   geom_text()


ggplot(data=tdn, aes(x=logFC, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


###-----------------------------------------------------------------------###

library(pheatmap)
library(RColorBrewer)
library(viridis)


probe <- read.table("C:/Export/Betsy/TDN/hope_indeed.csv", header=TRUE, sep=",")

cols <- c(1, 14:23)
probe2 <- probe[,cols]


#DEG list from volcano plot
DEG <- tdn$Gene[tdn$diffexpressed != "NO"]


query <- subset(probe2, Gene %in% DEG)

rownames(query) <- query$Gene

cols2 <- c(2:11)

query2 <- query[,cols2]

query2$tot <- rowSums(query2[,c(1:10)])
query3 <- subset(query2, tot >= 5)

query3[query3==1]<-1.5
query3[query3==0]<-1


ord <- c(2,3,5,7,1,8,9,4,6,10)

query4 <- query3[,ord]

colnames(query4) <-substr(colnames(query4), 55, 57)

# as.numeric(query3)
# 
# test <-as.matrix(query3)
# 
# 
# 
# dat <- data.frame(values = as.numeric(test))
# ggplot(dat, aes(values)) + geom_density(bw = "SJ")
# 
#------------------------------------------------------##
# col_groups <- substr(colnames(query4), 1, 1)
# table(col_groups)
# # Data frame with column annotations.
# mat_col <- data.frame(group = col_groups)
# rownames(mat_col) <- colnames(mat)
# 
# # List with colors for each annotation.
# mat_colors <- list(group = brewer.pal(3, "Set1"))
# names(mat_colors$group) <- unique(col_groups)

my_sample_col <- data.frame(Group = rep(c("G1", "G2", "G3"), c(4,3,3)))
row.names(my_sample_col) <- colnames(query4)


#------------------------------------------------------##

pheatmap(
  mat               = log10(query4),
  color             = inferno(20),
  # breaks            = mat_breaks,
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  annotation_col    = my_sample_col,
  annotation_colors = NA,
  drop_levels       = TRUE,
  fontsize          = 14,
  cluster_cols =  FALSE,
  main              = "UPCC35313 TDN log10 Transformed NanoStringDiff Expression Data"
  
)

###---------------------------------##
library(dplyr)
library(data.table)
library(supclust)
data(leukemia, package="supclust")

probe <- read.table("C:/Export/Betsy/TDN/hope_indeed.csv", header=TRUE, sep=",")

cols <- c(1, 14:23)
probe2 <- probe[,cols]


# ord <- c(2,3,5,7,1,8,9,4,6,10)

colnames(probe2) <-substr(colnames(probe2), 55, 57)

# probe2 <- probe2[,ord]

q <- select(probe2,"202","204","207","211","205","209","217")

# probes_transpose <- as.data.frame(t(as.matrix(probe2)))

iris <- transpose(q)

iris_mat <- as.matrix(iris)

# colnames(iris) <- c(1:776)


design <- c(0,0,0,0,1,1,1)



fit  <- wilma(iris_mat, design, noc = 2, trace = 1)

summary(fit)
plot(fit)
fitted(fit)


##--------------------------------------------------------###

q2 <- select(probe2,"202","204","207","211","201","213","216")

# probes_transpose <- as.data.frame(t(as.matrix(probe2)))

iris2 <- transpose(q2)

iris_mat2 <- as.matrix(iris2)

# colnames(iris) <- c(1:776)


design <- c(0,0,0,0,1,1,1)



fit  <- wilma(iris_mat2, design, noc = 2, trace = 1)

###----------------VIOLIN---------------------------##

vio <- read.table("C:/Export/Betsy/TDN/supclust_vio.csv", header=TRUE, sep=",")

vio_gene <- vio[vio$Gene == "PTGER2",]

p <- ggplot(vio_gene, aes(x=Group, y=Expr)) + geom_violin()

fin <- p + stat_summary(fun=median, geom="point", size=2, color="red") + geom_boxplot(width=0.1) + ggtitle("PTGER2")
