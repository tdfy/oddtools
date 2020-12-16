setwd("c:/Export/Tang/TCR_Seq/local_gamma/")
library(ggplot2)
library(GenomicFeatures)
library(ggbio)
library(rtracklayer)
library(GenomicRanges)
library(magick)
library(dplyr)

 
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# ###======================================================================#####


try <- "sort_YHT-1_1.sam.depth.txt"
depth <- read.table(try,header=F,sep="\t")
depth[,1] = "chr7"
depth[,2] = depth[,2] + 38240024

try2 <- "sort_YHT-3_3.sam.depth.txt"
depth2 <- read.table(try2,header=F,sep="\t")
depth2[,1] = "chr7"
depth2[,2] = depth2[,2] + 38240024


try3 <- "sort_YHT-4_4.sam.depth.txt"
depth3 <- read.table(try3,header=F,sep="\t")
depth3[,1] = "chr7"
depth3[,2] = depth3[,2] + 38240024

dep <- ggplot() + 
  geom_line(data = depth, aes(x=V2, y=V3), color = "blue",alpha=0.5) +
  geom_line(data = depth2, aes(x=V2, y=V3), color = "red",alpha=0.5) +
  geom_line(data = depth3, aes(x=V2, y=V3), color = "green",alpha=0.5) + theme_bw()
  

# dep <- ggplot(depth, aes(x=V2, y=V3)) + geom_bar()+theme_bw()+scale_y_continuous(position = "left")+ylab(NULL)
        
# ###======================================================================#####

rep <- "repELM.bed"
rep_elm <- read.table(rep,header=F,sep="\t")

rep_elm[,1] = "chr7"  

avs1.granges <- GRanges(seqnames=rep_elm[,1],
                        ranges=IRanges(start=rep_elm[,2],
                                       end=rep_elm[,3]),
                        strand="*")


mark <- autoplot(avs1.granges,geom_rect())+ theme_bw()

        
# ###======================================================================#####

# hip <- "trans_annotate.csv"
# anno <- read.table(hip,header=T,sep=",")
# 
# 
# alignID <- list(anno$alignID)
# names(alignID) <- list(anno$Gene)


##__________________________________________________________________________##

locus <- GRanges(seqnames='chr7',
                 ranges=IRanges(start=38240024,
                                end=38368055),
                 strand="*")


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene




transcript <- autoplot(txdb, which=locus,names.expr ="tx_name") + theme_bw()

# gr <- subset(transcripts(txdb), tx_name == "uc062uto.1")

# tx <- makeTxDbFromUCSC(
#   genome="hg38",
#   tablename="knownGene",
#   transcript_ids="uc062uto.1",
#   url="http://genome.ucsc.edu/cgi-bin/",
#   goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath" )
# 
# p <- autoplot(tx,
#               which = gr,
#               geom = "full",
#               coord = "genome",
#               layout = "linear" )

###---------------------------------------------------------------##
# geneid <- mapIds(org.Hs.eg.db, "ENTREZID", "SYMBOL")
# 
# hg38.genes <- genes(keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene))
# geneid %in% hg38.genes$gene_id
# 
# ###======================================================================#####
# 
sample2 <- tracks("Coverage"=dep,"Rep\nElements"=mark,"Transcript"=transcript,heights=c(0.1,0.01,0.1),label.text.cex=c(1,1),title = "TCR Gamma YHT 1-1",label.bg.fill=('white'))

     
# sample2 <- tracks("Coverage"=dep,"TRANS"=transcript,heights=c(0.1,0.1),label.text.cex=c(1,1),title = "TCR Gamma YHT",label.bg.fill=('white'))
