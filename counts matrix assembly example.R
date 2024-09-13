
setwd('')

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(PoiClaClu)
library(ggrepel)
library(BiocParallel)
library(glmpca)
register(SerialParam())
library(DESeq2)


dir <- file.path('/deduped bam files/')
list.files(dir)
csvfile <- file.path(dir, "sample_table_pigs.csv")
coldata <- read.csv(csvfile, stringsAsFactors=FALSE)
coldata
coldata$seqID[1:10]
coldata$seqID = paste0('P', coldata$Sample_ID, '_deduped.bam')
row.names(coldata) = coldata$seqID
coldata$names <- coldata$seqID
coldata$files <- file.path(dir,  paste0(coldata$seqID))
file.exists(coldata$files)
bamfiles <- BamFileList(coldata$files)
seqinfo(bamfiles[1])

#Include your file name for gtf
gtffile <- file.path(dir, "Sscrofa11.1_genomic.gff")
#change to gtf
txdb <- makeTxDbFromGFF(gtffile, format = "gff", circ_seqs = character())
txdb

ebg <- exonsBy(txdb, by="gene")
ebg

seqinfo(bamfiles)
#This will take some time
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=FALSE,
                        fragments=FALSE)

counts_matrix = as.data.frame(assay(se))

colData(se)$group = coldata$Tissue
colData(se)

dds =  DESeqDataSet(se, design = ~ group)
#filtering by number of counts in individuals
keep <- rowSums(counts(dds) >= 10) >= 2

#filtering rowsums
#keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
nrow(dds)

gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$group <- dds$group

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = group, label=row.names(gpca.dat))) +
  geom_point(size =2.5)  + ggtitle("glmpca - Generalized PCA") + theme_minimal() + geom_label_repel( label.size = NA,  
                                                                                                     label.padding=.1, 
                                                                                                     na.rm=TRUE,
                                                                                                     fill = alpha(c("white"),0.5), size = 2, max.overlaps = Inf)



poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$group, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd, main = 'Poisson clustering Distances')

new_counts = as.data.frame(assay(dds))
new_counts$gene_symbol = row.names(new_counts)
nn1 = new_counts %>% dplyr::select(gene_symbol, everything())
write.csv(nn1, file = 'counts filteres 10 in min2 rows.csv', row.names = F)
sample_table = coldata
sample_table$sample_name = row.names(sample_table)
ss1 = sample_table %>%  dplyr::select(sample_name, everything())
write.csv(ss1, file ='sample table for counts matrix.csv', row.names = F)





