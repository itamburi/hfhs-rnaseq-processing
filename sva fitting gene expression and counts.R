## Marcus Seldin, Ian Tamburini, Hosung Bae

# Process raw RNAseq counts from two batches
# Filter out samples, lowly expressed genes, and use combat() from SVA to apply batch correction


setwd("/hfhs-rnaseq-processing")

library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(WGCNA)
library(qvalue)
library(pheatmap)
library(qgraph)
#library(ggdendro)
library(RColorBrewer)
library(sva)
library(scales)

#### 1.0 - Format batch 1 and batch 2 rnaseq counts ####

# batch 1 rnaseq counts - 33 samples from normal chow animals
gene_counts = read.csv("./data/raw counts/counts filtered 10 in min2 rows.csv")
expr_sample_table = read.csv("./data/metadata/sample table for counts matrix.csv")

melted_counts = gene_counts
row.names(melted_counts) = melted_counts$gene_symbol
melted_counts$gene_symbol=NULL
melted_counts = gene_counts
row.names(melted_counts) = melted_counts$gene_symbol
melted_counts$gene_symbol=NULL
melted_counts[melted_counts==0] = NA
melted_counts$zer = rowSums(is.na(melted_counts))
melted_counts = melted_counts[melted_counts$zer < 6,]
melted_counts$zer=NULL
melted_counts = melt(as.matrix(melted_counts))
colnames(melted_counts) = c('variable_ID', 'Pig_ID', 'value')
melted_counts$data_type = paste0('raw RNA-Seq')
melted_counts$batch = paste0('counts1')
full_counts_table = melted_counts

# batch 2 rnaseq counts - 96 samples from high-fat-high-sucrose and normal chow animals
gene_counts = read.csv('.data/raw counts/full counts matrix HF 96 files.csv')
melted_counts = gene_counts
row.names(melted_counts) = melted_counts$X
melted_counts$X=NULL
melted_counts[melted_counts==0] = NA
melted_counts$zer = rowSums(is.na(melted_counts))
hist(melted_counts$zer)
melted_counts = melted_counts[melted_counts$zer < 85,]
melted_counts$zer=NULL
melted_counts = melt(as.matrix(melted_counts))
colnames(melted_counts) = c('variable_ID', 'Pig_ID', 'value')
melted_counts$data_type = paste0('raw RNA-Seq')
melted_counts$batch = paste0('counts2')
colnames(full_counts_table)

# combine batch1 and batch2 objects
full_counts_table = as.data.frame(rbind(melted_counts, full_counts_table))
summary(full_counts_table$value)
length(unique(full_counts_table$Pig_ID))
full_counts_table$value = ifelse(is.na(full_counts_table$value), as.numeric('0'), paste0(full_counts_table$value))
full_counts_table$value = as.numeric(full_counts_table$value)
full_counts_table$Pig_ID
# compute lower and upper whiskers
# scale y limits based on ylim1
#ggplot(full_counts_table, aes(x=Pig_ID, y=value, fill=Pig_ID)) + geom_boxplot(outlier.shape = NA)+ ylim(c(0,100))
#ggplot(full_counts_table, aes(x=Pig_ID, y=value, fill=Pig_ID)) + geom_boxplot(outlier.shape = NA)+ ylim(c(0,10)) +  geom_hline(yintercept = 0.5, col='red') + theme(legend.position="none")


#### 2.0 - Filter genes/samples based on expression values #### 

# Samples with mean count > 20
filtered_data = full_counts_table %>%  dplyr::select(Pig_ID, value) %>% dplyr::group_by(Pig_ID) %>%  dplyr::summarise(mean=mean(value))
hist(log10(filtered_data$mean))
ff1 = filtered_data[filtered_data$mean>20,]

# genes with mean count > 3
filtered_data = full_counts_table %>%  dplyr::select(variable_ID, value) %>% dplyr::group_by(variable_ID) %>%  dplyr::summarise(mean=mean(value))
hist(log10(filtered_data$mean))
ff2 = filtered_data[filtered_data$mean>3,]

mm1 = full_counts_table[full_counts_table$Pig_ID %in% ff1$Pig_ID,]
mm1 = mm1[mm1$variable_ID %in% ff2$variable_ID,]
bb1 = as.data.frame(table(mm1$Pig_ID))
bb1$batch = mm1$batch[match(bb1$Var1, mm1$Pig_ID)]
table(bb1$batch)
second_batch = bb1[bb1$batch=='counts2',]
mm1$nvalue = mm1$value+1.5
mm1$lvalue = log10(mm1$nvalue)
summary(mm1$lvalue)


#### 3.0 - SVA fitting of batch 1 and batch 2 using sva::ComBat() ####
mm2 = dcast(mm1, variable_ID ~ Pig_ID, value.var = 'lvalue', fun.aggregate = mean, na.rm=T)
row.names(mm2) = mm2$variable_ID
mm2$variable_ID=NULL

sva_rdy= mm2
batch_set = colnames(sva_rdy)
batch1 = ifelse(batch_set %in% second_batch$Var1, 'counts2', 'counts1')
table(batch1)
sva_rdy[1:10,1:10]
sva_rdy1 = sva_rdy %>%  mutate_all(~replace(., is.na(.), 0))
#sva_rdy1[is.na(sva_rdy1)]=0
combat_mydata= ComBat(dat=sva_rdy1, batch=batch1)
adj_data = as.data.frame(combat_mydata)
adj_data[is.na(adj_data)] =0


#### 4.0 - Rename columns to sample IDs ####
sample_meta_data = as.data.frame(colnames(adj_data))
colnames(sample_meta_data) = 'RNAseq_ID'

sample_meta_data$newID1 = expr_sample_table$Flux.sample.name[match(sample_meta_data$RNAseq_ID, expr_sample_table$sample_name)]
sample_meta_data$newID1

chow_flux = read.csv('./data/metadata/NC 96 pig samples sample table.csv')
hf_flux = read.csv('./data/metadata/HF 96 pig samples sample table.csv')

chow_flux$Sample.ID
sample_meta_data$newID1
summary(match(chow_flux$Sample.ID, sample_meta_data$newID1))

expr_sample_table$Flux.sample.name[1:10]

tt1 = melt(as.matrix(adj_data))

summary(tt1$value)
length(row.names(tt1[tt1$value<0.001,]))/length(row.names(tt1))
hc <- hclust(dist(t(adj_data)), "ave")
ggdendrogram(hc, rotate = FALSE, size = 2)


subset1 = sample_meta_data[grepl('bam_filesnR222.L4.G3.', sample_meta_data$RNAseq_ID),]
subset1$newID2 = gsub('bam_filesnR222.L4.G3.', '', subset1$RNAseq_ID, fixed = T) 
subset1$newID2 = gsub('_deduped.bam', '', subset1$newID2, fixed = T) 
subset1$newID2 = substring(subset1$newID2, 6)
subset1$newID2 = gsub("\\..*","",subset1$newID2)
subset1$newID2
subset1$sampleID = hf_flux$Sample.ID[match(subset1$newID2, hf_flux$i7.adapter)]
subset1$sampleID2 = chow_flux$Sample.ID[match(subset1$newID2, chow_flux$i7.adapter)]
subset1$newID = ifelse(is.na(subset1$sampleID), paste0(subset1$sampleID2), paste0(subset1$sampleID))
sample_meta_data$newID2 = subset1$newID[match(sample_meta_data$RNAseq_ID, subset1$RNAseq_ID)]
expr_sample_table$pig_ID = paste0(expr_sample_table$tissue_name_match, ' ', expr_sample_table$Pig)
sample_meta_data$newID3 = expr_sample_table$pig_ID[match(sample_meta_data$RNAseq_ID, expr_sample_table$seqID)]
sample_meta_data$final_ID = ifelse(is.na(sample_meta_data$newID2), paste0(sample_meta_data$newID3), paste0(sample_meta_data$newID2))

final_counts_matrix = adj_data
colnames(final_counts_matrix) = sample_meta_data$final_ID[match(colnames(final_counts_matrix), sample_meta_data$RNAseq_ID)]
# the colname format for high-fat-high-sucrose animals is "HF TISSUE animal#"
# the colname format for nomral chow-fed animals is "TISSUE animal#

write.csv(final_counts_matrix, './data/processed/fltered and sva adjusted counts - not integrated with flux.csv')


#### 5.0 - Combine the batch adjusted gene counts, and the flux values and perform combatseq() to fit gene counts to fluxes ####
hc <- hclust(dist(t(final_counts_matrix)), "ave")
col_scheme = hc$labels
pdf(file = 'reads-based clustering.pdf')
cols_list = ifelse(grepl('bam_files', col_scheme), 'darkorange2', 'dodgerblue2')
ggdendrogram(hc, rotate = F, size = 1)
dev.off()

mm1 = reshape2::melt(as.matrix(final_counts_matrix))
head(mm1)
colnames(mm1) = c('variable_ID', 'Tissue_pig', 'value')
mm1$data_type = paste0('RNA-seq') 
mm1 = mm1 %>% dplyr::select(Tissue_pig, everything())

hf1 = read.csv('./data/raw flux/HF flux cleaned_updated.csv', check.names = F)
row.names(hf1) = hf1$`Sample ID`
hf1$`Sample ID` = NULL
mm2 = reshape2::melt(as.matrix(hf1))
head(mm2)
colnames(mm2) = c('Tissue_pig', 'variable_ID', 'value')
mm2$data_type = paste0('flux')

hf1 = read.csv('./data/raw flux/chowflux cleaned_updated.csv', check.names = F)
row.names(hf1) = hf1$`Sample ID`
hf1$`Sample ID` = NULL
mm3 = reshape2::melt(as.matrix(hf1))
head(mm3)
colnames(mm3) = c('Tissue_pig', 'variable_ID', 'value')
mm3$data_type = paste0('flux')

full_data = as.data.frame(rbind(mm1, mm2, mm3))
nn1 = dcast(full_data, Tissue_pig ~ variable_ID, value.var = 'value', fun.aggregate = mean)
row.names(nn1) = nn1$Tissue_pig
nn1$Tissue_pig=NULL


sva_rdy1 = nn1 %>%  mutate_all(~replace(., is.na(.), 0))
#sva_rdy1[is.na(sva_rdy1)]=0
bb1 = colnames(sva_rdy1)
batch1 = full_data$data_type[match(bb1, full_data$variable_ID)]
table(batch1)
combat_mydata= ComBat(dat=sva_rdy1, batch=batch1)
adj_data = as.data.frame(combat_mydata)

nn1 = melt(as.matrix(adj_data))
head(nn1)
colnames(nn1) = c('Tissue_pig', 'Variable_ID', 'value')
nn1$data_type = full_data$data_type[match(nn1$Variable_ID, full_data$variable_ID)]
# write out the adjusted and fitted flux and gene counts
write.csv(nn1, file = './data/processed/batch-adjusted fitted and melted data ready for analysis.csv', row.names = F)

# write out the gene counts only
rr1 = nn1[nn1$data_type=='RNA-seq',]
rr2 = dcast(rr1, Variable_ID~Tissue_pig, value.var = 'value', fun.aggregate = mean)
write.csv(rr2, file = './data/processed/fitted filtered and batch adjusted counts combined.csv', row.names = F)

# write out the fluxes only
rr1 = nn1[!nn1$data_type=='RNA-seq',]
rr2 = dcast(rr1, Variable_ID~Tissue_pig, value.var = 'value', fun.aggregate = mean)
write.csv(rr2, file = './data/processed/fitted filtered and batch adjusted fluxes combined.csv', row.names = F)






