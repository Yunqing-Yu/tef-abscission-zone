####### load the libraries ######
library(DESeq2)
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(writexl)

###### import and format table for DESeq2 ######
file_directory <- "/Users/yunqingyu/Dropbox/postdoc/data/fonio_tef_ERAC/RNAseq/tef_AZ_202012/STAR_htseq_output"
list.files(file_directory)
#prepare sample table for DEseq analysis
Files <- grep('HTseq-count_star_dabbi.*.out$', list.files(file_directory), value=TRUE)
Files
Names <- gsub("HTseq-count_star_(.*)_[A|G|C|T]{6}.out", "\\1", Files)
Names
Conditions <- gsub("dabbi_(s[1-3])_([A-Z])([1-6])", "\\1-\\2-\\3", Names)
Conditions
Table <- data.frame(sampleName=Names, 
                          fileName=Files, 
                          condition=Conditions)
Table <- separate(Table, condition, into = c("stage", "tissue", "rep"), sep = "-")
Table$stage <- factor(Table$stage, levels = c("s1", "s2", "s3"))
Table$tissue <- factor(Table$tissue, levels = c("L", "A", "U"))
Table$rep <- factor(Table$rep, levels = c(1, 2, 3, 4, 5, 6))
#group the stage and tissue variable to reduce factors
Table$group <- factor(paste0(Table$stage, Table$tissue))
Table

####### input table into DESeq2 and DESeq2 analysis ######
HTseq <- DESeqDataSetFromHTSeqCount(sampleTable = Table, 
                                        directory = file_directory, 
                                        design = ~ rep + group)

####analysis with filtering at least six counts > 10 and at one average of each condition is larger than 10
HTseq1 <- estimateSizeFactors(HTseq)

## calculate average of each condition
HTseq_count <- counts(HTseq1, normalized=TRUE)
dim(HTseq_count)
HTseq_count_mean <- cbind(s1_A = apply(HTseq_count[,1:6], 1, mean), s1_L = apply(HTseq_count[,7:12], 1, mean), s1_U = apply(HTseq_count[,13:18], 1, mean),
                          s2_A = apply(HTseq_count[,19:23], 1, mean), s2_L = apply(HTseq_count[,24:29], 1, mean), s2_U = apply(HTseq_count[,30:32], 1, mean),
                          s3_A = apply(HTseq_count[,33:38], 1, mean), s3_L = apply(HTseq_count[,39:44], 1, mean), s3_U = apply(HTseq_count[,45:50], 1, mean))
dim(HTseq_count_mean)

HTseq2 <- HTseq1[(rowSums(counts(HTseq1, normalized=TRUE) >= 10) >= 6) & (rowSums(HTseq_count_mean >= 10) >= 1), ]


dds <- DESeq(HTseq2)
plotDispEsts(dds)

#get the normalized gene counts into a dataframe
counts <- as.data.frame(counts(dds, normalized = TRUE))
dim(counts)
names(counts)
counts <- counts[,c(7:12,1:6,13:18,24:29,19:23,30:32,39:44,33:38,45:50)]

#add a new column of gene names
row.names(counts) <- gsub("(.*),(maker.*,maker.*|snap.*,snap.*)", "\\1", rownames(counts))
counts$Et_geneID <- rownames(counts)
head(counts)
write.table(counts, "dabbi_counts.txt", sep = "\t", row.names = T)

#extracting transformed values
vst <- vst(dds, blind = FALSE)
vst_df <- as.data.frame(assay(vst))[,c(7:12,1:6,13:18,24:29,19:23,30:32,39:44,33:38,45:50)]
head(vst_df)
write.table(vst_df, "dabbi_vst.txt", sep = "\t", row.names = T)

#extract experiment design table
design_table <- as.data.frame(colData(dds))
write.table(design_table, "dabbi_design_table.txt", sep = "\t", row.names = T)

#PCA graph
pcaData <- plotPCA(vst, intgroup = c("stage", "tissue", "rep"), returnData = TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))
pdf("dabbi_PCA.pdf", width = 4, height = 4)
ggplot(pcaData, aes(x = PC1, y = PC2, color = stage, shape = tissue)) + geom_point() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))
dev.off()

#sample distances
sampleDists <- dist(t(assay(vst)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf("sampledistance.pdf", width = 9, height = 9)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()
 

###### differential expression analysis ######

#### tissue effect
##res_s1_A vs. L
res_s1_AvsL <- results(dds,contrast = c("group", "s1A", "s1L"), 
                               lfcThreshold = 1,
                               alpha = 0.05,
                               altHypothesis = "greaterAbs")
summary(res_s1_AvsL)

##res_s1_A vs. U
res_s1_AvsU <- results(dds,contrast = c("group", "s1A", "s1U"), 
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s1_AvsU)

##res_s2_A vs. L
res_s2_AvsL <- results(dds,contrast = c("group", "s2A", "s2L"), 
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s2_AvsL)

##res_s2_A vs. U
res_s2_AvsU <- results(dds,contrast = c("group", "s2A", "s2U"), 
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s2_AvsU)

##res_s3_A vs. L
res_s3_AvsL <- results(dds,contrast = c("group", "s3A", "s3L"), 
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s3_AvsL)

##res_s3_A vs. U
res_s3_AvsU <- results(dds,contrast = c("group", "s3A", "s3U"), 
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s3_AvsU)



#### stage effect
## res_s2 vs. s1_L
res_s2vs1_L <- results(dds, contrast = c("group", "s2L", "s1L"),
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s2vs1_L)

## res_s3 vs. s2_L
res_s3vs2_L <- results(dds, contrast = c("group", "s3L", "s2L"),
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s3vs2_L)

## res_s3 vs. s1_L
res_s3vs1_L <- results(dds, contrast = c("group", "s3L", "s1L"),
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s3vs1_L)

## res_s2 vs. s1_A
res_s2vs1_A <- results(dds, contrast = c("group", "s2A", "s1A"),
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s2vs1_A)

## res_s3 vs. s2_A
res_s3vs2_A <- results(dds, contrast = c("group", "s3A", "s2A"),
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s3vs2_A)

## res_s3 vs. s1_A
res_s3vs1_A <- results(dds, contrast = c("group", "s3A", "s1A"),
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s3vs1_A)

## res_s2 vs. s1_U
res_s2vs1_U <- results(dds, contrast = c("group", "s2U", "s1U"),
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s2vs1_U)

## res_s3 vs. s2_U
res_s3vs2_U <- results(dds, contrast = c("group", "s3U", "s2U"),
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s3vs2_U)

## res_s3 vs. s1_U
res_s3vs1_U <- results(dds, contrast = c("group", "s3U", "s1U"),
                       lfcThreshold = 1,
                       alpha = 0.05,
                       altHypothesis = "greaterAbs")
summary(res_s3vs1_U)


names(res_s1_AvsL) <- paste0(names(res_s1_AvsL), "_s1_AvsL")
names(res_s1_AvsU) <- paste0(names(res_s1_AvsU), "_s1_AvsU")
names(res_s2_AvsL) <- paste0(names(res_s2_AvsL), "_s2_AvsL")
names(res_s2_AvsU) <- paste0(names(res_s2_AvsU), "_s2_AvsU")
names(res_s3_AvsL) <- paste0(names(res_s3_AvsL), "_s3_AvsL")
names(res_s3_AvsU) <- paste0(names(res_s3_AvsU), "_s3_AvsU")
names(res_s2vs1_L) <- paste0(names(res_s2vs1_L), "_s2vs1_L")
names(res_s3vs2_L) <- paste0(names(res_s3vs2_L), "_s3vs2_L")
names(res_s3vs1_L) <- paste0(names(res_s3vs1_L), "_s3vs1_L")
names(res_s2vs1_A) <- paste0(names(res_s2vs1_A), "_s2vs1_A")
names(res_s3vs2_A) <- paste0(names(res_s3vs2_A), "_s3vs2_A")
names(res_s3vs1_A) <- paste0(names(res_s3vs1_A), "_s3vs1_A")
names(res_s2vs1_U) <- paste0(names(res_s2vs1_U), "_s2vs1_U")
names(res_s3vs2_U) <- paste0(names(res_s3vs2_U), "_s3vs2_U")
names(res_s3vs1_U) <- paste0(names(res_s3vs1_U), "_s3vs1_U")

# Add average to the result table
counts_mean <- HTseq_count_mean[rownames(counts),c(2,1,3,5,4,6,8,7,9)]
dim(counts_mean)

res <- cbind(res_s1_AvsL[,c(1,2,6)], res_s1_AvsU[,c(2,6)], 
             res_s2_AvsL[,c(2,6)], res_s2_AvsU[,c(2,6)], 
             res_s3_AvsL[,c(2,6)], res_s3_AvsU[,c(2,6)],
             res_s2vs1_L[,c(2,6)], res_s3vs2_L[,c(2,6)], res_s3vs1_L[,c(2,6)], 
             res_s2vs1_A[,c(2,6)], res_s3vs2_A[,c(2,6)], res_s3vs1_A[,c(2,6)],
             res_s2vs1_U[,c(2,6)], res_s3vs2_U[,c(2,6)], res_s3vs1_U[,c(2,6)], counts, counts_mean)[,c(1:81,83:91,82)]

dim(res)
names(res)
head(res)

##input Etef annotation
Etef_anno <- read.table("/Users/yunqingyu/Dropbox/postdoc/data/LCM/rna seq/reference genome/Eragrostis_tef_2019/ortho_finder/orthologs_Et_Os_anno.txt",
                        sep = "\t", header = T)
head(Etef_anno)
dim(Etef_anno)


##join the result table and annotation together
res_anno <- left_join(as.data.frame(res), Etef_anno, by = "Et_geneID")
dim(res_anno)
write.table(res_anno, "res_anno_dabbi_DEseq2.txt", row.names = F, sep = "\t")



########## make a table with significant differences in all of the comparisons ##########
#get the up and down regulated genes by 2 fold
gene_s1_AvsL_up <- rownames(subset(res_s1_AvsL, padj_s1_AvsL < 0.05 & log2FoldChange_s1_AvsL > 1))
gene_s1_AvsL_down <- rownames(subset(res_s1_AvsL, padj_s1_AvsL < 0.05 & log2FoldChange_s1_AvsL < -1))
gene_s1_AvsU_up <- rownames(subset(res_s1_AvsU, padj_s1_AvsU < 0.05 & log2FoldChange_s1_AvsU > 1))
gene_s1_AvsU_down <- rownames(subset(res_s1_AvsU, padj_s1_AvsU < 0.05 & log2FoldChange_s1_AvsU < -1))

gene_s2_AvsL_up <- rownames(subset(res_s2_AvsL, padj_s2_AvsL < 0.05 & log2FoldChange_s2_AvsL > 1))
gene_s2_AvsL_down <- rownames(subset(res_s2_AvsL, padj_s2_AvsL < 0.05 & log2FoldChange_s2_AvsL < -1))
gene_s2_AvsU_up <- rownames(subset(res_s2_AvsU, padj_s2_AvsU < 0.05 & log2FoldChange_s2_AvsU > 1))
gene_s2_AvsU_down <- rownames(subset(res_s2_AvsU, padj_s2_AvsU < 0.05 & log2FoldChange_s2_AvsU < -1))

gene_s3_AvsL_up <- rownames(subset(res_s3_AvsL, padj_s3_AvsL < 0.05 & log2FoldChange_s3_AvsL > 1))
gene_s3_AvsL_down <- rownames(subset(res_s3_AvsL, padj_s3_AvsL < 0.05 & log2FoldChange_s3_AvsL < -1))
gene_s3_AvsU_up <- rownames(subset(res_s3_AvsU, padj_s3_AvsU < 0.05 & log2FoldChange_s3_AvsU > 1))
gene_s3_AvsU_down <- rownames(subset(res_s3_AvsU, padj_s3_AvsU < 0.05 & log2FoldChange_s3_AvsU < -1))

gene_res_s2vs1_L_up <- rownames(subset(res_s2vs1_L, padj_s2vs1_L < 0.05 & log2FoldChange_s2vs1_L > 1))
gene_res_s2vs1_L_down <- rownames(subset(res_s2vs1_L, padj_s2vs1_L < 0.05 & log2FoldChange_s2vs1_L < -1))
gene_res_s3vs1_L_up <- rownames(subset(res_s3vs1_L, padj_s3vs1_L < 0.05 & log2FoldChange_s3vs1_L > 1))
gene_res_s3vs1_L_down <- rownames(subset(res_s3vs1_L, padj_s3vs1_L < 0.05 & log2FoldChange_s3vs1_L < -1))
gene_res_s3vs2_L_up <- rownames(subset(res_s3vs2_L, padj_s3vs2_L < 0.05 & log2FoldChange_s3vs2_L > 1))
gene_res_s3vs2_L_down <- rownames(subset(res_s3vs2_L, padj_s3vs2_L < 0.05 & log2FoldChange_s3vs2_L < -1))

gene_res_s2vs1_A_up <- rownames(subset(res_s2vs1_A, padj_s2vs1_A < 0.05 & log2FoldChange_s2vs1_A > 1))
gene_res_s2vs1_A_down <- rownames(subset(res_s2vs1_A, padj_s2vs1_A < 0.05 & log2FoldChange_s2vs1_A < -1))
gene_res_s3vs1_A_up <- rownames(subset(res_s3vs1_A, padj_s3vs1_A < 0.05 & log2FoldChange_s3vs1_A > 1))
gene_res_s3vs1_A_down <- rownames(subset(res_s3vs1_A, padj_s3vs1_A < 0.05 & log2FoldChange_s3vs1_A < -1))
gene_res_s3vs2_A_up <- rownames(subset(res_s3vs2_A, padj_s3vs2_A < 0.05 & log2FoldChange_s3vs2_A > 1))
gene_res_s3vs2_A_down <- rownames(subset(res_s3vs2_A, padj_s3vs2_A < 0.05 & log2FoldChange_s3vs2_A < -1))

gene_res_s2vs1_U_up <- rownames(subset(res_s2vs1_U, padj_s2vs1_U < 0.05 & log2FoldChange_s2vs1_U > 1))
gene_res_s2vs1_U_down <- rownames(subset(res_s2vs1_U, padj_s2vs1_U < 0.05 & log2FoldChange_s2vs1_U < -1))
gene_res_s3vs1_U_up <- rownames(subset(res_s3vs1_U, padj_s3vs1_U < 0.05 & log2FoldChange_s3vs1_U > 1))
gene_res_s3vs1_U_down <- rownames(subset(res_s3vs1_U, padj_s3vs1_U < 0.05 & log2FoldChange_s3vs1_U < -1))
gene_res_s3vs2_U_up <- rownames(subset(res_s3vs2_U, padj_s3vs2_U < 0.05 & log2FoldChange_s3vs2_U > 1))
gene_res_s3vs2_U_down <- rownames(subset(res_s3vs2_U, padj_s3vs2_U < 0.05 & log2FoldChange_s3vs2_U < -1))


#merge the genes from A vs. L and A vs. U and stage comparisons
gene_DEG_full <- unique(c(gene_s1_AvsL_up, gene_s1_AvsL_down, gene_s1_AvsU_up, gene_s1_AvsU_down, 
                          gene_s2_AvsL_up, gene_s2_AvsL_down, gene_s2_AvsU_up, gene_s2_AvsU_down, 
                          gene_s3_AvsL_up, gene_s3_AvsL_down, gene_s3_AvsU_up, gene_s3_AvsU_down, 
                          gene_res_s2vs1_L_up, gene_res_s2vs1_L_down, gene_res_s3vs1_L_up, gene_res_s3vs1_L_down, gene_res_s3vs2_L_up, gene_res_s3vs2_L_down, 
                          gene_res_s2vs1_A_up, gene_res_s2vs1_A_down, gene_res_s3vs1_A_up, gene_res_s3vs1_A_down, gene_res_s3vs2_A_up, gene_res_s3vs2_A_down, 
                          gene_res_s2vs1_U_up, gene_res_s2vs1_U_down, gene_res_s3vs1_U_up, gene_res_s3vs1_U_down, gene_res_s3vs2_U_up, gene_res_s3vs2_U_down))
length(gene_DEG_full)

res_DEG_full <- res_anno[res_anno$Et_geneID %in% gene_DEG_full,]
dim(res_DEG_full)
write.table(res_DEG_full, "res_dabbi_DEG_full.txt", row.names = FALSE, sep = "\t")

#get DEG genes in each comparisons
gene_DEG_s1_AvsLU <- unique(c(gene_s1_AvsL_up, gene_s1_AvsL_down, gene_s1_AvsU_up, gene_s1_AvsU_down))
gene_DEG_s1_Aup <- (intersect(gene_s1_AvsL_up, gene_s1_AvsU_up))
length(gene_DEG_s1_Aup)
gene_DEG_s1_Adown <- (intersect(gene_s1_AvsL_down, gene_s1_AvsU_down))
length(gene_DEG_s1_Adown)

gene_DEG_s2_AvsLU <- unique(c(gene_s2_AvsL_up, gene_s2_AvsL_down, gene_s2_AvsU_up, gene_s2_AvsU_down))
gene_DEG_s2_Aup <- (intersect(gene_s2_AvsL_up, gene_s2_AvsU_up))
length(gene_DEG_s2_Aup)
gene_DEG_s2_Adown <- (intersect(gene_s2_AvsL_down, gene_s2_AvsU_down))
length(gene_DEG_s2_Adown)

gene_DEG_s3_AvsLU <- unique(c(gene_s3_AvsL_up, gene_s3_AvsL_down, gene_s3_AvsU_up, gene_s3_AvsU_down))
gene_DEG_s3_Aup <- (intersect(gene_s3_AvsL_up, gene_s3_AvsU_up))
length(gene_DEG_s3_Aup)
gene_DEG_s3_Adown <- (intersect(gene_s3_AvsL_down, gene_s3_AvsU_down))
length(gene_DEG_s3_Adown)

gene_DEG_s2vs1_A <- unique(c(gene_res_s2vs1_A_up, gene_res_s2vs1_A_down))
length(gene_res_s2vs1_A_up)
length(gene_res_s2vs1_A_down)

gene_DEG_s3vs1_A <- unique(c(gene_res_s3vs1_A_up, gene_res_s3vs1_A_down))

gene_DEG_s3vs2_A <- unique(c(gene_res_s3vs2_A_up, gene_res_s3vs2_A_down))
length(gene_res_s3vs2_A_up)
length(gene_res_s3vs2_A_down)

gene_DEG_part <- unique(c(gene_DEG_s1_AvsLU, gene_DEG_s2_AvsLU, gene_DEG_s3_AvsLU, gene_DEG_s2vs1_A, gene_DEG_s3vs1_A, gene_DEG_s3vs2_A))
gene_DEG_AvsLU <- unique(c(gene_DEG_s1_AvsLU, gene_DEG_s2_AvsLU, gene_DEG_s3_AvsLU))
length(gene_DEG_part)      

res_DEG_list <- list(s1_AvsLU = res_anno[res_anno$Et_geneID %in% gene_DEG_s1_AvsLU,],
                     s2_AvsLU = res_anno[res_anno$Et_geneID %in% gene_DEG_s2_AvsLU,],
                     s3_AvsLU = res_anno[res_anno$Et_geneID %in% gene_DEG_s3_AvsLU,],
                     s2vs1_A = res_anno[res_anno$Et_geneID %in% gene_DEG_s2vs1_A,],
                     s3vs1_A = res_anno[res_anno$Et_geneID %in% gene_DEG_s3vs1_A,],
                     s3vs2_A = res_anno[res_anno$Et_geneID %in% gene_DEG_s3vs2_A,])

write_xlsx(res_DEG_list, "res_dabbi_DEG_AvsLU_A.xlsx")

res_DEG_AZ_list <- list(s1_Aup = res_anno[res_anno$Et_geneID %in% gene_DEG_s1_Aup,],
                     s1_Adown = res_anno[res_anno$Et_geneID %in% gene_DEG_s1_Adown,],
                     s2_Aup = res_anno[res_anno$Et_geneID %in% gene_DEG_s2_Aup,],
                     s2_Adown = res_anno[res_anno$Et_geneID %in% gene_DEG_s2_Adown,],
                     s3_Aup = res_anno[res_anno$Et_geneID %in% gene_DEG_s3_Aup,],
                     s3_Adown = res_anno[res_anno$Et_geneID %in% gene_DEG_s3_Adown,],
                     s2vs1_Aup = res_anno[res_anno$Et_geneID %in% gene_res_s2vs1_A_up,],
                     s2vs1_Adown = res_anno[res_anno$Et_geneID %in% gene_res_s2vs1_A_down,],
                     s3vs2_Aup = res_anno[res_anno$Et_geneID %in% gene_res_s3vs2_A_up,],
                     s3vs2_Adown = res_anno[res_anno$Et_geneID %in% gene_res_s3vs2_A_down,])

write_xlsx(res_DEG_AZ_list, "res_dabbi_DEG_AZ.xlsx")


############# make venn diagram to show up and down regulated genes #############
#make the table for running venn diagram funciton in limma
library(limma)
df <- data.frame(geneID = counts$Et_geneID, s1_AvsL = 0, s1_AvsU = 0, s2_AvsL = 0, s2_AvsU = 0, s3_AvsL = 0, s3_AvsU = 0,
                                            s1_Aup = 0, s1_Adown = 0, s2_Aup = 0, s2_Adown = 0, s3_Aup = 0, s3_Adown = 0,
                                            s2vs1_Aup = 0, s2vs1_Adown = 0, s3vs2_Aup = 0, s3vs2_Adown = 0, s3vs1_Aup = 0, s3vs1_Adown = 0)


#assign upregulated genes as 1 and downregulated genes as -1
df$s1_AvsL[df$geneID %in% gene_s1_AvsL_up] <- 1
df$s1_AvsL[df$geneID %in% gene_s1_AvsL_down] <- -1
df$s1_AvsU[df$geneID %in% gene_s1_AvsU_up] <- 1
df$s1_AvsU[df$geneID %in% gene_s1_AvsU_down] <- -1

df$s2_AvsL[df$geneID %in% gene_s2_AvsL_up] <- 1
df$s2_AvsL[df$geneID %in% gene_s2_AvsL_down] <- -1
df$s2_AvsU[df$geneID %in% gene_s2_AvsU_up] <- 1
df$s2_AvsU[df$geneID %in% gene_s2_AvsU_down] <- -1

df$s3_AvsL[df$geneID %in% gene_s3_AvsL_up] <- 1
df$s3_AvsL[df$geneID %in% gene_s3_AvsL_down] <- -1
df$s3_AvsU[df$geneID %in% gene_s3_AvsU_up] <- 1
df$s3_AvsU[df$geneID %in% gene_s3_AvsU_down] <- -1

df$s1_Aup[df$geneID %in% gene_DEG_s1_Aup] <- 1
df$s1_Adown[df$geneID %in% gene_DEG_s1_Adown] <- -1
df$s2_Aup[df$geneID %in% gene_DEG_s2_Aup] <- 1
df$s2_Adown[df$geneID %in% gene_DEG_s2_Adown] <- -1
df$s3_Aup[df$geneID %in% gene_DEG_s3_Aup] <- 1
df$s3_Adown[df$geneID %in% gene_DEG_s3_Adown] <- -1

df$s2vs1_Aup[df$geneID %in% gene_res_s2vs1_A_up] <- 1
df$s2vs1_Adown[df$geneID %in% gene_res_s2vs1_A_down] <- -1
df$s3vs1_Aup[df$geneID %in% gene_res_s3vs1_A_up] <- 1
df$s3vs1_Adown[df$geneID %in% gene_res_s3vs1_A_down] <- -1
df$s3vs2_Aup[df$geneID %in% gene_res_s3vs2_A_up] <- 1
df$s3vs2_Adown[df$geneID %in% gene_res_s3vs2_A_down] <- -1


head(df)

#run vennDiagram in limma package
pdf("Venn/Venn_s1_AvsLU.pdf", width = 3, height = 3)
vennDiagram(df[,2:3], include = c("up", "down"), mar=rep(0.1, 4), 
            circle.col = c("turquoise", "salmon"), cex = c(0.9, 0.5, 0.3))
#title(main = "s1_AvsLU", cex.main = 1)
dev.off()

pdf("Venn/Venn_s2_AvsLU.pdf", width = 3, height = 3)
vennDiagram(df[,4:5], include = c("up", "down"), mar=rep(0.1, 4), 
            circle.col = c("turquoise", "salmon"), cex = c(0.9, 0.5, 0.3))
#title(main = "s2_AvsLU", cex.main = 1)
dev.off()

pdf("Venn/Venn_s3_AvsLU.pdf", width = 3, height = 3)
vennDiagram(df[,6:7], include = c("up", "down"), mar=rep(0.1, 4), 
            circle.col = c("turquoise", "salmon"), cex = c(0.9, 0.5, 0.3))
#title(main = "s3_AvsLU", cex.main = 1)
dev.off()


pdf("Venn/Venn_Aup.pdf", width = 3, height = 3)
vennDiagram(df[,c(8, 10, 12)], mar=rep(0.1, 4), 
            circle.col = c("blue", "green", "red"), cex = c(0.9, 0.5, 0.3))
title(main = "upregulated genes in AZ", cex.main = 0.8)
dev.off()

pdf("Venn/Venn_Adown.pdf", width = 3, height = 3)
vennDiagram(df[,c(9, 11, 13)], mar=rep(0.1, 4), 
            circle.col = c("blue", "green", "red"), cex = c(0.9, 0.5, 0.3))
title(main = "downregulated genes in AZ", cex.main = 0.8)
dev.off()


pdf("Venn/Venn_AZupdown.pdf", width = 4, height = 4)
vennDiagram(df[,c(14, 16, 15, 17)], mar=rep(0.1, 4), 
            circle.col = c("red", "gold", "blue", "green"), cex = c(0.9, 0.9, 0.9))
title(main = "temporal regulated genes in AZ", cex.main = 0.8)
dev.off()





############# heatmap of known genes #############
vst_df <- read.table("dabbi_vst.txt", sep = "\t", header = T)
### YABBY
## Et_1B_011799 is not expressed not include in the heatmap
EtYAB <- c("Et_10A_000629", "Et_10B_002778", #YAB2L1
           #"Et_4B_037484", "Et_4A_033303", #SH1
           "Et_1B_013298", "Et_1A_008482", "Et_7A_052455", "Et_7B_055093", "Et_1A_006975", #FIL/YAB1
           "Et_4A_034697", "Et_4B_038843", #CRC
           "Et_5A_042135", "Et_5B_044811", #YAB3
           "Et_5B_044217", "Et_5A_041575") #YAB2L2


pdf("heatmap_dabbi_EtYAB.pdf", width = 8, height = 2.5)
pheatmap_out <- pheatmap(vst_df[EtYAB,], show_rownames = TRUE, cluster_cols = FALSE, show_colnames = FALSE,
                         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                         scale = "row",  #if use scale, the expression values are normalized by gene
                         breaks = seq(-3,3,length.out = 101),
                         annotation_col = as.data.frame(colData(dds)[,c("stage", "tissue")]),
                         annotation_colors = list(stage = c(s1 = "limegreen", s2 = "green4", s3 = "goldenrod"), 
                                                  tissue = c(L = "#8DD3C7", A = "#FFED6F", U ="#BEBADA")),
                         main = "EtYAB")
dev.off()

### Known AZ genes

EtAZ_known <- c("Et_4B_037484", "Et_4A_033303", #SH1
                "Et_3A_025418", "Et_3B_031395", #qSH1
                "Et_9A_061962", "Et_9B_064531", #SH5
                "Et_7B_053674", "Et_7A_051053", #SH4
                "Et_7B_054496", "Et_7A_052909", #SHAT1
                "Et_4A_034179", "Et_4B_038326", #Q
                "Et_3A_026212", "Et_3B_030583") #LES1

vst_df[EtAZ_known,]

pdf("heatmap_dabbi_knownAZgenes.pdf", width = 8, height = 2.5)
pheatmap_out <- pheatmap(vst_df[EtAZ_known,], show_rownames = TRUE, cluster_cols = FALSE, show_colnames = FALSE,
                         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                         scale = "row",  #if use scale, the expression values are normalized by gene
                         breaks = seq(-3,3,length.out = 101),
                         annotation_col = as.data.frame(colData(dds)[,c("stage", "tissue")]),
                         annotation_colors = list(stage = c(s1 = "limegreen", s2 = "green4", s3 = "goldenrod"), 
                                                  tissue = c(L = "#8DD3C7", A = "#FFED6F", U ="#BEBADA")),
                         main = "known AZ genes")
dev.off()



res_AZknown <- res_anno[res_anno$Et_geneID %in% c(EtYAB, EtAZ_known),]
write.table(res_AZknown, "res_dabbi_knownAZgenes.txt", row.names = FALSE, sep = "\t")



### Genes in GO floral organ abscission
Et_GO_AZ <- c("Et_3A_027120", "Et_3B_030237", "Et_6A_047003",  "Et_6B_049757", "Et_10A_001133", #BOP2
              "Et_2A_017313", #ATH1
              "Et_10B_003967")  #ARF

vst_df[Et_GO_AZ,]

pdf("heatmap_dabbi_GO_AZ.pdf", width = 8, height = 1.6)
pheatmap_out <- pheatmap(vst_df[Et_GO_AZ,], show_rownames = TRUE, cluster_cols = FALSE, show_colnames = FALSE,
                         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                         scale = "row",  #if use scale, the expression values are normalized by gene
                         #breaks = seq(-3,3,length.out = 101),
                         annotation_col = as.data.frame(colData(dds)[,c("stage", "tissue")]),
                         annotation_colors = list(stage = c(s1 = "limegreen", s2 = "green4", s3 = "goldenrod"), 
                                                  tissue = c(L = "#8DD3C7", A = "#FFED6F", U ="#BEBADA")),
                         main = "GO:floral abscission")
dev.off()



Et_celldeath_blue <- c("Et_9A_061611", "Et_4A_033777", "Et_2B_021147", "Et_1A_007251", "Et_3B_031211", "Et_9A_063238", "Et_4B_038061", "Et_2B_021099", 
                       "Et_7B_054335", "Et_1B_011453", "Et_2B_019537", "Et_7A_051649", "Et_8A_056928", "Et_1A_007909", "Et_1A_007232", "Et_3B_031043", 
                       "Et_4A_034873", "Et_3B_028127", "Et_2B_018998", "Et_2B_020845", "Et_5A_042191", "Et_9B_064873", "Et_1B_012704", "Et_5A_040664", 
                       "Et_1B_012043", "Et_1B_013036", "Et_3A_027264", "Et_5B_045119", "Et_4B_037930", "Et_9A_062318", "Et_6A_047240", "Et_2A_016600", 
                       "Et_1A_005085", "Et_9B_064802", "Et_5A_042834", "Et_9B_064312", "Et_4B_038027", "Et_8B_059106", "Et_4A_035047", "Et_1B_010652", 
                       "Et_9B_065111", "Et_4A_035691", "Et_10B_003111", "Et_3B_031597", "Et_3B_029882", "Et_5A_042990", "Et_9A_062502", "Et_8A_057241", 
                       "Et_6A_047279", "Et_10B_003812", "Et_3A_024548", "Et_4B_038062", "Et_9A_062538", "Et_1A_009251", "Et_7A_051832", "Et_8B_060059", 
                       "Et_5B_044315", "Et_4A_034699", "Et_2B_019519", "Et_3B_028790", "Et_3A_024380", "Et_9A_061635", "Et_6B_049043", "Et_8B_059601", 
                       "Et_2B_021851", "Et_4A_032967", "Et_7A_050419", "Et_5A_041184", "Et_1B_013844", "Et_7B_055501", "Et_4B_039831", "Et_4B_037132", "Et_5B_044372", "Et_2A_015248")
res_GOcelldeath_blue <- res_anno[res_anno$Et_geneID %in% Et_celldeath_blue,]
write.table(res_GOcelldeath_blue, "res_dabbi_GOcelldeath_blue.txt", row.names = FALSE, sep = "\t")

Et_celldeath_violet <- c("Et_1A_008280", "Et_2A_018642", "Et_2B_021942", "Et_2A_017387", "Et_4B_039920", "Et_1B_010626", "Et_2A_016784", "Et_4B_037721", 
                         "Et_4B_037315", "Et_2B_021372", "Et_1B_012887", "Et_2A_017914", "Et_3A_024938", "Et_2B_020974", "Et_1A_006624", "Et_4B_038321", 
                         "Et_1A_005728", "Et_1B_013708", "Et_2A_018302", "Et_3A_024972", "Et_5B_044899", "Et_10A_002029", "Et_2A_016944", "Et_3B_030288", 
                         "Et_7B_054416", "Et_1A_008817", "Et_8A_058063", "Et_6A_047434", "Et_9A_061731", "Et_4B_038267", "Et_2B_020427", "Et_6B_049357", 
                         "Et_4B_039184", "Et_4A_032825", "Et_4A_033572", "Et_4A_034137", "Et_9A_063502", "Et_3A_024176", "Et_1A_009233", "Et_4A_033874", 
                         "Et_4B_038847", "Et_1B_013903", "Et_4B_037026", "Et_2B_020304", "Et_3B_029365", "Et_1B_011471", "Et_8B_060375", "Et_6B_049227", "Et_8B_060466")
res_GOcelldeath_violet <- res_anno[res_anno$Et_geneID %in% Et_celldeath_blue,]
write.table(res_GOcelldeath_violet, "res_dabbi_GOcelldeath_violet.txt", row.names = FALSE, sep = "\t")

Et_aging_blue <- c("Et_1A_008002", "Et_2A_016778", "Et_5A_041990", "Et_4B_039323", "Et_1B_010966", "Et_7A_052829", "Et_3B_030142", "Et_8A_057709", 
                   "Et_2B_020141", "Et_4A_034816", "Et_2A_016027", "Et_4A_034468", "Et_4A_032649", "Et_2A_017503", "Et_9A_063414", "Et_8A_057002", 
                   "Et_8A_056436", "Et_1B_012889", "Et_2A_017107", "Et_5A_040341", "Et_8B_059272", "Et_4B_036833", "Et_1A_008764", "Et_1A_004824", 
                   "Et_4A_032554", "Et_2B_020372", "Et_3A_025041", "Et_3A_025775", "Et_1A_004906", "Et_1A_008300", "Et_3A_023512", "Et_3A_025830", 
                   "Et_3B_030431", "Et_1B_010518", "Et_4A_035047", "Et_8A_057314", "Et_1A_006615", "Et_2A_016272", "Et_2A_016170", "Et_4B_038675", 
                   "Et_9B_065467", "Et_7A_050807", "Et_7A_052761", "Et_2B_022848", "Et_2B_020012", "Et_3B_030362", "Et_2A_018184", "Et_8B_060257", 
                   "Et_1B_011458", "Et_3B_031408", "Et_3A_024410", "Et_3A_026006", "Et_9A_062929", "Et_2B_018949", "Et_2B_019838", "Et_7A_050808", 
                   "Et_4B_037935", "Et_10A_002251", "Et_9B_065090", "Et_3A_027321", "Et_4B_038880", "Et_3B_028790", "Et_9B_064413", "Et_9A_063413", 
                   "Et_5A_042838", "Et_3A_026086", "Et_4A_032496", "Et_1B_014381", "Et_1A_009558", "Et_4B_037445", "Et_8B_059675", "Et_4B_036683", 
                   "Et_5A_042939", "Et_3B_030451", "Et_9B_064167", "Et_3A_024891", "Et_5B_043726", "Et_3B_029275")
res_GOage_blue <- res_anno[res_anno$Et_geneID %in% Et_aging_blue,]
write.table(res_GOage_blue, "res_dabbi_GOage_blue.txt", row.names = FALSE, sep = "\t")

Et_aging_violet <- c("Et_2B_019726", "Et_4B_038971", "Et_1A_006743", "Et_9A_061739", "Et_1B_013383", "Et_1A_008693", "Et_8B_059909", "Et_6B_049643", 
                     "Et_8B_060837", "Et_1A_006114", "Et_5A_042723", "Et_2A_018302", "Et_2A_018303", "Et_5B_045368", "Et_5B_044749", "Et_4A_034340", 
                     "Et_10A_000515", "Et_6A_047434", "Et_3A_026165", "Et_4B_039184", "Et_3B_030009", "Et_4B_036680", "Et_1A_007447", "Et_1B_012266", 
                     "Et_4A_032495", "Et_9A_062264", "Et_3B_030197", "Et_2B_020304", "Et_2A_017879", "Et_4B_039098", "Et_4B_038481", "Et_5B_045668", 
                     "Et_1B_012993", "Et_8A_057935", "Et_5B_044801", "Et_2B_020897", "Et_4A_033275", "Et_3A_025653", "Et_8B_060466")
res_GOage_violet <- res_anno[res_anno$Et_geneID %in% Et_aging_violet,]
write.table(res_GOage_violet, "res_dabbi_GOage_violet.txt", row.names = FALSE, sep = "\t")

library(readxl)
hub_blue <- read_xlsx("/Users/yunqingyu/Dropbox/postdoc/data/fonio_tef_ERAC/RNAseq/tef_AZ_202012/WGCNA/normalized_counts/hub_genes/dabbi_AZmodules_hub_kIM_KME_consensus.xlsx", sheet = "blue")$Et_geneID
hub_violet <- read_xlsx("/Users/yunqingyu/Dropbox/postdoc/data/fonio_tef_ERAC/RNAseq/tef_AZ_202012/WGCNA/normalized_counts/hub_genes/dabbi_AZmodules_hub_kIM_KME_consensus.xlsx", sheet = "violet")$Et_geneID

intersect(Et_celldeath_blue, hub_blue)
intersect(Et_celldeath_violet, hub_violet)
intersect(Et_aging_blue, hub_blue)
intersect(Et_aging_violet, hub_violet)

res_GOcelldeath_hub_blue <- res_anno[res_anno$Et_geneID %in% intersect(Et_celldeath_blue, hub_blue),]
write.table(res_GOcelldeath_hub_blue, "res_dabbi_GOcelldeath_hub_blue.txt", row.names = FALSE, sep = "\t")

vst_df[intersect(Et_celldeath_blue, hub_blue),]

pdf("heatmap_dabbi_GOcelldeath_hub_blue.pdf", width = 8, height = 2)
pheatmap_out <- pheatmap(vst_df[intersect(Et_celldeath_blue, hub_blue),], show_rownames = TRUE, cluster_cols = FALSE, show_colnames = FALSE,
                         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                         scale = "row",  #if use scale, the expression values are normalized by gene
                         breaks = seq(-3,3,length.out = 101),
                         annotation_col = as.data.frame(colData(dds)[,c("stage", "tissue")]),
                         annotation_colors = list(stage = c(s1 = "limegreen", s2 = "green4", s3 = "goldenrod"), 
                                                  tissue = c(L = "#8DD3C7", A = "#FFED6F", U ="#BEBADA")),
                         main = "GO:cell death blue hub")
dev.off()

res_GOcelldeath_hub_violet <- res_anno[res_anno$Et_geneID %in% intersect(Et_celldeath_violet, hub_violet),]
write.table(res_GOcelldeath_hub_violet, "res_dabbi_GOcelldeath_hub_violet.txt", row.names = FALSE, sep = "\t")

vst_df[intersect(Et_celldeath_violet, hub_violet),]

pdf("heatmap_dabbi_GOcelldeath_hub_violet.pdf", width = 8, height = 1.5)
pheatmap_out <- pheatmap(vst_df[intersect(Et_celldeath_violet, hub_violet),], show_rownames = TRUE, cluster_cols = FALSE, show_colnames = FALSE,
                         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                         scale = "row",  #if use scale, the expression values are normalized by gene
                         breaks = seq(-3,3,length.out = 101),
                         annotation_col = as.data.frame(colData(dds)[,c("stage", "tissue")]),
                         annotation_colors = list(stage = c(s1 = "limegreen", s2 = "green4", s3 = "goldenrod"), 
                                                  tissue = c(L = "#8DD3C7", A = "#FFED6F", U ="#BEBADA")),
                         main = "GO:cell death violet hub")
dev.off()

res_GOaging_hub_blue <- res_anno[res_anno$Et_geneID %in% intersect(Et_aging_blue, hub_blue),]
write.table(res_GOaging_hub_blue, "res_dabbi_GOaging_hub_blue.txt", row.names = FALSE, sep = "\t")

res_GOaging_hub_violet <- res_anno[res_anno$Et_geneID %in% intersect(Et_aging_violet, hub_violet),]
write.table(res_GOaging_hub_violet, "res_dabbi_GOaging_hub_violet.txt", row.names = FALSE, sep = "\t")





