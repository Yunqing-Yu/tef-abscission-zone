################################################################################
####################        WGCNA module detection        ######################
################################################################################
# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
allowWGCNAThreads()


######Step 1###### Data input and cleaning
# Read in the data set. Filter the genes with variance smaller than 0.05
gene_expressed_counts = read.table("/Users/yunqingyu/Dropbox/postdoc/data/fonio_tef_ERAC/RNAseq/tef_AZ_202012/DEseq2/dabbi_counts.txt", header = TRUE)
gene_expressed_vst = read.table("/Users/yunqingyu/Dropbox/postdoc/data/fonio_tef_ERAC/RNAseq/tef_AZ_202012/DEseq2/dabbi_vst.txt", header = TRUE)
sum(gene_expressed_counts$Et_geneID == rownames(gene_expressed_vst))
gene_expressed = gene_expressed_counts[,1:50][apply(gene_expressed_vst[,1:50], 1, var) > 0.05,] 

# Take a quick look at what is in the data set:
dim(gene_expressed);
names(gene_expressed);

# Data exrtaction  and transformation
datExpr = as.data.frame(t(gene_expressed));
dim(datExpr) 
rownames(datExpr)
colnames(datExpr)

#Check for genes and samples with too many missing values:
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
# If not returns TRUE, need to remove the offending genes and samples from the data, see Tutorials.

#Cluster the samples to identify any obvious outliers, and plot the  tree
sampleTree = hclust(dist(datExpr), method = "average");
sizeGrWindow(12,9) # Adjust dimensions if the window is too large or too small.

pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

######Step 2###### Step-by-step network construction and module detection
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


# Plot the results:
pdf(file = "Scale-free topology.pdf", width = 12, height = 9)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()



# Check Scale-free topology
softPower = 9
k <- softConnectivity(datExpr,power=softPower) 
sizeGrWindow(10, 5)
pdf(file = "hist-Scale-free topology.pdf", width = 12, height = 9)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main=paste("Check Scale-free topology\n Power=",softPower))
dev.off()

# Calculate the adjacencies, using the soft thresholding power
adjacency = adjacency(datExpr, power = softPower); #this step could take some time. 

# Transform the adjacency into Topological Overlap Matrix (TOM), and calculate the corresponding dissimilarity
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency); #this step could take very long time (several hours). 
dissTOM = 1-TOM #this step also takes very long time. 

# Clustering using TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# Set the minimum module size #My settings
minModuleSize = 60;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Plot the module assignment under the gene dendrogram
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

pdf("Et_MEtree.pdf", width = 10, height = 4)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#Choose a height cut #Value of 0.25, corresponding to correlation of 0.75 
MEDissThres = 0.1 #My settings

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
# Plot the gene dendrogram again, with the original and merged module colors underneath 

pdf("Et_Dendrogram.pdf", width = 5, height = 3)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
table(mergedColors)
dim(table(mergedColors))

# Rename to moduleColors
moduleColors = mergedColors


# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Save module colors and labels for use in subsequent parts if neccessary
save(MEs, moduleLabels, moduleColors, geneTree, file = "Etef_AZ.RData")

# displaying module heatmap and the eigengene barplot
Colors = mergedColors
umc = unique(mergedColors)
lumc = length(umc)

for (i in c(1:lumc)){
  ME=MEs[, paste("ME",umc[i], sep="")]
  pdf_file_out= paste("dabbi-counts_var0.05-60-0.1-",umc[i],".pdf",sep="")
  pdf(file = pdf_file_out, wi = 9, he = 6)
  par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
  plotMat(t(scale(datExpr[,Colors==umc[i]])),nrgcols=30,rlabels=F,clabels=rownames(datExpr),rcols=umc[i], main=umc[i], cex.main=3)
  par(mar=c(5, 4.2, 0, 0.7))
  barplot(ME, col=umc[i], main="", cex.main=2,ylab="eigengene expression",xlab="array sample")
  dev.off()
}

library(pheatmap)
library(RColorBrewer)
for (i in c(1:lumc)){
  ME=MEs[, paste("ME",umc[i], sep="")]
  pdf_file_out= paste("dabbi-counts_var0.05-60-0.1-pheatmap-",umc[i],".pdf",sep="")
  pdf(file = pdf_file_out, wi = 9, he = 3)
  pheatmap(t(scale(datExpr[,Colors==umc[i]])), show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, 
           color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100), border_color = NA, main = umc[i], width = 8, height = 3)
  dev.off()
}

######Step 3###### Exporting network data to network visualization software

# Get outputs by colors
for (module in unique(mergedColors))
{
  
  # Select module probes
  probes = names(datExpr)
  inModule = is.finite(match(moduleColors, module));
  modProbes = probes[inModule];
  
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("dabbi-CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("dabbi-CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule]);
  
}



################################################################################
####################    detect hub genes in each module    ######################
################################################################################
library(tidyr)
library(dplyr)
library(writexl)


### Load DEseq2 annotation data
res_anno_dabbi <- read.table("/Users/yunqingyu/Dropbox/postdoc/data/fonio_tef_ERAC/RNAseq/tef_AZ_202012/DEseq2/res_anno_dabbi_DEseq2.txt", header = T, sep = "\t")

############## Detect hub genes using intramodular connectivitiy ###############
### Calculate intramodular connectivitiy
kIM <- intramodularConnectivity(adjMat = adjacency, colors = moduleColors, scaleByMax = TRUE) 
kIM_color <- mutate(kIM, module_color = moduleColors)%>%select("kTotal", "kWithin", "module_color")  # Add module color to each gene
head(kIM_color)

### Look at the known genes
# YAB
kIM_color["Et_10A_000629",]    #YAB2L1a, this is a hub gene (top 10% most connected gene).
kIM_color["Et_10B_002778",]
kIM_color["Et_4A_034697",]
kIM_color["Et_4B_038843",]


### Select modules of interest and find the genes with the highest connectivity
# module colors of interest
color_AZ <- c("grey60", "brown4", "darkred", "ivory", "lightyellow", "blue", "violet")
hub_AZ_kIM <- list()

for (c in color_AZ){
  kIM_module <- filter(kIM_color, module_color == c)%>%arrange(desc(kWithin))
  kIM_hub <- head(kIM_module, n = round(nrow(kIM_module)/10))%>%mutate(Et_geneID = row.names(.))  # Select the top 10% of genes
  
  # Add expression and annotation
  hub_AZ_kIM[[paste0("kIM_hub_anno_", c)]] <- left_join(kIM_hub, res_anno_dabbi, by = "Et_geneID")
  
}

# Check results
lapply(hub_AZ_kIM, head)

write_xlsx(hub_AZ_kIM, "dabbi_AZmodules_hub_kIM.xlsx")


############## Detect hub genes using kME ###############
#Signed eigengene-based connectivity of a gene in a module is defined as the correlation of the gene with the corresponding module eigengene. 
# This method is a lot faster compared with kIM
KME <- signedKME(datExpr = datExpr, datME = MEs, outputColumnName = "")
KME_color <- mutate(KME, module_color = moduleColors)  # Add module color to each gene

write.table(KME_color, "dabbi_KME_color.txt", sep = "\t")

hub_AZ_KME <- list()
for (c in color_AZ){
  KME_module <- KME_color%>%filter(module_color == c)
  n <- which(colnames(KME_module) == c)   # The column number corespond to the module color
  KME_module$kME_abs <- abs(KME_module[,n])  # 0 value of kME means no similarity to the eigengene. So values close to -1 and 1 are both significant.
  KME_module <- KME_module%>%
    mutate(Et_geneID = row.names(.))%>%
    arrange(desc(kME_abs))%>%
    select(c, kME_abs, module_color, Et_geneID)
  KME_hub <- head(KME_module, n = round(nrow(KME_module)/10))   # Select the top 10% of genes
  
  # Add expression and annotation
  hub_AZ_KME[[paste0("KME_hub_anno_", c)]] <- left_join(KME_hub, res_anno_dabbi, by = "Et_geneID")
  
}

# Check results
lapply(hub_AZ_KME, head)

write_xlsx(hub_AZ_KME, "dabbi_AZmodules_hub_KME.xlsx")
