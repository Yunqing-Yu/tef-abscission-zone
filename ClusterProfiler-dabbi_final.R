#Install SvOrgDb from local for the 1st run 
install.packages("/Users/yunqingyu/Dropbox/postdoc/data/fonio_tef_ERAC/RNAseq/tef_AZ_202012/clusterProfiler/org.Etef.eg.db",repos = NULL, type="source")

#Load clusterProfiler for GO analysis
library(clusterProfiler)
library(org.Etef.eg.db)
library(GOSemSim)
library(enrichplot)
library(writexl)

######### Prepare input for GO analysis #########
# Load gene list from GO_input
files <- list.files("/Users/yunqingyu/Dropbox/postdoc/data/fonio_tef_ERAC/RNAseq/tef_AZ_202012/clusterProfiler/GO_input/", pattern = ".txt")
files
gene_list <- list()
for (i in files){
  f <- read.table(paste0("/Users/yunqingyu/Dropbox/postdoc/data/fonio_tef_ERAC/RNAseq/tef_AZ_202012/clusterProfiler/GO_input/", i), header = F)[,1]
  module <- gsub("(.*)_Aup_(.*).txt", "\\1_\\2", i)
  gene_list[[module]] <- f
}
names(gene_list)
lapply(gene_list, length)

# Make GOsemsimDATA object used for enrichment analysis
SemData_tef <- godata(OrgDb = "org.Etef.eg.db", keytype = "GID", ont = "BP", computeIC = TRUE)

######### test parameters for enrichGO #########
# Use s1_grey60 to test method
ego_s1 <- enrichGO(gene = gene_list[["s1_grey60"]], 
                OrgDb = "org.Etef.eg.db", keyType="GID", ont="BP", 
                pAdjustMethod = "BH", qvalueCutoff = 0.05, 
                minGSSize = 10, maxGSSize = 1000)
ego_s1_Rel <- simplify(ego_s1, cutoff = 0.6, measure = "Rel", semData = SemData_tef)
ego_s1_Rel@result$Description

ego_s1_Jiang <- simplify(ego_s1, cutoff = 0.6, measure = "Jiang", semData = SemData_tef)
ego_s1_Jiang@result$Description

ego_s1_Wang <- simplify(ego_s1, cutoff = 0.6, measure = "Wang", semData = SemData_tef)
ego_s1_Wang@result$Description

ego_s1_Lin <- simplify(ego_s1, cutoff = 0.6, measure = "Lin", semData = SemData_tef)
ego_s1_Lin@result$Description

ego_s1_Resnik <- simplify(ego_s1, cutoff = 0.6, measure = "Resnik", semData = SemData_tef)
ego_s1_Resnik@result$Description
# These methods are all similar but not exactly the same. I decide to use "Rel".



######### Perform enrichGO #########
# Use enricher for visualizing a single gene list
ego_list <- list()  #store enrichGO output
ego_list_Rel <- list()   #store enrichGO output after simplify with "Rel" method
GO_list <- list()    #store result table
GO_list_Rel <- list()   #store result table after simplify with "Rel" method
for (module in names(gene_list)){
  GENE_ID <- gene_list[[module]]
  ego <- enrichGO(gene = GENE_ID, 
                  OrgDb = "org.Etef.eg.db", keyType="GID", ont="BP", 
                  pAdjustMethod = "BH", qvalueCutoff = 0.05, 
                  minGSSize = 10, maxGSSize = 1000)
  ego_list[[module]] <- ego
  GO_list[[module]] <- ego@result[ego@result$p.adjust < 0.05,]
  ego_Rel <- simplify(ego, cutoff = 0.6, measure = "Rel", semData = SemData_tef)
  ego_list_Rel[[module]] <- ego_Rel
  GO_list_Rel[[module]] <- ego_Rel@result
}

#Checking the number of GO terms before and after simplify and also check the identical terms
names(GO_list)
dim(GO_list[["s1_grey60"]])
dim(GO_list_Rel[["s1_grey60"]])
GO_dup_s1_grey60 <- GO_list_Rel[["s1_grey60"]]$ID[duplicated(GO_list_Rel[["s1_grey60"]]$geneID)]
GO_dup_s1_grey60 <- c(GO_dup_s1_grey60[c(1:7, 9:12)], "GO:0052158")
GO_list_Rel[["s1_grey60"]]$Description[GO_list_Rel[["s1_grey60"]]$ID %in% GO_dup_s1_grey60]

dim(GO_list[["s2_brown4"]])
dim(GO_list_Rel[["s2_brown4"]])
GO_dup_s2_brown4 <- GO_list_Rel[["s2_brown4"]]$ID[duplicated(GO_list_Rel[["s2_brown4"]]$geneID)]
GO_list_Rel[["s2_brown4"]]$Description[GO_list_Rel[["s2_brown4"]]$ID %in% GO_dup_s2_brown4]

dim(GO_list[["s2_darkred"]])
dim(GO_list_Rel[["s2_darkred"]])
GO_darkred_diff <- GO_list[["s2_darkred"]]$ID[!(GO_list[["s2_darkred"]]$ID %in% GO_list_Rel[["s2_darkred"]]$ID)]
length(GO_darkred_diff)
GO_dup_s2_darkred1 <- GO_list_Rel[["s2_darkred"]]$ID[duplicated(GO_list_Rel[["s2_darkred"]]$geneID)]
GO_list_Rel[["s2_darkred"]]$Description[GO_list_Rel[["s2_darkred"]]$ID %in% GO_dup_s2_darkred1]
GO_dup_s2_darkred <- c(GO_dup_s2_darkred1, GO_darkred_diff[2:185], "GO:0019288")
length(GO_dup_s2_darkred)

dim(GO_list[["s2_ivory"]])
dim(GO_list_Rel[["s2_ivory"]])
GO_dup_s2_ivory <- GO_list_Rel[["s2_ivory"]]$ID[duplicated(GO_list_Rel[["s2_ivory"]]$geneID)]
GO_dup_s2_ivory <- c(GO_dup_s2_ivory[2:6], "GO:0019288")
GO_list_Rel[["s2_ivory"]]$Description[GO_list_Rel[["s2_ivory"]]$ID %in% GO_dup_s2_ivory]

dim(GO_list[["s2_lightyellow"]])
dim(GO_list_Rel[["s2_lightyellow"]])
GO_dup_s2_lightyellow <- GO_list_Rel[["s2_lightyellow"]]$ID[duplicated(GO_list_Rel[["s2_lightyellow"]]$geneID)]
GO_dup_s2_lightyellow <- c("GO:0010116", GO_dup_s2_lightyellow[2:6])
GO_list_Rel[["s2_lightyellow"]]$Description[GO_list_Rel[["s2_lightyellow"]]$ID %in% GO_dup_s2_lightyellow]

dim(GO_list[["s3_blue"]])
dim(GO_list_Rel[["s3_blue"]])
GO_dup_s3_blue <- GO_list_Rel[["s3_blue"]]$ID[duplicated(GO_list_Rel[["s3_blue"]]$geneID)]
GO_list_Rel[["s3_blue"]]$Description[GO_list_Rel[["s3_blue"]]$ID %in% GO_dup_s3_blue]

dim(GO_list[["s3_violet"]])
dim(GO_list_Rel[["s3_violet"]])
GO_dup_s3_violet <- GO_list_Rel[["s3_violet"]]$ID[duplicated(GO_list_Rel[["s3_violet"]]$geneID)]
GO_list_Rel[["s3_violet"]]$Description[GO_list_Rel[["s3_violet"]]$ID %in% GO_dup_s3_violet]


# Remove redundant terms and manually selected terms
ego_Rel_unique <- list()
ego_Rel_unique[["s1_grey60"]] <- dropGO(ego_list_Rel[["s1_grey60"]], term = GO_dup_s1_grey60)
ego_Rel_unique[["s2_brown4"]] <- ego_list_Rel[["s2_brown4"]]
ego_Rel_unique[["s2_darkred"]] <- dropGO(ego_list[["s2_darkred"]], term = c(GO_dup_s2_darkred, "GO:0051147", "GO:0043489", "GO:1902456"))
ego_Rel_unique[["s2_ivory"]] <- dropGO(ego_list_Rel[["s2_ivory"]], term = GO_dup_s2_ivory)
ego_Rel_unique[["s2_lightyellow"]] <- dropGO(ego_list_Rel[["s2_lightyellow"]], term = c(GO_dup_s2_lightyellow, "GO:0000741", "GO:0043100"))
ego_Rel_unique[["s3_blue"]] <- dropGO(ego_list_Rel[["s3_blue"]], term = c(GO_dup_s3_blue, "GO:0120255", "GO:0120251", "GO:0120252", "GO:0007610", "GO:0120254", "GO:0051216"))
ego_Rel_unique[["s3_violet"]] <- dropGO(ego_list_Rel[["s3_violet"]], term = c(GO_dup_s3_violet, "GO:0015919", "GO:0009804", "GO:0002213", "GO:0120255", "GO:0006591"))

GO_Rel_unique <- list(s1_grey60 = ego_Rel_unique[["s1_grey60"]]@result, 
                      s2_brown4 = ego_Rel_unique[["s2_brown4"]]@result, 
                      s2_darkred = ego_Rel_unique[["s2_darkred"]]@result,
                      s2_ivory = ego_Rel_unique[["s2_ivory"]]@result, 
                      s2_lightyellow = ego_Rel_unique[["s2_lightyellow"]]@result, 
                      s3_blue = ego_Rel_unique[["s3_blue"]]@result, 
                      s3_violet = ego_Rel_unique[["s3_violet"]]@result)

write_xlsx(GO_list, "GO_clusterProfiler.xlsx")
write_xlsx(GO_Rel_unique, "GO_clusterProfiler_Rel0.6.xlsx")

# Check length after removing redundant terms
lapply(ego_Rel_unique, function(x){dim(x@result)})

# Make dotplot

dot_dabbi <- lapply(names(ego_Rel_unique), function(x){dotplot(ego_Rel_unique[[x]], x = "count", showCategory=100, title = x, font.size =11, label_format = 30)})
pdf("dotplot/dotplot_dabbi_s1_grey60.pdf", width = 6.3, height = 8)
dot_dabbi[[1]]
dev.off()
pdf("dotplot/dotplot_dabbi_s2_brown4.pdf", width = 5.5, height = 1.5)
dot_dabbi[[2]]
dev.off()
pdf("dotplot/dotplot_dabbi_s2_darkred.pdf", width = 6.5, height = 8)
dot_dabbi[[3]]
dev.off()
pdf("dotplot/dotplot_dabbi_s2_ivory.pdf", width = 5.5, height = 5)
dot_dabbi[[4]]
dev.off()
pdf("dotplot/dotplot_dabbi_s2_lightyellow.pdf", width = 6.3, height = 6)
dot_dabbi[[5]]
dev.off()
pdf("dotplot/dotplot_dabbi_s3_blue.pdf", width = 8.8, height = 11)
dot_dabbi[[6]]
dev.off()
pdf("dotplot/dotplot_dabbi_s3_violet.pdf", width = 8, height = 11)
dot_dabbi[[7]]
dev.off()




######################## Analyze different modules together ####################3#####3
#Use compareCluster for comparing between gene lists
Cluster_AZ <- compareCluster(gene_list, fun = "enrichGO",
                             OrgDb = "org.Etef.eg.db", keyType="GID", ont="BP", 
                             pAdjustMethod = "BH", qvalueCutoff = 0.05, 
                             minGSSize = 10, maxGSSize = 1000)
# Remove the redundant term with long name "GO:0019288"
Cluster_AZ1 <- dropGO(Cluster_AZ, term = "GO:0019288")

#Reduce GO redundancy
Cluster_AZ_s <- simplify(Cluster_AZ1, cutoff=0.6, by="p.adjust", measure = "Rel", semData = SemData_tef)


# Drop redundant GO terms and terms I manually selected and removed
GO_dup <- c(GO_dup_s1_grey60, GO_dup_s2_brown4, GO_dup_s2_darkred1, GO_dup_s2_ivory, GO_dup_s2_lightyellow, GO_dup_s3_blue, GO_dup_s3_violet, 
             "GO:0051147", "GO:0043489", "GO:1902456", "GO:0000741", "GO:0043100", "GO:0120255", "GO:0120251", "GO:0120252", "GO:0007610", 
             "GO:0120254", "GO:0051216", "GO:0015919", "GO:0009804", "GO:0002213", "GO:0120255", "GO:0006591")  

Cluster_AZ_s_unique <- dropGO(Cluster_AZ_s, term = GO_dup) # Essentially, this creates the same results as the individual analysis. 

pdf("dotplot/dotplot_dabbi_AZ_GO_top10.pdf", width = 8, height = 12)
dotplot(Cluster_AZ_s_unique, showCategory = 10, font.size =11, label_format = 30)
dev.off()






