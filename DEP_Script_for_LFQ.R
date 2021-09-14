library(DEP)
library(dplyr)

#load your data
data<- read.delim("enter proteingroups.txt here")
dim(data)
colnames(data)

#filter data based on valid values
data <- filter(data, Reverse != "+", Potential.contaminant != "+", Only.identified.by.site!= "+")
dim(data)

#are there duplicated gene names? if so how many?
data$Gene.names %>% duplicated() %>% any()

#Replace the duplicated gene names by their proteinID 
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
data$name %>% duplicated() %>% any()

#Load experimental design containing sample names and replicates and work with the LFQ intensities
experimental_design <- read.delim('Load your experimental design here')

LFQ_columns <- grep("LFQ.", colnames(data_unique))
data_se <- make_se(data_unique, LFQ_columns, experimental_design)
data_se

#Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 1)
plot_numbers(data_filt)

#normalize data and log2 transform
data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm)

#Inspect missing values with wrProteo
library(knitr)
library(wrMisc)
library(wrGraph)
library(wrProteo)

dataMQ <- readMaxQuantFile('enter proteingroups.txt here', specPref=NULL, normalizeMeth="median")
dim(dataMQ$quant)
summary(dataMQ$quant[,1:27]) 

#write your conditions (each= 3 for triplicates, change if more replicates are used)
conditions<-c('control_D1', '50mg/ml_D1','100mg/ml_D1','control_D2', '50mg/ml_D2', '100mg/ml_D2', 'control_D5', '50mg/ml_D5','100mg/ml_D5')
grp9<- paste0(rep(conditions, each=3))

matrixNAinspect(dataMQ$quant, gr=grp9, tit="Inspection of NA values") 
dataMQimp <- matrixNAneighbourImpute(dataMQ$quant, gr=grp9, imputMethod = 'informed', tit="NA-replacement") 

matrixNAinspect(dataMQ$quant, gr=grp9, tit="Inspection of NA values") 
dataMQimp1 <- matrixNAneighbourImpute(dataMQ$quant, gr=grp9, tit="NA-replacement") 

plotPCAw(dataMQimp, sampleGrp=grp9, tit="PCA on MaxQuant (NAs imputed)", rowTyName="proteins", useSymb2=0)

#Plot missing values
plot_missval(data_filt)
plot_detect(data_filt)

data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)
plot_imputation(data_norm, data_imp)
plot_imputation(data_norm, data_imp_knn)

#Diff. enrichment analysis
data_diff <- test_diff(data_imp, type = "all")
dep <- add_rejections(data_diff, alpha = 0.01, lfc = log2(1.5))

#PCA
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Spectral")

#Heatmap
plot_heatmap(dep, type = "centered", kmeans = TRUE, k = 6, col_limit = 4, show_row_names = FALSE, indicate = c("condition", "replicate"))

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, k = 6, col_limit = 10, show_row_names = FALSE)

#Volcano plot comp 1, the conditions given in the experimental design can be compared using the contrast argument
volcano_plot<-plot_volcano(dep, contrast = "condition1_vs_condition2", label_size = 2, add_names = TRUE)


#Change the colors of the volcano plot
library(ggplot2)
library(ggsci)
volcano_plot + scale_color_lancet()

#get results
#the p.val and p.adj columns contain the raw and adjusted p values, respectively
#the ratio columns contain the average log2 fold changes
#The significant columns indicate whether the protein is differentially enriched/expressed, as defined by the chosen cutoffs.
#The centered columns contain the average log2 fold changes scaled by protein-wise centering.
data_results <- get_results(dep)
data_results %>% filter(significant) %>% nrow()
significant_proteins <- data_results[data_results$significant,]
dim(data_results)
dim(significant_proteins)

#Select up and down regulated proteins from a comparison and GO analysis
up<-significant_proteins[significant_proteins$condition1_vs_condition2_ratio>0 ,c('ID','condition1_vs_condtion2_ratio', 'codition1_vs_condition2_significant')]
up_reg<-up[up$condition1_vs_condition2_significant==TRUE, c('ID')]

dn<-significant_proteins[significant_proteins$condition1_vs_condition2_ratio<0 ,c('ID','condition1_vs_condtion2_ratio', 'codition1_vs_condition2_significant')]
dn_reg<-up[dn$condition1_vs_condition2_significant==TRUE, c('ID')]


#Convert Uniprot to Entrez or KEGG and viseversa
#Download the correct ordganism-database (here yeast)
library(clusterProfiler)
Entrez<-bitr(significant_proteins$ID, 'UNIPROT', 'ENTREZID', org.Sc.sgd.db, drop=TRUE)
kegg<-bitr_kegg(significant_proteins$ID, 'uniprot', 'kegg', 'sce', drop = TRUE)
  
enrichBP<-enrichGO(Entrez$ENTREZID, org.Sc.sgd.db, ont = 'BP', pvalueCutoff = 0.05, pAdjustMethod = 'BH', qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500)
enrichMF<-enrichGO(Entrez$ENTREZID, org.Sc.sgd.db, ont = 'MF', pvalueCutoff = 0.05, pAdjustMethod = 'BH', qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500)
enrichCC<-enrichGO(Entrez$ENTREZID, org.Sc.sgd.db, ont = 'BP', pvalueCutoff = 0.05, pAdjustMethod = 'BH', qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500)
enrichkegg<-enrichKEGG(kegg$kegg, organism = 'sce', keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = 'BH', qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500)

plot_BP_ratio<-clusterProfiler::dotplot(enrichBP, showCategory=30)
plot_BP_count<-clusterProfiler::dotplot(enrichBP, x='count', showCategory=30)

plot_MF_ratio<-clusterProfiler::dotplot(enrichMF, showCategory=30)
plot_MF_count<-clusterProfiler::dotplot(enrichMF, x='count', showCategory=30)

plot_CC_ratio<-clusterProfiler::dotplot(enrichCC, showCategory=30)
plot_CC_count<-clusterProfiler::dotplot(enrichCC, x='count', showCategory=30)

plot_kegg_ratio<-clusterProfiler::dotplot(enrichkegg, showCategory=30)
plot_kegg_count<-clusterProfiler::dotplot(enrichkegg, x='count', showCategory=30)

#Change colors of the dotplots containing GO annotations
plot_BP_ratio+theme_grey(base_size = 11)+scale_color_gsea()+ggtitle('Biological process')
plot_BP_count+theme_grey(base_size = 11)+scale_color_gsea()+ggtitle('Biological process')

plot_MF_ratio+theme_grey(base_size = 11)+scale_color_gsea()+ggtitle('Molecular function')
plot_MF_count+theme_grey(base_size = 11)+scale_color_gsea()+ggtitle('Molecular function')

plot_CC_ratio+theme_grey(base_size = 11)+scale_color_gsea()+ggtitle('Cellular component')
plot_CC_count+theme_grey(base_size = 11)+scale_color_gsea()+ggtitle('Cellular component')

plot_kegg_ratio+theme_grey(base_size = 11)+scale_color_gsea()+ggtitle('KEGG pathway')
plot_kegg_count+theme_grey(base_size = 11)+scale_color_gsea()+ggtitle('KEGG pathway')

#EXTRA (needs more work)
library(fgsea)
reactome<-reactomePathways(entez$entrezID)

#Get Gene ontolgy Information 
library(UniprotR)
library(writexl)
protein_functions<-GetProteinFunction(significant_proteins$ID, directorypath = NULL)
Enrichment.BP(significant_proteins$ID,OS="scerevisiae",p_value=0.05,directorypath=NULL)
Enrichment.KEGG(significant_proteins$ID,OS="scerevisiae",p_value=0.05,directorypath=NULL)
Enrichment.CC(significant_proteins$ID,OS="scerevisiae",p_value=0.05,directorypath=NULL)
Enrichment.MF(significant_proteins$ID,OS="scerevisiae",p_value=0.05,directorypath=NULL)
GO_info<-GetProteinGOInfo(significant_proteins$ID)
write_xlsx(protein_functions,"/Users/pedarf/Desktop/Julie/protein_functions.xls")
write_xlsx(GO_info,"/Users/pedarf/Desktop/Julie/GO_info.xls")

#Plot Gene ontology Information
Plot.GOMolecular(GO_info, Top = 20)
Plot.GOSubCellular(GO_info)
PlotGoInfo(GO_info)
#Enrichment analysis using KEGG, Reactome of protein list
Pathway.Enr(significant_proteins$ID, OS="scerevisiae",p_value=0.05,directorypath=NULL)

#Get diseases associated with protein list
Pathology <- GetPathology_Biotech(GO_info)
Diseases <- Get.diseases(Pathology)

#Protein- Protein interaction using STRING
GetproteinNetwork(GO_info)
# Generate a wide data.frame
df_wide <- get_df_wide(dep)
# Generate a long data.frame
df_long <- get_df_long(dep)

#Clean data frame with protein names, p.adj, log2 FC, Significance (good when comparing two conditions only)
results_table<-select(data_results,'name',"ID","significant")
results_table

# Save analyzed data
save(data_se, data_norm, data_imp, data_diff, dep, file = "data.RData")
# These data can be loaded in future R sessions using this command
load("data.RData")

