###########################
#Libraries
###########################
library("scater", quietly = TRUE)
library("Seurat")
library("gridExtra")
library("ggpubr")
library("biomaRt")
library("scran")


library('KEGG.db')
library('ReactomePA')
library('pathview')
library('GSEABase')
library('org.Hs.eg.db')
library('annotate')
library('GOstats')
library('AnnotationDbi')
library("xlsx")
library("gtools")

library("minet")

set.seed(1234567)
options(stringsAsFactors = FALSE)

#DLLs controls and changed for >100 DLLs :)
Sys.getenv()
length(getLoadedDLLs()) 

###########################
#LOAD DATA
###########################
#Control & Infarto
#READ BOTH FILES
setwd("C:/Users/Analyzer/Desktop/seurat_Melanoma_timeSeries")
data_pre = read10xResults("C:/Users/Analyzer/Desktop/seurat_Melanoma_timeSeries/P1020Pre/GRCh38")
data_post = read10xResults("C:/Users/Analyzer/Desktop/seurat_Melanoma_timeSeries/P1020Post/GRCh38")
data_prog1 = read10xResults("C:/Users/Analyzer/Desktop/seurat_Melanoma_timeSeries/P1020Prog1/GRCh38")
data_prog2 = read10xResults("C:/Users/Analyzer/Desktop/seurat_Melanoma_timeSeries/P1020Prog2/GRCh38")
data_response = read10xResults("C:/Users/Analyzer/Desktop/seurat_Melanoma_timeSeries/P1020Response/GRCh38")

#file_list = c("data_pre","data_post","data_prog1","data_prog2","data_response")
#list_all = c(data_pre=data_pre,data_post=data_post,data_prog1=data_prog1,data_prog2=data_prog2, data_response=data_response)
file_list = c("data_pre","data_response")
list_all = c(data_pre=data_pre, data_response=data_response)

list_all_Seurat = c()
rm(list = file_list)
###########################
#START
###########################

####
#PRE-PROCESSING -> FILTERING EMPTY INSTANCES / ERCC / MT AND COMPUTING QC
####
#######################################################
###Rename cols to evade downstream errors
for (n in 1:length(file_list)){
  colnames(list_all[[n]]) = paste0(file_list[n], "_", colnames(list_all[[n]]))
}

#######################################################
#Remove Feature(genes,rna,etc) not expressed in any cell
for (n in 1:length(file_list)){
  keep_feature = rowSums(as.matrix(counts(list_all[[n]])) > 0) > 0
  list_all[[n]] = list_all[[n]][keep_feature, ]
}
rm(keep_feature)

#######################################################
#Change ENSEMBLE rownames by gene symbols
for (n in 1:length(file_list)){
  rownames(list_all[[n]]) = make.unique(as.character(list_all[[n]]@rowRanges@elementMetadata@listData$symbol), sep = "_")
}

#######################################################
#Define with features(genes,rna,etc) are the ERCC spike-ins and mitochondrial genes
for (n in 1:length(file_list)){
  isSpike(list_all[[n]], "ERCC") = grepl("^ERCC", rownames(list_all[[n]]))
  isSpike(list_all[[n]], "MT") = grepl("^MT", rownames(list_all[[n]]))
}

#######################################################
#Calculate de quality metrics
for (n in 1:length(file_list)){
  list_all[[n]]  = calculateQCMetrics(list_all[[n]],
                                      feature_controls = list(
                                        ERCC = isSpike(list_all[[n]], "ERCC"), 
                                        MT = isSpike(list_all[[n]], "MT")
                                        )
                                      )
  print (sprintf("QC added to sample -> %s", file_list[n]))
}

####
#PRE-PROCESSING -> FILTERING BAD QUALITY INSTANCES AND VARIABLES
####
#######################################################
### total_counts / total_features
lista_dataframes_counts = filter_instances_variables("total_counts", "Filtering_total_counts_plots_Histograms.pdf", "Filtering_total_counts_plots_Tables.pdf")
lista_dataframes_features = filter_instances_variables("total_features", "Filtering_total_features_plots_Histograms.pdf", "Filtering_total_features_plots_Tables.pdf")

#######################################################
###ERCC/ MT 
lista_dataframes_ERCC = filter_ERCC_MT("is_spike_MT")
lista_dataframes_MT = filter_ERCC_MT("is_spike_ERCC")

#######################################################
###Apply all previous created filters and filter umi based on those filters
for (n in 1:length(file_list)){
  list_all[[n]]$use = (unlist(lista_dataframes_features[n]) & unlist(lista_dataframes_counts[n]))
  print (table(list_all[[n]]$use))
}

#######################################################
#Gene filtering after cell filtering !!!
#Low quality genes + MT + ERCC
lista_filter_genes = list()
for (n in 1:length(file_list)){
  #
  filter_genes = apply(
    counts(list_all[[n]][ , colData(list_all[[n]])$use]), 
    1, 
    function(x) length(x[x > 1]) >= 2
  )
  lista_filter_genes = append(lista_filter_genes, list(filter_genes))
  #Filter rows (genes)
  rowData(list_all[[n]])$use = (unlist(lista_filter_genes[n]) & !unlist(lista_dataframes_MT[n]) & !unlist(lista_dataframes_ERCC[n]))
  print(table(rowData(list_all[[n]])$use))
}

rm(lista_dataframes_counts, lista_dataframes_features, lista_dataframes_ERCC, lista_dataframes_MT, lista_filter_genes, filter_genes)
####
#PROCESSING -> COMBINE DATASETS SEURAT CCA 
####
#######################################################
#Creating Seurat objects
for (n in 1:length(file_list)){
  list_all_Seurat = append(list_all_Seurat, CreateSeuratObject(
    raw.data = counts(list_all[[n]][rowData(list_all[[n]])$use , colData(list_all[[n]])$use])
    )) 
}
rm(list_all)
#######################################################
#Visualize data for filtering propouse
seurat_data_plot("Seurat_data_filtered.pdf")
rm(blankPlot)

#######################################################
#Normalization and log2 transformation
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = NormalizeData(
    object = list_all_Seurat[[n]], 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
  )
}

#######################################################
#Scale data
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = ScaleData(object = list_all_Seurat[[n]], 
                                   vars.to.regress = c("nUMI")
  )
}

#######################################################
#Highly variable genes obtain for downstream analysis
pdf(file="Seurat_data_filtered_MOST_VARIABLE_GENES_SELECTED.pdf")
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = FindVariableGenes(
    object = list_all_Seurat[[n]],
    mean.function = ExpMean, 
    dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, 
    x.high.cutoff = 3, 
    y.cutoff = 0.5
  )
}
dev.off()

#######################################################
#Number of selected genes
lista_var_genes = c()
for (n in 1:length(file_list)){
  lista_var_genes = append(lista_var_genes ,length(x = list_all_Seurat[[n]]@var.genes))
}

#######################################################
#Use the top min(lista_var_genes) most variable genes in each dataset for alignment 
lista_var_genes_all = list()
for (n in 1:length(file_list)){
  lista_var_genes_all = append(lista_var_genes_all, list(rownames(head(list_all_Seurat[[n]]@hvg.info, min(lista_var_genes)))))
}
hvg_union = as.data.frame(table(unlist(lista_var_genes_all)))
hvg_union = hvg_union[hvg_union$Freq > 1,]
genes.use  = as.vector(unlist(hvg_union['Var1']))

#######################################################
#Set the protocol and name for easy identification
for (n in 1:length(file_list)){
  list_all_Seurat[[n]]@meta.data[, "protocol"] = file_list[n]
}

rm(lista_var_genes, lista_var_genes_all)
#######################################################
#Run CCA
gc()
if (length(file_list) >2){
  for (n in 1:length(file_list)) {
    genes.use = genes.use[genes.use %in% rownames(list_all_Seurat[[n]]@scale.data)]
  }
  multiCCA = RunMultiCCA(list_all_Seurat, genes.use = genes.use, num.ccs = 2)
}else{
  multiCCA <- RunCCA(list_all_Seurat[[1]], list_all_Seurat[[2]], genes.use = genes.use)
}

#######################################################
#Scale data based in batches
multiCCA = ScaleData(object = multiCCA, vars.to.regress = c("protocol", "nUMI"))

#######################################################
#PC selection
multiCCA = RunPCA(
  object = multiCCA, 
  pc.genes = genes.use, 
  do.print = FALSE
)

pdf(file="PC_selection.pdf")
PCElbowPlot(multiCCA, num.pc = 20)
dev.off()

#######################################################
#Run CCA
pdf(file="MultiCCA_cc1_cc2_Before_Align.pdf", width=12, height=6)
#plot CC1 versus CC2 and look at a violin plot
p1 = DimPlot(multiCCA, reduction.use = "cca", group.by = "protocol", pt.size = 0.5,
              do.return = T)
p2 = VlnPlot(multiCCA, features.plot = "CC1", group.by = "protocol", do.return = T)
plot_grid(p1, p2)
dev.off()

#######################################################
#Before align search cells whose profile cannot be well_explained by low_dimensional CCA compared to low-dimensional PCA
multiCCA = CalcVarExpRatio(multiCCA, reduction.type = "pca", grouping.var = "protocol", dims.use = 1:2)
multiCCA.all = multiCCA

#######################################################
#Those are in both samples
multiCCA = SubsetData(multiCCA, subset.name = "var.ratio.pca", accept.low = 0.5)

#Discard ones becouse are specific of one sample
multiCCA.discard = SubsetData(multiCCA.all, subset.name = "var.ratio.pca", accept.high = 0.5)

#Align subspaces for obtaining CCA_ALIGNED
multiCCA = AlignSubspace(multiCCA, reduction.type = "cca", grouping.var = "protocol", dims.align = 1:2)

pdf(file="MultiCCA_cc1_cc2_After_Align.pdf", width=12, height=6)
p1 = DimPlot(multiCCA, reduction.use = "cca", group.by = "protocol", pt.size = 0.5,
              do.return = T)
p2 = VlnPlot(multiCCA, features.plot = "ACC1", group.by = "protocol", do.return = T)
plot_grid(p1, p2)
dev.off()

#######################################################
#TSNE plot generation
multiCCA = RunTSNE(multiCCA, reduction.use = "cca.aligned", dims.use = 1:2, do.fast = T, perplexity = 100, check_duplicates = FALSE)
multiCCA = FindClusters(multiCCA, reduction.type = "cca.aligned", dims.use = 1:2, save.SNN = T)

pdf(file="TSNE.pdf", width=12, height=6)
p1 = TSNEPlot(multiCCA, group.by = "protocol", do.return = T, pt.size = 0.5)
p2 = TSNEPlot(multiCCA, do.return = T, pt.size = 0.5, do.label = T)
plot_grid(p1, p2)
dev.off()

rownames(multiCCA@dr$tsne@cell.embeddings)
ls(pattern="data_pre_")

#, cells.use=rowns
#rowns = rownames(multiCCA@dr$tsne@cell.embeddings)[grepl("data_pre_", rownames(multiCCA@dr$tsne@cell.embeddings)) | grepl("data_response_", rownames(multiCCA@dr$tsne@cell.embeddings))]

####
#PROCESSING -> DE ANALYSIS
####
#######################################################
#Normalized and log transformed data for DE analysis
data_matrix=data.frame(as.matrix(multiCCA@data))

cluster_ident=data.frame(multiCCA@ident)
cluster_ident['Cells'] = rownames(cluster_ident)
cluster_ident = cluster_ident[with(cluster_ident, order(pbmc2.ident)), ]

#######################################################
#FINDMARKERS all clusters agains rest of cells
multiCCA.markers = FindAllMarkers(object = multiCCA, only.pos = TRUE, min.pct = 0.25)
print(x = head(x = multiCCA.markers, n = 5))
multiCCA.markers = (subset(multiCCA.markers,(multiCCA.markers['p_val_adj'] <= 0.01 )))

#######################################################
#Save
save.image("Data_SAVE.RData")

################################################################################################################
################################################################################################################
###
#FUNCTIONAL ANALYSIS -> FOR EACH DE CLUSTER
####
#######################################################
gsc.GO = create_GO_universe(org.Hs.egGO, "Homo sapiens")
gsc.KEGG = create_KEGG_universe(org.Hs.egPATH, "Homo sapiens")

##############
#GOstats
##############
#universe & geneList
#######################################################
universe = Lkeys(org.Hs.egGO)

####
#GO BP / CC / MF
#######################################################
sectors = c('BP', 'CC' , 'MF')
for (n in 0: (max(as.integer(multiCCA.markers$cluster))-1)) {
  #GO
  cluster_markers = multiCCA.markers$gene[multiCCA.markers$cluster == n]
  annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', cluster_markers)
  entrezgene = annotation$entrezgene[match(cluster_markers, annotation$external_gene_name)]
  annotations = as.character(entrezgene)
  for (i in 1:length(sectors)){
    tmp = gsea_HyperG_GO("GOstats", gsc.GO, annotations, sectors[i], 1)
    assign(paste0("GOresults", sectors[i], "_", "cluster", "_", n), tmp)
  }
  #KEGG
  tmp = gsea_HyperG_KEGG("KEGG", gsc.KEGG, annotations, 1)
  assign(paste0("KEGGresults", "_", "cluster", "_", n), tmp)
}

sectors = c("GOresultsBP_cluster_", "GOresultsCC_cluster_", "GOresultsMF_cluster_")
for (i in 1:length(sectors)){
  list = ls(pattern=sectors[i])
  assign(paste0(sectors[i], "list"), mixedsort(list))
}
KEGGresults_cluster_list = mixedsort(ls(pattern="KEGGresults_cluster_"))

#######################################################
#Output -> one cluster for each .xlsx // content (BP/CC/MF/KEGG)
count = -1
for (i in 1:(length(GOresultsBP_cluster_list)-1)){
  count = count + 1
  ###
  book=createWorkbook(type="xlsx")
  
  sheet1 = createSheet(book, sheetName = "GOresultsBP")
  sheet2 = createSheet(book, sheetName = "GOresultsCC")
  sheet3 = createSheet(book, sheetName = "GOresultsMF")
  sheet4 = createSheet(book, sheetName = "KEEGresult")
  
  addDataFrame(get(GOresultsBP_cluster_list[i]), sheet = sheet1, startRow = 1, startColumn = 1)
  addDataFrame(get(GOresultsCC_cluster_list[i]), sheet = sheet2, startRow = 1, startColumn = 1)
  addDataFrame(get(GOresultsMF_cluster_list[i]), sheet = sheet3, startRow = 1, startColumn = 1)
  addDataFrame(get(KEGGresults_cluster_list[i]), sheet = sheet4, startRow = 1, startColumn = 1)
  
  saveWorkbook(book, paste0("Cluster_", count, ".xlsx"))
}







#clusters 4
#genes_to_be_used = subset(multiCCA.markers,(multiCCA.markers['cluster'] == c(4)))
#genes_to_be_used = as.character(unlist(genes_to_be_used['Genes']))
#full_table_marker_DE = cc[rownames(cc) %in% genes_to_be_used,]
#write.table(full_table_marker_DE,"full_table_marker_DE.csv", sep=",")

################################################################################################################
################################################################################################################
###
#NETWORK ANALYSIS -> FOR EACH DE CLUSTER
####
#######################################################
#Filter for selecting rows and column from count matrix
aracne_use_cluster = create_selection(0,0)

#######################################################
#MINET
aracne_res_cluster = minet(t(aracne_use_cluster), method="aracne",
                             estimator="spearman", disc="equalfreq",
                             nbins=sqrt(NROW(t(aracne_use_cluster))))

#Clean matrix for Cytoscape
indices = data.frame(which(aracne_res_cluster!=0,arr.ind = T))
indices = indices[order(indices$row),]

tabla = matrix(, nrow = length(indices[,1]), ncol = 3)

for (i in 1:(length(indices[,1]))) {
  tabla[i,1] = rownames(aracne_res_cluster)[indices['row'][i,1]]
  tabla[i,2] = colnames(aracne_res_cluster)[indices['col'][i,1]]
  tabla[i,3] = aracne_res_cluster[(indices['row'][i,1]),(indices['col'][i,1])]
}

aracne_res_cluster_final = tabla

aracne_res_cluster_final_filt = data.frame(aracne_res_cluster_final)[!duplicated(data.frame(aracne_res_cluster_final)[3]),]

plot(density(as.numeric(tabla[,3])))
abline(v = c(0.13), col = "red")

aracne_res_cluster_final_filt = (subset(aracne_res_cluster_final_filt,(aracne_res_cluster_final_filt['X3'] >= 0.13 )))

write.table(aracne_res_cluster_final_filt,"aracne_res_cluster_final.txt", sep="\t",
            quote=F,row.names = F,col.names = F)

#########
#########
#Compare 2 clusters to have de same edges to go Cytoscape
both_4 = read.table("stroke_4.txt")
both_4_to_8 = read.table("stroke_4_to_8.txt")

tabla = matrix(, nrow = 200, ncol = 4)

for (i in 1:(length(both_4[,1]))){
  print ("no")
  for (j in 1:(length(both_4_to_8[,1]))){ 
    print ("no")
    if (both_4[i,1] == both_4_to_8[j,1]){
      print ("si")
      if (both_4[i,2] == both_4_to_8[j,2]){
        print ("si")
        tabla[i,1] = both_4[i,1]
        tabla[i,2] = both_4[i,2]
        tabla[i,3] = both_4[i,3]
        tabla[i,4] = both_4_to_8[i,3]
      }
    }
    if (both_4[i,1] == both_4_to_8[j,2]){ 
      print ("si")
      if (both_4[i,2] == both_4_to_8[j,1]){
        tabla[i,1] = both_4[i,1]
        tabla[i,2] = both_4[i,2]
        tabla[i,3] = both_4[i,3]
        tabla[i,4] = both_4_to_8[i,3] 
        print ("si")
      }
    }
  } 
}

tabla = tabla[complete.cases(tabla),]

write.table(tabla,"tabla_no_especificos_stroke.txt", sep="\t",
            quote=F,row.names = F,col.names = F)







