#Importing Packages
  install.packages('patchwork')
  install.packages('Seurat')
  install.packages('dplyr')
  install.packages('ggplot2')

#Initializing Packages
  library(patchwork)
  library(Seurat)
  library(dplyr)
  library(ggplot2)

#Setting Working Directory
  setwd("C:/Raghav/Coding/R/UCI Research/Pathogenic Activity Project")

#Importing Dataset
  imp.data = Read10X(data.dir = "C:/Raghav/Coding/R/UCI Research/Pathogenic Activity Project/data/filtered_feature_bc_matrix-20210630T115618Z-001/filtered_feature_bc_matrix")
  imp.data

#Creating matr Seurat Object
  matr1 = CreateSeuratObject(counts = imp.data, project = "filtered_feature_bc_matrix", min.cells = 3, min.features = 4)
  matr1[["percent.mt"]] <- PercentageFeatureSet(matr1, pattern = "^mt-")

#Visualize QC metric as a violin plot
  VlnPlot(matr1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

  matr1 <- subset(matr1, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
  matr1

  matr1 <- NormalizeData(matr1, normalization.method = "LogNormalize", scale.factor = 10000)
  matr1 <- FindVariableFeatures(matr1, selection.method = "vst", nfeatures = 2000)

  matr = matr1[-(36)]
  matr


# Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(matr), 10)

# plot variable features with and without labels
  plot1 <- VariableFeaturePlot(matr)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2

  all.genes <- rownames(matr)
  matr <- ScaleData(matr, features = all.genes)

#Calculate PCA Score
  matr <- RunPCA(matr)

# Examine and visualize PCA results a few different ways
  print(matr[["pca"]], dims = 1:5, nfeatures = 5)

  VizDimLoadings(matr, dims = 1:2, reduction = "pca")

  DimPlot(matr, reduction = "pca")

  DimHeatmap(matr, dims = 1, cells = 500, balanced = TRUE)

  DimHeatmap(matr, dims = 1:2, cells = 500, balanced = TRUE)

# Determining Dimensionality
  matr <- JackStraw(matr, num.replicate = 150)
  matr <- ScoreJackStraw(matr, dims = 1:20)

  JackStrawPlot(matr, dims = 1:20)

  ElbowPlot(matr)

#Cluster cells
  matr <- FindNeighbors(matr, dims = 1:10)
  matr <- FindClusters(matr, resolution = 0.5)

#Run Non-Linear Dimensional Reduction (Create UMAP)  
  head(Idents(matr), 5)
  matr <- RunUMAP(matr, dims=1:10)
  DimPlot(matr, reduction = "umap", group.by = "Shh_2")
  
saveRDS(matr, file ="../figures/UMAP.rds")
  
# find markers for every cluster compared to all remaining cells, report only the positive ones
  matr.markers <- FindAllMarkers(matr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  matr.markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC)
  
# Determining top expressed genes for each cluster
  cluster0.markers <- FindMarkers(matr, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  cluster1.markers <- FindMarkers(matr, ident.1 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  cluster2.markers <- FindMarkers(matr, ident.1 = 2, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  cluster3.markers <- FindMarkers(matr, ident.1 = 3, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  cluster4.markers <- FindMarkers(matr, ident.1 = 4, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  cluster5.markers <- FindMarkers(matr, ident.1 = 5, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  cluster6.markers <- FindMarkers(matr, ident.1 = 6, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  cluster7.markers <- FindMarkers(matr, ident.1 = 7, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  cluster8.markers <- FindMarkers(matr, ident.1 = 8, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  cluster9.markers <- FindMarkers(matr, ident.1 = 9, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  cluster10.markers <- FindMarkers(matr, ident.1 = 10, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  cluster11.markers <- FindMarkers(matr, ident.1 = 11, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  
#Generate Expression Heatmap
  matr.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  DoHeatmap(matr, features = top10$gene) + NoLegend()

#Assigning Cell Type Identities to Clusters
  new.cluster.ids <- c("Mesenchyme", "Anterior/Posterior Mesenchyme", "Mesenchyme", "Proximal limb mesenchyme", "Chondrocyte precursors", "Erythocytes", "Mesenchhyme",
                       "Proximal limb mesenchyme", "Vascular cells", "Muscle precursors", "Ectoderm", "Non-erythroid blood")
  names(new.cluster.ids) <- levels(matr)
  matr <- RenameIdents(matr, new.cluster.ids)
  DimPlot(matr, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  
saveRDS(matr, 'matr.rds')  

#Isolate Mesenchymal-like Cells
mese = subset(x = matr, idents= c(0,1,2,3,4,6,7)) 

#Reducing the number of mCherry counted
  count_matrix = mese[["RNA"]]@counts

  #find row index of mCherry
    index_mCherry = which(row.names(count_matrix)=="mCherry")
    count_matrix[index_mCherry,]>=2

  #find mCherry having >= 2 counts
    which(count_matrix[index_mCherry,] >=2)
    mese[['mCherry_2']] = rep('mCherry_count<2', dim(mese)[2])
    mese$mCherry_2[which(count_matrix[index_mCherry,]>=2)] = 'mCherry_count>=2'


#Reducing the number of Shh counted
  count_matrix = mese[["RNA"]]@counts

  #find row index of Shh
    index_Shh2 = which(row.names(count_matrix)=="Shh")
    count_matrix[index_Shh2,]>=2
    
  #find Shh having >= 2 counts
    which(count_matrix[index_Shh2,] >=2)
    mese[['Shh_2']] = rep('Shh_count<2', dim(mese)[2])
    mese$Shh_2[which(count_matrix[index_Shh2,]>=2)] = 'Shh_count>=2'

#Import Transcription Factor List
  list_TF = c('Eno1', 'Ybx1','Sox11','Hmga1','Hoxa10','Msx1','Trp53','Sox9','Ewsr1','Tcf4','Fubp1','Tcf3','Hmga2','Tcf7l2','Prrx1','Smarcc1','Xbp1','Sox4','Mecom','Shox2','Ets2','Atf4','Pitx1','Tcf12','E2f4','Lhx9','Hbp1','Mycn','Nfyc','Prrx2','Rarg','Hoxa9','Jund','Pbx1','Tgif1','Hif1a','Srebf2','Atf2','Hoxd12','Sp3','Foxp1','Hes1','Snai2','Osr1','Cdc5l','Smad2','Usf1','Meox2','Nrf1','E2f5','Osr2','Zkscan3','Nfatc4','Ebf1','Nr2f2','Nfyb','Rest','Nr2f6','Gabpb1','Meis1','Maz','Sp1','Hoxd11','Hey1','Tgif2','Hoxd13','Hoxc10','Lef1','Cux1','Zfp410','Six1','Etv5','Smad1','Yy1','Gabpa','Nr1h2','Zfp423','Nfe2l2','Foxm1','Zfp281','Tbx15','Rara','Arnt','Creb1','Atf1','Msx2','Ctcf','Thap1','Arid5b','Gli3','Jun','Rela','Creb3','Mlx','Elk3','Hoxa11','Tbp','Cxxc1','E2f6','Irf3','Rbpj','Cebpz','Smad4','Ets1','Myc','Meis2','Cebpg','Tfap4','Tbx4','Pknox1','Pitx2','Zfp105','Elf2','Nfya','Stat3','Pbx3','Hsf2','E2f7','Hsf1','Plagl1','Myog','Nfib','Hoxa13','Nfil3','Lhx2','Mef2a','Zfp263','Nfe2l1','Mybl2','Srebf1','Hoxd10','Rxrb','Zfp691','Zfp740','Creb3l2','Sox12','Gsc','Zfp148','Usf2','Tfap2b','Nfia','Zbtb14','Foxj3','Etv4','Fli1','Myod1','Rfx3','Pou2f1','Max','Zscan26','Sox18','Msc','Glis2','Isl1','Bptf','Foxp2','Smad3','Brca1','Hoxd9','Zfx','Irx3','Nfic','E2f3','Nfat5','Tead1','Klf13','Glis1','Zeb1','Tbx5','Etv6','Alx3','Hmbox1','Trp63','Erf','Mecp2','Irf6','Dlx5','Atf7','Pml','Pbx2','Hltf','Stat5b','Mga','Thra','Bbx','Tef','Crem','Arnt2','Nr2f1','Sp4','Zfp282','Irx5','Mef2d','Klf6','Six4','Zbtb7a','E2f8','Ddit3','Gli2','Mef2c','Lmx1b','Rarb','Homez','Meox1','Bach1','Prdm4','Emx2','Zbtb40','Mtf1','Hoxc6','Zic3','Zfp143','Rfx5','Runx2','Tbx2','Foxc1','Pknox2','Hoxa5','Mbd2','Hic2','Nr2c2','Hoxd8','Klf12','Zfp524','Zfp652','Zic2','Hlx','Rfx7','Klf1','Nr2c1','Sox17','Nfatc1','Rreb1','Stat2','Hic1','Zbtb33','E2f1','Batf3','Elk1','Arid3a','Cebpd','Irf1','Clock','Rxra','Junb','Id4','Sox7','Rfx1','Elf3','Hinfp','Foxo4','Elf1','Sox5','Dlx2','Klf5','Sp2','Epas1','Sox8','Etv1','Pax1','Six2','Mafb','Nr4a2','Hivep1','Tead3','Hoxa7','Hoxc11','Tfap2a','Hoxc8','Bcl6b','Stat6','Foxo1','Dlx3','Klf7','Zfp219','Irf2','Zfp128','Meis3','Zbtb6','Zfp384','Stat1','Nfe2','Hoxb2','Mafk','Tfap2c','Pitx3','Zbtb3','Nfkb1','Prdm1','Sox13','Klf16','Eno2')

#Find most hihgly expressed Transcription Facotrs
  tf_list <- FindMarkers(matr, ident.1 = c(0,1,2,3,4,6,7), ident.2 = c(5,8,9,10,11), features = list_TF, logfc.threshold = .25, test.use = "roc", only.pos = TRUE)
  
  #Find most highly expressed transcription factors in Cluster 1
    Cluster1_TFs <- FindMarkers(matr, ident.1 = 1, ident.2 = c(0,2,3,4,5,6,7,8,9,10,11), features = list_TF, logfc.threshold = .25, test.use = "roc", only.pos = TRUE)
  
  #Find most highly expressed transcription factors in cells that exhibit expression of Shh 
    Shh_markers = FindMarkers(mese, ident.1 = "Shh_count>=2", group.by = "Shh_2", features = list_TF, logfc.threshold = .5, test.use = "roc", only.pos = TRUE)
  
  #Find most highly expressed transcription factors in cells that exhibit expression of mCherry
    mCherry_markers = FindMarkers(mese,ident.1 = "mCherry_count>=2", group.by = "mCherry_2", features = list_TF, logfc.threshold = .5, test.use = "roc", only.pos = TRUE)


#Isolating mCherry cells that exhibit pathogenic misexpression
  Idents(mese) = 'mCherry_2'
  pathogenicmCherry = subset(x = mese, idents = 'mCherry_count>=2')

  #Determine positions of cells that exhibit Hoxd13, Hand2, and Shh
    count_matrix1 = pathogenicmCherry[["RNA"]]@counts
    index_nohoxd13 = which(row.names(count_matrix1) == "Hoxd13")
    index_nohand2 = which(row.names(count_matrix1) == "Hand2")
    index_noshh = which(row.names(count_matrix1) == "Shh")

    count_matrix1[index_nohoxd13,]==0
    count_matrix1[index_nohand2,]==0
    count_matrix1[index_noshh,]==0

  #find mCherry cells that do not express Hoxd13,Hand2, and Shh
    which(count_matrix1[index_nohoxd13,] ==0)
    which(count_matrix1[index_nohand2,] ==0)
    which(count_matrix1[index_noshh,] ==0)
    
  #Combining Lists
    pathogenicmCherrycells = intersect(which(count_matrix1[index_nohoxd13,] ==0), which(count_matrix1[index_nohand2,] ==0))
    finalpathogenicmCherrycells = intersect(pathogenicmCherrycells, which(count_matrix1[index_noshh,] ==0))

    pathogenicmCherry[['mutatedmCherrycells']] = rep('pathomCherrycells', dim(pathogenicmCherry)[2])
    pathogenicmCherry$mutatedmCherrycells[finalpathogenicmCherrycells] = 'nonexpressedcells'
    
  #Find pathogenicmCherry_markers
    pathogenicmCherry_markers = FindMarkers(pathogenicmCherry,  ident.1 = "nonexpressedcells", group.by = "mutatedmCherrycells", features = list_TF, logfc.threshold = .2, test.use = "roc", only.pos = TRUE)