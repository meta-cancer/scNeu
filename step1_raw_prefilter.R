cat("<<<scNeu: A pipeline for identify neutrophils from single-cell data>>>")
library(BiocParallel)
library(Seurat)
library(SeuratObject)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(igraph)
library(dplyr)
library(qs)

args=commandArgs(T)
coi <- args[1]
input_path <- args[2]
output_path <- args[3]
config_path <- args[4]

features.IR <- read.table(file.path(config_path, "gene_name_type.TR_IG.txt"), stringsAsFactors = FALSE)
features.prot <- read.table(file.path(config_path, "features.prot.tsv"), stringsAsFactors = FALSE)
features <- c(features.prot$V2,features.IR$V2) # 20234

############ input file #####
data.dir <- paste0(input_path,coi,"/outs/raw_feature_bc_matrix/")
# print(args)
print(coi)
print(input_path)
print(output_path)
# print(.libPaths())

########## output file #####
{
  # if(!dir.exists(paste0(output_path,coi))) dir.create(paste0(output_path,coi))
  # if(!dir.exists(paste0(output_path,"00log"))) dir.create(paste0(output_path,"00log"))
  vlnplot0 <- paste0(output_path,'/',coi,'_plot_vln0.pdf') # before 
  vlnplot1 <- paste0(output_path,'/',coi,'_plot_vln1.pdf') # after soft filter
  vlnplot2 <- paste0(output_path,'/',coi,'_plot_vln2.pdf') # after soft filter
  
  dimplot0 <- paste0(output_path,'/',coi,'_plot_cluster.pdf') # before 
  dimplot1 <- paste0(output_path,'/',coi,'_plot_immune.pdf') # after soft filter
  dimplot2 <- paste0(output_path,'/',coi,'_plot_stromal.pdf') # after soft filter
  
  output_rds  <- paste0(output_path,'/',coi,'.prefilter.rds') 
  output_neu_rds  <- paste0(output_path,'/',coi,'.neu.rds')
  output_qs  <- paste0(output_path,'/',coi,'.prefilter.qs') 
  output_neu_qs  <- paste0(output_path,'/',coi,'.neu.qs') 
  output_meta  <- paste0(output_path,'/',coi,'.anno.meta.txt') 
  output_cluster  <- paste0(output_path,'/',coi,'.cluster.meta.txt') 
  
  output_log  <- paste0(output_path,'/',coi,'.log.txt') 
  
  dimplot1_neu <- paste0(output_path,'/',coi,'_plot_immune_neu.pdf') # after soft filter
  vlnplot2_neu <- paste0(output_path,'/',coi,'_plot_vln2_neu.pdf') # after soft filter
  
  dimplot1_neu2 <- paste0(output_path,'/',coi,'_plot_immune_neu2.pdf') # after soft filter
  vlnplot2_neu2 <- paste0(output_path,'/',coi,'_plot_vln2_neu2.pdf') # after soft filter
}


############### Step 0: read 10X file ============

## read matrix 
seu_count <- Read10X(data.dir)
## create seurat object seu_data
seu_data <- CreateSeuratObject(counts = seu_count)
##  pre-filter, rm empty droplet and these high-quality droplet with cells
seu_data
# seu_data <- subset(seu_data, nCount_RNA > 100 & nFeature_RNA > 50); seu_data
seu_data <- subset(seu_data, nCount_RNA > 200 & nFeature_RNA > 100); seu_data
seu_data <- subset(seu_data,  nCount_RNA < 10000 & nFeature_RNA < 3000 ); seu_data
## print info
print(coi)
seu_data
seu_data@meta.data[1:5,]

## format seu_data
seu_data$Sample <- coi
seu_data <- RenameCells(seu_data,new.names = paste0(coi,'_',colnames(seu_data),'-0'))
# seu_data <- seu_data[rownames(seu_data) %in% features,]  # 这步骤基因先不过滤吧
seu_data
seu_data@meta.data[1:5,]

seu_data$prefilter <- 'N'
seu_data$celltypes <- 'N'
seu_data$celltypes_detail <- 'N'

####### Step 1: soft filter ==========
## Percentage
GeneMT <- read.table(file.path(config_path, "gene_name_type.MT.txt"), stringsAsFactors = F)
GeneHSP_data <- read.table(file.path(config_path, "gene_name_type.HSP.txt"), stringsAsFactors = FALSE)
GeneHSP <- GeneHSP_data[GeneHSP_data$V2 %in% features, "V2"]
GeneRP_data <- read.table(file.path(config_path, "gene_name_type.RP.txt"), stringsAsFactors = FALSE)
GeneRP <- GeneRP_data[GeneRP_data$V2 %in% features, "V2"]
GeneSC_data <- read.table(file.path(config_path, "gene_name_type.sc_dissociation_induced_gene.txt"), stringsAsFactors = FALSE)
GeneSC <- GeneSC_data[GeneSC_data$V2 %in% features, "V2"]
GeneRBC_data <- read.table(file.path(config_path, "gene_list.RBC.txt"), stringsAsFactors = FALSE)
GeneRBC <- GeneRBC_data[GeneRBC_data$V1 %in% features, ]

seu_data[["percent.mt"]] <- PercentageFeatureSet(seu_data, pattern = "^MT-")
seu_data[["percent.RP"]] <- PercentageFeatureSet(seu_data, features = GeneRP[GeneRP %in% rownames(seu_data)])
seu_data[["percent.HSP"]] <- PercentageFeatureSet(seu_data, features = GeneHSP[GeneHSP %in% rownames(seu_data)])
seu_data[["percent.SC"]] <- PercentageFeatureSet(seu_data, features = GeneSC[GeneSC %in% rownames(seu_data)])
seu_data[["percent.RBC"]] <- PercentageFeatureSet(seu_data, features = GeneRBC[GeneRBC %in% rownames(seu_data)])

VlnPlot(seu_data, features = c('nCount_RNA','nFeature_RNA',"percent.mt","percent.HSP","percent.RP","percent.SC","percent.RBC"), pt.size = 0, ncol = 7)
ggsave(vlnplot0, width = 10,height = 3)

#  cutoff 3*MAD
cutoff_mt <-  median( seu_data$percent.mt ) + 3* mad( seu_data$percent.mt ) ## !!!!!!!
if(cutoff_mt < 20) cutoff_mt <- 20

cutoff_RP <-  median( seu_data$percent.RP ) + 3* mad( seu_data$percent.RP ) ## !!!!!!!
if(cutoff_RP < 50) cutoff_RP <- 50

cutoff_HSP <- 20
cutoff_SC <- 20
cutoff_RBC <- 20

seu_data$prefilter <- 'N'
for( cant in c('mt','RP','HSP','SC','RBC')){
  cant_true <- which(seu_data@meta.data[,paste0('percent.',cant)] > get(paste0('cutoff_',cant)) )
  if( length(cant_true) == 0 ) next
  seu_data@meta.data[cant_true,]$prefilter <- paste0('Cant_',cant)
}

table(seu_data$prefilter)
VlnPlot(seu_data[,grep('Cant',seu_data$prefilter,invert = T)],
        features = c('nCount_RNA','nFeature_RNA',"percent.mt","percent.HSP","percent.RP","percent.SC","percent.RBC"), pt.size = 0, ncol = 7)
ggsave(vlnplot1, width = 10,height = 3)

####### Step 2: Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
#### 2.1 first cluster #########
{
  set.seed(12)
  seu_raw <- seu_data
  seu_data <- NormalizeData(seu_data)
  seu_data <- FindVariableFeatures(seu_data, selection.method = "vst", nfeatures = 2000)
  seu_data <- ScaleData(seu_data, features = VariableFeatures(seu_data))
  seu_data <- RunPCA(seu_data)
  seu_data <- RunUMAP(seu_data, dims = 1:10)
  
  seu_data <- FindNeighbors(seu_data, dims = 1:10)
  seu_data <- FindClusters(seu_data, resolution = 1)
}
DimPlot(seu_data, label = T) + theme(aspect.ratio = 1)
ggsave(dimplot0, width = 5,height = 4)
FeaturePlot(seu_data, features = c('PTPRC','CSF3R','CD3D','S100A8','CD68','CD79A'), label = T, ncol = 3)
ggsave(dimplot1, width = 10,height = 6)
FeaturePlot(seu_data, features = c('EPCAM','KRT10','KRT1','VWF','ACTA2','MPG'), label = T, ncol = 3)
ggsave(dimplot2, width = 10,height = 6)
VlnPlot(seu_data[,grep('Cant',seu_data$celltypes,invert = T)],
        features = c('nCount_RNA','nFeature_RNA',"percent.mt","percent.HSP","percent.RP","percent.SC","percent.RBC"), pt.size = 0, ncol = 1)
ggsave(vlnplot2, width = 10,height = 15)

### filter 
Immune_score <- AddModuleScore(seu_data, features = c('PTPRC','CD3D','CD79A','MS4A1','CD68','CD14','FCGR3A'))$Cluster1
Stromal_score <- AddModuleScore(seu_data, features = c('KRT10','VWF','ACTA2','MPG'))$Cluster1
seu_data_meta <- seu_data@meta.data
seu_data_cluster <- data.frame(
  nCell = table(seu_data_meta$seurat_clusters),
  nCount_RNA = tapply(seu_data_meta$nCount_RNA, seu_data_meta$seurat_clusters, median),
  nFeature_RNA = tapply(seu_data_meta$nFeature_RNA, seu_data_meta$seurat_clusters, median),
  Immune_score = tapply(Immune_score, seu_data_meta$seurat_clusters, median),
  Stromal_score = tapply(Stromal_score, seu_data_meta$seurat_clusters, median),
  PTPRC.exp = DotPlot(seu_data, features = c('PTPRC'))$data[,c('avg.exp')],
  PTPRC.pct = DotPlot(seu_data, features = c('PTPRC'))$data[,c('pct.exp')],
  CSF3R.exp = DotPlot(seu_data, features = c('CSF3R'))$data[,c('avg.exp')],
  CSF3R.pct = DotPlot(seu_data, features = c('CSF3R'))$data[,c('pct.exp')],
  S100A8.exp = DotPlot(seu_data, features = c('S100A8'))$data[,c('avg.exp')],
  S100A8.pct = DotPlot(seu_data, features = c('S100A8'))$data[,c('pct.exp')]
)
seu_data_cluster$PTPRC_score <- seu_data_cluster$PTPRC.exp * seu_data_cluster$PTPRC.pct
seu_data_cluster$CSF3R_score <- seu_data_cluster$CSF3R.exp * seu_data_cluster$CSF3R.pct 
seu_data_cluster$celltypes <- 'N'
result <- try({
  # if(sum(seu_data_cluster$nCount_RNA < 500 & seu_data_cluster$nFeature_RNA < 500)) seu_data_cluster[seu_data_cluster$nCount_RNA < 500 & seu_data_cluster$nFeature_RNA < 500, ]$celltypes <- 'Cant'
  if(sum(seu_data_cluster$Stromal_score > 0.1) > 0) {
    seu_data_cluster[seu_data_cluster$Stromal_score > 0.1,]$celltypes <- 'Stromal'
  }
  # if(sum(seu_data_cluster$nCount_RNA > 25000 & seu_data_cluster$nFeature_RNA > 6000)) seu_data_cluster[seu_data_cluster$nCount_RNA > 25000 & seu_data_cluster$nFeature_RNA > 6000, ]$celltypes <- 'Cant'
  if(sum(seu_data_cluster$PTPRC_score > 500) > 0) {
    seu_data_cluster[seu_data_cluster$PTPRC_score > 500,]$celltypes <- 'Immune'
  }
  # if( sum(seu_data_cluster$CSF3R_score > 10) ) if( sum(seu_data_cluster$PTPRC_score > 200) ) seu_data_cluster[seu_data_cluster$CSF3R_score > 10 & seu_data_cluster$PTPRC_score > 200,]$celltypes <- 'Neu'
  if(sum(seu_data_cluster$CSF3R_score > 10 & seu_data_cluster$PTPRC_score > 200) > 0) {
    seu_data_cluster[seu_data_cluster$CSF3R_score > 10 & seu_data_cluster$PTPRC_score > 200,]$celltypes <- 'Neu'
  }
})
table(seu_data_cluster$celltypes)
seu_data_cluster

seu_data$celltypes <- seu_data_cluster[seu_data$seurat_clusters,]$celltypes
table(seu_data$celltypes, seu_data$prefilter)

DimPlot(seu_data, group.by = "celltypes", label = T) + theme(aspect.ratio = 1)

seu_data$celltypes %>% table()
seu_data_neu <- tryCatch({
    seu_data_subset <- seu_data[, seu_data$celltypes == 'Neu']
    if (ncol(seu_data_subset) == 0) stop("结果集为空")
    seu_data_subset
}, error = function(e) {
    message("没有找到属于 'Neu' 的细胞: ", e$message)
    NULL  # 可以返回NULL或者一个空的Seurat对象
})
seu_data_neu

# 映射celltype信息
seu_raw <- AddMetaData(object = seu_raw, metadata = seu_data$celltypes[match(colnames(seu_raw), colnames(seu_data))], col.name = "celltypes")
seu_raw_neu <- tryCatch({
    seu_data_subset <- seu_raw[, seu_raw$celltypes == 'Neu']
    if (ncol(seu_data_subset) == 0) stop("结果集为空")
    seu_data_subset
}, error = function(e) {
    message("没有找到属于 'Neu' 的细胞: ", e$message)
    NULL  # 可以返回NULL或者一个空的Seurat对象
})

if(is.null(seu_data_neu))
{
  qsave(seu_raw,output_qs)
}

qsave(seu_raw,output_qs)
qsave(seu_raw_neu,output_neu_qs)
# saveRDS(seu_data_neu, output_neu_rds)

#### 2.2 second cluster #########
{
  seu_data_neu <- NormalizeData(seu_data_neu)
  seu_data_neu <- FindVariableFeatures(seu_data_neu, selection.method = "vst", nfeatures = 2000)
  seu_data_neu <- ScaleData(seu_data_neu, features = VariableFeatures(seu_data_neu))
  seu_data_neu <- RunPCA(seu_data_neu)
  seu_data_neu <- RunUMAP(seu_data_neu, dims = 1:10)
  
  seu_data_neu <- FindNeighbors(seu_data_neu, dims = 1:10)
  seu_data_neu <- FindClusters(seu_data_neu, resolution = 1)
}

DimPlot(seu_data_neu, label = T) + theme(aspect.ratio = 1)
ggsave(dimplot0, width = 5,height = 4)

# FeaturePlot(seu_data_neu, features = c('CD14','FCGR3A','NAMPT','CXCL8','S100A11','CD68','S100A8','B2M','FTH1'), label = T, ncol = 3)
FeaturePlot(seu_data_neu, features = c('CSF3R','FCGR3B','NAMPT','CXCL8','S100A11','CD68','S100A8','B2M','FTH1'), label = T, ncol = 3)
ggsave(dimplot1_neu, width = 10,height = 9)

VlnPlot(seu_data_neu,
        features = c('nCount_RNA','nFeature_RNA',"percent.mt","percent.HSP","percent.RP","percent.SC"), pt.size = 0, ncol = 1)
ggsave(vlnplot2_neu, width = 10,height = 15)

table(seu_data_neu$seurat_clusters)
seu_data_neu_meta <- seu_data_neu@meta.data
seu_data_cluster <- data.frame(
  nCell = table(seu_data_neu_meta$seurat_clusters),
  nCount_RNA = tapply(seu_data_neu_meta$nCount_RNA, seu_data_neu_meta$seurat_clusters, median),
  nFeature_RNA = tapply(seu_data_neu_meta$nFeature_RNA, seu_data_neu_meta$seurat_clusters, median),
  PTPRC.exp = DotPlot(seu_data_neu, features = c('PTPRC'))$data[,c('avg.exp')],
  PTPRC.pct = DotPlot(seu_data_neu, features = c('PTPRC'))$data[,c('pct.exp')],
  CSF3R.exp = DotPlot(seu_data_neu, features = c('CSF3R'))$data[,c('avg.exp')],
  CSF3R.pct = DotPlot(seu_data_neu, features = c('CSF3R'))$data[,c('pct.exp')],
  S100A8.exp = DotPlot(seu_data_neu, features = c('S100A8'))$data[,c('avg.exp')],
  S100A8.pct = DotPlot(seu_data_neu, features = c('S100A8'))$data[,c('pct.exp')]
)
seu_data_cluster$PTPRC_score <- seu_data_cluster$PTPRC.exp * seu_data_cluster$PTPRC.pct
seu_data_cluster$CSF3R_score <- seu_data_cluster$CSF3R.exp * seu_data_cluster$CSF3R.pct 
seu_data_cluster

seu_data_cluster$celltypes <- 'Neu'
result <- try({
  # if(sum(seu_data_cluster$nCount_RNA > 1500 )) seu_data_cluster[seu_data_cluster$nCount_RNA > 1500,]$celltypes <- 'Cant'
  # if(sum(seu_data_cluster$PTPRC_score < 200 )) seu_data_cluster[seu_data_cluster$PTPRC_score < 200,]$celltypes <- 'Cant'
  # if(sum(seu_data_cluster$CSF3R_score < 50 )) seu_data_cluster[seu_data_cluster$CSF3R_score < 50,]$celltypes <- 'Cant'
  if(sum(seu_data_cluster$nCount_RNA > 1500) > 0) {
    seu_data_cluster[seu_data_cluster$nCount_RNA > 1500,]$celltypes <- 'Cant'
  }
  if(sum(seu_data_cluster$PTPRC_score < 200) > 0) {
    seu_data_cluster[seu_data_cluster$PTPRC_score < 200,]$celltypes <- 'Cant'
  }
  if(sum(seu_data_cluster$CSF3R_score < 50) > 0) {
    seu_data_cluster[seu_data_cluster$CSF3R_score < 50,]$celltypes <- 'Cant'
  }
})

seu_data_neu$celltypes <- seu_data_cluster[seu_data_neu$seurat_clusters,]$celltypes
table(seu_data_neu$celltypes, seu_data_neu$prefilter)

DimPlot(seu_data_neu, group.by = "celltypes", label = T) + theme(aspect.ratio = 1)

# 保存原始数据
seu_raw_neu <- AddMetaData(object = seu_raw_neu, metadata = seu_data_neu$celltypes[match(colnames(seu_raw_neu), colnames(seu_data_neu))], col.name = "celltypes")

seu_data_neu <- tryCatch({
    seu_data_subset <- seu_data_neu[, seu_data_neu$celltypes == 'Neu']
    if (ncol(seu_data_subset) == 0) stop("结果集为空")
    seu_data_subset
}, error = function(e) {
    message("没有找到属于 'Neu' 的细胞: ", e$message)
    NULL  # 可以返回NULL或者一个空的Seurat对象
})
seu_data_neu

# this is the final neutrophil 
seu_data$celltypes <- gsub('Neu','N',seu_data$celltypes)
if (!is.null(seu_data_neu) && ncol(seu_data_neu) > 0) {
    seu_data@meta.data[colnames(seu_data_neu), "celltypes"] <- 'Neu'
} else {
    message("seu_data_neu 是 NULL 或没有列，无法更新 seu_data 的元数据")
}
seu_data$celltypes %>% table

seu_raw <- AddMetaData(object = seu_raw, metadata = seu_data$celltypes[match(colnames(seu_raw), colnames(seu_data))], col.name = "celltypes")

qsave(seu_raw,output_qs)
qsave(seu_raw_neu,output_neu_qs)

# saveRDS(seu_data_neu, output_neu_rds)

#### 2.3 third cluster for visualization#########
{
  seu_data_neu <- NormalizeData(seu_data_neu)
  seu_data_neu <- FindVariableFeatures(seu_data_neu, selection.method = "vst", nfeatures = 2000)
  seu_data_neu <- ScaleData(seu_data_neu, features = VariableFeatures(seu_data_neu))
  seu_data_neu <- RunPCA(seu_data_neu)
  seu_data_neu <- RunUMAP(seu_data_neu, dims = 1:10)
  
  seu_data_neu <- FindNeighbors(seu_data_neu, dims = 1:10)
  seu_data_neu <- FindClusters(seu_data_neu, resolution = 1)
}

DimPlot(seu_data_neu, label = T) + theme(aspect.ratio = 1)
ggsave(dimplot0, width = 5,height = 4)

FeaturePlot(seu_data_neu, features = c('PTPRC','CSF3R','IFITM2','S100A8','B2M','FTH1'), label = T, ncol = 3)
FeaturePlot(seu_data_neu, features = c('CSF3R','FCGR3B','IFITM2','NAMPT','CXCL8','S100A11'), label = T, ncol = 3)
ggsave(dimplot1_neu2, width = 10,height = 6)

VlnPlot(seu_data_neu,
        features = c('nCount_RNA','nFeature_RNA',"percent.mt","percent.HSP","percent.RP","percent.SC","percent.RBC"), pt.size = 0, ncol = 1)
ggsave(vlnplot2_neu2, width = 10,height = 15)

# ########## Step 4: singleR ##########################
# 
# 
# #进行singleR注释
# seu_data_for_SingleR <- GetAssayData(seu_data[,seu_data$prefilter == 'N' & seu_data$celltypes != "Cant"], slot="data") ##获取标准化矩阵
# seu_data.hesc <- SingleR(test = seu_data_for_SingleR, ref = hpca.se, labels = hpca.se$label.main) #
# seu_data.hesc
# 
# #seurat 和 singleR的table表
# # table(seu_data.hesc$labels,seu_data.hesc$seurat_clusters)
# 
# seu_data$singleR <- 'N'
# seu_data@meta.data[rownames(seu_data.hesc),]$singleR <- seu_data.hesc$labels 
# 
# table(seu_data$celltypes, seu_data$singleR)

############  Step 5: output --------------------------------------------------------------------------------------

write.table(seu_data_cluster, output_cluster, quote = F, sep = '\t')
seu_data_meta <- seu_data_neu@meta.data
write.table(seu_data_meta, output_meta, quote = F, sep = '\t')
# seu_data <- subset(seu_data, prefilter == 'N' & celltypes != "Cant")
# saveRDS(seu_data, output_rds)
##  log file
log.file <- data.frame('process'='seurat_prefilter')
log.file$name <- coi
log.file$totalCell <- ncol(seu_data)
log.file$Cant <- nrow(seu_data_meta[seu_data_meta$prefilter != 'N' | seu_data_meta$celltypes == 'Cant',])
log.file$HQcell <- nrow(seu_data_meta[seu_data_meta$prefilter == 'N' & seu_data_meta$celltypes != 'Cant',])
log.file$Neu <- nrow(seu_data_meta[seu_data_meta$prefilter == 'N' & seu_data_meta$celltypes == 'Neu',])
log.file$Neu.pct <- (log.file$Neu / log.file$totalCell * 100 )%>% round(digits = 2)
write.table(log.file, output_log, quote = F, sep = '\t')
saveRDS(seu_data_neu, output_neu_rds)
