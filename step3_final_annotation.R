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

custom_colors <- c(
  "Immune" = "#b2df8a",
  "Neu" = "#00BEC4",
  "N" = "#F8736A",
  "Cant" = "lightgrey"
)

args=commandArgs(T)
coi <- args[1]
input_path <- args[2]
output_path <- args[3]
pred_path <-  args[4]
final_path <- args[5]
config_path <- args[6]

############ input file #####
data.dir <- paste0(input_path,coi,"/outs/raw_feature_bc_matrix/")
print(coi)
print(input_path)
print(output_path)
print(pred_path)
print(final_path)

{
  vlnplot0 <- paste0(output_path,'/',coi,'_plot_vln0.pdf') # before 
  vlnplot1 <- paste0(output_path,'/',coi,'_plot_vln1.pdf') # after soft filter
  vlnplot2 <- paste0(output_path,'/',coi,'_plot_vln2.pdf') # after soft filter
  
  dimplot0 <- paste0(output_path,'/',coi,'_plot_cluster.pdf') # before 
  dimplot1 <- paste0(output_path,'/',coi,'_plot_immune.pdf') # after soft filter
  dimplot2 <- paste0(output_path,'/',coi,'_plot_stromal.pdf') # after soft filter
  
  output_rds  <- paste0(input_path,'/',coi,'.prefilter.rds') 
  output_neu_rds  <- paste0(input_path,'/',coi,'.neu.rds')
  output_qs  <- paste0(input_path,'/',coi,'.prefilter.qs') 
  output_neu_qs  <- paste0(input_path,'/',coi,'.neu.qs') 
  output_meta  <- paste0(input_path,'/',coi,'.anno.meta.txt') 
  output_cluster  <- paste0(input_path,'/',coi,'.cluster.meta.txt') 
  
  output_log  <- paste0(input_path,'/',coi,'.log.txt') 
  
  dimplot1_neu <- paste0(output_path,'/',coi,'_plot_immune_neu.pdf') # after soft filter
  vlnplot2_neu <- paste0(output_path,'/',coi,'_plot_vln2_neu.pdf') # after soft filter
  
  dimplot1_neu2 <- paste0(output_path,'/',coi,'_plot_immune_neu2.pdf') # after soft filter
  vlnplot2_neu2 <- paste0(output_path) # after soft filter
  
  pred_csv <- paste0(pred_path,'/','predictions.csv')
  labels_npy <- paste0(pred_path,'/','labels.npy')
  features_npy <- paste0(pred_path,'/','features.npy')
  SVMplot1 <- paste0(output_path,'/',coi,'_SVMpred_vs_Clustering.png')
  SVMplot2 <- paste0(output_path,'/',coi,'_SVMpred_vs_Clustering.pdf')
  finalplot1 <- paste0(output_path,'/',coi,'_Final.png')
  finalplot2 <- paste0(output_path,'/',coi,'_Final.pdf')
  final_qs  <- paste0(output_path,'/',coi,'.total.qs')
  final_rds <- paste0(final_path,'/',coi,'_neu.rds')
  neu_qs <- paste0(final_path,'/',coi,'_neu.qs')

  # delete
  fp1 <- paste0(final_path,'/',coi,'_Summary.png')
  fp2 <- paste0(final_path,'/',coi,'_Summary.pdf')
}
result <- try({
seu_data <- qread(output_qs)
seu_data_neu <- qread(output_neu_qs)
y_pred <- as.vector(read.csv(pred_csv, header = FALSE))$V1
seu_data$y_pred <- y_pred

{
  set.seed(123)
  seu_data <- NormalizeData(seu_data)
  seu_data <- FindVariableFeatures(seu_data, selection.method = "vst", nfeatures = 2000)
  seu_data <- ScaleData(seu_data, features = VariableFeatures(seu_data))
  seu_data <- RunPCA(seu_data)
  seu_data <- RunUMAP(seu_data, dims = 1:12)
  
  seu_data <- FindNeighbors(seu_data, dims = 1:12)
  seu_data <- FindClusters(seu_data, resolution = 1.5)
}

p1 <- DimPlot(seu_data, group.by = "y_pred")  + theme(aspect.ratio = 1) +ggtitle("SVM predict")
p2 <- DimPlot(seu_data, group.by = "celltypes")  + theme(aspect.ratio = 1) +ggtitle("final celltypes(second)")
combined_plot <- p1 + p2
combined_plot
ggsave(filename=SVMplot1, plot = combined_plot, width = 12, height = 6)
ggsave(filename=SVMplot2, plot = combined_plot, width = 12, height = 6)
ggsave(SVMplot1, width = 12, height = 6)
ggsave(SVMplot2, width = 12, height = 6)
})

result <- try({
seu_data$celltypes_step1 <- seu_data$celltypes
seu_data$celltypes_step3 <- seu_data$celltypes
### filter2
seu_data_meta <- seu_data@meta.data
seu_data_cluster <- data.frame(
  nCell = table(seu_data_meta$seurat_clusters),
  nCount_RNA = tapply(seu_data_meta$nCount_RNA, seu_data_meta$seurat_clusters, median),
  nFeature_RNA = tapply(seu_data_meta$nFeature_RNA, seu_data_meta$seurat_clusters, median),
  PTPRC.exp = DotPlot(seu_data, features = c('PTPRC'))$data[,c('avg.exp')],
  PTPRC.pct = DotPlot(seu_data, features = c('PTPRC'))$data[,c('pct.exp')],
  CSF3R.exp = DotPlot(seu_data, features = c('CSF3R'))$data[,c('avg.exp')],
  CSF3R.pct = DotPlot(seu_data, features = c('CSF3R'))$data[,c('pct.exp')],
  CD68.exp = DotPlot(seu_data, features = c('CD68'))$data[,c('avg.exp')],
  CD68.pct = DotPlot(seu_data, features = c('CD68'))$data[,c('pct.exp')]
)
seu_data_cluster$PTPRC_score <- seu_data_cluster$PTPRC.exp * seu_data_cluster$PTPRC.pct
seu_data_cluster$CSF3R_score <- seu_data_cluster$CSF3R.exp * seu_data_cluster$CSF3R.pct
seu_data_cluster$CD68_score <- seu_data_cluster$CD68.exp * seu_data_cluster$CD68.pct

seu_data_cluster$celltypes <- 'Neu'
if(sum(seu_data_cluster$nCount_RNA > 2750 | seu_data_cluster$CSF3R_score < 10 | seu_data_cluster$PTPRC_score < 100 | seu_data_cluster$CD68_score > 120) > 0) {
  seu_data_cluster[seu_data_cluster$nCount_RNA > 2750 | seu_data_cluster$CSF3R_score < 10 | seu_data_cluster$PTPRC_score < 100 | seu_data_cluster$CD68_score > 120,]$celltypes <- 'N'
}
table(seu_data_cluster$celltypes)

seu_data$celltypes_step3 <- seu_data_cluster[seu_data$seurat_clusters,]$celltypes
seu_data$celltypes <- ifelse(seu_data$celltypes_step1 == "Neu" & !is.na(seu_data$celltypes_step3) & seu_data$celltypes_step3 != "Neu", "N", seu_data$celltypes)

})

# 先保存一下确保有结果输出
# final_neu <- seu_data[,seu_data$celltypes == 'Neu' & seu_data$prefilter == 'N']
final_neu <- tryCatch({
  subset_result <- seu_data[, seu_data$celltypes == 'Neu' & seu_data$prefilter == 'N']
  if (ncol(subset_result) == 0) {
    # 如果子集为空，返回 NULL
    NULL
  } else {
    # 如果子集非空，返回子集
    subset_result
  }
}, error = function(e) {
  # 捕获错误并返回 NULL
  cat("No cells found, returning NULL.\n")
  NULL
})
saveRDS(final_neu,final_rds)
qsave(final_neu,neu_qs)
final_neu

# 使用png函数保存图像
# png(paste0('00_SVMpred.png'), width = 12, height = 6, units = "in", res = 300)
# print(combined_plot)
# dev.off()

result <- try({
# 开始进行处理
# 有争议的neu为模型判断是但是聚类不是的细胞 或者 某个seurat_clusters预测是细胞的概率达5% 
clusters <- unique(seu_data$seurat_clusters)
seu_data$marked <- FALSE
# 遍历每一个群集，计算y_pred为1的细胞占比
for (cluster in clusters) {
  # 获取当前群集的所有细胞
  cells_in_cluster <- WhichCells(seu_data, idents = cluster)
  # 计算y_pred为1的细胞数
  y_pred_1_count <- sum(seu_data$y_pred[cells_in_cluster] == 1)
  # 计算总细胞数
  total_count <- length(cells_in_cluster)
  # 计算y_pred为1的细胞占比
  y_pred_1_ratio <- y_pred_1_count / total_count
  # 如果占比超过一定值，则将这些细胞标记出来
  if (y_pred_1_ratio > 0.10) {
    seu_data$marked[cells_in_cluster] <- TRUE
  }
}
cells_to_mark <- WhichCells(seu_data, expression = celltypes != "Neu" & y_pred == 1)
seu_data$marked[cells_to_mark] <- TRUE
table(seu_data$marked)

print(table(seu_data$celltypes))
p1 <- DimPlot(seu_data, group.by = "celltypes",cols = custom_colors)  + theme(aspect.ratio = 1) +ggtitle("step1: Cluster Annotation")
})

########################## 加法步骤 #########################
result <- try({
Controversial_neu <- tryCatch({
  # 提取子集
  subset_result <- subset(seu_data, subset = celltypes != "Neu" & marked == TRUE)
  # 检查子集细胞数量
  if (ncol(subset_result) < 20) {
    # 细胞太少的话直接按照每个细胞来看，如果基因PTPRC和CSF3R表达量达标，则标记为Neu
    cells_to_check <- colnames(subset_result)
    # 对每个细胞进行检查
    for (cell in cells_to_check) {
      # 获取该细胞的基因表达量
      ptprc_expr <- subset_result@assays$RNA@data["PTPRC", cell]
      csf3r_expr <- subset_result@assays$RNA@data["CSF3R", cell]
      # 如果表达量超过阈值，则将其标记为 'Neu'
      if (ptprc_expr > 2 & csf3r_expr > 1) {
        seu_data$celltypes[cell] <- "Neu"
      }
    }
    NULL
  } else {
    subset_result
  }
}, error = function(e) {
  cat("Too few Cells, Skipping.\n")
  NULL
})

if (!is.null(Controversial_neu)) {
{
  Controversial_neu <- NormalizeData(Controversial_neu)
  Controversial_neu <- FindVariableFeatures(Controversial_neu, selection.method = "vst", nfeatures = 2000)
  Controversial_neu <- ScaleData(Controversial_neu, features = VariableFeatures(Controversial_neu))
  Controversial_neu <- RunPCA(Controversial_neu)
  Controversial_neu <- RunUMAP(Controversial_neu, dims = 1:10)
  
  Controversial_neu <- FindNeighbors(Controversial_neu, dims = 1:10)
  Controversial_neu <- FindClusters(Controversial_neu, resolution = 1)
}

DimPlot(Controversial_neu, label = T) + theme(aspect.ratio = 1)
ggsave(dimplot0, width = 5,height = 4)

# FeaturePlot(Controversial_neu, features = c('CD14','FCGR3A','NAMPT','CXCL8','S100A11','CD68','S100A8','B2M','FTH1'), label = T, ncol = 3)
FeaturePlot(Controversial_neu, features = c('CSF3R','FCGR3B','NAMPT','CXCL8','S100A11','CD68','S100A8','B2M','FTH1'), label = T, ncol = 3)
ggsave(dimplot1_neu, width = 10,height = 9)

VlnPlot(Controversial_neu,
        features = c('nCount_RNA','nFeature_RNA',"percent.mt","percent.HSP","percent.RP","percent.SC"), pt.size = 0, ncol = 1)
ggsave(vlnplot2_neu, width = 10,height = 15)

table(Controversial_neu$seurat_clusters)
seu_data_neu <- Controversial_neu
seu_data_neu_meta <- Controversial_neu@meta.data
seu_data_cluster <- data.frame(
  nCell = table(seu_data_neu_meta$seurat_clusters),
  nCount_RNA = tapply(seu_data_neu_meta$nCount_RNA, seu_data_neu_meta$seurat_clusters, median),
  nFeature_RNA = tapply(seu_data_neu_meta$nFeature_RNA, seu_data_neu_meta$seurat_clusters, median),
  PTPRC.exp = DotPlot(seu_data_neu, features = c('PTPRC'))$data[,c('avg.exp')],
  PTPRC.pct = DotPlot(seu_data_neu, features = c('PTPRC'))$data[,c('pct.exp')],
  CSF3R.exp = DotPlot(seu_data_neu, features = c('CSF3R'))$data[,c('avg.exp')],
  CSF3R.pct = DotPlot(seu_data_neu, features = c('CSF3R'))$data[,c('pct.exp')],
  S100A8.exp = DotPlot(seu_data_neu, features = c('S100A8'))$data[,c('avg.exp')],
  S100A8.pct = DotPlot(seu_data_neu, features = c('S100A8'))$data[,c('pct.exp')],
  CD68.exp = DotPlot(seu_data_neu, features = c('CD68'))$data[,c('avg.exp')],
  CD68.pct = DotPlot(seu_data_neu, features = c('CD68'))$data[,c('pct.exp')]
)
seu_data_cluster$PTPRC_score <- seu_data_cluster$PTPRC.exp * seu_data_cluster$PTPRC.pct
seu_data_cluster$CSF3R_score <- seu_data_cluster$CSF3R.exp * seu_data_cluster$CSF3R.pct 
seu_data_cluster$CD68_score <- seu_data_cluster$CD68.exp * seu_data_cluster$CD68.pct 
seu_data_cluster

seu_data_cluster$celltypes <- 'Neu'
# if(sum(seu_data_cluster$nCount_RNA > 2500 )) seu_data_cluster[seu_data_cluster$nCount_RNA > 2500,]$celltypes <- 'Cant'
# if(sum(seu_data_cluster$PTPRC_score < 200 )) seu_data_cluster[seu_data_cluster$PTPRC_score < 200,]$celltypes <- 'Cant'
# if(sum(seu_data_cluster$CSF3R_score < 50 )) seu_data_cluster[seu_data_cluster$CSF3R_score < 50,]$celltypes <- 'Cant'
# if(sum(seu_data_cluster$CD68_score > 150 )) seu_data_cluster[seu_data_cluster$CD68_score < 150,]$celltypes <- 'Cant'
if (sum(seu_data_cluster$nCount_RNA > 2500) > 0) {
  matching_cells <- seu_data_cluster$nCount_RNA > 2500
  if (any(matching_cells)) {
    seu_data_cluster[matching_cells, ]$celltypes <- 'Cant'
  }
}
if (sum(seu_data_cluster$PTPRC_score < 20) > 0) {
  matching_cells <- seu_data_cluster$PTPRC_score < 20
  if (any(matching_cells)) {
    seu_data_cluster[matching_cells, ]$celltypes <- 'Cant'
  }
}
if (sum(seu_data_cluster$CSF3R_score < 30) > 0) {
  matching_cells <- seu_data_cluster$CSF3R_score < 30
  if (any(matching_cells)) {
    seu_data_cluster[matching_cells, ]$celltypes <- 'Cant'
  }
}
if (sum(seu_data_cluster$CD68_score > 80) > 0) {
  matching_cells <- seu_data_cluster$CD68_score > 80
  if (any(matching_cells)) {
    seu_data_cluster[matching_cells, ]$celltypes <- 'Cant'
  }
}

table(seu_data_cluster$celltypes)

Controversial_neu$celltypes <- seu_data_cluster[Controversial_neu$seurat_clusters,]$celltypes
table(Controversial_neu$celltypes, Controversial_neu$prefilter)

# 映射回去
common_cells <- intersect(Cells(seu_data), Cells(Controversial_neu))
controversial_celltypes <- Controversial_neu$celltypes[common_cells]
# 将 celltypes 信息映射回 seu_data
seu_data$celltypes[common_cells] <- controversial_celltypes
print(table(seu_data$celltypes))
p2 <- DimPlot(seu_data, group.by = "celltypes")  + theme(aspect.ratio = 1) +ggtitle("final Annotation")
combined_plot <- p1 + p2
combined_plot
ggsave(filename=finalplot1, plot = combined_plot, width = 12, height = 6)
ggsave(filename=finalplot2, plot = combined_plot, width = 12, height = 6)
ggsave(finalplot1, width = 12, height = 6)
ggsave(finalplot2, width = 12, height = 6)

final_neu <- tryCatch({
  subset_result <- seu_data[, seu_data$celltypes == 'Neu' & seu_data$prefilter == 'N']
  if (ncol(subset_result) == 0) {
    NULL
  } else {
    subset_result
  }
}, error = function(e) {
  # 捕获错误并返回 NULL
  cat("No cells found, returning NULL.\n")
  NULL
})
saveRDS(final_neu,final_rds)
qsave(final_neu,neu_qs)
}
})

############### 最后一轮聚类，减法步骤 #################
# 首先检查 final_neu 中的细胞数量
result <- try({
if (!is.null(final_neu) && ncol(final_neu) > 30) {
  final_neu <- NormalizeData(final_neu)
  final_neu <- FindVariableFeatures(final_neu, selection.method = "vst", nfeatures = 2000)
  final_neu <- ScaleData(final_neu, features = VariableFeatures(final_neu))
  final_neu <- RunPCA(final_neu)
  final_neu <- RunUMAP(final_neu, dims = 1:15)
  final_neu <- FindNeighbors(final_neu, dims = 1:15)
  final_neu <- FindClusters(final_neu, resolution = 1.5)
  table(final_neu$seurat_clusters)
  seu_data_neu <- final_neu
  seu_data_neu_meta <- final_neu@meta.data
  seu_data_cluster <- data.frame(
    nCell = table(seu_data_neu_meta$seurat_clusters),
    nCount_RNA = tapply(seu_data_neu_meta$nCount_RNA, seu_data_neu_meta$seurat_clusters, median),
    nFeature_RNA = tapply(seu_data_neu_meta$nFeature_RNA, seu_data_neu_meta$seurat_clusters, median),
    PTPRC.exp = DotPlot(seu_data_neu, features = c('PTPRC'))$data[,c('avg.exp')],
    PTPRC.pct = DotPlot(seu_data_neu, features = c('PTPRC'))$data[,c('pct.exp')],
    CSF3R.exp = DotPlot(seu_data_neu, features = c('CSF3R'))$data[,c('avg.exp')],
    CSF3R.pct = DotPlot(seu_data_neu, features = c('CSF3R'))$data[,c('pct.exp')],
    S100A8.exp = DotPlot(seu_data_neu, features = c('S100A8'))$data[,c('avg.exp')],
    S100A8.pct = DotPlot(seu_data_neu, features = c('S100A8'))$data[,c('pct.exp')],
    CD68.exp = DotPlot(seu_data_neu, features = c('CD68'))$data[,c('avg.exp')],
    CD68.pct = DotPlot(seu_data_neu, features = c('CD68'))$data[,c('pct.exp')]
  )
  seu_data_cluster$PTPRC_score <- seu_data_cluster$PTPRC.exp * seu_data_cluster$PTPRC.pct
  seu_data_cluster$CSF3R_score <- seu_data_cluster$CSF3R.exp * seu_data_cluster$CSF3R.pct 
  seu_data_cluster$CD68_score <- seu_data_cluster$CD68.exp * seu_data_cluster$CD68.pct 
  seu_data_cluster
  seu_data_cluster$celltypes <- 'Neu'
  if (sum(seu_data_cluster$nCount_RNA > 2600) > 0) {
    matching_cells <- seu_data_cluster$nCount_RNA > 2600
    if (any(matching_cells)) {
      seu_data_cluster[matching_cells, ]$celltypes <- 'Cant'
    }
  }
  if (sum(seu_data_cluster$PTPRC_score < 20) > 0) {
    matching_cells <- seu_data_cluster$PTPRC_score < 20
    if (any(matching_cells)) {
      seu_data_cluster[matching_cells, ]$celltypes <- 'Cant'
    }
  }
  if (sum(seu_data_cluster$CSF3R_score < 30) > 0) {
    matching_cells <- seu_data_cluster$CSF3R_score < 30
    if (any(matching_cells)) {
      seu_data_cluster[matching_cells, ]$celltypes <- 'Cant'
    }
  }
  if (sum(seu_data_cluster$CD68_score > 100) > 0) {
    matching_cells <- seu_data_cluster$CD68_score > 100
    if (any(matching_cells)) {
      seu_data_cluster[matching_cells, ]$celltypes <- 'Cant'
    }
  }
  table(seu_data_cluster$celltypes)
  final_neu$celltypes <- seu_data_cluster[final_neu$seurat_clusters,]$celltypes
  table(final_neu$celltypes, final_neu$prefilter)
  # 映射回去
  common_cells <- intersect(Cells(seu_data), Cells(final_neu))
  controversial_celltypes <- final_neu$celltypes[common_cells]
  # 将 celltypes 信息映射回 seu_data
  seu_data$celltypes[common_cells] <- controversial_celltypes
  print(table(seu_data$celltypes))

  final_neu <- tryCatch({
    subset_result <- final_neu[, final_neu$celltypes == 'Neu' & final_neu$prefilter == 'N']
    if (ncol(subset_result) == 0) {
      NULL
    } else {
      subset_result
    }
  }, error = function(e) {
    # 捕获错误并返回 NULL
    cat("No cells found, returning NULL.\n")
    NULL
  })
  saveRDS(final_neu,final_rds)
  qsave(final_neu,neu_qs)
} else {
  cat("The number of cells in final_neu too few, no analysis was performed.\n")
}

# 更新 seu_data$celltypes，对于 seu_data$prefilter 标记不为 'N' 的细胞
# seu_data$celltypes <- ifelse(seu_data$prefilter != 'N', seu_data$prefilter, seu_data$celltypes)
seu_data$celltypes <- ifelse(seu_data$prefilter != 'N', "Cant", seu_data$celltypes)
qsave(seu_data,final_qs)
})

############### 和大群注释所得的正常pipeline做对比展示（可删除） #################
filtered_df <- tryCatch({
  final_df <- read.csv(file.path(config_path, "reference_Normal_Annotation.csv"), header = TRUE)
  filtered <- final_df[final_df$SampleID == coi, ]
  
  if (nrow(filtered) == 0) {
    cat("No cells found in filtered_df for SampleID =", coi, ". Skipping this step.\n")
    NULL
  } else {
    filtered
  }
}, error = function(e) {
  cat("Error loading or filtering reference_Normal_Annotation.csv:", e$message, "\n")
  NULL
})

# 如果 filtered_df 为 NULL，则跳过后续步骤
if (!is.null(filtered_df)) {
  print(nrow(filtered_df))
  
  filtered_values <- filtered_df$X
  cell_names <- colnames(seu_data)
  seu_data$Neu_mark <- "N"
  seu_data$Neu_mark[cell_names %in% filtered_values] <- "Neu"
  print(table(seu_data$Neu_mark))
  
  final_neu <- tryCatch({
    subset_result <- seu_data[, seu_data$celltypes == 'Neu' & seu_data$prefilter == 'N']
    if (ncol(subset_result) == 0) {
      cat("No cells in final_neu subset, skipping final_neu processing.\n")
      NULL
    } else {
      subset_result
    }
  }, error = function(e) {
    cat("Error during Neu subset filtering:", e$message, "\n")
    NULL
  })
  
  if (!is.null(final_neu)) {
    print(table(final_neu$Neu_mark))
  }
  
  # DimPlot 绘图
  p2 <- DimPlot(seu_data, group.by = "y_pred") + theme(aspect.ratio = 1) + ggtitle("step2: SVM Annotation")
  p3 <- DimPlot(seu_data, group.by = "Neu_mark") + theme(aspect.ratio = 1) + ggtitle("ref: General Annotation")
  p4 <- DimPlot(seu_data, group.by = "celltypes", cols = custom_colors) + theme(aspect.ratio = 1) + ggtitle("scNeu: Final Annotation")
  
  combined_plot <- p1 + p2 + p3 + p4
  combined_plot
  
  # 保存图像
  ggsave(filename = fp1, plot = combined_plot, width = 12, height = 12)
  ggsave(filename = fp2, plot = combined_plot, width = 12, height = 12)
  ggsave(fp1, width = 12, height = 12)
  ggsave(fp2, width = 12, height = 12)
} else {
  cat("Filtered dataframe is NULL, skipping annotation and plotting steps.\n")
}
