
# data_loading ------------------------------------------------------------

source("config.R")

# 读包 ----------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

# 数据处理 --------------------------------------------------------------------

#  读取原始文件 
counts <- readRDS("data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
meta   <- read.delim("data/GSE131907_Lung_Cancer_cell_annotation.txt.gz", 
                     stringsAsFactors = FALSE, check.names = FALSE)

######读取数据后查看矩阵文件和注释文件的内容，确认数据是否正确加载

#  筛选肿瘤样本
tumor_samples <- grep("^LUNG_T", unique(meta$Sample), value = TRUE)
meta <- meta %>% filter(Sample %in% tumor_samples) 

#  重新排序细胞
counts <- counts[, meta$Index]

#  设置细胞ID
cell_id <- paste0(meta$Sample, "_", meta$Index)
colnames(counts) <- cell_id
rownames(meta)   <- cell_id
 
#  创建Seurat对象
obj <- CreateSeuratObject(counts = counts, meta.data = meta, project = "LUAD_SLC7A11")

#  计算线粒体基因比例
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
#  质量控制，设置阈值
VlnPlot(obj, features = c("nFeature_RNA", "percent.mt"), ncol = 2)



# 参考文献“Single-cell best practices”，
# 线粒体阈值，基因数阈值可以根据数据集的特点进行调整，
# 过高的线粒体比例可能表示细胞质量较差，过少或过多的基因数可能表示细胞状态异常或双细胞等问题，
# 根据实际情况设置合理的过滤标准，可以提高后续分析的准确性和可靠性，
# 过滤细胞，去除线粒体基因比例过高的细胞。
# 设置过滤阈值
limit_mt_high <- 15   # 线粒体不能超过 15%
limit_gene_low <- 200 # 基因数最少 200
limit_gene_high <- 6000 # 基因数最多 6000

#  过滤细胞
sc_obj <- subset(obj, subset = nFeature_RNA > limit_gene_low & 
                   nFeature_RNA < limit_gene_high & 
                   percent.mt < limit_mt_high)

# 释放一下内存
rm(obj,counts,meta)
gc()


# 归一化
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 2000)

# 寻找高变基因 (HVGs)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)

# 数据缩放
sc_obj <- ScaleData(sc_obj, features = VariableFeatures(sc_obj))
# 线性降维 (PCA)
sc_obj <- RunPCA(sc_obj, features = VariableFeatures(object = sc_obj), verbose = FALSE)


# 可视化QC -------------------------------------------------------------------

# PCA 
p1 <- DimPlot(sc_obj, reduction = "pca", group.by = "Sample") + 
  ggtitle("PCA Plot by Sample")

# 维度贡献基因图 (Loadings Plot) 
p2 <- VizDimLoadings(sc_obj, dims = 1:2, reduction = "pca")

# 维度热图 (Heatmap)
p3 <- DimHeatmap(sc_obj, dims = 1:6, cells = 500, balanced = TRUE, fast = FALSE)

# 碎石图 (Elbow Plot) 
p4 <- ElbowPlot(sc_obj, ndims = 40)

# 显示这一张最重要的图
print(p3)

# 保存所有 PCA 诊断图 
save_dual_plots(plot = p1, filename = "01_PCA_Scatter", w = 6, h = 5)
save_dual_plots(plot = p2, filename = "01_PCA_Loadings", w = 8, h = 6)
save_dual_plots(plot = p4, filename = "01_PCA_Elbow", w = 6, h = 4)


# UMAP 可视化与聚类 -------------------------------------------------------------

# UMAP (非线性降维)
sc_obj <- RunUMAP(sc_obj, dims = 1:20)

# 构建细胞间的关系网
sc_obj <- FindNeighbors(sc_obj, dims = 1:20)

# 寻找聚类 (FindClusters)
sc_obj <- FindClusters(sc_obj, resolution = 0.5)

#  UMAP 图
p_umap <- DimPlot(sc_obj, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP Clustering (Res 0.5)")

# print(p_umap)

# 保存图
save_dual_plots(p_umap, "01_Final_UMAP")

save_data(sc_obj, "module01_seurat.rds", subdir = "temp")

