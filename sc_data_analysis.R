
# sc_data_analysis --------------------------------------------------------


library(Seurat)
library(dplyr)

# 细胞纯度
# 预期结果: 只能看到 "Epithelial cells" 这一种类型
# 如果看到其他类型 (如 T cells)，说明 subset 切割不干净
table(sc_data$Cell_type)

# [样本细胞量
# 预期结果: 列出所有病人 (Sample) ID
# 重点观察: 是否有样本的细胞数极少 (例如 < 10 个)?
# 如果太少，计算该病人的平均表达量会产生巨大误差
table(sc_data$Sample)

#  目标基因 存在性与分布
# 预期结果:
# 1. Min/Max 不全是 0
# 2. Mean (均值) > 0
# 如果全是 0，说明这个基因在这个亚群里完全没表达，没法做分组
summary(FetchData(sc_data, vars = "SLC7A11"))

#  直观分布图
# 预期结果: 看到小提琴图有“肚子” (有宽度)，而不是一条死线
# 这能最直观地告诉你这个基因的数据质量
VlnPlot(sc_data, features = "SLC7A11", pt.size = 0.1)



#  宏观结构:主要部件 (Slots)
# 预期: assays, meta.data, reductions, graphs 等
slotNames(epi_obj)

#  元数据结构: "临床记录表"
# 预期: 每一行是一个细胞，每一列是一个变量 (如 Sample, nFeature_RNA)
print(dim(epi_obj@meta.data)) 
print(colnames(epi_obj@meta.data)) 

#  表达矩阵结构: 看看最核心的数据矩阵
# 预期: 行是基因，列是细胞
# 只看左上角 5x5 的小角落，否则屏幕会炸
# (dgCMatrix 是稀疏矩阵的存储格式)
print(epi_obj[["RNA"]]$counts[1:5, 1:5])

#  查看看系统当前是按什么给细胞分类的
# 预期: 应该是 "Epithelial cells" 或者 样本名
Idents(epi_obj) <- "Cell_type" # 强制设为细胞类型
print(levels(epi_obj))