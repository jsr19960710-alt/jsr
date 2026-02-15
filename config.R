
# 配置文件 --------------------------------------------------------------------


# 1. 根目录设置 ----------------------------------------------------------------
PROJECT_ROOT <- "D:/LUAD_sc_SLC7A11" 
setwd(PROJECT_ROOT)

# 2. 核心生物学参数 ------------------------------------------------------------
PARAMS <- list(
  target_gene = "SLC7A11",  # 目标基因
  cancer_type = "LUAD",     # 疾病类型
  species     = "Human",    # 物种
  p_cutoff    = 0.05,       
  fdr_cutoff  = 0.25        
)

# 3. 全局绘图标准 --------------------------------------------------------------
PLOT_STYLE <- list(
  width  = 8,    # 单细胞推荐宽度
  height = 6,    # 单细胞推荐高度
  dpi    = 300,  
  font   = "sans"
)

# 4. 定义目录结构 --------------------------------------------------------------
DIR <- list(
  root    = PROJECT_ROOT,
  data    = file.path(PROJECT_ROOT, "data"),           # 原始数据
  output  = file.path(PROJECT_ROOT, "output"),         # 总输出
  figures = file.path(PROJECT_ROOT, "output/figures"), # 最终图表
  tables  = file.path(PROJECT_ROOT, "output/tables"),  # 最终表格
  temp    = file.path(PROJECT_ROOT, "temp")            # 中间文件
)

# 5. 工具函数 ------------------------------------------------------------------

# RDS保存 (含文件大小监控)
save_data <- function(obj, filename, subdir = "temp") {
#  路径
filepath <- file.path(DIR[[subdir]], filename)
#  保存
saveRDS(obj, filepath)
#  反馈 
cat(sprintf("-> Saved: %s (%s)\n", filename, format(object.size(obj), units = "Mb")))
}

# RDS 加载
load_data <- function(filename, subdir = "temp") {
# 路径
filepath <- file.path(DIR[[subdir]], filename)
# 2. 读取 (
obj <- readRDS(filepath)
# 3. 反馈
cat(sprintf("-> Loaded: %s\n", filename))
return(obj)
}
# 双格式绘图 (PDF + TIFF)
save_dual_plots <- function(plot, filename, w = PLOT_STYLE$width, h = PLOT_STYLE$height) {
  # Path
  f_pdf <- file.path(DIR$figures, paste0(filename, ".pdf"))
  f_tif <- file.path(DIR$figures, paste0(filename, ".tiff"))
  # PDF
  cairo_pdf(f_pdf, width = w, height = h)
  print(plot)
  dev.off()
  # TIFF
  tiff(f_tif, width = w, height = h, units = "in", res = PLOT_STYLE$dpi, compression = "lzw")
  print(plot)
  dev.off()
}


# CSV 保存 
save_csv <- function(data, filename, subdir = "tables") {
  # 路径 
  filepath <- file.path(DIR[[subdir]], paste0(filename, ".csv"))
  # 写入文件 
  write.csv(data, file = filepath, row.names = FALSE)
  # 反馈
  cat(sprintf("-> Saved CSV: %s\n", filename))
}







# 空值合并操作符
`%||%` <- function(a, b) if (is.null(a)) b else a
