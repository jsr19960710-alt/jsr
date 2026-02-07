
# 全局配置--------------------------------------------------------------------

# 根目录设置 -------------------------------------------------------------------

PROJECT_ROOT <- "D:/LUAD_sc_SLC7A11" 

# 设置工作目录 ------------------------------------------------------------------

setwd(PROJECT_ROOT)


# 目录创建 --------------------------------------------------------------------
# source("config.R")        # 加载配置 (确保路径正确)
# 
# for (path in DIR) {
#   if (!dir.exists(path)) dir.create(path, recursive = TRUE)
# }




#  核心生物学参数
# ------------------------------------------------------------------------------
PARAMS <- list(
  target_gene = "SLC7A11",  # 目标基因
  cancer_type = "LUAD",     # 疾病类型
  species     = "Human",    # 物种 (Human/Mouse)
  p_cutoff    = 0.05,       # 显著性阈值
  fdr_cutoff  = 0.25        # FDR 阈值
)

#  全局绘图标准 
# ------------------------------------------------------------------------------
PLOT_STYLE <- list(
  width  = 6,              # 默认宽度 (inch) - 单细胞项目建议改为 8
  height = 5,              # 默认高度 (inch) - 单细胞项目建议改为 6
  dpi    = 300,            # 出版级分辨率
  font   = "sans"
)

# 定义目录结构 ------------------------------------------------------------------


DIR <- list(
  root    = PROJECT_ROOT,
  data    = file.path(PROJECT_ROOT, "data"),           # 原始数据 (10X/FASTQ)
  output  = file.path(PROJECT_ROOT, "output"),         # 总输出
  figures = file.path(PROJECT_ROOT, "output/figures"), # 最终图表
  tables  = file.path(PROJECT_ROOT, "output/tables"),  # 最终表格
  temp    = file.path(PROJECT_ROOT, "temp")            # 中间文件 (RDS)
)


# RDS文件保存设置 ---------------------------------------------------------------
# 先保留
save_data <- function(obj, filename, subdir = "temp") {
  # 检查 subdir 是否有效
  if (!subdir %in% names(DIR)) stop("Invalid subdir! Use: ", paste(names(DIR), collapse=", "))
  
  filepath <- file.path(DIR[[subdir]], filename)
  saveRDS(obj, filepath)
  
  # 打印保存信息 
  file_size <- format(object.size(obj), units = "Mb")
  cat(sprintf("-> Saved: %s (%s) -> %s\n", basename(filename), file_size, subdir))
  
  invisible(filepath)
}


# 加载RDS文件 -----------------------------------------------------------------

load_data <- function(filename, subdir = "temp") {
  filepath <- file.path(DIR[[subdir]], filename)
  if (!file.exists(filepath)) stop("❌ File not found: ", filepath)
  
  obj <- readRDS(filepath)
  cat(sprintf("-> Loaded: %s\n", basename(filename)))
  return(obj)
}


# 保存文件 --------------------------------------------------------------------
save_dual_plots <- function(plot, filename, w = PLOT_STYLE$width, h = PLOT_STYLE$height) {
  # 构建路径
  f_pdf <- file.path(DIR$figures, paste0(filename, ".pdf"))
  f_tif <- file.path(DIR$figures, paste0(filename, ".tiff"))
  
  # 1. Save PDF (Vector)
  cairo_pdf(f_pdf, width = w, height = h)
  print(plot)
  dev.off()
  
  # 2. Save TIFF (Raster)
  tiff(f_tif, width = w, height = h, units = "in", res = PLOT_STYLE$dpi, compression = "lzw")
  print(plot)
  dev.off()}
  

# 空值合并操作符 -----------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

