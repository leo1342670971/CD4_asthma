# ============================================================
# 5组韦恩图：更居中、更规整版本
# CD4 / DICE / OneK1K / eQTLgen / GTEX
# ============================================================

# 如果没安装
# install.packages(c("data.table", "VennDiagram"))

library(data.table)
library(VennDiagram)
library(grid)

# 1. 路径
work_dir <- ""
setwd(work_dir)

# 如果你的文件名是 评分(1).csv，就改这里
input_file <- file.path(work_dir, "评分.csv")

# 2. 读取数据
df <- fread(input_file, encoding = "UTF-8")

# 3. 提取并去重函数
extract_genes <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  
  # 支持逗号、分号、竖线、空格分隔
  genes <- unlist(strsplit(x, "[,;，；\\|\\s]+"))
  genes <- trimws(genes)
  
  genes <- genes[
    genes != "" &
      !is.na(genes) &
      !toupper(genes) %in% c("NA", "NAN", "NULL", "NONE", ".", "-")
  ]
  
  unique(genes)
}

# 4. 固定集合顺序
# 这个顺序对应的布局通常是：
# 1 = 上, 2 = 左, 3 = 左下, 4 = 右下, 5 = 右
gene_sets <- list(
  CD4      = extract_genes(df$CD4),
  DICE     = extract_genes(df$DICE),
  OneK1K   = extract_genes(df$OneK1K),
  eQTLgen  = extract_genes(df$eQTLgen),
  GTEX     = extract_genes(df$GTEX)
)

# 5. 计算每组基因数和百分比
all_genes <- unique(unlist(gene_sets))
total_n <- length(all_genes)

set_n <- sapply(gene_sets, length)
set_pct <- round(set_n / total_n * 100)

label_names <- paste0(names(gene_sets), "\n(", set_pct, "%)")

# 输出计数表
count_table <- data.table(
  Group = names(gene_sets),
  Unique_gene_number = set_n,
  Percent_of_union = paste0(set_pct, "%")
)

print(count_table)
fwrite(count_table, file.path(work_dir, "Venn_gene_count.csv"))

# 6. 绘图
venn_plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,
  category.names = label_names,
  
  fill = c(
    "#F4CFCF",  # CD4
    "#D8CDE6",  # DICE
    "#D5E2F0",  # OneK1K
    "#D8E4D8",  # eQTLgen
    "#E8DDB7"   # GTEX
  ),
  alpha = 0.45,
  
  col = "white",
  lwd = 2,
  
  cex = 1.2,
  fontface = "plain",
  fontfamily = "sans",
  
  cat.cex = 1.35,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  
  # 手动调整标签位置（很关键）
  # 顺序对应：CD4, DICE, OneK1K, eQTLgen, GTEX
  cat.pos  = c(0, 180, 235, 305, 0),
  cat.dist = c(0.06, 0.08, 0.08, 0.08, 0.08),
  
  # 整体轻微旋转；如果还歪，可以改成 -3, 3, 5 试试
  rotation.degree = -5,
  
  margin = 0.03
)

# 7. 保存为正方形图片（比长方形更不容易“歪”）
png(
  filename = file.path(work_dir, "Venn_CD4_DICE_OneK1K_eQTLgen_GTEX_square.png"),
  width = 3200,
  height = 3200,
  res = 300,
  bg = "white"
)

grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.95, height = 0.95))
grid.draw(venn_plot)
popViewport()

dev.off()

# 8. 保存PDF
pdf(
  file = file.path(work_dir, "Venn_CD4_DICE_OneK1K_eQTLgen_GTEX_square.pdf"),
  width = 10,
  height = 10
)

grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.95, height = 0.95))
grid.draw(venn_plot)
popViewport()

dev.off()

cat("绘图完成：\n")
cat(file.path(work_dir, "Venn_CD4_DICE_OneK1K_eQTLgen_GTEX_square.png"), "\n")
cat(file.path(work_dir, "Venn_CD4_DICE_OneK1K_eQTLgen_GTEX_square.pdf"), "\n")