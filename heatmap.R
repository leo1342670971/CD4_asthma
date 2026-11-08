# =========================
# 0. 加载包
# =========================
library(ComplexHeatmap)
library(grid)

# =========================
# 1. 设置路径并读取数据
# =========================
setwd("")

df <- read.csv("评分结果.csv",
               check.names = FALSE,
               stringsAsFactors = FALSE)

print(head(df))
print(colnames(df))

# =========================
# 2. 整理数据
# 原始文件列名仍然是：
# Gene | CD4 | DICE | OneK1K | eQTLgen | GTEX
# 图中把 CD4 显示为 Dynamic eQTLs
# =========================

gene_names <- df$Gene

# 原始数据中的列名
row_order <- c("CD4", "DICE", "OneK1K", "eQTLgen", "GTEX")

# 图中显示的行名
row_labels <- c("Dynamic eQTLs", "DICE", "OneK1K", "eQTLgen", "GTEX")

score_df <- df[, row_order]

# 统一 Yes/No
score_df[] <- lapply(score_df, function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- "No"
  ifelse(toupper(x) == "YES", "Yes", "No")
})

# 转置
mat <- t(as.matrix(score_df))
colnames(mat) <- gene_names
rownames(mat) <- row_labels

print(dim(mat))
print(table(mat))

# =========================
# 3. 设置颜色
# =========================

yes_col <- c(
  "Dynamic eQTLs" = "#F2A6B8",
  "DICE"          = "#B8D760",
  "OneK1K"        = "#BBA7F2",
  "eQTLgen"       = "#F6B36B",
  "GTEX"          = "#8FC6E8"
)

no_col <- "#D9D9D9"

# =========================
# 4. 自定义 Q版圆角正方形
# =========================

round_cell_fun <- function(j, i, x, y, width, height, fill) {
  
  value_now <- mat[i, j]
  row_now   <- rownames(mat)[i]
  
  fill_now <- if (value_now == "Yes") {
    yes_col[row_now]
  } else {
    no_col
  }
  
  w_cm <- convertWidth(width, "cm", valueOnly = TRUE)
  h_cm <- convertHeight(height, "cm", valueOnly = TRUE)
  
  # 控制方块大小
  # 0.80 越大方块越大，间隔越小
  # 0.70 越小方块越小，间隔越大
  side <- min(w_cm, h_cm) * 0.80
  
  grid.roundrect(
    x = x,
    y = y,
    width = unit(side, "cm"),
    height = unit(side, "cm"),
    r = unit(0.08, "cm"),
    gp = gpar(
      fill = fill_now,
      col = "white",
      lwd = 1
    )
  )
}

# =========================
# 5. 设置整体格子大小
# =========================

cell_size <- 0.90

# =========================
# 6. 构建热图主图
# =========================

ht <- Heatmap(
  mat,
  name = "Status",
  
  col = c("Yes" = "white", "No" = "white"),
  rect_gp = gpar(col = NA),
  cell_fun = round_cell_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 18, fontface = "bold"),
  column_names_gp = gpar(fontsize = 11),
  column_names_rot = 45,
  
  show_heatmap_legend = FALSE,
  border = FALSE,
  
  width  = ncol(mat) * unit(cell_size, "cm"),
  height = nrow(mat) * unit(cell_size, "cm")
)

# =========================
# 7. 构建单独图例
# No 显示为 No evidence
# CD4 显示为 Dynamic eQTLs
# =========================

lgd <- Legend(
  title = "Annotation",
  labels = c(
    "No evidence",
    "Dynamic eQTLs",
    "DICE",
    "OneK1K",
    "eQTLgen",
    "GTEX"
  ),
  legend_gp = gpar(
    fill = c(
      no_col,
      yes_col["Dynamic eQTLs"],
      yes_col["DICE"],
      yes_col["OneK1K"],
      yes_col["eQTLgen"],
      yes_col["GTEX"]
    ),
    col = "white"
  ),
  grid_width = unit(7, "mm"),
  grid_height = unit(7, "mm"),
  labels_gp = gpar(fontsize = 13),
  title_gp = gpar(fontsize = 14, fontface = "bold"),
  nrow = 1
)

# =========================
# 8. 保存主图，不带图例
# =========================

pdf("Q版热图_主图_final.pdf",
    width = 22,
    height = 7)

draw(
  ht,
  padding = unit(c(10, 10, 10, 10), "mm")
)

dev.off()

png("Q版热图_主图_final.png",
    width = 9000,
    height = 2400,
    res = 300)

draw(
  ht,
  padding = unit(c(10, 10, 10, 10), "mm")
)

dev.off()

# =========================
# 9. 保存单独图例
# =========================

pdf("Q版热图_图例_final.pdf",
    width = 12,
    height = 2)

grid.newpage()

draw(
  lgd,
  x = unit(0.5, "npc"),
  y = unit(0.5, "npc"),
  just = c("center", "center")
)

dev.off()

png("Q版热图_图例_final.png",
    width = 3800,
    height = 550,
    res = 300)

grid.newpage()

draw(
  lgd,
  x = unit(0.5, "npc"),
  y = unit(0.5, "npc"),
  just = c("center", "center")
)

dev.off()

# =========================
# 10. 在当前绘图设备显示主图
# =========================

draw(
  ht,
  padding = unit(c(10, 10, 10, 10), "mm")
)