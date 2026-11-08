# =========================
# 0. 加载包
# =========================
library(data.table)
library(ggplot2)
library(ggrepel)

# =========================
# 1. 设置路径并读取数据
# =========================
setwd("")

# 按你的实际文件名修改
df <- fread("data.csv")
# 如果是 data(3).csv 就改成：
# df <- fread("data(3).csv")

print(head(df))
print(colnames(df))

# =========================
# 2. 检查列名
# =========================
if (!"Gene" %in% colnames(df)) {
  stop("文件中没有 Gene 列")
}

if (!"Z" %in% colnames(df)) {
  stop("文件中没有 Z 列")
}

# 自动识别 P 值列
p_col_candidates <- c("P_Z", "P", "p", "P.value", "p.value", "pval", "Pval")
p_col <- p_col_candidates[p_col_candidates %in% colnames(df)][1]

if (is.na(p_col)) {
  stop("没有找到 P 值列，请确认列名是否为 P_Z、P、p、P.value、pval 等")
}

cat("使用的 P 值列为：", p_col, "\n")

# 转换数据类型
df[, Gene := as.character(Gene)]
df[, Z := as.numeric(Z)]
df[, P_value := as.numeric(get(p_col))]

# 去掉缺失
df <- df[!is.na(Gene) & !is.na(Z) & !is.na(P_value)]

# =========================
# 3. 计算 -log10(P-value)
# =========================
p_floor <- 1e-300
df[P_value <= 0, P_value := p_floor]

df[, neg_log10_p := -log10(P_value)]

# =========================
# 4. 设置显著性
# 显著点统一黄色，不显著点黑色
# =========================
p_cutoff <- 0.05
y_cutoff <- -log10(p_cutoff)

df[, group := "Not significant"]
df[P_value < p_cutoff, group := "Significant"]

# =========================
# 5. 选择需要标注的基因
# 这里只标注显著基因里最显著的前若干个
# =========================
label_n <- 10

label_df <- df[P_value < p_cutoff]

# 同一个 Gene 如果有多行，只保留最显著的一条
label_df <- label_df[
  order(Gene, -neg_log10_p, -abs(Z))
][
  , .SD[1], by = Gene
]

# 取最显著的前 label_n 个
label_df <- label_df[
  order(-neg_log10_p, -abs(Z))
][
  1:min(.N, label_n)
]

# =========================
# 6. 按左右分开标签
# 左侧标签往左上方分散
# 右侧标签往右上方分散
# 中间附近的点也单独处理
# =========================
label_left   <- label_df[Z < -0.8]
label_mid    <- label_df[Z >= -0.8 & Z <= 0.8]
label_right  <- label_df[Z > 0.8]

# =========================
# 7. 绘图
# =========================
p <- ggplot(df, aes(x = Z, y = neg_log10_p)) +
  
  geom_point(
    aes(color = group),
    size = 3.8,
    alpha = 0.9
  ) +
  
  geom_hline(
    yintercept = y_cutoff,
    color = "red",
    linewidth = 0.8
  ) +
  
  annotate(
    "text",
    x = max(df$Z, na.rm = TRUE) * 0.58,
    y = y_cutoff + 0.12,
    label = "-log10(0.05)",
    size = 5,
    fontface = "bold",
    hjust = 0
  ) +
  
  # 左侧：往左上方分散
  geom_text_repel(
    data = label_left,
    aes(label = Gene),
    size = 6,
    fontface = "bold",
    color = "black",
    direction = "both",
    nudge_x = -0.45,
    nudge_y = 0.35,
    hjust = 1,
    box.padding = 0.7,
    point.padding = 0.4,
    force = 3.5,
    force_pull = 0.4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.3,
    seed = 123
  ) +
  
  # 中间：主要向上分散
  geom_text_repel(
    data = label_mid,
    aes(label = Gene),
    size = 6,
    fontface = "bold",
    color = "black",
    direction = "y",
    nudge_y = 0.55,
    box.padding = 0.75,
    point.padding = 0.45,
    force = 4,
    force_pull = 0.35,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.3,
    seed = 123
  ) +
  
  # 右侧：往右上方分散
  geom_text_repel(
    data = label_right,
    aes(label = Gene),
    size = 6,
    fontface = "bold",
    color = "black",
    direction = "both",
    nudge_x = 0.45,
    nudge_y = 0.35,
    hjust = 0,
    box.padding = 0.7,
    point.padding = 0.4,
    force = 3.5,
    force_pull = 0.4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.3,
    seed = 123
  ) +
  
  scale_color_manual(
    values = c(
      "Not significant" = "black",
      "Significant" = "#F2C94C"
    )
  ) +
  
  labs(
    x = "Z-score",
    y = "-log10(P-value)"
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    legend.position = "none",
    
    axis.title = element_text(size = 28, face = "bold"),
    axis.text = element_text(size = 20, color = "black"),
    
    axis.line = element_line(linewidth = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    
    plot.margin = margin(20, 55, 20, 25)
  )

print(p)

# =========================
# 8. 保存图片
# =========================
ggsave(
  filename = "PairWise_Zscore_plot_repel_no_overlap.pdf",
  plot = p,
  width = 8.5,
  height = 7,
  dpi = 300
)

ggsave(
  filename = "PairWise_Zscore_plot_repel_no_overlap.png",
  plot = p,
  width = 8.5,
  height = 7,
  dpi = 300
)