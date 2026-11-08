
library(data.table)
library(ggplot2)
library(ggrepel)

setwd("")
df <- fread("data.csv")

print(head(df))
print(colnames(df))

if (!"Gene" %in% colnames(df)) {
  stop("missing data")
}

if (!"Z" %in% colnames(df)) {
  stop("missing data")
}

p_col_candidates <- c("P_Z", "P", "p", "P.value", "p.value", "pval", "Pval")
p_col <- p_col_candidates[p_col_candidates %in% colnames(df)][1]

df[, Gene := as.character(Gene)]
df[, Z := as.numeric(Z)]
df[, P_value := as.numeric(get(p_col))]

df <- df[!is.na(Gene) & !is.na(Z) & !is.na(P_value)]

p_floor <- 1e-300
df[P_value <= 0, P_value := p_floor]

df[, neg_log10_p := -log10(P_value)]


p_cutoff <- 0.05
y_cutoff <- -log10(p_cutoff)

df[, group := "Not significant"]
df[P_value < p_cutoff, group := "Significant"]

label_n <- 10

label_df <- df[P_value < p_cutoff]


label_df <- label_df[
  order(Gene, -neg_log10_p, -abs(Z))
][
  , .SD[1], by = Gene
]

label_df <- label_df[
  order(-neg_log10_p, -abs(Z))
][
  1:min(.N, label_n)
]


label_left   <- label_df[Z < -0.8]
label_mid    <- label_df[Z >= -0.8 & Z <= 0.8]
label_right  <- label_df[Z > 0.8]

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
