library(ComplexHeatmap)
library(grid)

setwd("")

df <- read.csv("122.csv",
               check.names = FALSE,
               stringsAsFactors = FALSE)
print(head(df))
print(colnames(df))


gene_names <- df$Gene

row_order <- c("CD4", "DICE", "OneK1K", "eQTLgen", "GTEX")

row_labels <- c("Dynamic eQTLs", "DICE", "OneK1K", "eQTLgen", "GTEX")

score_df <- df[, row_order]

score_df[] <- lapply(score_df, function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- "No"
  ifelse(toupper(x) == "YES", "Yes", "No")
})


mat <- t(as.matrix(score_df))
colnames(mat) <- gene_names
rownames(mat) <- row_labels

print(dim(mat))
print(table(mat))


yes_col <- c(
  "Dynamic eQTLs" = "#F2A6B8",
  "DICE"          = "#B8D760",
  "OneK1K"        = "#BBA7F2",
  "eQTLgen"       = "#F6B36B",
  "GTEX"          = "#8FC6E8"
)

no_col <- "#D9D9D9"

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

cell_size <- 0.90

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

pdf("final.pdf",
    width = 22,
    height = 7)

draw(
  ht,
  padding = unit(c(10, 10, 10, 10), "mm")
)

dev.off()

png("final.png",
    width = 9000,
    height = 2400,
    res = 300)

draw(
  ht,
  padding = unit(c(10, 10, 10, 10), "mm")
)

dev.off()

pdf("final1.pdf",
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

png("final1.png",
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

draw(
  ht,
  padding = unit(c(10, 10, 10, 10), "mm")
)
