library(data.table)
library(stringr)
library(ggplot2)
library(ggrepel)

work_dir <- ""
setwd(work_dir)

plot_df <- fread(
  file.path(work_dir, "05_plot_data_eGenes_vs_significant_signals.csv")
)

plot_df[, time := fifelse(str_detect(cell_type, "_0h$"),  "0h",
                          fifelse(str_detect(cell_type, "_16h$"), "16h",
                                  fifelse(str_detect(cell_type, "_40h$"), "40h",
                                          fifelse(str_detect(cell_type, "_5d$"),  "5d",
                                                  fifelse(str_detect(cell_type, "_LA$"),  "LA", "Other")))))]


plot_df[, label := cell_type]

plot_df[, label := str_replace(label, "_(0h|16h|40h|5d|LA)$", "")]
plot_df[, label := str_replace(label, "_stim$", "")]
plot_df[, label := str_replace(label, "_uns$", "")]


plot_df[, time := factor(time, levels = c("0h", "16h", "40h", "5d", "LA", "Other"))]

spearman_res <- cor.test(
  plot_df$N_eGenes,
  plot_df$N_sig,
  method = "spearman",
  exact = FALSE
)

spearman_summary <- data.table(
  method = "Spearman",
  rho = as.numeric(spearman_res$estimate),
  p_value = spearman_res$p.value
)

fwrite(
  spearman_summary,
  file.path(work_dir, "Figure4B_spearman_result.csv")
)

print(spearman_summary)

p <- ggplot(
  plot_df,
  aes(
    x = N_sig,
    y = N_eGenes,
    color = time,
    label = label
  )
) +
  geom_point(size = 3.2, alpha = 0.85) +
  geom_text_repel(
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    segment.color = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "0h"    = "#66C2A5",
      "16h"   = "#FC8D62",
      "40h"   = "#8DA0CB",
      "5d"    = "#E78AC3",
      "LA"    = "#A6D854",
      "Other" = "grey60"
    ),
    name = "Time\npoint"
  ) +
  labs(
    x = "N_MR_signals",
    y = "N_eGenes",
    title = paste0(
      "Spearman rho = ",
      round(as.numeric(spearman_res$estimate), 2),
      ", P = ",
      signif(spearman_res$p.value, 3)
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

print(p)

ggsave(
  filename = file.path(work_dir, "Figure4B_style_scatterplot.pdf"),
  plot = p,
  width = 9,
  height = 8
)

ggsave(
  filename = file.path(work_dir, "Figure4B_style_scatterplot.png"),
  plot = p,
  width = 9,
  height = 8,
  dpi = 300
)
