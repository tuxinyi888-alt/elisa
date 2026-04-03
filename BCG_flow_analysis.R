# ============================================================
# BCG Growth Assay - Flow Cytometry Analysis (v2)
# Two conditions: Media only (files 1-3) vs +Carprofen (files 4-6)
# EC and IC plotted separately, each MOI as its own panel
# ============================================================

library(readxl)
library(tidyverse)

# ---- Configuration ----
# setwd("path/to/your/data")

# Group A: 7d = media only (no carprofen re-added after 4h wash)
# Group B: 7d = carprofen re-added after 4h wash
group_files <- list(
  "Media only" = c(D1 = "AIF007-1_Table.xls",
                   D2 = "AIF007-2_Table.xls",
                   D3 = "AIF007-3_Table.xls"),
  "+Carprofen" = c(D4 = "AIF007-4_Table.xls",
                   D5 = "AIF007-5_Table.xls",
                   D6 = "AIF007-6_Table.xls")
)

# ---- Read and combine data ----
read_donor <- function(file, donor_id, condition) {
  df <- read_xls(file, sheet = 1)
  colnames(df) <- c("filename", "beads_freq", "bead_count",
                     "neg_beads_freq", "neg_beads_count",
                     "bcg_freq", "bcg_events", "geomean_PE_A",
                     "debris_freq")
  df <- df %>%
    filter(!filename %in% c("Mean", "SD")) %>%
    mutate(
      donor = donor_id,
      condition = condition,
      timepoint = str_extract(filename, "^[^_]+"),
      type = str_extract(filename, "(?<=_)[A-Z]+(?=_)"),
      carprofen = str_extract(filename, "C\\d+"),
      bcg_moi = str_extract(filename, "B\\d+"),
      bead_count = as.numeric(bead_count),
      bcg_events = as.numeric(bcg_events),
      geomean_PE_A = as.numeric(geomean_PE_A)
    )
  return(df)
}

all_data <- bind_rows(
  lapply(names(group_files), function(cond) {
    files <- group_files[[cond]]
    bind_rows(lapply(names(files), function(d) {
      read_donor(files[d], d, cond)
    }))
  })
)

# ---- Calculate Integrated Fluorescence and Fold Change ----
all_data <- all_data %>%
  mutate(integrated_fluorescence = bcg_events * geomean_PE_A / bead_count)

fold_change <- all_data %>%
  select(donor, condition, timepoint, type, carprofen, bcg_moi, integrated_fluorescence) %>%
  pivot_wider(names_from = timepoint, values_from = integrated_fluorescence, names_prefix = "IF_") %>%
  mutate(fold_change = IF_7dpi / IF_4hpi)

# Set factor levels
fold_change$carprofen <- factor(fold_change$carprofen, levels = c("C0", "C1", "C10", "C100"))
fold_change$bcg_moi <- factor(fold_change$bcg_moi, levels = c("B0", "B1", "B10", "B100"))
fold_change$condition <- factor(fold_change$condition, levels = c("Media only", "+Carprofen"))
fold_change$type <- factor(fold_change$type, levels = c("EC", "IC"))

# ---- Summary statistics ----
summary_stats <- fold_change %>%
  group_by(carprofen, bcg_moi, type, condition) %>%
  summarise(
    mean_fc = mean(fold_change, na.rm = TRUE),
    sd_fc = sd(fold_change, na.rm = TRUE),
    sem_fc = sd(fold_change, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# ---- Generate plots: EC and IC separately, 4 MOI panels each ----
for (st in c("EC", "IC")) {
  plots <- list()

  for (moi in levels(fold_change$bcg_moi)) {
    sum_sub <- summary_stats %>% filter(type == st, bcg_moi == moi)
    fc_sub  <- fold_change %>% filter(type == st, bcg_moi == moi)

    p <- ggplot() +
      geom_col(data = sum_sub,
               aes(x = carprofen, y = mean_fc, fill = condition),
               position = position_dodge(width = 0.7),
               width = 0.6, alpha = 0.7, color = "black") +
      geom_errorbar(data = sum_sub,
                    aes(x = carprofen,
                        ymin = mean_fc - sem_fc,
                        ymax = mean_fc + sem_fc,
                        group = condition),
                    position = position_dodge(width = 0.7),
                    width = 0.2, linewidth = 0.5) +
      geom_point(data = fc_sub,
                 aes(x = carprofen, y = fold_change, group = condition),
                 position = position_dodge(width = 0.7),
                 size = 2, alpha = 0.8) +
      labs(title = moi, x = "Carprofen", y = NULL, fill = "7d Treatment") +
      scale_fill_manual(values = c("Media only" = "#74a9cf", "+Carprofen" = "#ef6548")) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = if (moi == "B0") "top" else "none"
      )

    plots[[moi]] <- p
  }

  # Combine 4 panels
  library(patchwork)
  combined <- (plots[["B0"]] | plots[["B1"]] | plots[["B10"]] | plots[["B100"]]) +
    plot_annotation(
      title = paste0("BCG Growth Assay \u2014 ", st),
      subtitle = "Fold Change (7d / 4h) | Mean \u00b1 SEM | n = 3 donors per group",
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 11)
      )
    )

  # Add shared y-axis label
  combined <- combined & ylab("Fold Change (Integrated Fluorescence)")

  ggsave(paste0("BCG_fold_change_", st, "_v2.pdf"), combined, width = 18, height = 5)
  ggsave(paste0("BCG_fold_change_", st, "_v2.png"), combined, width = 18, height = 5, dpi = 300)
  cat(paste0("Saved: BCG_fold_change_", st, "_v2.pdf / .png\n"))
}

# ---- Print summary ----
cat("\n=== Fold Change Summary ===\n")
print(as.data.frame(summary_stats %>% arrange(type, bcg_moi, condition, carprofen)))
