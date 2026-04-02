# ============================================================
# BCG Growth Assay - Flow Cytometry Analysis
# Integrated Fluorescence & Fold Change Calculation
# Each BCG MOI produces a separate plot with EC and IC grouped
# ============================================================

# Load required packages
library(readxl)
library(tidyverse)

# ---- Configuration ----
# Set your working directory to the folder containing the data files
# setwd("path/to/your/data")

# Donor file mapping
donor_files <- c(
  "D1" = "AIF007-1_Table.xls",
  "D2" = "AIF007-2_Table.xls",
  "D3" = "AIF007-3_Table.xls",
  "D4" = "AIF007-4_Table.xls"
)

# ---- Read and combine data from all donors ----
read_donor <- function(file, donor_id) {
  df <- read_xls(file, sheet = 1)

  colnames(df) <- c("filename", "beads_freq", "bead_count",
                     "neg_beads_freq", "neg_beads_count",
                     "bcg_freq", "bcg_events", "geomean_PE_A",
                     "debris_freq")

  # Remove Mean and SD rows
  df <- df %>% filter(!filename %in% c("Mean", "SD"))

  # Parse filename to extract experimental variables
  df <- df %>%
    mutate(
      donor = donor_id,
      timepoint = str_extract(filename, "^[^_]+"),        # 4hpi or 7dpi
      type = str_extract(filename, "(?<=_)[A-Z]+(?=_)"),  # EC or IC
      carprofen = str_extract(filename, "C\\d+"),          # C0, C1, C10, C100
      bcg_moi = str_extract(filename, "B\\d+"),            # B0, B1, B10, B100
      bead_count = as.numeric(bead_count),
      bcg_events = as.numeric(bcg_events),
      geomean_PE_A = as.numeric(geomean_PE_A)
    )

  return(df)
}

# Read all donors
all_data <- bind_rows(
  lapply(names(donor_files), function(d) {
    read_donor(donor_files[d], d)
  })
)

# ---- Calculate Integrated Fluorescence ----
# IF = BCG events x Geometric Mean (PE-A) / Bead count
# (Total beads added is constant and cancels out in fold change)
all_data <- all_data %>%
  mutate(
    integrated_fluorescence = bcg_events * geomean_PE_A / bead_count
  )

# ---- Calculate Fold Change (7d / 4h) ----
fold_change <- all_data %>%
  select(donor, timepoint, type, carprofen, bcg_moi, integrated_fluorescence) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = integrated_fluorescence,
    names_prefix = "IF_"
  ) %>%
  mutate(
    fold_change = IF_7dpi / IF_4hpi
  )

# Set factor levels for correct ordering
fold_change$carprofen <- factor(fold_change$carprofen,
                                levels = c("C0", "C1", "C10", "C100"))
fold_change$bcg_moi <- factor(fold_change$bcg_moi,
                              levels = c("B0", "B1", "B10", "B100"))
fold_change$type <- factor(fold_change$type, levels = c("EC", "IC"))

# ---- Summary statistics ----
summary_stats <- fold_change %>%
  group_by(carprofen, bcg_moi, type) %>%
  summarise(
    mean_fc = mean(fold_change, na.rm = TRUE),
    sd_fc = sd(fold_change, na.rm = TRUE),
    sem_fc = sd(fold_change, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# ---- Generate one plot per BCG MOI ----
moi_levels <- levels(fold_change$bcg_moi)

for (moi in moi_levels) {

  sum_sub <- summary_stats %>% filter(bcg_moi == moi)
  fc_sub  <- fold_change %>% filter(bcg_moi == moi)

  p <- ggplot() +
    # Grouped bar chart: EC vs IC side by side
    geom_col(data = sum_sub,
             aes(x = carprofen, y = mean_fc, fill = type),
             position = position_dodge(width = 0.7),
             width = 0.6, alpha = 0.7, color = "black") +
    # Error bars (mean +/- SEM)
    geom_errorbar(data = sum_sub,
                  aes(x = carprofen,
                      ymin = mean_fc - sem_fc,
                      ymax = mean_fc + sem_fc,
                      group = type),
                  position = position_dodge(width = 0.7),
                  width = 0.2, linewidth = 0.5) +
    # Individual donor data points
    geom_point(data = fc_sub,
               aes(x = carprofen, y = fold_change, group = type),
               position = position_dodge(width = 0.7),
               size = 2, alpha = 0.8) +
    # Labels and theme
    labs(
      title = paste0("BCG Growth Assay - ", moi),
      subtitle = "Fold Change (7d / 4h) | Mean \u00b1 SEM | n = 4 donors",
      x = "Carprofen Concentration",
      y = "Fold Change (Integrated Fluorescence)",
      fill = "Sample Type"
    ) +
    scale_fill_manual(values = c("EC" = "#4292c6", "IC" = "#ef6548")) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "top"
    )

  print(p)

  # Save each plot
  ggsave(paste0("BCG_fold_change_", moi, ".pdf"), p, width = 7, height = 5)
  ggsave(paste0("BCG_fold_change_", moi, ".png"), p, width = 7, height = 5, dpi = 300)

  cat(paste0("Saved: BCG_fold_change_", moi, ".pdf / .png\n"))
}

# ---- Print summary table ----
cat("\n=== Fold Change Summary ===\n")
print(as.data.frame(summary_stats))

cat("\n=== Individual Donor Fold Changes ===\n")
print(as.data.frame(fold_change %>%
  select(donor, type, carprofen, bcg_moi, fold_change) %>%
  arrange(bcg_moi, type, carprofen, donor)))
