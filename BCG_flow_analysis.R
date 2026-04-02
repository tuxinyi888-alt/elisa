# ============================================================
# BCG Growth Assay - Flow Cytometry Analysis
# Integrated Fluorescence & Fold Change Calculation
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

# Choose sample type: "EC" (extracellular) or "IC" (intracellular)
sample_type <- "EC"

# ---- Read and combine data from all donors ----
read_donor <- function(file, donor_id) {
  df <- read_xls(file, sheet = 1)

  # Rename columns for easier handling
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
      # Clean numeric columns (remove % signs)
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

# ---- Filter and Calculate ----
# Filter for chosen sample type and non-zero BCG MOI
data_filtered <- all_data %>%
  filter(type == sample_type, bcg_moi != "B0")

# Calculate Integrated Fluorescence (IF)
# IF = BCG events × Geometric Mean (PE-A) × (Total beads added / Bead count)
# Since Total beads added is constant and cancels in fold change,
# we use: IF_norm = BCG events × Geometric Mean (PE-A) / Bead count
data_filtered <- data_filtered %>%
  mutate(
    integrated_fluorescence = bcg_events * geomean_PE_A / bead_count
  )

# ---- Calculate Fold Change (7d / 4h) ----
# Pivot to get 4h and 7d side by side
fold_change <- data_filtered %>%
  select(donor, timepoint, carprofen, bcg_moi, integrated_fluorescence) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = integrated_fluorescence,
    names_prefix = "IF_"
  ) %>%
  mutate(
    fold_change = IF_7dpi / IF_4hpi
  )

# ---- Set factor levels for correct ordering ----
fold_change$carprofen <- factor(fold_change$carprofen,
                                levels = c("C0", "C1", "C10", "C100"))
fold_change$bcg_moi <- factor(fold_change$bcg_moi,
                              levels = c("B1", "B10", "B100"))

# ---- Summary statistics for plotting ----
summary_stats <- fold_change %>%
  group_by(carprofen, bcg_moi) %>%
  summarise(
    mean_fc = mean(fold_change, na.rm = TRUE),
    sd_fc = sd(fold_change, na.rm = TRUE),
    sem_fc = sd(fold_change, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# ---- Plot: Bar chart with SEM error bars + individual data points ----
p <- ggplot() +
  # Bar chart (mean)
  geom_col(data = summary_stats,
           aes(x = carprofen, y = mean_fc, fill = carprofen),
           width = 0.7, alpha = 0.7, color = "black") +
  # Error bars (mean ± SEM)
  geom_errorbar(data = summary_stats,
                aes(x = carprofen, ymin = mean_fc - sem_fc, ymax = mean_fc + sem_fc),
                width = 0.2, linewidth = 0.5) +
  # Individual donor data points
  geom_jitter(data = fold_change,
              aes(x = carprofen, y = fold_change),
              width = 0.1, size = 2, alpha = 0.8) +
  # Facet by BCG MOI
  facet_wrap(~ bcg_moi, scales = "free_y",
             labeller = labeller(bcg_moi = c("B1" = "BCG MOI: B1",
                                              "B10" = "BCG MOI: B10",
                                              "B100" = "BCG MOI: B100"))) +
  # Labels and theme
  labs(
    title = "BCG Growth Assay - Fold Change (7d / 4h)",
    subtitle = paste0("Sample type: ", sample_type, " | Mean ± SEM | n = 4 donors"),
    x = "Carprofen Concentration",
    y = "Fold Change (Integrated Fluorescence)",
    fill = "Carprofen"
  ) +
  scale_fill_brewer(palette = "Blues") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

print(p)

# Save plot
ggsave("BCG_fold_change_barplot.pdf", p, width = 10, height = 5)
ggsave("BCG_fold_change_barplot.png", p, width = 10, height = 5, dpi = 300)

cat("\nPlot saved as BCG_fold_change_barplot.pdf and BCG_fold_change_barplot.png\n")

# ---- Print summary table ----
cat("\n=== Fold Change Summary ===\n")
print(as.data.frame(summary_stats))

cat("\n=== Individual Donor Fold Changes ===\n")
print(as.data.frame(fold_change %>% select(donor, carprofen, bcg_moi, fold_change)))
