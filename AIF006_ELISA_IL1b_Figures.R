###############################################
# IL-1β ELISA Figures – LPS & BCG
# Adapted from SandwichELISA-Automated script
###############################################

# Load packages (install if missing)
packages <- c("ggplot2", "dplyr", "cowplot")
installed <- rownames(installed.packages())
to_install <- packages[!packages %in% installed]
if (length(to_install) > 0) install.packages(to_install)

library(ggplot2)
library(dplyr)
library(cowplot)

#####
# USER INPUT
base_dir <- "."  # current directory; change if needed

#####
# Import data
lps <- read.csv(file.path(base_dir, "AIF006_ELISA_IL1b(LPS).csv"))
bcg <- read.csv(file.path(base_dir, "AIF006_ELISA_IL1b(BCG).csv"))

# Combine
all_data <- bind_rows(lps, bcg)

# Convert Carprofen and Stim_Conc to factors for plotting
all_data <- all_data %>%
  mutate(
    Carprofen = factor(Carprofen, levels = c(0, 1, 10, 100)),
    Stim_Conc = factor(Stim_Conc, levels = c(0, 1, 10, 100))
  )

#####
# Summary statistics (mean +/- SEM)
summary_data <- all_data %>%
  group_by(Stimulus, Carprofen, Stim_Conc) %>%
  summarise(
    mean_conc = mean(IL1b_pgml, na.rm = TRUE),
    sd_conc   = sd(IL1b_pgml, na.rm = TRUE),
    n         = sum(!is.na(IL1b_pgml)),
    sem_conc  = sd_conc / sqrt(n),
    .groups   = "drop"
  )

#####
# Shared theme matching the original script style
elisa_theme <- theme_minimal(base_size = 16) +
  theme(
    aspect.ratio = 0.8,
    axis.title = element_text(size = 18),
    axis.text  = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length = unit(2, "mm"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 14)
  )

#####
# Plotting function for each stimulus
make_il1b_plot <- function(stimulus_name) {

  plot_title <- paste0("IL-1\u03B2 (", stimulus_name, ")")

  df_summary <- summary_data %>% filter(Stimulus == stimulus_name)
  df_points  <- all_data %>% filter(Stimulus == stimulus_name)

  # Drop Stim_Conc == 0 rows (unstimulated baseline, all NA)
  df_summary <- df_summary %>% filter(Stim_Conc != 0)
  df_points  <- df_points %>% filter(Stim_Conc != 0)

  p <- ggplot(df_summary, aes(x = Stim_Conc, y = mean_conc, fill = Carprofen)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
    geom_errorbar(
      aes(ymin = pmax(mean_conc - sem_conc, 0), ymax = mean_conc + sem_conc),
      position = position_dodge(width = 0.8),
      width = 0.25
    ) +
    geom_point(
      data = df_points,
      aes(x = Stim_Conc, y = IL1b_pgml, group = Carprofen),
      position = position_dodge(width = 0.8),
      size = 2, shape = 21, fill = "white", color = "black"
    ) +
    scale_fill_manual(
      values = c("0" = "#4DAF4A", "1" = "#377EB8", "10" = "#FF7F00", "100" = "#E41A1C"),
      name = expression("Carprofen (" * mu * "M)")
    ) +
    labs(
      title = plot_title,
      x = paste0(stimulus_name, " concentration (ng/mL)"),
      y = expression("IL-1" * beta * " (pg/mL)")
    ) +
    elisa_theme

  return(p)
}

#####
# Generate figures
p_lps <- make_il1b_plot("LPS")
p_bcg <- make_il1b_plot("BCG")

print(p_lps)
print(p_bcg)

#####
# Combined panel
p_combined <- plot_grid(p_lps, p_bcg, ncol = 2, labels = c("A", "B"), label_size = 20)
print(p_combined)

#####
# Export figures
export_dir <- file.path(base_dir, "AIF006_IL1b_Figures")
dir.create(export_dir, showWarnings = FALSE)

fig_width  <- 18
fig_height <- 14

# LPS
ggsave(file.path(export_dir, "AIF006_IL1b_LPS.png"),
       plot = p_lps, width = fig_width, height = fig_height, units = "cm", dpi = 300)
ggsave(file.path(export_dir, "AIF006_IL1b_LPS.pdf"),
       plot = p_lps, width = fig_width, height = fig_height, units = "cm")

# BCG
ggsave(file.path(export_dir, "AIF006_IL1b_BCG.png"),
       plot = p_bcg, width = fig_width, height = fig_height, units = "cm", dpi = 300)
ggsave(file.path(export_dir, "AIF006_IL1b_BCG.pdf"),
       plot = p_bcg, width = fig_width, height = fig_height, units = "cm")

# Combined
ggsave(file.path(export_dir, "AIF006_IL1b_Combined.png"),
       plot = p_combined, width = 34, height = 16, units = "cm", dpi = 300)
ggsave(file.path(export_dir, "AIF006_IL1b_Combined.pdf"),
       plot = p_combined, width = 34, height = 16, units = "cm")

cat("Figures exported to:", export_dir, "\n")
