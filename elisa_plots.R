# ELISA Cytokine Plots
# Run this script in RStudio with the ELISA.xlsx file in the same directory

# Install packages if needed
if (!require("readxl")) install.packages("readxl")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")

library(readxl)
library(ggplot2)
library(dplyr)

# Set working directory to where the ELISA.xlsx file is
# setwd("path/to/your/elisa/folder")

# Create output directory for plots
dir.create("plots", showWarnings = FALSE)

# Define sheet info
sheets <- list(
  list(name = "IL1b", col = "IL1b_pgml", label = "IL-1\u03b2 (pg/ml)"),
  list(name = "IL6",  col = "IL6_pgml",  label = "IL-6 (pg/ml)"),
  list(name = "TNFa", col = "TNFa_pgml", label = "TNF-\u03b1 (pg/ml)")
)

# Carprofen order on x-axis
carp_levels <- c("C0", "C1", "C10", "C100")
carp_labels <- c("0", "1", "10", "100")

# Stimulus concentration levels
lps_levels <- c("L0", "L1", "L10", "L100")
lps_labels <- c("0", "1", "10", "100")
bcg_levels <- c("B0", "B1", "B10", "B100")
bcg_labels <- c("0", "1", "10", "100")

# Color palettes
lps_colors <- c("L0" = "#1b9e77", "L1" = "#d95f02", "L10" = "#7570b3", "L100" = "#e7298a")
bcg_colors <- c("B0" = "#1b9e77", "B1" = "#d95f02", "B10" = "#7570b3", "B100" = "#e7298a")

# Function to create a single plot
make_plot <- function(data, value_col, y_label, stim_type, title_text) {
  if (stim_type == "LPS") {
    data$Stim_Conc <- factor(data$Stim_Conc, levels = lps_levels, labels = paste0("LPS ", lps_labels, " ng/ml"))
    color_vals <- setNames(
      c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
      paste0("LPS ", lps_labels, " ng/ml")
    )
  } else {
    data$Stim_Conc <- factor(data$Stim_Conc, levels = bcg_levels, labels = paste0("BCG ", bcg_labels, " MOI"))
    color_vals <- setNames(
      c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
      paste0("BCG ", bcg_labels, " MOI")
    )
  }

  data$Carprofen <- factor(data$Carprofen, levels = carp_levels, labels = carp_labels)

  # Replace NA with 0 for plotting (below detection = 0)
  data[[value_col]][is.na(data[[value_col]])] <- 0

  p <- ggplot(data, aes(x = Carprofen, y = .data[[value_col]],
                         color = Stim_Conc, group = Stim_Conc)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(values = color_vals) +
    labs(
      title = title_text,
      x = expression("Carprofen (" * mu * "M)"),
      y = y_label,
      color = if (stim_type == "LPS") "LPS (ng/ml)" else "BCG (MOI)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  return(p)
}

# Generate all plots
donors <- c("D1", "D2", "D3", "D4")
stimuli <- c("LPS", "BCG")

for (s in sheets) {
  # Read the sheet
  df <- read_excel("ELISA.xlsx", sheet = s$name)
  colnames(df) <- c("Donor", "Carprofen", "Stimulus", "Stim_Conc", s$col)

  for (stim in stimuli) {
    stim_data <- df %>% filter(Stimulus == stim)

    # Per-donor plots
    for (d in donors) {
      donor_data <- stim_data %>% filter(Donor == d)

      title_text <- paste0(s$label, " - ", d, " (", stim, " + Carprofen)")
      p <- make_plot(donor_data, s$col, s$label, stim, title_text)

      filename <- paste0("plots/", s$name, "_", stim, "_", d, ".png")
      ggsave(filename, p, width = 6, height = 4, dpi = 300)
      cat("Saved:", filename, "\n")
    }

    # Average plot across 4 donors
    avg_data <- stim_data %>%
      group_by(Carprofen, Stimulus, Stim_Conc) %>%
      summarise(
        mean_val = mean(.data[[s$col]], na.rm = TRUE),
        sd_val = sd(.data[[s$col]], na.rm = TRUE),
        .groups = "drop"
      )
    # Replace NaN (all NA group) with 0
    avg_data$mean_val[is.nan(avg_data$mean_val)] <- 0
    avg_data$sd_val[is.na(avg_data$sd_val)] <- 0

    # Rename for plotting
    colnames(avg_data)[colnames(avg_data) == "mean_val"] <- s$col

    title_text <- paste0(s$label, " - Average of 4 Donors (", stim, " + Carprofen)")

    # Custom plot with error bars for average
    if (stim == "LPS") {
      avg_data$Stim_Conc <- factor(avg_data$Stim_Conc, levels = lps_levels, labels = paste0("LPS ", lps_labels, " ng/ml"))
      color_vals <- setNames(
        c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
        paste0("LPS ", lps_labels, " ng/ml")
      )
    } else {
      avg_data$Stim_Conc <- factor(avg_data$Stim_Conc, levels = bcg_levels, labels = paste0("BCG ", bcg_labels, " MOI"))
      color_vals <- setNames(
        c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
        paste0("BCG ", bcg_labels, " MOI")
      )
    }
    avg_data$Carprofen <- factor(avg_data$Carprofen, levels = carp_levels, labels = carp_labels)

    p <- ggplot(avg_data, aes(x = Carprofen, y = .data[[s$col]],
                               color = Stim_Conc, group = Stim_Conc)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = .data[[s$col]] - sd_val,
                        ymax = .data[[s$col]] + sd_val),
                    width = 0.15, linewidth = 0.5) +
      scale_color_manual(values = color_vals) +
      labs(
        title = title_text,
        x = expression("Carprofen (" * mu * "M)"),
        y = s$label,
        color = if (stim == "LPS") "LPS (ng/ml)" else "BCG (MOI)"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.position = "right",
        panel.grid.minor = element_blank()
      )

    filename <- paste0("plots/", s$name, "_", stim, "_Average.png")
    ggsave(filename, p, width = 6, height = 4, dpi = 300)
    cat("Saved:", filename, "\n")
  }
}

cat("\nAll plots saved in the 'plots/' directory!\n")
