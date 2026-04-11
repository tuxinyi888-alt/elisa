# ==============================================================================
# IL-1β Dose-Response Figures — LPS and BCG (separate)
# Stratified by Carprofen Dose
# ==============================================================================
#
# Data files:
#   AIF006_ELISA_IL1b(LPS).csv   (64 rows, 27 NAs)
#   AIF006_ELISA_IL1b(BCG).csv   (64 rows, 51 NAs)
#
# Columns (both files):
#   Donor      – donor ID  (D1–D4)
#   Carprofen  – carprofen dose in µM  (0, 1, 10, 100)
#   Stimulus   – stimulation group     (LPS or BCG)
#   Stim_Conc  – stimulus concentration in ng/mL  (0, 1, 10, 100)
#   IL1b_pgml  – IL-1β concentration in pg/mL (or NA)
#
# Output:
#   IL1b_LPS.png / .pdf
#   IL1b_BCG.png / .pdf
# ==============================================================================

# Set working directory to the folder containing this script (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(ggplot2)

# ---------- 1. Read data ----------
# Files are CSVs with "NA" as missing-value string (read.csv handles this)
lps <- read.csv("AIF006_ELISA_IL1b(LPS).csv")
bcg <- read.csv("AIF006_ELISA_IL1b(BCG).csv")

# ---------- 2. Prepare factors ----------
# Convert Stim_Conc and Carprofen to ordered factors for clean axis labels
prep <- function(df) {
  df %>%
    mutate(
      Stim_Conc  = factor(Stim_Conc,  levels = c(0, 1, 10, 100)),
      Carprofen  = factor(Carprofen,  levels = c(0, 1, 10, 100))
    )
}

lps <- prep(lps)
bcg <- prep(bcg)

# ---------- 3. Summary statistics (mean ± SEM) ----------
summarise_il1b <- function(df) {
  df %>%
    group_by(Stim_Conc, Carprofen) %>%
    summarise(
      mean  = mean(IL1b_pgml, na.rm = TRUE),
      sd    = sd(IL1b_pgml, na.rm = TRUE),
      n     = sum(!is.na(IL1b_pgml)),
      sem   = sd / sqrt(n),
      .groups = "drop"
    )
}

lps_summary <- summarise_il1b(lps)
bcg_summary <- summarise_il1b(bcg)

# ---------- 4. Plotting function ----------
# Reusable function so both figures share the same style
plot_il1b <- function(df_summary, df_raw, stim_label) {

  # Axis label includes the stimulus name
  x_label <- paste0(stim_label, " dose (ng/mL)")

  ggplot(df_summary,
         aes(x = Stim_Conc, y = mean,
             colour = Carprofen, group = Carprofen)) +

    # Error bars (mean ± SEM)
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                  width = 0.15, linewidth = 0.5) +

    # Lines connecting group means
    geom_line(linewidth = 0.8) +

    # Mean points
    geom_point(size = 2.5) +

    # Individual donor points (semi-transparent, slightly jittered)
    geom_point(data = df_raw,
               aes(x = Stim_Conc, y = IL1b_pgml, colour = Carprofen),
               alpha = 0.35, size = 1.5,
               position = position_dodge(width = 0.3),
               show.legend = FALSE) +

    # Colour palette (colour-blind friendly)
    scale_colour_brewer(palette = "Dark2",
                        name = "Carprofen (\u00b5M)") +

    labs(
      x = x_label,
      y = "IL-1\u03b2 (pg/mL)",
      title = paste0("IL-1\u03b2 response to ", stim_label,
                      " — stratified by carprofen dose")
    ) +

    theme_bw(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      legend.position  = "bottom",
      panel.grid.minor = element_blank()
    )
}

# ---------- 5. Build and save figures ----------

# --- LPS figure ---
p_lps <- plot_il1b(lps_summary, lps, "LPS")

ggsave("IL1b_LPS.png", p_lps, width = 6, height = 5, dpi = 300)
ggsave("IL1b_LPS.pdf", p_lps, width = 6, height = 5)
cat("Saved IL1b_LPS.png and IL1b_LPS.pdf\n")

# --- BCG figure ---
p_bcg <- plot_il1b(bcg_summary, bcg, "BCG")

ggsave("IL1b_BCG.png", p_bcg, width = 6, height = 5, dpi = 300)
ggsave("IL1b_BCG.pdf", p_bcg, width = 6, height = 5)
cat("Saved IL1b_BCG.png and IL1b_BCG.pdf\n")

cat("Done!\n")
