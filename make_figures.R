# ─────────────────────────────────────────────────────────────────────────────
# make_figures.R
# Generates Figure 2, 3, 4 for carprofen / MDM report
#
# Requirements: ggplot2, dplyr, tidyr, readxl, patchwork, stringr
# Install if needed:
#   install.packages(c("tidyverse","readxl","patchwork"))
# ─────────────────────────────────────────────────────────────────────────────

library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(patchwork)
library(stringr)

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_report <- function() {
  theme_classic(base_size = 8, base_family = "sans") +
    theme(
      axis.line        = element_line(linewidth = 0.5),
      axis.ticks       = element_line(linewidth = 0.4),
      axis.text        = element_text(size = 7.5),
      axis.title       = element_text(size = 8),
      plot.title       = element_text(size = 9, hjust = 0.5),
      legend.title     = element_text(size = 7.5),
      legend.text      = element_text(size = 7.5),
      legend.key.size  = unit(0.4, "cm"),
      plot.tag         = element_text(size = 11, face = "bold"),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3)
    )
}

# Carprofen colour palette (light → dark blue)
carp_colors <- c("0"   = "#AECDE8",
                 "1"   = "#5BA4CF",
                 "10"  = "#2171B5",
                 "100" = "#084594")
carp_labels <- c("0" = "0 µg/mL", "1" = "1 µg/mL",
                 "10" = "10 µg/mL", "100" = "100 µg/mL")

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Load & reshape ELISA data
# ─────────────────────────────────────────────────────────────────────────────
read_elisa_sheet <- function(path, sheet, cytokine_col) {
  read_excel(path, sheet = sheet) |>
    select(Donor, Carprofen, Stimulus, Stim_Conc, value = all_of(cytokine_col)) |>
    mutate(
      cytokine = cytokine_col,
      value    = suppressWarnings(as.numeric(value))
    )
}

elisa_path <- "ELISA.xlsx"

elisa_raw <- bind_rows(
  read_elisa_sheet(elisa_path, "IL1b", "IL1b_pgml"),
  read_elisa_sheet(elisa_path, "IL6",  "IL6_pgml"),
  read_elisa_sheet(elisa_path, "TNFa", "TNFa_pgml")
)

# Numeric carprofen & stimulus concentration from labels
elisa <- elisa_raw |>
  mutate(
    carp_num  = as.numeric(str_remove(Carprofen, "C")),
    stim_num  = as.numeric(str_remove(Stim_Conc, "[LB]")),
    carp_fac  = factor(carp_num, levels = c(0, 1, 10, 100)),
    cytokine  = factor(cytokine,
                       levels = c("IL1b_pgml", "IL6_pgml", "TNFa_pgml"),
                       labels = c("IL-1\u03b2 (pg/mL)", "IL-6 (pg/mL)", "TNF-\u03b1 (pg/mL)"))
  )

# Summary: mean ± SEM per (carp, stim)
elisa_summary <- function(df) {
  df |>
    group_by(cytokine, Carprofen, carp_num, carp_fac, stim_num) |>
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      sem_val  = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
      n        = sum(!is.na(value)),
      .groups  = "drop"
    ) |>
    filter(!is.na(mean_val))
}

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Build a single ELISA panel (line plot)
# ─────────────────────────────────────────────────────────────────────────────
elisa_panel <- function(raw_df, summ_df, cyt_label, show_legend = FALSE) {
  raw_sub  <- filter(raw_df,  cytokine == cyt_label, !is.na(value))
  summ_sub <- filter(summ_df, cytokine == cyt_label)

  p <- ggplot(summ_sub, aes(x = stim_num, colour = carp_fac, group = carp_fac)) +
    # Individual donor dots (small, transparent, same colour)
    geom_jitter(data = raw_sub,
                aes(x = stim_num, y = value, colour = carp_fac),
                width = 0.08, size = 1.4, alpha = 0.45,
                inherit.aes = FALSE) +
    # Mean line
    geom_line(aes(y = mean_val), linewidth = 0.9) +
    # Mean points
    geom_point(aes(y = mean_val), size = 2.8, shape = 16) +
    # SEM error bars
    geom_errorbar(aes(ymin = mean_val - sem_val,
                      ymax = mean_val + sem_val),
                  width = 0.15, linewidth = 0.7) +
    scale_colour_manual(values = carp_colors, labels = carp_labels,
                        name = "Carprofen") +
    scale_x_continuous(breaks = unique(raw_sub$stim_num)) +
    expand_limits(y = 0) +
    labs(title = str_remove(cyt_label, " \\(pg/mL\\)"),
         y = cyt_label, x = NULL) +
    theme_report()

  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — ELISA LPS
# ─────────────────────────────────────────────────────────────────────────────
lps_raw  <- filter(elisa, Stimulus == "LPS")
lps_summ <- elisa_summary(lps_raw)

p2a <- elisa_panel(lps_raw, lps_summ, "IL-1\u03b2 (pg/mL)") +
  labs(x = "LPS concentration (ng/mL)")
p2b <- elisa_panel(lps_raw, lps_summ, "IL-6 (pg/mL)") +
  labs(x = "LPS concentration (ng/mL)")
p2c <- elisa_panel(lps_raw, lps_summ, "TNF-\u03b1 (pg/mL)", show_legend = TRUE) +
  labs(x = "LPS concentration (ng/mL)")

fig2 <- (p2a + p2b + p2c) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

ggsave("Figure2_ELISA_LPS.pdf", fig2, width = 7.2, height = 2.8, units = "in")
ggsave("Figure2_ELISA_LPS.png", fig2, width = 7.2, height = 2.8, units = "in", dpi = 300)
message("Saved Figure 2")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — ELISA BCG
# ─────────────────────────────────────────────────────────────────────────────
bcg_raw  <- filter(elisa, Stimulus == "BCG")
bcg_summ <- elisa_summary(bcg_raw)

p3a <- elisa_panel(bcg_raw, bcg_summ, "IL-1\u03b2 (pg/mL)") +
  labs(x = "BCG (MOI)")
p3b <- elisa_panel(bcg_raw, bcg_summ, "IL-6 (pg/mL)") +
  labs(x = "BCG (MOI)")
p3c <- elisa_panel(bcg_raw, bcg_summ, "TNF-\u03b1 (pg/mL)", show_legend = TRUE) +
  labs(x = "BCG (MOI)")

fig3 <- (p3a + p3b + p3c) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

ggsave("Figure3_ELISA_BCG.pdf", fig3, width = 7.2, height = 2.8, units = "in")
ggsave("Figure3_ELISA_BCG.png", fig3, width = 7.2, height = 2.8, units = "in", dpi = 300)
message("Saved Figure 3")

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Load BCG flow-cytometry data (AIF007 xls files)
#     Files 1–3 = Media only (transient)  |  Files 4–6 = +Carprofen (sustained)
# ─────────────────────────────────────────────────────────────────────────────
read_aif <- function(file_num) {
  condition <- if (file_num <= 3) "Media only" else "+Carprofen"
  donor     <- paste0("D", ((file_num - 1) %% 3) + 1)
  path      <- sprintf("AIF007-%d_Table.xls", file_num)

  df <- read_excel(path, col_names = FALSE) |>
    setNames(c("sample", "bead_pct", "beads", "nobead_pct", "nobeads",
               "bcg_pct", "bcg_n", "bcg_gmean", "debris_pct"))

  df |>
    filter(str_detect(sample, "\\.fcs$")) |>
    mutate(
      sample    = str_remove(sample, "\\.fcs$"),
      condition = condition,
      donor     = donor,
      beads     = as.numeric(beads),
      bcg_n     = as.numeric(bcg_n)
    ) |>
    separate(sample, into = c("timepoint", "fraction", "carp_moi"), sep = "_") |>
    mutate(
      timepoint = as.integer(str_extract(timepoint, "\\d+")),
      carprofen = as.integer(str_extract(carp_moi, "(?<=C)\\d+")),
      moi       = as.integer(str_extract(carp_moi, "(?<=B)\\d+"))
    ) |>
    select(condition, donor, timepoint, fraction, carprofen, moi, beads, bcg_n)
}

bcg_raw_all <- bind_rows(lapply(1:6, read_aif))

# Fold change = (7dpi BCG/beads) / (4hpi BCG/beads)
bcg_fc <- bcg_raw_all |>
  mutate(norm = bcg_n / beads) |>
  select(condition, donor, fraction, carprofen, moi, timepoint, norm) |>
  pivot_wider(names_from = timepoint, values_from = norm,
              names_prefix = "t") |>
  mutate(fold_change = t7 / t4) |>
  filter(!is.na(fold_change))

# Summary
bcg_summ <- bcg_fc |>
  group_by(condition, fraction, carprofen, moi) |>
  summarise(
    mean_fc = mean(fold_change),
    sem_fc  = sd(fold_change) / sqrt(n()),
    vals    = list(fold_change),
    .groups = "drop"
  )

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 — BCG growth assay (2 rows × 4 MOI columns)
# ─────────────────────────────────────────────────────────────────────────────
# Condition colours
cond_colors <- c("Media only" = "#5B9BD5", "+Carprofen" = "#E06060")

# Factorise
bcg_fc_plot <- bcg_fc |>
  mutate(
    condition  = factor(condition, levels = c("Media only", "+Carprofen")),
    fraction   = factor(fraction,  levels = c("EC", "IC")),
    moi_label  = factor(paste("MOI", moi), levels = paste("MOI", c(0, 1, 10, 100))),
    carp_fac   = factor(carprofen, levels = c(0, 1, 10, 100))
  )

bcg_summ_plot <- bcg_summ |>
  mutate(
    condition  = factor(condition, levels = c("Media only", "+Carprofen")),
    fraction   = factor(fraction,  levels = c("EC", "IC")),
    moi_label  = factor(paste("MOI", moi), levels = paste("MOI", c(0, 1, 10, 100))),
    carp_fac   = factor(carprofen, levels = c(0, 1, 10, 100))
  )

# Unnest individual points for dot overlay
bcg_dots <- bcg_summ_plot |>
  select(condition, fraction, moi_label, carp_fac, vals) |>
  unnest(vals)

fig4 <- ggplot(bcg_summ_plot,
               aes(x = carp_fac, y = mean_fc, fill = condition)) +
  # Bars
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65, colour = NA) +
  # Error bars
  geom_errorbar(aes(ymin = mean_fc - sem_fc,
                    ymax = mean_fc + sem_fc),
                position = position_dodge(width = 0.75),
                width = 0.25, linewidth = 0.6) +
  # Individual donor dots
  geom_point(data = bcg_dots,
             aes(x = carp_fac, y = vals, fill = condition),
             position = position_jitterdodge(jitter.width = 0.12,
                                             dodge.width  = 0.75),
             shape = 21, size = 1.4, colour = "black",
             stroke = 0.3, alpha = 0.85) +
  # Reference line at fold change = 1
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey60", linewidth = 0.5) +
  scale_fill_manual(values = cond_colors, name = NULL) +
  scale_x_discrete(labels = c("0","1","10","100")) +
  # Independent y-scale per facet
  facet_grid(fraction ~ moi_label, scales = "free_y") +
  labs(x = "Carprofen (µg/mL)",
       y = "Fold change (normalised to 4 hpi)") +
  theme_report() +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 8.5, face = "plain"),
    legend.position  = "bottom",
    panel.spacing    = unit(0.5, "lines")
  )

ggsave("Figure4_BCG_growth_assay.pdf", fig4,
       width = 7.5, height = 4.8, units = "in")
ggsave("Figure4_BCG_growth_assay.png", fig4,
       width = 7.5, height = 4.8, units = "in", dpi = 300)
message("Saved Figure 4")

message("\nDone — all 3 figures saved.")
