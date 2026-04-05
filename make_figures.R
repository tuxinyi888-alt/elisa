# ─────────────────────────────────────────────────────────────────────────────
# make_figures.R
# Figures 2, 3, 4 for carprofen / MDM report
#
# Install once:  install.packages(c("tidyverse","readxl","patchwork"))
# Run:           source("make_figures.R")    (working dir = repo root)
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(patchwork)
  library(scales)
})

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_report <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line          = element_line(linewidth = 0.45),
      axis.ticks         = element_line(linewidth = 0.35),
      axis.text          = element_text(size = base_size - 0.5),
      axis.title         = element_text(size = base_size),
      plot.title         = element_text(size = base_size + 0.5, hjust = 0.5,
                                        face = "plain"),
      plot.tag           = element_text(size = 10, face = "bold"),
      legend.title       = element_text(size = base_size - 0.5, face = "bold"),
      legend.text        = element_text(size = base_size - 0.5),
      legend.key.size    = unit(0.35, "cm"),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
      plot.margin        = margin(4, 6, 2, 4)
    )
}

# ── Colour palettes ──────────────────────────────────────────────────────────
# LPS doses:  warm colours (light → dark red)
lps_pal  <- c("0"="#D9D9D9", "1"="#FDBB84", "10"="#E34A33", "100"="#7F0000")
# BCG doses:  cool colours (light → dark blue)
bcg_pal  <- c("0"="#D9D9D9", "1"="#9ECAE1", "10"="#3182BD", "100"="#08306B")
# BCG MOI for Figure 4
moi_pal  <- c("0"="#D9D9D9", "1"="#9ECAE1", "10"="#3182BD", "100"="#08306B")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load & tidy ELISA data
# ─────────────────────────────────────────────────────────────────────────────
read_elisa_sheet <- function(path, sheet, col_name) {
  read_excel(path, sheet = sheet) |>
    select(1:5) |>
    setNames(c("Donor","Carprofen","Stimulus","Stim_Conc","value")) |>
    mutate(
      cytokine  = col_name,
      value     = suppressWarnings(as.numeric(value)),
      carp_num  = as.numeric(str_remove(Carprofen, "C")),
      stim_num  = as.numeric(str_remove(Stim_Conc, "[LBb]")),
      # 0 µg/mL carprofen → 0.5 for log-scale positioning
      carp_x    = if_else(carp_num == 0, 0.5, as.numeric(carp_num)),
      stim_fac  = factor(stim_num)
    )
}

elisa <- bind_rows(
  read_elisa_sheet("ELISA.xlsx", "IL1b", "IL-1\u03b2 (pg/mL)"),
  read_elisa_sheet("ELISA.xlsx", "IL6",  "IL-6 (pg/mL)"),
  read_elisa_sheet("ELISA.xlsx", "TNFa", "TNF-\u03b1 (pg/mL)")
) |>
  mutate(cytokine = factor(cytokine,
                           levels = c("IL-1\u03b2 (pg/mL)",
                                      "IL-6 (pg/mL)",
                                      "TNF-\u03b1 (pg/mL)")))

# ─────────────────────────────────────────────────────────────────────────────
# 2. Two-way ANOVA helper (returns formatted string for annotation)
# ─────────────────────────────────────────────────────────────────────────────
run_anova <- function(df, cytokine_name, stimulus) {
  d <- df |>
    filter(cytokine == cytokine_name, Stimulus == stimulus, !is.na(value))
  if (nrow(d) < 4) return("")
  tryCatch({
    m  <- aov(value ~ factor(carp_num) * factor(stim_num), data = d)
    s  <- summary(m)[[1]]
    rn <- rownames(s)
    fmt <- function(p) {
      if (is.na(p))    return("n.s.")
      if (p < 0.0001) return("p < 0.0001")
      if (p < 0.001)  return("p < 0.001")
      sprintf("p = %.3f", p)
    }
    p_carp <- fmt(s[grep("carp", rn)[1], "Pr(>F)"])
    p_stim <- fmt(s[grep("stim", rn)[1], "Pr(>F)"])
    p_int  <- fmt(s[grep(":",    rn)[1], "Pr(>F)"])
    paste0("2-way ANOVA\n",
           "Carprofen: ",   p_carp, "\n",
           "Stimulus: ",    p_stim, "\n",
           "Interaction: ", p_int)
  }, error = function(e) "")
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. Single ELISA panel (one donor or average)
# ─────────────────────────────────────────────────────────────────────────────
elisa_panel <- function(df, cytokine_name, stimulus,
                        donor       = NULL,   # NULL = all donors (average)
                        dose_pal    = lps_pal,
                        dose_title  = "LPS (ng/mL)",
                        show_legend = FALSE,
                        show_ylabel = TRUE,
                        show_xlabel = TRUE,
                        anova_text  = "") {

  d <- df |>
    filter(cytokine == cytokine_name, Stimulus == stimulus)
  if (!is.null(donor)) d <- filter(d, Donor == donor)

  # Drop rows with NA or non-positive value (can't log-transform)
  d <- filter(d, !is.na(value), value > 0)

  # Title
  ttl <- if (is.null(donor)) "Average" else donor

  # x / y limits
  x_breaks <- c(0.5, 1, 10, 100)
  x_labels <- c("0", "1", "10", "100")

  p <- ggplot(d, aes(x = carp_x, y = value,
                     colour = stim_fac, fill = stim_fac)) +
    # Individual data points
    geom_point(size = 1.0, alpha = 0.55, shape = 16) +
    # Smooth fit (lm in log-log space = power-law dose-response)
    geom_smooth(method  = "lm",
                formula = y ~ x,
                se      = TRUE,
                alpha   = 0.12,
                linewidth = 0.75) +
    scale_x_log10(breaks = x_breaks, labels = x_labels,
                  limits = c(0.4, 130)) +
    scale_y_log10(labels = label_comma()) +
    scale_colour_manual(values = dose_pal, name = dose_title) +
    scale_fill_manual(  values = dose_pal, name = dose_title, guide = "none") +
    labs(title = ttl,
         y = if (show_ylabel) cytokine_name else NULL,
         x = if (show_xlabel) "Carprofen (\u00b5g/mL)" else NULL) +
    theme_report()

  # 2-way ANOVA annotation (average panel only)
  if (nchar(anova_text) > 0) {
    p <- p + annotate("text",
                      x = 0.45, y = Inf,
                      label = anova_text,
                      hjust = 0, vjust = 1.25,
                      size  = 2.0, colour = "grey30",
                      fontface = "plain")
  }

  if (!show_legend) p <- p + theme(legend.position = "none")

  p
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. Build full ELISA figure  (3 cytokines × 5 panels)
# ─────────────────────────────────────────────────────────────────────────────
make_elisa_figure <- function(stimulus,
                               dose_pal,
                               dose_title,
                               fig_title,
                               filename_stem) {

  cytokines_vec <- levels(elisa$cytokine)
  donors_vec    <- c("D1","D2","D3","D4")
  col_names     <- c(donors_vec, "Average")

  # Precompute ANOVA strings (one per cytokine)
  anova_texts <- setNames(
    lapply(cytokines_vec, run_anova,
           df = elisa, stimulus = stimulus),
    cytokines_vec
  )

  # Build 3×5 panel matrix
  panel_list <- list()
  for (ri in seq_along(cytokines_vec)) {
    cyt <- cytokines_vec[ri]
    for (ci in seq_along(col_names)) {
      donor      <- if (col_names[ci] == "Average") NULL else col_names[ci]
      is_avg     <- is.null(donor)
      show_y     <- (ci == 1)
      show_x     <- (ri == length(cytokines_vec))
      show_leg   <- (ci == length(col_names) && ri == 1)
      anova_ann  <- if (is_avg) anova_texts[[cyt]] else ""

      panel_list[[length(panel_list) + 1]] <-
        elisa_panel(
          df          = elisa,
          cytokine_name = cyt,
          stimulus    = stimulus,
          donor       = donor,
          dose_pal    = dose_pal,
          dose_title  = dose_title,
          show_legend = show_leg,
          show_ylabel = show_y,
          show_xlabel = show_x,
          anova_text  = anova_ann
        )
    }
  }

  # Arrange as 3 rows × 5 cols
  fig <- wrap_plots(panel_list, nrow = 3, ncol = 5,
                    guides = "collect") +
    plot_annotation(
      title     = fig_title,
      tag_levels = list(
        c("A","B","C","D","E",
          "F","G","H","I","J",
          "K","L","M","N","O")
      )
    ) &
    theme(plot.tag = element_text(size = 8, face = "bold"),
          legend.position = "right")

  ggsave(paste0(filename_stem, ".pdf"),
         fig, width = 13, height = 7.5, units = "in")
  ggsave(paste0(filename_stem, ".png"),
         fig, width = 13, height = 7.5, units = "in", dpi = 300)
  message("Saved ", filename_stem)
  invisible(fig)
}

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — ELISA LPS
# ─────────────────────────────────────────────────────────────────────────────
make_elisa_figure(
  stimulus      = "LPS",
  dose_pal      = lps_pal,
  dose_title    = "LPS (ng/mL)",
  fig_title     = "Figure 2 | ELISA: Effect of carprofen on LPS-induced cytokine production",
  filename_stem = "Figure2_ELISA_LPS"
)

# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — ELISA BCG
# ─────────────────────────────────────────────────────────────────────────────
make_elisa_figure(
  stimulus      = "BCG",
  dose_pal      = bcg_pal,
  dose_title    = "BCG (MOI)",
  fig_title     = "Figure 3 | ELISA: Effect of carprofen on BCG-induced cytokine production",
  filename_stem = "Figure3_ELISA_BCG"
)

# ─────────────────────────────────────────────────────────────────────────────
# 5. Load BCG flow-cytometry data (AIF007 xls files)
#    Files 1–3 = media only after 4 h (transient)
#    Files 4–6 = +carprofen continued (sustained)
# ─────────────────────────────────────────────────────────────────────────────
read_aif <- function(fnum) {
  cond  <- if (fnum <= 3) "Media only" else "+Carprofen"
  donor <- paste0("D", ((fnum - 1) %% 3) + 1)
  df <- read_excel(sprintf("AIF007-%d_Table.xls", fnum),
                   col_names = FALSE) |>
    setNames(c("sample","bead_pct","beads","nb_pct","nb_n",
               "bcg_pct","bcg_n","bcg_gmean","debris_pct"))
  df |>
    filter(str_detect(sample, "\\.fcs$")) |>
    transmute(
      sample    = str_remove(sample, "\\.fcs$"),
      condition = cond,
      donor     = donor,
      beads     = as.numeric(beads),
      bcg_n     = as.numeric(bcg_n)
    ) |>
    separate(sample, into = c("tp","fraction","cm"), sep = "_") |>
    mutate(
      timepoint = as.integer(str_extract(tp, "\\d+")),
      carprofen = as.numeric(str_extract(cm, "(?<=C)\\d+")),
      moi       = as.numeric(str_extract(cm, "(?<=B)\\d+")),
      carp_x    = if_else(carprofen == 0, 0.5, carprofen)
    ) |>
    select(condition, donor, timepoint, fraction, carprofen, carp_x, moi, beads, bcg_n)
}

fc_raw <- bind_rows(lapply(1:6, read_aif))

# Fold change = (7 dpi BCG/beads) / (4 hpi BCG/beads)
fc_df <- fc_raw |>
  mutate(norm = bcg_n / beads) |>
  pivot_wider(id_cols = c(condition, donor, fraction, carprofen, carp_x, moi),
              names_from = timepoint, names_prefix = "t",
              values_from = norm) |>
  mutate(fold_change = t7 / t4) |>
  filter(!is.na(fold_change)) |>
  mutate(
    condition = factor(condition, levels = c("Media only","+Carprofen")),
    fraction  = factor(fraction,  levels = c("EC","IC")),
    moi_fac   = factor(moi, levels = c(0,1,10,100))
  )

# ─────────────────────────────────────────────────────────────────────────────
# 6. Figure 4 — BCG growth assay
#    2 rows (EC, IC) × 2 cols (Media only, +Carprofen)
#    x = carprofen (log), y = fold change (linear), lines = MOI
# ─────────────────────────────────────────────────────────────────────────────
flow_panel <- function(frac, cond, show_legend = FALSE,
                       show_ylabel = TRUE, show_xlabel = TRUE) {
  d <- filter(fc_df, fraction == frac, condition == cond,
              !is.na(fold_change))

  p <- ggplot(d, aes(x = carp_x, y = fold_change,
                     colour = moi_fac, fill = moi_fac)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey50", linewidth = 0.5) +
    geom_point(size = 1.2, alpha = 0.55, shape = 16) +
    geom_smooth(method = "lm", formula = y ~ x,
                se = TRUE, alpha = 0.12, linewidth = 0.75) +
    scale_x_log10(breaks = c(0.5,1,10,100),
                  labels = c("0","1","10","100"),
                  limits = c(0.4, 130)) +
    scale_colour_manual(values = moi_pal, name = "BCG (MOI)") +
    scale_fill_manual(  values = moi_pal, name = "BCG (MOI)", guide = "none") +
    labs(
      title = as.character(cond),
      y = if (show_ylabel) paste0(frac, " fold change\n(normalised to 4 hpi)") else NULL,
      x = if (show_xlabel) "Carprofen (\u00b5g/mL)" else NULL
    ) +
    theme_report()

  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

fig4 <-
  (flow_panel("EC", "Media only",  show_xlabel = FALSE) +
   flow_panel("EC", "+Carprofen",  show_xlabel = FALSE, show_ylabel = FALSE,
              show_legend = TRUE)) /
  (flow_panel("IC", "Media only") +
   flow_panel("IC", "+Carprofen",  show_ylabel = FALSE)) +
  plot_annotation(
    title      = "Figure 4 | BCG growth assay: fold change in extracellular and intracellular BCG",
    tag_levels = "A"
  ) &
  theme(plot.tag       = element_text(size = 10, face = "bold"),
        legend.position = "right")

ggsave("Figure4_BCG_growth_assay.pdf",
       fig4, width = 7.5, height = 6.0, units = "in")
ggsave("Figure4_BCG_growth_assay.png",
       fig4, width = 7.5, height = 6.0, units = "in", dpi = 300)
message("Saved Figure 4")

message("\nDone — Figures 2, 3, 4 saved.")
