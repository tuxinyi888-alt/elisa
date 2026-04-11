# ==============================================================================
# LPS Dose-Response Figure — 4-Panel Cytokine Plot
# Stratified by Carprofen Dose
# ==============================================================================
#
# Data source : ELISA.xlsx  (sheets: IL1b, IL6, TNFa)
# Columns     : Donor, Carprofen (C0/C1/C10/C100),
#               Stimulus (LPS/BCG), Stim_Conc (L0/L1/L10/L100),
#               <cytokine>_pgml
#
# NOTE — PGE2 concentrations have not yet been calculated from the raw
#   plate-reader OD files.  Once you have processed PGE2 data, add a "PGE2"
#   sheet to ELISA.xlsx with the same column layout
#   (Donor, Carprofen, Stimulus, Stim_Conc, PGE2_pgml)
#   and uncomment the four lines marked "## PGE2" below.
# ==============================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)

# ---------- 1. Read each cytokine sheet ----------
il1b_raw <- read_excel("ELISA.xlsx", sheet = "IL1b")
il6_raw  <- read_excel("ELISA.xlsx", sheet = "IL6")
tnfa_raw <- read_excel("ELISA.xlsx", sheet = "TNFa")
## PGE2: pge2_raw <- read_excel("ELISA.xlsx", sheet = "PGE2")

# ---------- 2. Helper: decode coded factor levels to numeric ----------
decode_carprofen <- function(x) as.numeric(gsub("C", "", x))
decode_lps_dose  <- function(x) as.numeric(gsub("L", "", x))

# ---------- 3. Filter to LPS, decode factors, standardise columns ----------
clean_lps <- function(df, value_col, cytokine_label) {
  df %>%
    filter(Stimulus == "LPS") %>%
    mutate(
      Carprofen_uM = decode_carprofen(Carprofen),
      LPS_ngml     = decode_lps_dose(Stim_Conc),
      Cytokine     = cytokine_label,
      Conc_pgml    = .data[[value_col]]
    ) %>%
    select(Donor, Carprofen_uM, LPS_ngml, Cytokine, Conc_pgml)
}

il1b <- clean_lps(il1b_raw, "IL1b_pgml", "IL-1\u03b2")
il6  <- clean_lps(il6_raw,  "IL6_pgml",  "IL-6")
tnfa <- clean_lps(tnfa_raw, "TNFa_pgml", "TNF\u03b1")
## PGE2: pge2 <- clean_lps(pge2_raw, "PGE2_pgml", "PGE2")

# ---------- 4. Combine into one long data frame ----------
dat <- bind_rows(il1b, il6, tnfa)
## PGE2: dat <- bind_rows(il1b, il6, tnfa, pge2)

# Ordered factors for clean axis labels
dat <- dat %>%
  mutate(
    LPS_ngml     = factor(LPS_ngml,     levels = c(0, 1, 10, 100)),
    Carprofen_uM = factor(Carprofen_uM, levels = c(0, 1, 10, 100)),
    Cytokine     = factor(Cytokine,     levels = c("IL-1\u03b2", "IL-6",
                                                    "TNF\u03b1", "PGE2"))
  )

# ---------- 5. Summary statistics (mean ± SEM) per group ----------
dat_summary <- dat %>%
  group_by(Cytokine, LPS_ngml, Carprofen_uM) %>%
  summarise(
    mean  = mean(Conc_pgml, na.rm = TRUE),
    sem   = sd(Conc_pgml, na.rm = TRUE) / sqrt(sum(!is.na(Conc_pgml))),
    n     = sum(!is.na(Conc_pgml)),
    .groups = "drop"
  )

# ---------- 6. Build the figure ----------
p <- ggplot(dat_summary,
            aes(x = LPS_ngml, y = mean,
                colour = Carprofen_uM, group = Carprofen_uM)) +

  # Mean ± SEM lines
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0.15, linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +

  # Individual donor points (semi-transparent, slightly jittered)
  geom_point(data = dat,
             aes(x = LPS_ngml, y = Conc_pgml, colour = Carprofen_uM),
             alpha = 0.3, size = 1.2,
             position = position_dodge(width = 0.3),
             show.legend = FALSE) +

  # One panel per cytokine, free y-axis so scales suit each analyte
  facet_wrap(~ Cytokine, scales = "free_y", ncol = 2) +

  # Colour-blind-friendly palette (Dark2)
  scale_colour_brewer(palette = "Dark2", name = "Carprofen (\u00b5M)") +

  labs(x = "LPS dose (ng/mL)",
       y = "Concentration (pg/mL)") +

  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 11),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# ---------- 7. Save ----------
ggsave("LPS_cytokine_4panel.png", p, width = 8, height = 7, dpi = 300)
ggsave("LPS_cytokine_4panel.pdf", p, width = 8, height = 7)

cat("Done — saved LPS_cytokine_4panel.png and .pdf\n")
