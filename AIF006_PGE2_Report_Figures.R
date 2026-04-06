###############################################
# PGE2 Competitive ELISA - Report Figures
# Fixed column offset, clean analysis
###############################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(cowplot)

base_dir <- "/home/user/elisa"
export_dir <- file.path(base_dir, "AIF006_PGE2_Results")
dir.create(export_dir, showWarnings = FALSE)

#####
# 1. Parse Layout - FIXED column offset
# Layout CSV structure per data row:
#   V1=row_letter, V2=col1_welltype, V3=col2_stdconc, V4=col3_stdconc,
#   V5=col4_sample, V6=col5_sample, ...
# So plate column j => CSV column (j+1) => block[i, j+1]
#####
layout_raw <- read.csv(file.path(base_dir, "AIF006_ELISA_PGE2_Analysis(Layout).csv"),
                       header = FALSE, stringsAsFactors = FALSE)

parse_layout_block <- function(layout_raw, start_row, plate_name) {
  block <- layout_raw[start_row:(start_row + 7), ]
  result <- data.frame()
  
  for (i in 1:8) {
    row_letter <- LETTERS[i]
    col1_type <- trimws(as.character(block[i, 2]))  # V2 = plate col 1 well type
    
    for (j in 1:12) {
      well <- paste0(row_letter, j)
      # FIXED: plate col j content is at CSV column (j+1)
      content <- trimws(as.character(block[i, j + 1]))
      
      if (j == 1) {
        role <- col1_type
        std_conc <- NA
      } else if (j %in% c(2, 3)) {
        role <- "Standard"
        std_conc <- suppressWarnings(as.numeric(content))
      } else {
        role <- "Sample"
        std_conc <- NA
      }
      
      result <- rbind(result, data.frame(
        Plate = plate_name, Row = row_letter, Col = j, Well = well,
        Role = role, Content = content, StdConc = std_conc,
        stringsAsFactors = FALSE
      ))
    }
  }
  return(result)
}

layout_p2 <- parse_layout_block(layout_raw, 22, "Plate2")
layout_p3 <- parse_layout_block(layout_raw, 32, "Plate3")
all_layout <- rbind(layout_p2, layout_p3)

#####
# 2. Parse OD data (405nm only)
#####
parse_od_file <- function(filepath, plate_name) {
  raw <- read.csv(filepath, header = FALSE, stringsAsFactors = FALSE)
  result <- data.frame()
  for (letter in LETTERS[1:8]) {
    row_idx <- which(raw[, 2] == letter)
    if (length(row_idx) == 0) next
    row_405 <- row_idx[1]
    od405 <- suppressWarnings(as.numeric(raw[row_405, 3:14]))
    for (j in 1:12) {
      result <- rbind(result, data.frame(
        Plate = plate_name, Row = letter, Col = j,
        Well = paste0(letter, j), OD405 = od405[j],
        stringsAsFactors = FALSE
      ))
    }
  }
  return(result)
}

od_p2 <- parse_od_file(file.path(base_dir, "20260401_AIF006_ELISA_PGE2_Plate2(Plate 1 - Sheet1).csv"), "Plate2")
od_p3 <- parse_od_file(file.path(base_dir, "20260401_AIF006_ELISA_PGE2_Plate3(Plate 2 - Sheet1).csv"), "Plate3")
all_od <- rbind(od_p2, od_p3)

merged <- merge(all_layout, all_od, by = c("Plate", "Row", "Col", "Well"))

# Verify layout-OD alignment
cat("=== Verification ===\n")
cat("Standard concentrations at Col 2:\n")
print(merged %>% filter(Plate == "Plate2", Col == 2) %>% select(Well, Role, Content, StdConc, OD405))
cat("\nStandard concentrations at Col 3:\n")
print(merged %>% filter(Plate == "Plate2", Col == 3) %>% select(Well, Role, Content, StdConc, OD405))
cat("\nFirst few samples (Col 4):\n")
print(merged %>% filter(Plate == "Plate2", Col == 4) %>% select(Well, Role, Content, OD405))

#####
# 3. 4PL model functions
#####
fourPL <- function(x, a, d, c_val, b) {
  d + (a - d) / (1 + (x / c_val)^b)
}

inverse_4pl <- function(y, a, d, c_val, b) {
  if (y >= a || y <= d) return(NA)
  ratio <- (a - d) / (y - d) - 1
  if (ratio <= 0) return(NA)
  return(c_val * ratio^(1/b))
}

#####
# 4. Process each plate
#####
process_plate <- function(plate_data, plate_name) {
  cat("\n========================================\n")
  cat("Processing:", plate_name, "\n")
  cat("========================================\n")
  
  blk_od <- plate_data$OD405[plate_data$Role == "Blk"]
  nsb_od <- plate_data$OD405[plate_data$Role == "NSB"]
  b0_od  <- plate_data$OD405[plate_data$Role == "B0"]
  
  blk_mean <- mean(blk_od, na.rm = TRUE)
  nsb_mean <- mean(nsb_od, na.rm = TRUE)
  b0_mean  <- mean(b0_od, na.rm = TRUE)
  
  cat("Blank mean:", round(blk_mean, 4), "\n")
  cat("NSB mean:", round(nsb_mean, 4), "\n")
  cat("B0 mean:", round(b0_mean, 4), "\n")
  
  plate_data$pct_B_B0 <- ((plate_data$OD405 - nsb_mean) / (b0_mean - nsb_mean)) * 100
  
  stds <- plate_data %>% filter(Role == "Standard", !is.na(StdConc), StdConc > 0)
  
  std_summary <- stds %>%
    group_by(StdConc) %>%
    summarise(
      mean_OD = mean(OD405, na.rm = TRUE),
      sd_OD = sd(OD405, na.rm = TRUE),
      mean_pctBB0 = mean(pct_B_B0, na.rm = TRUE),
      sd_pctBB0 = sd(pct_B_B0, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(StdConc)
  
  cat("\nStandards (n per conc):\n")
  print(as.data.frame(std_summary))
  
  std_fit <- std_summary %>% filter(mean_pctBB0 > 0, mean_pctBB0 < 110)
  
  model <- tryCatch({
    nls(mean_pctBB0 ~ d + (a - d) / (1 + (StdConc / c_val)^b),
        data = std_fit,
        start = list(a = 100, d = 5, c_val = 80, b = 1.0),
        algorithm = "port",
        lower = c(a = 60, d = -20, c_val = 5, b = 0.1),
        upper = c(a = 150, d = 40, c_val = 2000, b = 5),
        control = nls.control(maxiter = 1000))
  }, error = function(e) {
    cat("NLS error:", e$message, "\n")
    NULL
  })
  
  if (is.null(model)) return(NULL)
  
  coefs <- coef(model)
  a <- coefs["a"]; d <- coefs["d"]; c_val <- coefs["c_val"]; b <- coefs["b"]
  
  ss_res <- sum(residuals(model)^2)
  ss_tot <- sum((std_fit$mean_pctBB0 - mean(std_fit$mean_pctBB0))^2)
  r_squared <- 1 - ss_res / ss_tot
  
  cat("\n4PL: a=", round(a,2), "d=", round(d,2), "c=", round(c_val,2), "b=", round(b,2), "R2=", round(r_squared,4), "\n")
  
  # Predict sample concentrations
  samples <- plate_data %>% filter(Role == "Sample")
  samples$predicted_conc <- sapply(samples$pct_B_B0, function(y) inverse_4pl(y, a, d, c_val, b))
  
  min_std <- min(std_fit$StdConc)
  max_std <- max(std_fit$StdConc)
  samples$result_flag <- case_when(
    is.na(samples$predicted_conc) & samples$pct_B_B0 >= a ~ paste("<", round(min_std, 2)),
    is.na(samples$predicted_conc) & samples$pct_B_B0 <= d ~ paste(">", round(max_std, 2)),
    TRUE ~ "OK"
  )
  
  # Standard curve for plotting
  conc_range <- exp(seq(log(min(std_fit$StdConc) * 0.5), log(max(std_fit$StdConc) * 2), length.out = 300))
  curve_data <- data.frame(Concentration = conc_range, pctBB0 = fourPL(conc_range, a, d, c_val, b))
  
  p_std <- ggplot() +
    geom_line(data = curve_data, aes(x = Concentration, y = pctBB0), linewidth = 0.9) +
    geom_point(data = std_summary, aes(x = StdConc, y = mean_pctBB0), size = 3) +
    geom_errorbar(data = std_summary %>% filter(!is.na(sd_pctBB0)),
                  aes(x = StdConc, ymin = mean_pctBB0 - sd_pctBB0, ymax = mean_pctBB0 + sd_pctBB0),
                  width = 0.05) +
    scale_x_log10() +
    labs(title = paste("PGE2 Standard Curve -", plate_name),
         subtitle = bquote(paste("4PL | ", R^2, " = ", .(round(r_squared, 4)))),
         x = "PGE2 (pg/mL)", y = "%B/B0") +
    theme_minimal(base_size = 14) +
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          panel.grid = element_blank(),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(face = "bold"))
  
  cat("\n=== Sample Results ===\n")
  print(as.data.frame(samples %>% select(Well, Content, OD405, pct_B_B0, predicted_conc, result_flag) %>% arrange(Content)))
  
  return(list(model = model, coefs = coefs, r_squared = r_squared,
              std_summary = std_summary, samples = samples, p_std = p_std,
              nsb_mean = nsb_mean, b0_mean = b0_mean))
}

res_p2 <- process_plate(merged %>% filter(Plate == "Plate2"), "Plate2")
res_p3 <- process_plate(merged %>% filter(Plate == "Plate3"), "Plate3")

# Save standard curves
ggsave(file.path(export_dir, "PGE2_Plate2_std_curve_FIXED.png"), res_p2$p_std, width = 15, height = 12.5, units = "cm", dpi = 300)
ggsave(file.path(export_dir, "PGE2_Plate3_std_curve_FIXED.png"), res_p3$p_std, width = 15, height = 12.5, units = "cm", dpi = 300)

#####
# 5. Combine samples from Plate2 (Exp1/2) and Plate3 (Exp3/4)
#####
combine_samples <- function(res, plate_name) {
  s <- res$samples %>%
    filter(result_flag == "OK") %>%
    select(Well, Content, OD405, pct_B_B0, predicted_conc) %>%
    mutate(Plate = plate_name)
  
  # Parse sample info: "Exp1 C0 B0 (1/1)" or "Exp3 DMSO 0 (1/1)"
  s <- s %>%
    mutate(
      Experiment = gsub("^(Exp[0-9]+)\\s.*", "\\1", Content),
      Drug = ifelse(grepl("DMSO", Content), "DMSO",
                    gsub(".*\\s(C[0-9]+)\\s.*", "\\1", Content)),
      Stimulus_type = case_when(
        grepl("DMSO", Content) ~ "DMSO",
        grepl("\\sB[0-9]", Content) ~ "BCG",
        grepl("\\sL[0-9]", Content) ~ "LPS",
        TRUE ~ "Unknown"
      ),
      Stimulus_dose = case_when(
        grepl("DMSO", Content) ~ gsub(".*DMSO\\s+([0-9]+).*", "\\1", Content),
        grepl("\\sB([0-9]+)", Content) ~ gsub(".*\\s[BL]([0-9]+).*", "\\1", Content),
        grepl("\\sL([0-9]+)", Content) ~ gsub(".*\\s[BL]([0-9]+).*", "\\1", Content),
        TRUE ~ NA_character_
      ),
      Stimulus_dose_num = as.numeric(Stimulus_dose),
      Carprofen_conc = case_when(
        Drug == "DMSO" ~ 0,
        Drug == "C0" ~ 0,
        Drug == "C1" ~ 1,
        Drug == "C10" ~ 10,
        Drug == "C100" ~ 100,
        TRUE ~ NA_real_
      )
    )
  return(s)
}

all_s <- rbind(combine_samples(res_p2, "Plate2"), combine_samples(res_p3, "Plate3"))

cat("\n\n========================================\n")
cat("COMBINED DATA (Plate2 + Plate3)\n")
cat("========================================\n")

# Show LPS conditions
cat("\n--- LPS Stimulation ---\n")
lps_data <- all_s %>% filter(Stimulus_type == "LPS") %>%
  arrange(Experiment, Carprofen_conc, Stimulus_dose_num)
print(as.data.frame(lps_data %>% select(Experiment, Drug, Stimulus_dose, predicted_conc)))

# Show BCG conditions
cat("\n--- BCG Stimulation ---\n")
bcg_data <- all_s %>% filter(Stimulus_type == "BCG") %>%
  arrange(Experiment, Carprofen_conc, Stimulus_dose_num)
print(as.data.frame(bcg_data %>% select(Experiment, Drug, Stimulus_dose, predicted_conc)))

# Show DMSO controls
cat("\n--- DMSO Controls ---\n")
dmso_data <- all_s %>% filter(Stimulus_type == "DMSO")
print(as.data.frame(dmso_data %>% select(Experiment, Drug, Stimulus_dose, predicted_conc)))

write.csv(all_s, file.path(export_dir, "PGE2_all_samples_FIXED.csv"), row.names = FALSE)

#####
# 6. Report Figure: PGE2 in LPS-stimulated MDMs
#####
lps_data$Carprofen_label <- factor(paste0("C", lps_data$Carprofen_conc),
                                    levels = c("C0", "C1", "C10", "C100"))
lps_data$LPS_dose <- factor(paste0("L", lps_data$Stimulus_dose),
                             levels = c("L0", "L1", "L10", "L100"))

# Panel A: Individual donors
lps_summary <- lps_data %>%
  group_by(Carprofen_label, LPS_dose) %>%
  summarise(
    mean_conc = mean(predicted_conc, na.rm = TRUE),
    sd_conc = sd(predicted_conc, na.rm = TRUE),
    sem_conc = sd(predicted_conc, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

cat("\n--- LPS Summary (mean ± SEM) ---\n")
print(as.data.frame(lps_summary))

p_lps <- ggplot(lps_summary, aes(x = LPS_dose, y = mean_conc, fill = Carprofen_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(0, mean_conc - sem_conc), ymax = mean_conc + sem_conc),
                position = position_dodge(width = 0.8), width = 0.25) +
  geom_point(data = lps_data, aes(x = LPS_dose, y = predicted_conc, group = Carprofen_label),
             position = position_dodge(width = 0.8), size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = c("C0" = "#4DAF4A", "C1" = "#377EB8", "C10" = "#FF7F00", "C100" = "#E41A1C"),
                    labels = c("0", "1", "10", "100")) +
  labs(title = "PGE2 - LPS Stimulation",
       x = "LPS (ng/mL)", y = expression("PGE"[2]~"(pg/mL)"),
       fill = "Carprofen\n(\u00b5g/mL)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(2, "mm"),
    plot.title = element_text(face = "bold", size = 16),
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

#####
# 7. Report Figure: PGE2 in BCG-stimulated MDMs
#####
bcg_data$Carprofen_label <- factor(paste0("C", bcg_data$Carprofen_conc),
                                    levels = c("C0", "C1", "C10", "C100"))
bcg_data$BCG_MOI <- factor(paste0("B", bcg_data$Stimulus_dose),
                            levels = c("B0", "B1", "B10", "B100"))

bcg_summary <- bcg_data %>%
  group_by(Carprofen_label, BCG_MOI) %>%
  summarise(
    mean_conc = mean(predicted_conc, na.rm = TRUE),
    sd_conc = sd(predicted_conc, na.rm = TRUE),
    sem_conc = sd(predicted_conc, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

cat("\n--- BCG Summary (mean ± SEM) ---\n")
print(as.data.frame(bcg_summary))

p_bcg <- ggplot(bcg_summary, aes(x = BCG_MOI, y = mean_conc, fill = Carprofen_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(0, mean_conc - sem_conc), ymax = mean_conc + sem_conc),
                position = position_dodge(width = 0.8), width = 0.25) +
  geom_point(data = bcg_data, aes(x = BCG_MOI, y = predicted_conc, group = Carprofen_label),
             position = position_dodge(width = 0.8), size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = c("C0" = "#4DAF4A", "C1" = "#377EB8", "C10" = "#FF7F00", "C100" = "#E41A1C"),
                    labels = c("0", "1", "10", "100")) +
  labs(title = "PGE2 - BCG Stimulation",
       x = "BCG (MOI)", y = expression("PGE"[2]~"(pg/mL)"),
       fill = "Carprofen\n(\u00b5g/mL)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(2, "mm"),
    plot.title = element_text(face = "bold", size = 16),
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

#####
# 8. Combined figure
#####
p_combined <- plot_grid(p_lps, p_bcg, ncol = 2, labels = c("A", "B"), label_size = 18, align = "h")

ggsave(file.path(export_dir, "Figure5_PGE2_LPS_BCG.png"),
       plot = p_combined, width = 30, height = 14, units = "cm", dpi = 300)
ggsave(file.path(export_dir, "Figure5_PGE2_LPS.png"),
       plot = p_lps, width = 15, height = 12.5, units = "cm", dpi = 300)
ggsave(file.path(export_dir, "Figure5_PGE2_BCG.png"),
       plot = p_bcg, width = 15, height = 12.5, units = "cm", dpi = 300)

#####
# 9. Two-way ANOVA
#####
cat("\n\n========================================\n")
cat("STATISTICAL ANALYSIS\n")
cat("========================================\n")

cat("\n--- Two-way ANOVA: PGE2 ~ Carprofen * LPS ---\n")
lps_anova_data <- lps_data %>%
  mutate(Carprofen_f = factor(Carprofen_conc),
         LPS_f = factor(Stimulus_dose_num))
lps_aov <- aov(predicted_conc ~ Carprofen_f * LPS_f, data = lps_anova_data)
print(summary(lps_aov))

cat("\n--- Two-way ANOVA: PGE2 ~ Carprofen * BCG ---\n")
bcg_anova_data <- bcg_data %>%
  mutate(Carprofen_f = factor(Carprofen_conc),
         BCG_f = factor(Stimulus_dose_num))
bcg_aov <- aov(predicted_conc ~ Carprofen_f * BCG_f, data = bcg_anova_data)
print(summary(bcg_aov))

# Save ANOVA results
sink(file.path(export_dir, "PGE2_ANOVA_results.txt"))
cat("Two-way ANOVA: PGE2 ~ Carprofen * LPS\n")
cat("========================================\n")
print(summary(lps_aov))
cat("\n\nTwo-way ANOVA: PGE2 ~ Carprofen * BCG\n")
cat("========================================\n")
print(summary(bcg_aov))
sink()

# Save summary tables
write.csv(lps_summary, file.path(export_dir, "PGE2_LPS_summary.csv"), row.names = FALSE)
write.csv(bcg_summary, file.path(export_dir, "PGE2_BCG_summary.csv"), row.names = FALSE)

cat("\n\nDone! Report figures saved.\n")
cat(list.files(export_dir, pattern = "FIXED|Figure5|ANOVA|summary"), sep = "\n")
