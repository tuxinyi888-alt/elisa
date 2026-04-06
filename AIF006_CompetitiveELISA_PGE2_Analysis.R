###############################################
# Competitive ELISA Analysis - PGE2 - 4PL
# AIF006 - All Plates
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
# 1. Parse Layout CSV
#####
# Layout structure per block (8 rows x 14 cols):
#   Col 1: Row letter (A-H)
#   Col 2: Well type for column 1 ONLY (Blk, NSB, B0, TA)
#   Col 3-4: Standard concentrations (duplicates) for columns 2-3
#   Col 5-14: Sample IDs for columns 4-12
# Well types apply ONLY to column 1:
#   A1,B1 = Blank; C1,D1 = NSB; E1,F1,G1 = B0; H1 = TA
# Columns 2-3 = Standards (all rows)
# Columns 4-12 = Samples

layout_raw <- read.csv(file.path(base_dir, "AIF006_ELISA_PGE2_Analysis(Layout).csv"),
                       header = FALSE, stringsAsFactors = FALSE)

parse_layout_block <- function(layout_raw, start_row, plate_name) {
  block <- layout_raw[start_row:(start_row + 7), ]
  result <- data.frame()
  
  for (i in 1:8) {
    row_letter <- LETTERS[i]
    col1_type <- trimws(as.character(block[i, 2]))
    
    for (j in 1:12) {
      well <- paste0(row_letter, j)
      content <- trimws(as.character(block[i, j + 2]))
      
      # Determine role of each well
      if (j == 1) {
        role <- col1_type  # Blk, NSB, B0, or TA
        std_conc <- NA
      } else if (j %in% c(2, 3)) {
        role <- "Standard"
        std_conc <- suppressWarnings(as.numeric(content))
      } else {
        role <- "Sample"
        std_conc <- NA
      }
      
      result <- rbind(result, data.frame(
        Plate = plate_name,
        Row = row_letter,
        Col = j,
        Well = well,
        Role = role,
        Content = content,
        StdConc = std_conc,
        stringsAsFactors = FALSE
      ))
    }
  }
  return(result)
}

# Parse 4 layout blocks
layout_p1_0306 <- parse_layout_block(layout_raw, 2, "Plate1_0306")
layout_p1_0308 <- parse_layout_block(layout_raw, 12, "Plate1_0308")
layout_p2 <- parse_layout_block(layout_raw, 22, "Plate2")
layout_p3 <- parse_layout_block(layout_raw, 32, "Plate3")
all_layout <- rbind(layout_p1_0306, layout_p1_0308, layout_p2, layout_p3)

#####
# 2. Parse OD data files
#####
parse_od_file <- function(filepath, plate_name) {
  raw <- read.csv(filepath, header = FALSE, stringsAsFactors = FALSE)
  
  # Count wavelengths
  wl_lines <- grep("Wavelengths", raw[,2])
  wl_text <- raw[wl_lines[1], 2]
  n_wl <- length(gregexpr("[0-9]+", wl_text)[[1]])
  
  result <- data.frame()
  for (letter in LETTERS[1:8]) {
    # Find the row with this letter
    row_idx <- which(raw[, 2] == letter)
    if (length(row_idx) == 0) next
    row_405 <- row_idx[1]
    
    od405 <- suppressWarnings(as.numeric(raw[row_405, 3:14]))
    
    for (j in 1:12) {
      result <- rbind(result, data.frame(
        Plate = plate_name,
        Row = letter,
        Col = j,
        Well = paste0(letter, j),
        OD405 = od405[j],
        stringsAsFactors = FALSE
      ))
    }
  }
  return(result)
}

od_p1_0306 <- parse_od_file(
  file.path(base_dir, "20260306_AIF006_ELISA_PGE2_Plate1(Plate 1 - Sheet1).csv"), "Plate1_0306")
od_p1_0308 <- parse_od_file(
  file.path(base_dir, "20260308_AIF006_ELISA_PGE2_Plate1(Plate 1 - Sheet1).csv"), "Plate1_0308")
od_p2 <- parse_od_file(
  file.path(base_dir, "20260401_AIF006_ELISA_PGE2_Plate2(Plate 1 - Sheet1).csv"), "Plate2")
od_p3 <- parse_od_file(
  file.path(base_dir, "20260401_AIF006_ELISA_PGE2_Plate3(Plate 2 - Sheet1).csv"), "Plate3")

all_od <- rbind(od_p1_0306, od_p1_0308, od_p2, od_p3)

# Merge
merged <- merge(all_layout, all_od, by = c("Plate", "Row", "Col", "Well"))

cat("Merged data:", nrow(merged), "wells\n")
cat("Roles:", paste(unique(merged$Role), collapse = ", "), "\n")

#####
# 3. Process each plate
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

process_plate <- function(plate_data, plate_name) {
  cat("\n========================================\n")
  cat("Processing:", plate_name, "\n")
  cat("========================================\n")
  
  # Control wells (column 1 only)
  blk_od <- plate_data$OD405[plate_data$Role == "Blk"]
  nsb_od <- plate_data$OD405[plate_data$Role == "NSB"]
  b0_od  <- plate_data$OD405[plate_data$Role == "B0"]
  ta_od  <- plate_data$OD405[plate_data$Role == "TA"]
  
  blk_mean <- mean(blk_od, na.rm = TRUE)
  nsb_mean <- mean(nsb_od, na.rm = TRUE)
  b0_mean  <- mean(b0_od, na.rm = TRUE)
  
  cat("Blank wells:", blk_od, "-> mean:", round(blk_mean, 4), "\n")
  cat("NSB wells:", nsb_od, "-> mean:", round(nsb_mean, 4), "\n")
  cat("B0 wells:", b0_od, "-> mean:", round(b0_mean, 4), "\n")
  if (length(ta_od) > 0) cat("TA wells:", ta_od, "\n")
  
  # %B/B0 = ((OD - NSB) / (B0 - NSB)) * 100
  plate_data$pct_B_B0 <- ((plate_data$OD405 - nsb_mean) / (b0_mean - nsb_mean)) * 100
  
  # Standards
  stds <- plate_data %>% filter(Role == "Standard", !is.na(StdConc), StdConc > 0)
  
  std_summary <- stds %>%
    group_by(StdConc) %>%
    summarise(
      mean_OD = mean(OD405, na.rm = TRUE),
      sd_OD = sd(OD405, na.rm = TRUE),
      mean_pctBB0 = mean(pct_B_B0, na.rm = TRUE),
      sd_pctBB0 = sd(pct_B_B0, na.rm = TRUE),
      n = n(),
      cv_pct = ifelse(n() > 1, (sd(OD405, na.rm = TRUE) / mean(OD405, na.rm = TRUE)) * 100, NA),
      .groups = "drop"
    ) %>%
    arrange(StdConc)
  
  cat("\nStandard Curve Data:\n")
  print(as.data.frame(std_summary))
  
  # Filter for fitting: remove if CV > 30% (for duplicates)
  std_fit <- std_summary %>%
    filter(!is.na(mean_pctBB0), mean_pctBB0 > 0, mean_pctBB0 < 110)
  
  if (any(!is.na(std_summary$cv_pct) & std_summary$cv_pct > 30)) {
    high_cv <- std_summary %>% filter(!is.na(cv_pct) & cv_pct > 30)
    cat("\nWARNING: Standards with CV > 30%:\n")
    print(as.data.frame(high_cv))
  }
  
  cat("\nStandards used for fitting:\n")
  print(as.data.frame(std_fit))
  
  # Fit 4PL: competitive ELISA -> a = upper (~100), d = lower (~0), b > 0
  a_start <- max(std_fit$mean_pctBB0)
  d_start <- min(std_fit$mean_pctBB0)
  c_start <- exp(mean(log(std_fit$StdConc)))
  b_start <- 1
  
  model <- tryCatch({
    nls(mean_pctBB0 ~ d + (a - d) / (1 + (StdConc / c_val)^b),
        data = std_fit,
        start = list(a = a_start, d = d_start, c_val = c_start, b = b_start),
        control = nls.control(maxiter = 1000, tol = 1e-8))
  }, error = function(e) {
    cat("First attempt failed:", e$message, "\n")
    cat("Trying with adjusted starting values...\n")
    tryCatch({
      nls(mean_pctBB0 ~ d + (a - d) / (1 + (StdConc / c_val)^b),
          data = std_fit,
          start = list(a = 100, d = 0, c_val = c_start, b = 1.5),
          control = nls.control(maxiter = 1000, tol = 1e-6))
    }, error = function(e2) {
      cat("ERROR fitting 4PL:", e2$message, "\n")
      return(NULL)
    })
  })
  
  if (is.null(model)) return(NULL)
  
  coefs <- coef(model)
  a <- coefs["a"]; d <- coefs["d"]; c_val <- coefs["c_val"]; b <- coefs["b"]
  
  ss_res <- sum(residuals(model)^2)
  ss_tot <- sum((std_fit$mean_pctBB0 - mean(std_fit$mean_pctBB0))^2)
  r_squared <- 1 - ss_res / ss_tot
  
  cat("\n4PL Model Parameters:\n")
  cat("  Upper asymptote (a):", round(a, 4), "\n")
  cat("  Lower asymptote (d):", round(d, 4), "\n")
  cat("  EC50 (c):", round(c_val, 4), "pg/mL\n")
  cat("  Slope (b):", round(b, 4), "\n")
  cat("  R-squared:", round(r_squared, 6), "\n")
  print(summary(model))
  
  # Smooth curve
  conc_range <- exp(seq(log(min(std_fit$StdConc) * 0.5),
                        log(max(std_fit$StdConc) * 2),
                        length.out = 300))
  curve_data <- data.frame(
    Concentration = conc_range,
    pctBB0 = fourPL(conc_range, a, d, c_val, b)
  )
  
  # ---- Standard Curve Plot ----
  std_summary$used_in_fit <- std_summary$StdConc %in% std_fit$StdConc
  
  p1 <- ggplot() +
    geom_line(data = curve_data, aes(x = Concentration, y = pctBB0),
              linewidth = 0.9, color = "black") +
    geom_point(data = std_summary,
               aes(x = StdConc, y = mean_pctBB0, color = used_in_fit),
               size = 3) +
    geom_errorbar(data = std_summary %>% filter(!is.na(sd_pctBB0)),
                  aes(x = StdConc,
                      ymin = mean_pctBB0 - sd_pctBB0,
                      ymax = mean_pctBB0 + sd_pctBB0),
                  width = 0.05, color = "black") +
    scale_x_log10() +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "#D81B60"),
                       labels = c("TRUE" = "Used", "FALSE" = "Excluded"),
                       name = "Fit Status") +
    labs(
      title = paste("PGE2 Competitive ELISA -", plate_name),
      subtitle = bquote(paste("4PL Fit | ", R^2, " = ", .(round(r_squared, 4)))),
      x = expression("PGE2 concentration (pg/mL)"),
      y = "%B/B0"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      aspect.ratio = 1,
      axis.title = element_text(size = 18),
      axis.text = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.ticks.length = unit(2, "mm"),
      axis.ticks = element_line(color = "black")
    )
  
  # ---- Predict sample concentrations ----
  samples <- plate_data %>% filter(Role == "Sample")
  
  if (nrow(samples) > 0) {
    samples$predicted_conc <- sapply(samples$pct_B_B0, function(y) {
      inverse_4pl(y, a, d, c_val, b)
    })
    
    min_std <- min(std_fit$StdConc)
    max_std <- max(std_fit$StdConc)
    
    samples$result_flag <- case_when(
      is.na(samples$predicted_conc) & samples$pct_B_B0 >= a ~ paste("<", round(min_std, 2)),
      is.na(samples$predicted_conc) & samples$pct_B_B0 <= d ~ paste(">", round(max_std, 2)),
      TRUE ~ "OK"
    )
    
    cat("\n=== Sample Results ===\n")
    print(as.data.frame(samples %>%
      select(Well, Content, OD405, pct_B_B0, predicted_conc, result_flag) %>%
      arrange(Content)))
    
    # Standard curve + samples plot
    samples_valid <- samples %>% filter(!is.na(predicted_conc))
    
    p2 <- p1 +
      geom_point(data = samples_valid,
                 aes(x = predicted_conc, y = pct_B_B0),
                 color = "orange2", size = 3, alpha = 0.7) +
      labs(title = paste("PGE2 -", plate_name, "- Standards + Samples"))
    
  } else {
    samples <- data.frame()
    p2 <- p1
  }
  
  # Save outputs
  ggsave(file.path(export_dir, paste0("PGE2_", plate_name, "_standard_curve.png")),
         plot = p1, width = 15, height = 12.5, units = "cm", dpi = 300)
  ggsave(file.path(export_dir, paste0("PGE2_", plate_name, "_std_curve_samples.png")),
         plot = p2, width = 15, height = 12.5, units = "cm", dpi = 300)
  
  write.csv(std_summary,
            file = file.path(export_dir, paste0("PGE2_", plate_name, "_standards.csv")),
            row.names = FALSE)
  
  if (nrow(samples) > 0) {
    write.csv(samples %>% select(Well, Content, OD405, pct_B_B0, predicted_conc, result_flag),
              file = file.path(export_dir, paste0("PGE2_", plate_name, "_samples.csv")),
              row.names = FALSE)
  }
  
  model_text <- c(
    paste("Plate:", plate_name),
    paste("Blank mean:", round(blk_mean, 4)),
    paste("NSB mean:", round(nsb_mean, 4)),
    paste("B0 mean:", round(b0_mean, 4)),
    "",
    "4PL Parameters: y = d + (a-d) / (1 + (x/c)^b)",
    paste("  a (upper):", round(a, 4)),
    paste("  d (lower):", round(d, 4)),
    paste("  c (EC50):", round(c_val, 4), "pg/mL"),
    paste("  b (slope):", round(b, 4)),
    paste("  R-squared:", round(r_squared, 6)),
    "",
    capture.output(summary(model))
  )
  writeLines(model_text, file.path(export_dir, paste0("PGE2_", plate_name, "_model.txt")))
  
  return(list(model = model, coefs = coefs, r_squared = r_squared,
              std_summary = std_summary, samples = samples,
              p1 = p1, p2 = p2, nsb_mean = nsb_mean, b0_mean = b0_mean))
}

#####
# Process plates
# Plate1_0306: first reading (20260306) - note: standards had preparation mistake
# Plate1_0308: second reading (20260308) - same plate, better developed
# Plate2: Exp1 & Exp2 samples (20260401)
# Plate3: Exp3 & Exp4 samples (20260401)
#####

# Use Plate1_0308 (second reading gives better signal development)
# Skip Plate1_0306 (first day reading, lower signal)
plates_to_process <- c("Plate1_0308", "Plate2", "Plate3")
results <- list()

for (plate in plates_to_process) {
  plate_data <- merged %>% filter(Plate == plate)
  results[[plate]] <- process_plate(plate_data, plate)
}

#####
# 4. Combined sample results
#####
cat("\n\n========================================\n")
cat("COMBINED RESULTS ACROSS ALL PLATES\n")
cat("========================================\n")

all_samples <- data.frame()
for (plate in plates_to_process) {
  if (!is.null(results[[plate]]) && nrow(results[[plate]]$samples) > 0) {
    s <- results[[plate]]$samples %>%
      select(Well, Content, OD405, pct_B_B0, predicted_conc, result_flag) %>%
      mutate(Plate = plate)
    all_samples <- rbind(all_samples, s)
  }
}

if (nrow(all_samples) > 0) {
  # Parse sample identifiers from Content
  # Format examples: "Exp3 C0 B0 (1/1)", "Exp3 DMSO 0 (1/1)"
  all_samples <- all_samples %>%
    mutate(
      Experiment = gsub("^(Exp[0-9]+)\\s.*", "\\1", Content),
      # Carprofen/treatment: C0, C1, C10, C100, DMSO
      Drug = ifelse(grepl("DMSO", Content), "DMSO",
                    gsub(".*\\s(C[0-9]+)\\s.*", "\\1", Content)),
      # Infection: B0, B1, B10, B100, L0, L1, L10, L100, or dose after DMSO
      Infection = ifelse(grepl("DMSO", Content),
                         paste0("DMSO_", gsub(".*DMSO\\s+([0-9]+).*", "\\1", Content)),
                         gsub(".*\\s([BL][0-9]+)\\s.*", "\\1", Content))
    )
  
  cat("\nAll sample results:\n")
  print(as.data.frame(all_samples %>%
    select(Plate, Experiment, Drug, Infection, OD405, pct_B_B0, predicted_conc, result_flag) %>%
    arrange(Experiment, Drug, Infection)))
  
  write.csv(all_samples,
            file = file.path(export_dir, "PGE2_all_samples.csv"),
            row.names = FALSE)
  
  # Summary statistics
  valid_samples <- all_samples %>% filter(result_flag == "OK")
  
  if (nrow(valid_samples) > 0) {
    sample_summary <- valid_samples %>%
      group_by(Experiment, Drug, Infection) %>%
      summarise(
        mean_conc = mean(predicted_conc, na.rm = TRUE),
        sd_conc = sd(predicted_conc, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
    
    cat("\nSample summary (mean +/- SD):\n")
    print(as.data.frame(sample_summary))
    
    write.csv(sample_summary,
              file = file.path(export_dir, "PGE2_sample_summary.csv"),
              row.names = FALSE)
    
    # Bar plot by treatment
    p_bar <- ggplot(sample_summary,
                    aes(x = Infection, y = mean_conc, fill = Drug)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_errorbar(aes(ymin = pmax(0, mean_conc - sd_conc),
                        ymax = mean_conc + sd_conc),
                    position = position_dodge(width = 0.8), width = 0.3) +
      facet_wrap(~ Experiment, scales = "free_y") +
      labs(
        title = "PGE2 Concentration by Treatment Condition",
        x = "Infection Condition",
        y = "PGE2 (pg/mL)",
        fill = "Drug"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        plot.title = element_text(face = "bold", size = 16),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.ticks = element_line(color = "black"),
        strip.text = element_text(face = "bold", size = 13)
      )
    
    ggsave(file.path(export_dir, "PGE2_barplot_by_treatment.png"),
           plot = p_bar, width = 30, height = 18, units = "cm", dpi = 300)
    
    # Dot plot of individual samples
    p_dot <- ggplot(valid_samples,
                    aes(x = Infection, y = predicted_conc, color = Drug)) +
      geom_jitter(size = 3, width = 0.15, alpha = 0.8) +
      facet_wrap(~ Experiment, scales = "free_y") +
      labs(
        title = "PGE2 - Individual Sample Concentrations",
        x = "Infection Condition",
        y = "PGE2 (pg/mL)",
        color = "Drug"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        plot.title = element_text(face = "bold", size = 16),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.ticks = element_line(color = "black"),
        strip.text = element_text(face = "bold", size = 13)
      )
    
    ggsave(file.path(export_dir, "PGE2_dotplot_individual.png"),
           plot = p_dot, width = 30, height = 18, units = "cm", dpi = 300)
  }
}

cat("\n\nAnalysis complete!\n")
cat("Results saved to:", export_dir, "\n")
cat("\nFiles generated:\n")
cat(paste(list.files(export_dir), collapse = "\n"), "\n")
