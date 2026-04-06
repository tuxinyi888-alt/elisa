###############################################
# Sandwich ELISA Analysis – 4PL - Automated
###############################################

# Load packages (install if missing)
packages <- c("drc", "ggplot2", "dplyr", "readxl", "purrr", "cowplot")
installed <- rownames(installed.packages())
to_install <- packages[!packages %in% installed]
if (length(to_install) > 0) install.packages(to_install)

library(drc)
library(ggplot2)
library(dplyr)
library(readxl)
library(purrr)
library(cowplot)

#####
# USER INPUT — CHANGE ONLY THESE TWO
sheet_name <- "Plate1"     # e.g. "Plate1", "Plate2", "Plate7"
cytokine_name <- "TNFa"    # e.g. "IL1b", "TNFa", "IL6"

#####
# Helper function: convert IL1b → IL‑1β for plot titles
convert_to_greek <- function(name) {
  name <- gsub("a$", "α", name)
  name <- gsub("b$", "β", name)
  name <- gsub("([A-Za-z]+)([0-9]+)", "\\1-\\2", name)
  return(name)
}

cytokine_title <- convert_to_greek(cytokine_name)

#####
# Derived variables
plate_number <- gsub("[^0-9]", "", sheet_name)
plot_title   <- paste(cytokine_title, "Plate", plate_number)

base_dir <- "/Users/anafernandes/Desktop/UCL Div Infection Immunity/Noursadeghi Team/2026 Project Student Xinyi Tu"

export_dir <- file.path(base_dir, paste0("AIF006_", cytokine_name, "_", sheet_name))
dir.create(export_dir, showWarnings = FALSE)

#####
# Import data
file_path <- file.path(base_dir, "AIF006_ELISA_TNFa.xlsx")  # if file changes, update here
data <- read_excel(file_path, sheet = sheet_name)

#####
# Standard curve QC
standards <- data %>%
  filter(Type == "Standard") %>%
  group_by(Concentration) %>%
  summarise(
    mean_od = mean(OD, na.rm = TRUE),
    sd_od   = sd(OD, na.rm = TRUE),
    n       = n(),
    cv_percent = (sd_od / mean_od) * 100
  ) %>%
  arrange(Concentration)

print(standards)

# Remove standards with CV > 30%
standards_fit <- standards %>%
  filter(Concentration > 0, cv_percent <= 30)

#####
# Fit 4PL model
model <- drm(
  mean_od ~ Concentration,
  data = standards_fit,
  fct = LL.4(names = c("Slope", "Lower", "Upper", "EC50"))
)

summary(model)
confint(model)

# Goodness of fit
pseudo_r2 <- 1 - deviance(model) / sum((standards_fit$mean_od - mean(standards_fit$mean_od))^2)
residual_sd <- sqrt(deviance(model) / df.residual(model))

cat("Pseudo R²:", pseudo_r2, "\n")
cat("Residual SD:", residual_sd, "\n")

# Smooth curve
newdata <- data.frame(
  Concentration = exp(seq(
    log(min(standards_fit$Concentration)),
    log(max(standards_fit$Concentration)),
    length.out = 300
  ))
)

newdata$predicted <- predict(model, newdata)

#####
# Standard curve plot
standards$excluded <- standards$cv_percent > 30

p1 <- ggplot(standards, aes(x = Concentration, y = mean_od)) +
  geom_point(aes(color = excluded), size = 3) +
  geom_errorbar(aes(ymin = mean_od - sd_od,
                    ymax = mean_od + sd_od,
                    color = excluded),
                width = 0.05) +
  geom_line(data = newdata,
            aes(x = Concentration, y = predicted),
            linewidth = 0.9,
            color = "black") +
  scale_x_log10(
    breaks = c(3.9, 7.8, 15.6, 31.2, 62.5, 125, 250, 500),
    labels = scales::label_number()
  ) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "#D81B60"),
    labels = c("<30%", ">30%")
  ) +
  labs(
    title = plot_title,
    x = expression(log[10]~"concentration (pg/mL)"),
    y = expression(OD[450]),
    color = "CV%"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    aspect.ratio = 1,
    axis.title = element_text(size = 18),
    axis.text  = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length = unit(2, "mm"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 14)
  )

print(p1)

#####
# Predict sample concentrations
samples <- data %>% filter(Type == "Sample")

max_std <- max(standards_fit$Concentration)
min_std <- min(standards_fit$Concentration)

min_od <- min(predict(model, data.frame(Concentration = standards_fit$Concentration)))
max_od <- max(predict(model, data.frame(Concentration = standards_fit$Concentration)))

predict_conc <- function(od_value) {
  if (od_value < min_od) return(NA)
  if (od_value > max_od) return(NA)
  
  tryCatch({
    uniroot(
      function(x) predict(model, data.frame(Concentration = x)) - od_value,
      lower = min_std,
      upper = max_std
    )$root
  },
  error = function(e) NA)
}

samples <- samples %>%
  mutate(
    predicted_conc = map_dbl(OD, predict_conc),
    result_flag = case_when(
      is.na(predicted_conc) & OD < min_od ~ paste("<", min_std),
      is.na(predicted_conc) & OD > max_od ~ paste(">", max_std),
      TRUE ~ "OK"
    )
  )

print(samples)

#####
# Standard curve + samples plot
y_min <- min(standards$mean_od - standards$sd_od)
y_max <- max(standards$mean_od + standards$sd_od)

p2 <- ggplot() +
  geom_point(data = standards,
             aes(x = Concentration, y = mean_od),
             size = 3, color = "black") +
  geom_errorbar(data = standards,
                aes(x = Concentration,
                    ymin = mean_od - sd_od,
                    ymax = mean_od + sd_od),
                width = 0.05, color = "black") +
  geom_line(data = newdata,
            aes(x = Concentration, y = predicted),
            linewidth = 0.9, color = "black") +
  geom_point(data = samples,
             aes(x = predicted_conc, y = OD),
             color = "orange2", size = 3) +
  scale_x_log10(
    breaks = c(3.9, 7.8, 15.6, 31.2, 62.5, 125, 250, 500),
    labels = scales::label_number()
  ) +
  ylim(y_min, y_max) +
  labs(
    title = plot_title,
    x = expression(log[10]~"concentration (pg/mL)"),
    y = expression(OD[450])
  ) +
  theme_minimal(base_size = 16) +
  theme(
    aspect.ratio = 1,
    axis.title = element_text(size = 18),
    axis.text  = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length = unit(2, "mm"),
    axis.ticks = element_line(color = "black")
  )

print(p2)

#####
# Final output table
final_samples <- samples %>%
  select(ID, OD, predicted_conc, result_flag)

View(final_samples)

#####
# Export figures
p1_fixed <- ggdraw(p1)
p2_fixed <- ggdraw(p2)

fig_width  <- 15
fig_height <- 12.5

ggsave(file.path(export_dir, paste0("AIF006_", cytokine_name, "_", sheet_name, "_standard_curve.png")),
       plot = p1_fixed, width = fig_width, height = fig_height, units = "cm", dpi = 300)

ggsave(file.path(export_dir, paste0("AIF006_", cytokine_name, "_", sheet_name, "_standard_curve.pdf")),
       plot = p1_fixed, width = fig_width, height = fig_height, units = "cm")

ggsave(file.path(export_dir, paste0("AIF006_", cytokine_name, "_", sheet_name, "_standard_curve_samples.png")),
       plot = p2_fixed, width = fig_width, height = fig_height, units = "cm", dpi = 300)

ggsave(file.path(export_dir, paste0("AIF006_", cytokine_name, "_", sheet_name, "_standard_curve_samples.pdf")),
       plot = p2_fixed, width = fig_width, height = fig_height, units = "cm")

#####
# Export tables and model statistics
write.csv(standards,
          file = file.path(export_dir, paste0("AIF006_", cytokine_name, "_", sheet_name, "_standards_table.csv")),
          row.names = FALSE)

write.csv(final_samples,
          file = file.path(export_dir, paste0("AIF006_", cytokine_name, "_", sheet_name, "_samples_table.csv")),
          row.names = FALSE)

model_summary <- capture.output(summary(model))
writeLines(model_summary,
           con = file.path(export_dir, paste0("AIF006_", cytokine_name, "_", sheet_name, "_model-summary.txt")))

model_confint <- capture.output(confint(model))
writeLines(model_confint,
           con = file.path(export_dir, paste0("AIF006_", cytokine_name, "_", sheet_name, "_confidence-intervals.txt")))

gof_stats <- data.frame(
  pseudo_R2 = pseudo_r2,
  residual_SD = residual_sd
)

write.csv(gof_stats,
          file = file.path(export_dir, paste0("AIF006_", cytokine_name, "_", sheet_name, "_goodness-fit.csv")),
          row.names = FALSE)
