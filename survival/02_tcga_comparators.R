#!/usr/bin/env Rscript
# SPDX-License-Identifier: MIT
#
# 02_tcga_comparators.R
#
# TCGA-HNSC overall survival by Gfx and cell-composition metrics. For each of 13 metrics the
# cohort is split by tertile and by top third vs bottom two-thirds, and the script fits
# univariable Cox models and log-rank tests. It writes the source tables for
#   - Supplementary Table 4: 26 univariable models (13 metrics x {binary, tertile});
#   - Supplementary Fig. 4: tertile cut-offs, Kaplan-Meier curves and survival differences
#     (dOS) between the high and low groups;
#   - Supplementary Fig. 2, TCGA row: Tumor:Lymphoid Gfx at 10, 20 and 40 um.
#
# Method:
#   - Overall survival: OS.time in days; OS = 1 death, OS = 0 censored, other codes dropped.
#     A missing OS.time is set to 7000 days. Follow-up is censored at 1,460 days.
#   - Tertiles: type-7 quantiles at 1/3 and 2/3 of the non-missing values;
#     Bottom <= q(1/3) < Middle <= q(2/3) < Top. Binary: Top Third = value > q(2/3),
#     Bottom Two-Thirds otherwise.
#   - Binary models: 2-group log-rank P and Cox HR. Tertile models: 3-level Cox model with the
#     bottom tertile as reference (Wald P) and Tukey all-pairs comparisons (multcomp::glht).
#     AIC and the cox.zph proportional-hazards P are reported for both.
#   - dOS: Kaplan-Meier survival of the high group minus that of the low group (tertile: top
#     minus bottom tertile; binary: top third minus bottom two-thirds).
#
# Usage:
#   Rscript survival/02_tcga_comparators.R [INPUT_CSV] [OUT_DIR] [--dedup-patients] [--no-plots]
#     INPUT_CSV         input table; default data/tcga_analysis_table.csv
#     OUT_DIR           output folder; default outputs/02_tcga_comparators
#                       (both defaults are relative to the folder above this script's folder)
#     --dedup-patients  keep one row per case (see DEDUP_PATIENTS below)
#     --no-plots        skip the PDF plots (then only survival and multcomp are needed)
#
# Input: CSV with columns case (patient identifier), OS, OS.time and the 13 metrics in
#   columns_of_interest below: Gfx at 10, 20 and 40 um for Tumor:Lymphoid, Tumor:Tumor and
#   Lymphoid:Lymphoid (Gfx<r>_Tum_vs_Lym, Gfx<r>_Tum_vs_Tum, Gfx<r>_Lym_vs_Lym), the tumour and
#   lymphoid fractions of classified cells (Tumorpercent, Lymphoidpercent) and the cell counts
#   (Tumor_count, Lymphoid_count). Other columns are ignored.
#
# Outputs (in OUT_DIR):
#   binary_biomarker_survival_summary.csv    binary models (log-rank P, HR, AIC, PH P)
#   tertile_biomarker_survival_summary.csv   tertile models (Wald P, HR, AIC, PH P)
#   pairwise/<metric>_tertile_pairwise_comparisons.csv   Tukey comparisons of the tertiles
#   tertile_survival_differences.csv         survival (%) of the top and bottom tertiles at 300,
#                                            700, 1000 and 1300 days, and the difference
#   supp_table4.csv                          Supplementary Table 4, with n and events
#   km_panel_statistics.csv                  one row per Kaplan-Meier panel: group sizes, events,
#                                            log-rank and Wald P, significance stars
#   tertile_cutoffs.csv                      tertile cut-offs and group sizes
#   km_curves.csv                            Kaplan-Meier step functions
#   delta_os_daily.csv, delta_os_timepoints.csv   dOS on a daily grid (0-1460 days) and at the
#                                            four time points
#   supp_fig2_tcga_row.csv                   Supplementary Fig. 2, TCGA row
#   run_manifest.csv, session_info.txt       input md5, row and event counts, package versions
#   plots/*.pdf                              Kaplan-Meier plots and histograms (unless --no-plots)
#
# Dependencies: R with survival and multcomp; survminer, ggplot2 and ggpubr for the plots.
#   Tested with R 4.4.3, survival 3.8-3, multcomp 1.4-32, survminer 0.5.2, ggplot2 4.0.3 and
#   ggpubr 1.0.0.

# ---- arguments ----------------------------------------------------------------
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
stage_dir <- if (length(script_path) == 1) normalizePath(file.path(dirname(script_path), "..")) else getwd()
args <- commandArgs(trailingOnly = TRUE)
flags <- args[grepl("^--", args)]
pos <- args[!grepl("^--", args)]
IN_FILE <- if (length(pos) >= 1) pos[1] else file.path(stage_dir, "data", "tcga_analysis_table.csv")
OUT_DIR <- if (length(pos) >= 2) pos[2] else file.path(stage_dir, "outputs", "02_tcga_comparators")

# DEDUP_PATIENTS (--dedup-patients) --------------------------------------------
#   FALSE (default): all rows of the input table are analysed.
#   TRUE: one row per case is kept (the first row of each case, in input order).
DEDUP_PATIENTS <- "--dedup-patients" %in% flags
MAKE_PLOTS <- !("--no-plots" %in% flags)

PLOT_DIR <- file.path(OUT_DIR, "plots")
PAIR_DIR <- file.path(OUT_DIR, "pairwise")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PAIR_DIR, recursive = TRUE, showWarnings = FALSE)
if (MAKE_PLOTS) dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# multcomp's single-step (Tukey) P values use randomised quasi-Monte Carlo integration
# (mvtnorm). A seed is set here and again immediately before each glht() call below, because
# plotting also consumes random numbers. The estimates do not depend on the seed.
set.seed(20241205)

# Load necessary libraries
# (plotting packages are loaded only when plots are drawn)
library(survival)
library(multcomp)
if (MAKE_PLOTS) {
  library(survminer)  # for ggsurvplot
  library(ggplot2)    # for ggplot functions
  library(ggpubr)     # for ggarrange
}

# Read data
data <- read.csv(IN_FILE, check.names = FALSE, stringsAsFactors = FALSE)

# Convert data to data frame if it's not already
data <- as.data.frame(data)
rows_in <- nrow(data); cases_in <- length(unique(data$case))

# Check for duplicates
case_counts <- table(data$case)
duplicate_cases <- case_counts[case_counts > 1]

# Print cases with duplicates
cat("Cases with multiple entries:\n")
print(duplicate_cases)

# Remove duplicates keeping first occurrence
# Only with --dedup-patients (see DEDUP_PATIENTS above).
if (DEDUP_PATIENTS) {
  data <- data[!duplicated(data$case), ]
}

# Print final dimensions
cat("\nFinal dataset dimensions:", dim(data)[1], "rows by", dim(data)[2], "columns\n")

# Prepare the data
data$OS.time <- as.numeric(as.character(data$OS.time))
data$OS.time[is.na(data$OS.time)] <- 7000  # Handling missing values
data$event_indicator <- ifelse(data$OS == 1, 1, ifelse(data$OS == 0, 0, NA))
filtered_data <- data[!is.na(data$event_indicator), ]
filtered_data$adjusted_event_indicator <- ifelse(filtered_data$OS.time > 1460, 0, filtered_data$event_indicator)
filtered_data$adjusted_days_to_death <- pmin(filtered_data$OS.time, 1460)

# Columns of interest for survival analysis
columns_of_interest <- c("Gfx10_Tum_vs_Lym", "Gfx20_Tum_vs_Lym", "Gfx40_Tum_vs_Lym", "Gfx10_Tum_vs_Tum", "Gfx20_Tum_vs_Tum", "Gfx40_Tum_vs_Tum", "Gfx10_Lym_vs_Lym", "Gfx20_Lym_vs_Lym", "Gfx40_Lym_vs_Lym", "Tumorpercent", "Lymphoidpercent", "Tumor_count", "Lymphoid_count" )

# Initialize summary data frames to store p-values and hazard ratios
binary_summary_df <- data.frame(column = character(),
                                p_value = numeric(),
                                hazard_ratio = numeric(),
                                hr_lower = numeric(),
                                hr_upper = numeric(),
                                AIC = numeric(),
                                PH_p_value = numeric(),
                                stringsAsFactors = FALSE)

tertile_summary_df <- data.frame(column = character(),
                                 group_comparison = character(),
                                 p_value = numeric(),
                                 hazard_ratio = numeric(),
                                 hr_lower = numeric(),
                                 hr_upper = numeric(),
                                 AIC = numeric(),
                                 PH_p_value = numeric(),
                                 stringsAsFactors = FALSE)

# Collectors for model sizes, P values and Kaplan-Meier tables
added_model_info <- list()   # n, events, group sizes, log-rank P, Wald P
added_km <- list()           # Kaplan-Meier step functions

# Helper to turn a survfit object into a data frame (all event and censoring times)
km_table <- function(fit, metric, scheme) {
  s <- summary(fit, censored = TRUE)
  data.frame(metric = metric, scheme = scheme,
             group = sub("^[^=]*=", "", as.character(s$strata)),
             time_days = s$time, n_risk = s$n.risk, n_event = s$n.event, n_censor = s$n.censor,
             surv = s$surv, lower_95 = s$lower, upper_95 = s$upper, stringsAsFactors = FALSE)
}

# Perform survival analysis for each column
for (column in columns_of_interest) {
  filtered_data[[column]] <- as.numeric(filtered_data[[column]])

  #### Binary Grouping: Top Third vs. Bottom Two-Thirds ####
  # Define the cutoff for the top third
  top_third_cutoff <- quantile(filtered_data[[column]], probs = 2/3, na.rm = TRUE)

  # Create a binary grouping: Top third vs Bottom two-thirds
  filtered_data$binary_group <- ifelse(filtered_data[[column]] > top_third_cutoff, "Top Third", "Bottom Two-Thirds")

  # Perform the survival analysis
  surv_curves_binary <- survfit(Surv(filtered_data$adjusted_days_to_death, filtered_data$adjusted_event_indicator) ~ binary_group, data = filtered_data)

  # Log-rank test to compare the survival curves
  log_rank_test_binary <- survdiff(Surv(filtered_data$adjusted_days_to_death, filtered_data$adjusted_event_indicator) ~ binary_group, data = filtered_data)
  binary_p_value <- 1 - pchisq(log_rank_test_binary$chisq, df = length(log_rank_test_binary$n) - 1)

  # Cox proportional hazards model to calculate hazard ratios
  cox_model_binary <- coxph(Surv(filtered_data$adjusted_days_to_death, filtered_data$adjusted_event_indicator) ~ binary_group, data = filtered_data)
  cox_summary_binary <- summary(cox_model_binary)
  hr_binary <- cox_summary_binary$conf.int[,"exp(coef)"]
  hr_confint_lower_binary <- cox_summary_binary$conf.int[,"lower .95"]
  hr_confint_upper_binary <- cox_summary_binary$conf.int[,"upper .95"]

  # Calculate AIC for the binary model
  aic_binary <- AIC(cox_model_binary)

  # Test proportional hazards assumption
  ph_test_binary <- cox.zph(cox_model_binary)
  ph_p_value_binary <- ph_test_binary$table[1, "p"]

  # Append the results to the binary summary data frame
  binary_summary_df <- rbind(binary_summary_df, data.frame(column = column,
                                                           p_value = binary_p_value,
                                                           hazard_ratio = hr_binary,
                                                           hr_lower = hr_confint_lower_binary,
                                                           hr_upper = hr_confint_upper_binary,
                                                           AIC = aic_binary,
                                                           PH_p_value = ph_p_value_binary,
                                                           stringsAsFactors = FALSE))

  # Model size, group sizes and the Cox Wald P (p_value in the binary summary is the log-rank P)
  gb <- table(factor(filtered_data$binary_group, levels = c("Top Third", "Bottom Two-Thirds")))
  eb <- tapply(filtered_data$adjusted_event_indicator, factor(filtered_data$binary_group, levels = c("Top Third", "Bottom Two-Thirds")), sum)
  added_model_info[[length(added_model_info) + 1]] <- data.frame(
    column = column, scheme = "binary", cutoff_lower = NA_real_, cutoff_upper = unname(top_third_cutoff),
    n_model = cox_model_binary$n, events_model = cox_model_binary$nevent,
    n_high = unname(gb["Top Third"]), n_mid = NA_integer_, n_low = unname(gb["Bottom Two-Thirds"]),
    events_high = unname(eb["Top Third"]), events_mid = NA_real_, events_low = unname(eb["Bottom Two-Thirds"]),
    logrank_chisq = log_rank_test_binary$chisq, logrank_df = length(log_rank_test_binary$n) - 1,
    logrank_P = binary_p_value, cox_wald_P_high_vs_low = cox_summary_binary$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE)
  added_km[[length(added_km) + 1]] <- km_table(surv_curves_binary, column, "binary")

  # Generate the survival plot for binary grouping
  if (MAKE_PLOTS) {
  plot_title_binary <- paste(column, "Survival Analysis (Binary Grouping) - p-value:", signif(binary_p_value, 3))
  g_binary <- ggsurvplot(
    surv_curves_binary,
    data = filtered_data,
    pval = F,  # Include the log-rank test p-value in the plot automatically
    risk.table = TRUE,
    conf.int = TRUE,
    surv.median.line = "hv",  # Horizontal and vertical lines for median
    ggtheme = theme_bw(),
    palette = c('#DE3163', '#6495ED'),  # Adjust the palette if needed
    xlab = "Days to Death",
    ylab = "Survival Probability",
    title = plot_title_binary,
    linetype = "solid",  # Ensure solid line type
    size = 1.5  # Increase line thickness for better visibility
  )

  combined_plot_binary <- ggarrange(
    g_binary$plot,
    g_binary$table,
    ncol = 1,
    heights = c(2, 1)  # Adjust based on your preference
  )

  # Save the combined plot for binary grouping
  ggsave(file.path(PLOT_DIR, paste0(column, "_binary_survival_analysis.pdf")), combined_plot_binary, width = 6, height = 8)
  }

  #### Tertile Grouping: Bottom, Middle, and Top Thirds ####
  # Define the cutoffs for the tertiles
  lower_cutoff <- quantile(filtered_data[[column]], probs = 1/3, na.rm = TRUE)
  upper_cutoff <- quantile(filtered_data[[column]], probs = 2/3, na.rm = TRUE)

  # Create grouping variable based on tertiles
  filtered_data$tertile_group <- cut(filtered_data[[column]],
                                     breaks = c(-Inf, lower_cutoff, upper_cutoff, Inf),
                                     labels = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))

  # Ensure that group is a factor
  filtered_data$tertile_group <- factor(filtered_data$tertile_group, levels = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))

  # Remove any rows with NA in group
  filtered_data_non_na <- filtered_data[!is.na(filtered_data$tertile_group), ]

  # Perform the survival analysis
  surv_curves_tertile <- survfit(Surv(filtered_data_non_na$adjusted_days_to_death, filtered_data_non_na$adjusted_event_indicator) ~ tertile_group, data = filtered_data_non_na)

  # Log-rank test to compare the survival curves
  log_rank_test_tertile <- survdiff(Surv(filtered_data_non_na$adjusted_days_to_death, filtered_data_non_na$adjusted_event_indicator) ~ tertile_group, data = filtered_data_non_na)
  tertile_p_value <- 1 - pchisq(log_rank_test_tertile$chisq, df = length(log_rank_test_tertile$n) - 1)

  # Cox proportional hazards model to calculate hazard ratios
  cox_model_tertile <- coxph(Surv(filtered_data_non_na$adjusted_days_to_death, filtered_data_non_na$adjusted_event_indicator) ~ tertile_group, data = filtered_data_non_na)
  cox_summary_tertile <- summary(cox_model_tertile)

  # Extract hazard ratios and confidence intervals
  hr_tertile <- cox_summary_tertile$conf.int[,"exp(coef)"]
  hr_confint_lower_tertile <- cox_summary_tertile$conf.int[,"lower .95"]
  hr_confint_upper_tertile <- cox_summary_tertile$conf.int[,"upper .95"]

  # Calculate AIC for the tertile model
  aic_tertile <- AIC(cox_model_tertile)

  # Test proportional hazards assumption
  ph_test_tertile <- cox.zph(cox_model_tertile)
  ph_p_values_tertile <- ph_test_tertile$table[,"p"]

  # Append the results to the tertile summary data frame
  for (i in seq_along(hr_tertile)) {
    tertile_summary_df <- rbind(tertile_summary_df, data.frame(column = column,
                                                               group_comparison = rownames(cox_summary_tertile$coefficients)[i],
                                                               p_value = cox_summary_tertile$coefficients[i, "Pr(>|z|)"],
                                                               hazard_ratio = hr_tertile[i],
                                                               hr_lower = hr_confint_lower_tertile[i],
                                                               hr_upper = hr_confint_upper_tertile[i],
                                                               AIC = aic_tertile,
                                                               PH_p_value = ph_p_values_tertile[i],
                                                               stringsAsFactors = FALSE))
  }

  # Model size, group sizes and the 3-group log-rank P (also shown in the plot title)
  gt <- table(filtered_data_non_na$tertile_group)
  et <- tapply(filtered_data_non_na$adjusted_event_indicator, filtered_data_non_na$tertile_group, sum)
  added_model_info[[length(added_model_info) + 1]] <- data.frame(
    column = column, scheme = "tertile", cutoff_lower = unname(lower_cutoff), cutoff_upper = unname(upper_cutoff),
    n_model = cox_model_tertile$n, events_model = cox_model_tertile$nevent,
    n_high = unname(gt["Top Tertile"]), n_mid = unname(gt["Middle Tertile"]), n_low = unname(gt["Bottom Tertile"]),
    events_high = unname(et["Top Tertile"]), events_mid = unname(et["Middle Tertile"]), events_low = unname(et["Bottom Tertile"]),
    logrank_chisq = log_rank_test_tertile$chisq, logrank_df = length(log_rank_test_tertile$n) - 1,
    logrank_P = tertile_p_value, cox_wald_P_high_vs_low = cox_summary_tertile$coefficients[2, "Pr(>|z|)"],
    stringsAsFactors = FALSE)
  added_km[[length(added_km) + 1]] <- km_table(surv_curves_tertile, column, "tertile")

  # Perform pairwise comparisons using multcomp
  set.seed(20241205)  # fixes the Monte Carlo Tukey P values whether or not plots are drawn
  pairwise_comp <- summary(glht(cox_model_tertile, linfct = mcp(tertile_group = "Tukey")))

  # Extract comparison names
  comparison_names <- names(coef(pairwise_comp))

  # Create the data frame
  pairwise_results <- data.frame(column = column,
                                 comparison = comparison_names,
                                 estimate = pairwise_comp$test$coefficients,
                                 p_value = pairwise_comp$test$pvalues,
                                 stringsAsFactors = FALSE)

  # Save pairwise comparisons to a CSV file
  write.csv(pairwise_results, file = file.path(PAIR_DIR, paste0(column, "_tertile_pairwise_comparisons.csv")), row.names = FALSE)

  # Generate the survival plot for tertile grouping
  if (MAKE_PLOTS) {
  plot_title_tertile <- paste(column, "Survival Analysis (Tertile Grouping) - p-value:", signif(tertile_p_value, 3))
  g_tertile <- ggsurvplot(
    surv_curves_tertile,
    data = filtered_data_non_na,
    pval = F,  # Include the log-rank test p-value in the plot automatically
    risk.table = TRUE,
    conf.int = TRUE,
    surv.median.line = "hv",  # Horizontal and vertical lines for median
    ggtheme = theme_bw(),
    palette = c('#DE3163', '#64ed6b', '#6495ED'),  # Adjust the palette for three groups
    xlab = "Days to Death",
    ylab = "Survival Probability",
    title = plot_title_tertile,
    linetype = "solid",  # Ensure solid line type
    size = 1.5  # Increase line thickness for better visibility
  )

  combined_plot_tertile <- ggarrange(
    g_tertile$plot,
    g_tertile$table,
    ncol = 1,
    heights = c(2, 1)  # Adjust based on your preference
  )

  # Save the combined plot for tertile grouping
  ggsave(file.path(PLOT_DIR, paste0(column, "_tertile_survival_analysis.pdf")), combined_plot_tertile, width = 6, height = 8)
  }
}

# Sort the binary summary data frame by p-value in ascending order
binary_summary_df <- binary_summary_df[order(binary_summary_df$p_value),]

# Save the binary summary to a CSV file
write.csv(binary_summary_df, file = file.path(OUT_DIR, "binary_biomarker_survival_summary.csv"), row.names = FALSE)

# Print the binary summary data frame to view the results
print("Binary Grouping Survival Analysis Summary:")
print(binary_summary_df)

# Sort the tertile summary data frame by p-value in ascending order
tertile_summary_df <- tertile_summary_df[order(tertile_summary_df$p_value),]

# Save the tertile summary to a CSV file
write.csv(tertile_summary_df, file = file.path(OUT_DIR, "tertile_biomarker_survival_summary.csv"), row.names = FALSE)

# Print the tertile summary data frame to view the results
print("Tertile Grouping Survival Analysis Summary:")
print(tertile_summary_df)

# New section for histogram visualization
print("Creating distribution histograms for each variable...")

# HISTOGRAM of proportion of tertiles
# ----------------------------------------------------------------
added_cutoffs <- list()

# Create histograms for each column
for (column in columns_of_interest) {
  filtered_data[[column]] <- as.numeric(filtered_data[[column]])

  # Define the cutoffs for the tertiles
  lower_cutoff <- quantile(filtered_data[[column]], probs = 1/3, na.rm = TRUE)
  upper_cutoff <- quantile(filtered_data[[column]], probs = 2/3, na.rm = TRUE)

  # Create grouping variable based on tertiles
  filtered_data$tertile_group <- cut(filtered_data[[column]],
                                     breaks = c(-Inf, lower_cutoff, upper_cutoff, Inf),
                                     labels = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))

  # Create histogram with borders
  if (MAKE_PLOTS) {
  hist_plot <- ggplot(filtered_data, aes(x = .data[[column]], fill = tertile_group)) +
    geom_histogram(bins = 30,
                   color = "black",
                   size = 1,
                   alpha = 0.7) +
    scale_fill_manual(values = c('#DE3163', '#64ed6b', '#6495ED')) +
    theme_bw() +
    labs(title = paste("Distribution of", column),
         x = column,
         y = "Frequency",
         fill = "Tertile Group") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_line(color = "gray95"))

  # Save the histogram
  ggsave(file.path(PLOT_DIR, paste0(column, "_distribution_histogram.pdf")), hist_plot, width = 8, height = 6)
  }

  # Print the tertile cutoff values
  cat("\nCutoff values for", column, ":\n")
  cat("Lower tertile cutoff (33.33%):", lower_cutoff, "\n")
  cat("Upper tertile cutoff (66.67%):", upper_cutoff, "\n\n")

  # Cut-offs and group sizes (Supplementary Fig. 4a)
  tg <- table(filtered_data$tertile_group)
  added_cutoffs[[length(added_cutoffs) + 1]] <- data.frame(
    column = column, n_nonmissing = sum(!is.na(filtered_data[[column]])),
    lower_cutoff = unname(lower_cutoff), upper_cutoff = unname(upper_cutoff),
    n_bottom = unname(tg["Bottom Tertile"]), n_middle = unname(tg["Middle Tertile"]), n_top = unname(tg["Top Tertile"]),
    min = min(filtered_data[[column]], na.rm = TRUE), max = max(filtered_data[[column]], na.rm = TRUE),
    stringsAsFactors = FALSE)
}

print("Histogram generation completed.")


print("Beginning tertile difference analysis at specific timepoints...")

# Specify time points (in days)
time_points <- c(300, 700, 1000, 1300)

# Initialize list to store all results
all_results <- list()

# Analyze each column
for(column in columns_of_interest) {
  # Calculate survival probabilities for high and low tertiles
  filtered_data[[column]] <- as.numeric(filtered_data[[column]])

  # Define tertile cutoffs (reusing from earlier code)
  lower_cutoff <- quantile(filtered_data[[column]], probs = 1/3, na.rm = TRUE)
  upper_cutoff <- quantile(filtered_data[[column]], probs = 2/3, na.rm = TRUE)

  # Create tertile groups
  filtered_data$tertile_group <- cut(filtered_data[[column]],
                                     breaks = c(-Inf, lower_cutoff, upper_cutoff, Inf),
                                     labels = c("Low", "Middle", "High"))

  # Create survival object for high and low tertiles only
  data_high_low <- filtered_data[filtered_data$tertile_group %in% c("Low", "High"),]
  surv_obj <- Surv(data_high_low$adjusted_days_to_death, data_high_low$adjusted_event_indicator)
  fit <- survfit(surv_obj ~ tertile_group, data = data_high_low)

  # Initialize results for this column
  results <- data.frame(
    timepoint = time_points,
    high_survival = NA,
    low_survival = NA,
    difference = NA,
    high_n_risk = NA,
    low_n_risk = NA
  )

  # Calculate survival probabilities at each timepoint
  for(i in seq_along(time_points)) {
    t <- time_points[i]

    # Get survival probabilities
    surv_summary <- summary(fit, times = t)

    # Extract values for high and low tertiles
    high_idx <- which(surv_summary$strata == "tertile_group=High")
    low_idx <- which(surv_summary$strata == "tertile_group=Low")

    if(length(high_idx) > 0 && length(low_idx) > 0) {
      results$high_survival[i] <- surv_summary$surv[high_idx]
      results$low_survival[i] <- surv_summary$surv[low_idx]
      results$difference[i] <- surv_summary$surv[high_idx] - surv_summary$surv[low_idx]
      results$high_n_risk[i] <- surv_summary$n.risk[high_idx]
      results$low_n_risk[i] <- surv_summary$n.risk[low_idx]
    }
  }

  # Convert survival probabilities to percentages
  results$high_survival <- round(results$high_survival * 100, 1)
  results$low_survival <- round(results$low_survival * 100, 1)
  results$difference <- round(results$difference * 100, 1)

  # Store results
  all_results[[column]] <- results

  # Print results for this column
  cat("\nResults for", column, "\n")
  cat("Tertile cutoffs - Lower:", lower_cutoff, "Upper:", upper_cutoff, "\n")
  print(results)
  cat("\n-----------------------------------\n")
}

# Combine all results into a single dataframe
results_df <- do.call(rbind, lapply(names(all_results), function(name) {
  df <- all_results[[name]]
  data.frame(
    column = name,
    timepoint = df$timepoint,
    high_survival = df$high_survival,
    low_survival = df$low_survival,
    difference = df$difference,
    high_n_risk = df$high_n_risk,
    low_n_risk = df$low_n_risk,
    stringsAsFactors = FALSE
  )
}))

# Save to CSV
write.csv(results_df, file.path(OUT_DIR, "tertile_survival_differences.csv"), row.names = FALSE)

print("Tertile difference analysis completed. Results saved to 'tertile_survival_differences.csv'")

# =============================================================================
# Source tables for Supplementary Table 4 and Supplementary Figs 2 and 4.
# Everything below re-uses the objects created above; no model is refitted
# except the dOS curves, which re-apply the same groupings and KM estimator.
# =============================================================================
stars <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
display_name <- function(x) ifelse(x == "Tumorpercent", "Tumor_percent", ifelse(x == "Lymphoidpercent", "Lymphoid_percent", x))
info <- do.call(rbind, added_model_info)

# ---- Supplementary Table 4: 13 binary rows + 13 tertile (top vs bottom) rows ----
b <- binary_summary_df; b$Group <- "Binary"
b$contrast <- "Top third vs bottom two-thirds"; b$P_type <- "log-rank, 2 groups"
tt <- tertile_summary_df[tertile_summary_df$group_comparison == "tertile_groupTop Tertile", ]
tt$Group <- "Tertile"; tt$contrast <- "Top vs bottom tertile (3-level Cox model)"; tt$P_type <- "Cox Wald, top vs bottom"
keep <- c("column", "Group", "contrast", "P_type", "p_value", "hazard_ratio", "hr_lower", "hr_upper", "AIC", "PH_p_value")
st4 <- rbind(b[, keep], tt[, keep])
st4$Metric <- display_name(st4$column)
st4 <- st4[order(st4$Group, st4$Metric), ]   # row order: Binary block then Tertile block, metrics alphabetical
m <- info[match(paste(st4$column, tolower(st4$Group)), paste(info$column, info$scheme)), ]
st4$n_model <- m$n_model; st4$events_model <- m$events_model
st4$n_high <- m$n_high; st4$n_low <- m$n_low; st4$events_high <- m$events_high; st4$events_low <- m$events_low
st4$logrank_P <- m$logrank_P; st4$cox_wald_P_high_vs_low <- m$cox_wald_P_high_vs_low
st4_out <- data.frame(
  Metric = st4$Metric, Group = st4$Group, contrast = st4$contrast,
  P_value = st4$p_value, P_type = st4$P_type, HR = st4$hazard_ratio, HR_lower = st4$hr_lower, HR_upper = st4$hr_upper,
  AIC = st4$AIC, PH_P_value = st4$PH_p_value,
  n_model = st4$n_model, events_model = st4$events_model,
  n_high = st4$n_high, n_low = st4$n_low, events_high = st4$events_high, events_low = st4$events_low,
  logrank_P = st4$logrank_P, cox_wald_P_high_vs_low = st4$cox_wald_P_high_vs_low,
  P_3dp = round(st4$p_value, 3), HR_3dp = round(st4$hazard_ratio, 3), HR_lower_3dp = round(st4$hr_lower, 3),
  HR_upper_3dp = round(st4$hr_upper, 3), AIC_3dp = round(st4$AIC, 3), PH_P_3dp = round(st4$PH_p_value, 3),
  stringsAsFactors = FALSE)
write.csv(st4_out, file.path(OUT_DIR, "supp_table4.csv"), row.names = FALSE)

# ---- Supplementary Fig. 4 / Fig. 2 panel statistics (one row per KM panel) ----
info$metric_display <- display_name(info$column)
info$stars_logrank <- stars(info$logrank_P)
info$stars_cox_wald <- stars(info$cox_wald_P_high_vs_low)
info$plot_title_as_run <- paste(info$column, ifelse(info$scheme == "binary", "Survival Analysis (Binary Grouping) - p-value:",
                                                    "Survival Analysis (Tertile Grouping) - p-value:"), signif(info$logrank_P, 3))
write.csv(info, file.path(OUT_DIR, "km_panel_statistics.csv"), row.names = FALSE)
write.csv(do.call(rbind, added_cutoffs), file.path(OUT_DIR, "tertile_cutoffs.csv"), row.names = FALSE)

# ---- KM step functions (every event/censoring time; follow-up capped at 1,460 days) ----
write.csv(do.call(rbind, added_km), file.path(OUT_DIR, "km_curves.csv"), row.names = FALSE)

# ---- dOS: survival difference high minus low, daily grid 0-1460 days and the 4 fixed time points ----
# Tertile: top vs bottom tertile (Middle excluded), exactly as in the time-point loop above.
# Binary: top third vs bottom two-thirds (the "High vs Mid+Low" panels).
day_grid <- 0:1460
dos_rows <- list(); dos_fixed <- list()
for (column in columns_of_interest) {
  x <- filtered_data[[column]]
  lo <- quantile(x, probs = 1/3, na.rm = TRUE); hi <- quantile(x, probs = 2/3, na.rm = TRUE)
  grp_t <- as.character(cut(x, breaks = c(-Inf, lo, hi, Inf), labels = c("Low", "Middle", "High")))
  grp_b <- ifelse(x > hi, "High", "MidLow")
  for (scheme in c("tertile", "binary")) {
    g <- if (scheme == "tertile") grp_t else grp_b
    low_label <- if (scheme == "tertile") "Low" else "MidLow"
    sel <- !is.na(g) & g %in% c("High", low_label)
    d <- data.frame(t = filtered_data$adjusted_days_to_death[sel], e = filtered_data$adjusted_event_indicator[sel], g = g[sel])
    f_hi <- survfit(Surv(t, e) ~ 1, data = d[d$g == "High", ])
    f_lo <- survfit(Surv(t, e) ~ 1, data = d[d$g == low_label, ])
    s_hi <- summary(f_hi, times = day_grid, extend = TRUE); s_lo <- summary(f_lo, times = day_grid, extend = TRUE)
    dos_rows[[length(dos_rows) + 1]] <- data.frame(metric = column, scheme = scheme, comparison = paste("High minus", low_label),
      time_days = day_grid, time_months = day_grid / (365.25 / 12),
      surv_high = s_hi$surv, surv_low = s_lo$surv, delta_os_pct = 100 * (s_hi$surv - s_lo$surv),
      n_risk_high = s_hi$n.risk, n_risk_low = s_lo$n.risk, stringsAsFactors = FALSE)
    fx <- summary(f_hi, times = time_points, extend = TRUE); fy <- summary(f_lo, times = time_points, extend = TRUE)
    dos_fixed[[length(dos_fixed) + 1]] <- data.frame(metric = column, scheme = scheme, comparison = paste("High minus", low_label),
      timepoint = time_points, high_survival = round(fx$surv * 100, 1), low_survival = round(fy$surv * 100, 1),
      difference = round((fx$surv - fy$surv) * 100, 1), high_n_risk = fx$n.risk, low_n_risk = fy$n.risk, stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, dos_rows), file.path(OUT_DIR, "delta_os_daily.csv"), row.names = FALSE)
dos_fixed <- do.call(rbind, dos_fixed)
write.csv(dos_fixed, file.path(OUT_DIR, "delta_os_timepoints.csv"), row.names = FALSE)
# Consistency check: the tertile rows must equal tertile_survival_differences.csv.
chk <- merge(results_df, dos_fixed[dos_fixed$scheme == "tertile", ], by.x = c("column", "timepoint"), by.y = c("metric", "timepoint"))
stopifnot(nrow(chk) == nrow(results_df), isTRUE(all.equal(chk$difference.x, chk$difference.y)))

# ---- Supplementary Fig. 2, TCGA row (panels f, g): Tumor:Lymphoid at 10/20/40 um ----
s2 <- info[info$column %in% c("Gfx10_Tum_vs_Lym", "Gfx20_Tum_vs_Lym", "Gfx40_Tum_vs_Lym"), ]
hrs <- rbind(
  data.frame(column = binary_summary_df$column, scheme = "binary", HR = binary_summary_df$hazard_ratio,
             HR_lower = binary_summary_df$hr_lower, HR_upper = binary_summary_df$hr_upper),
  data.frame(column = tt$column, scheme = "tertile", HR = tt$hazard_ratio, HR_lower = tt$hr_lower, HR_upper = tt$hr_upper))
s2 <- merge(s2, hrs, by = c("column", "scheme"))
w <- reshape(dos_fixed[dos_fixed$metric %in% s2$column, c("metric", "scheme", "timepoint", "difference")],
             idvar = c("metric", "scheme"), timevar = "timepoint", direction = "wide")
names(w) <- sub("^difference\\.", "delta_os_pct_day", names(w))
s2 <- merge(s2, w, by.x = c("column", "scheme"), by.y = c("metric", "scheme"))
s2 <- s2[order(s2$scheme == "binary", s2$column), ]
write.csv(s2, file.path(OUT_DIR, "supp_fig2_tcga_row.csv"), row.names = FALSE)

# ---- run manifest and session info ----
manifest <- data.frame(
  item = c("input_file", "input_md5", "rows_in", "cases_in", "DEDUP_PATIENTS", "rows_analysed", "cases_analysed",
           "events_within_1460d", "follow_up_cap_days", "missing_OS.time_set_to", "MAKE_PLOTS", "run_time"),
  value = c(normalizePath(IN_FILE), unname(tools::md5sum(IN_FILE)), rows_in, cases_in, DEDUP_PATIENTS, nrow(filtered_data),
            length(unique(filtered_data$case)), sum(filtered_data$adjusted_event_indicator), 1460, 7000, MAKE_PLOTS,
            format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")),
  stringsAsFactors = FALSE)
manifest$value[1] <- basename(manifest$value[1])
write.csv(manifest, file.path(OUT_DIR, "run_manifest.csv"), row.names = FALSE)
# Host-specific lines are left out of session_info.txt: the BLAS/LAPACK library
# paths are reduced to file names, and the time-zone lines are dropped.
si <- capture.output(sessionInfo())
si <- sub("^(BLAS|LAPACK):(\\s+)\\S*/", "\\1:\\2", si)
si <- si[!grepl("^(time zone|tzcode source):", si)]
si <- si[!(si == "" & c(FALSE, head(si, -1) == ""))]
writeLines(si, file.path(OUT_DIR, "session_info.txt"))
cat("\nWrote outputs to", OUT_DIR, "\n")
