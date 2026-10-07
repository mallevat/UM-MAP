#!/usr/bin/env Rscript
# SPDX-License-Identifier: MIT
#
# 07_suppfig2_umich_comparators.R
#
# Gfx comparator survival analysis for the UMICH cohorts, or for any cohort table with an
# overall-survival time, a death indicator and one or more Gfx columns. The default mode writes
# the UMICH rows of Supplementary Fig. 2: Kaplan-Meier curves by Tumor:Lymphoid Gfx at 10, 20
# and 40 um, tertile grouping (top row) and binary grouping (bottom row), with significance
# stars, and the survival-difference (dOS) panels; UMICH1 (training cohort, panels b-d) and
# UMICH2 (testing cohort, panels h-j).
#
# Method (as in 02_tcga_comparators.R):
#   - Time in days. A missing time is set to 7000 days and so ends up censored at the cap.
#     Events: OS == 1 death, OS == 0 censored; other codes are dropped.
#   - Follow-up is censored at 1,460 days (4 years): an event after day 1,460 becomes censored,
#     and every time is truncated at 1,460.
#   - Tertiles: type-7 quantiles at 1/3 and 2/3 of the non-missing values of that metric;
#     groups (-Inf, q1], (q1, q2], (q2, Inf] = Bottom / Middle / Top tertile.
#   - Binary: top third = value > q(2/3) (strictly greater); everything else is the bottom
#     two-thirds.
#   - Records with a missing value of a metric are left out of that metric's models only.
#   - Binary models: 2-group log-rank P and Cox HR (top third vs bottom two-thirds).
#     Tertile models: 3-level Cox model with the bottom tertile as reference (Wald P and HR for
#     top vs bottom), and the 3-group log-rank P. AIC and the cox.zph proportional-hazards P
#     are reported.
#   - Stars: tertile panels use the Cox Wald P for top vs bottom tertile; binary panels use the
#     log-rank P. *** P < 0.001, ** P < 0.01, * P < 0.05, ns otherwise.
#   - dOS: Kaplan-Meier survival of the high group minus the low group, in percentage points
#     rounded to 0.1, at 300, 700, 1000 and 1300 days. Tertile: top minus bottom tertile
#     ("High vs Low"); binary: top third minus bottom two-thirds ("High vs Mid+Low").
#   - Tukey pairwise comparisons of the three tertiles: contrasts of the Cox coefficients with
#     single-step adjusted P values from the multivariate normal distribution (the test that
#     multcomp::glht reports for a Cox model), computed by numerical integration. They are not
#     used in any figure.
#   - No de-duplication by default. Generic mode offers --dedup-by COL (keeps the first row of
#     each COL value).
#   Only base R and the recommended package survival are used; figures use base graphics.
#
# Usage:
#   Rscript survival/07_suppfig2_umich_comparators.R [--umich1 PATH] [--umich2 PATH] [--out DIR] [--no-plots]
#   Generic mode, any cohort table:
#   Rscript survival/07_suppfig2_umich_comparators.R --input PATH --time COL --event COL --metrics C1,C2,...
#       [--time-scale 30.4375] [--subset COL=VALUE] [--dedup-by COL] [--label NAME] [--out DIR] [--no-plots]
#   --time-scale multiplies the time column to give days (default 1); --subset keeps the rows where
#   COL equals VALUE; --label names the run (default "cohort").
#   Defaults (relative to the folder above this script's folder): --umich1 data/umich1_recurrence_deid.csv,
#   --umich2 data/umich2_radius_deid.csv (or the environment variables GFX_UMICH1 / GFX_UMICH2),
#   --out outputs/07_suppfig2_umich.
#
# Inputs (the UMICH patient-level inputs are not distributed):
#   UMICH1  one row per patient: Gfx_10um, Gfx_20um, Gfx_40um; stime (months; days = stime x 30.4375);
#           deathstatus (1 death, 0 censored); site (used for the oral-cavity subset).
#   UMICH2  one row per record: gfx10, gfx20, gfx40 (gfx20 and gfx40 may be missing); os_time (days);
#           os_event (1 death, 0 censored).
#
# Runs written by the default mode:
#   UMICH1_all447          all UMICH1 patients
#   UMICH2_pooled295       all UMICH2 records (the Gfx20 and Gfx40 models use the records that have them)
#   UMICH2_linked268       sensitivity: the UMICH2 records that have all three radii
#   UMICH1_oralcavity223   sensitivity: UMICH1 oral-cavity patients only (site == "oral cavity")
#
# Outputs (in --out; aggregate only, no record-level rows; Kaplan-Meier tables on a monthly grid):
#   suppfig2_umich_rows.csv    one row per run, scheme and radius: n, events, cut-offs, group sizes,
#                              HR, P values, star, dOS at the four time points
#   model_table.csv            Supplementary Table 4 layout, per run
#   km_panel_statistics.csv, tertile_cutoffs.csv, tertile_survival_differences.csv,
#   delta_os_timepoints.csv, delta_os_monthly.csv, km_curves_monthly.csv, pairwise_tukey.csv,
#   binary_biomarker_survival_summary.csv, tertile_biomarker_survival_summary.csv,
#   suppfig2_<run>.pdf (unless --no-plots), run_manifest.csv, session_info.txt
#
# Dependencies: R with survival (base graphics only).
#   Tested with R 4.4.3 and survival 3.8-3.

suppressPackageStartupMessages(library(survival))

CAP_DAYS <- 1460
MISSING_TIME_DAYS <- 7000
TIME_POINTS <- c(300, 700, 1000, 1300)
DAYS_PER_MONTH <- 30.4375
MONTH_GRID_DAYS <- pmin((0:48) * DAYS_PER_MONTH, CAP_DAYS)
TL_COLUMNS <- c("Gfx10_Tum_vs_Lym", "Gfx20_Tum_vs_Lym", "Gfx40_Tum_vs_Lym")
PAL_BINARY <- c("#DE3163", "#6495ED")             # Bottom Two-Thirds, Top Third (as in 02_tcga_comparators.R)
PAL_TERTILE <- c("#DE3163", "#64ed6b", "#6495ED")  # Bottom, Middle, Top Tertile (as in 02_tcga_comparators.R)

stars <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

# ---- time/event preparation (as in 02_tcga_comparators.R) ---------------------------------------------------
prepare_os <- function(data) {
  data$OS.time <- as.numeric(as.character(data$OS.time))
  data$OS.time[is.na(data$OS.time)] <- MISSING_TIME_DAYS  # Handling missing values
  data$event_indicator <- ifelse(data$OS == 1, 1, ifelse(data$OS == 0, 0, NA))
  filtered_data <- data[!is.na(data$event_indicator), ]
  filtered_data$adjusted_event_indicator <- ifelse(filtered_data$OS.time > CAP_DAYS, 0, filtered_data$event_indicator)
  filtered_data$adjusted_days_to_death <- pmin(filtered_data$OS.time, CAP_DAYS)
  filtered_data
}

# survfit -> data frame (as km_table in 02_tcga_comparators.R)
km_table <- function(fit, metric, scheme) {
  s <- summary(fit, censored = TRUE)
  data.frame(metric = metric, scheme = scheme,
             group = sub("^[^=]*=", "", as.character(s$strata)),
             time_days = s$time, n_risk = s$n.risk, n_event = s$n.event, n_censor = s$n.censor,
             surv = s$surv, lower_95 = s$lower, upper_95 = s$upper, stringsAsFactors = FALSE)
}

# survfit -> data frame on a fixed time grid (used for UMICH outputs instead of event-time steps)
km_grid_table <- function(fit, metric, scheme, times) {
  s <- summary(fit, times = times, extend = TRUE)
  data.frame(metric = metric, scheme = scheme, group = sub("^[^=]*=", "", as.character(s$strata)),
             month = round(s$time / DAYS_PER_MONTH, 3), time_days = s$time, n_risk = s$n.risk,
             surv = s$surv, lower_95 = s$lower, upper_95 = s$upper, stringsAsFactors = FALSE)
}

# ---- Tukey all-pairs comparison of the tertile Cox coefficients (single-step test, as multcomp::glht) ------
# Contrasts of (b_Middle, b_Top), both relative to the bottom tertile:
#   Middle - Bottom = b_M, Top - Bottom = b_T, Top - Middle = b_T - b_M.
# Two-sided single-step adjusted P_i = P(max_j |Z_j| >= |z_i|), Z = standardised contrasts of U ~ N(0, V)
# (normal reference distribution, as multcomp uses for a Cox model). The three contrasts have rank 2, so the
# probability of the box is a one-dimensional integral over U1 with U2 | U1 normal.
tukey_single_step <- function(cox_fit) {
  b <- unname(coef(cox_fit)); V <- unname(vcov(cox_fit))
  stopifnot(length(b) == 2)
  K <- rbind(c(1, 0), c(0, 1), c(-1, 1))
  est <- drop(K %*% b); S <- K %*% V %*% t(K); se <- sqrt(diag(S)); z <- est / se
  s1 <- sqrt(V[1, 1]); m <- V[1, 2] / V[1, 1]; s21 <- sqrt(V[2, 2] - V[1, 2]^2 / V[1, 1])
  box_prob <- function(cc) {
    a1 <- cc * se[1]; a2 <- cc * se[2]; a3 <- cc * se[3]
    f <- function(u1) {
      lo <- pmax(-a2, u1 - a3); hi <- pmin(a2, u1 + a3)
      pr <- ifelse(hi > lo, pnorm((hi - m * u1) / s21) - pnorm((lo - m * u1) / s21), 0)
      dnorm(u1, sd = s1) * pr
    }
    brk <- sort(unique(c(-a1, a1, c(a3 - a2, a2 - a3, -a2 - a3, a2 + a3)[abs(c(a3 - a2, a2 - a3, -a2 - a3, a2 + a3)) < a1])))
    sum(vapply(seq_len(length(brk) - 1), function(k)
      integrate(f, brk[k], brk[k + 1], rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value, numeric(1)))
  }
  p <- vapply(abs(z), function(cc) max(0, 1 - box_prob(cc)), numeric(1))
  data.frame(comparison = c("Middle Tertile - Bottom Tertile", "Top Tertile - Bottom Tertile", "Top Tertile - Middle Tertile"),
             estimate = est, p_value = p, stringsAsFactors = FALSE)
}

# ---- the comparator analysis (model code as in 02_tcga_comparators.R) ---------------------------------------
# day_grid: time grid for the dOS curve (default 0:1460 days; cohort runs use MONTH_GRID_DAYS).
# km_steps = TRUE returns Kaplan-Meier tables at every event/censoring time; FALSE returns them on day_grid.
comparator_core <- function(filtered_data, columns_of_interest, fig_columns = TL_COLUMNS,
                            day_grid = 0:CAP_DAYS, km_steps = TRUE) {
  binary_summary_df <- data.frame(column = character(), p_value = numeric(), hazard_ratio = numeric(),
                                  hr_lower = numeric(), hr_upper = numeric(), AIC = numeric(),
                                  PH_p_value = numeric(), stringsAsFactors = FALSE)
  tertile_summary_df <- data.frame(column = character(), group_comparison = character(), p_value = numeric(),
                                   hazard_ratio = numeric(), hr_lower = numeric(), hr_upper = numeric(),
                                   AIC = numeric(), PH_p_value = numeric(), stringsAsFactors = FALSE)
  added_model_info <- list(); added_km <- list(); pairwise <- list(); fits <- list()

  for (column in columns_of_interest) {
    filtered_data[[column]] <- as.numeric(filtered_data[[column]])

    #### Binary Grouping: Top Third vs. Bottom Two-Thirds ####
    top_third_cutoff <- quantile(filtered_data[[column]], probs = 2/3, na.rm = TRUE)
    filtered_data$binary_group <- ifelse(filtered_data[[column]] > top_third_cutoff, "Top Third", "Bottom Two-Thirds")
    surv_curves_binary <- survfit(Surv(filtered_data$adjusted_days_to_death, filtered_data$adjusted_event_indicator) ~ binary_group, data = filtered_data)
    log_rank_test_binary <- survdiff(Surv(filtered_data$adjusted_days_to_death, filtered_data$adjusted_event_indicator) ~ binary_group, data = filtered_data)
    binary_p_value <- 1 - pchisq(log_rank_test_binary$chisq, df = length(log_rank_test_binary$n) - 1)
    cox_model_binary <- coxph(Surv(filtered_data$adjusted_days_to_death, filtered_data$adjusted_event_indicator) ~ binary_group, data = filtered_data)
    cox_summary_binary <- summary(cox_model_binary)
    hr_binary <- cox_summary_binary$conf.int[,"exp(coef)"]
    hr_confint_lower_binary <- cox_summary_binary$conf.int[,"lower .95"]
    hr_confint_upper_binary <- cox_summary_binary$conf.int[,"upper .95"]
    aic_binary <- AIC(cox_model_binary)
    ph_test_binary <- cox.zph(cox_model_binary)
    ph_p_value_binary <- ph_test_binary$table[1, "p"]
    binary_summary_df <- rbind(binary_summary_df, data.frame(column = column, p_value = binary_p_value,
                                                             hazard_ratio = hr_binary, hr_lower = hr_confint_lower_binary,
                                                             hr_upper = hr_confint_upper_binary, AIC = aic_binary,
                                                             PH_p_value = ph_p_value_binary, stringsAsFactors = FALSE))
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
    added_km[[length(added_km) + 1]] <- if (km_steps) km_table(surv_curves_binary, column, "binary") else km_grid_table(surv_curves_binary, column, "binary", day_grid)

    #### Tertile Grouping: Bottom, Middle, and Top Thirds ####
    lower_cutoff <- quantile(filtered_data[[column]], probs = 1/3, na.rm = TRUE)
    upper_cutoff <- quantile(filtered_data[[column]], probs = 2/3, na.rm = TRUE)
    filtered_data$tertile_group <- cut(filtered_data[[column]],
                                       breaks = c(-Inf, lower_cutoff, upper_cutoff, Inf),
                                       labels = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
    filtered_data$tertile_group <- factor(filtered_data$tertile_group, levels = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
    filtered_data_non_na <- filtered_data[!is.na(filtered_data$tertile_group), ]
    surv_curves_tertile <- survfit(Surv(filtered_data_non_na$adjusted_days_to_death, filtered_data_non_na$adjusted_event_indicator) ~ tertile_group, data = filtered_data_non_na)
    log_rank_test_tertile <- survdiff(Surv(filtered_data_non_na$adjusted_days_to_death, filtered_data_non_na$adjusted_event_indicator) ~ tertile_group, data = filtered_data_non_na)
    tertile_p_value <- 1 - pchisq(log_rank_test_tertile$chisq, df = length(log_rank_test_tertile$n) - 1)
    cox_model_tertile <- coxph(Surv(filtered_data_non_na$adjusted_days_to_death, filtered_data_non_na$adjusted_event_indicator) ~ tertile_group, data = filtered_data_non_na)
    cox_summary_tertile <- summary(cox_model_tertile)
    hr_tertile <- cox_summary_tertile$conf.int[,"exp(coef)"]
    hr_confint_lower_tertile <- cox_summary_tertile$conf.int[,"lower .95"]
    hr_confint_upper_tertile <- cox_summary_tertile$conf.int[,"upper .95"]
    aic_tertile <- AIC(cox_model_tertile)
    ph_test_tertile <- cox.zph(cox_model_tertile)
    ph_p_values_tertile <- ph_test_tertile$table[,"p"]
    for (i in seq_along(hr_tertile)) {
      tertile_summary_df <- rbind(tertile_summary_df, data.frame(column = column,
                                                                 group_comparison = rownames(cox_summary_tertile$coefficients)[i],
                                                                 p_value = cox_summary_tertile$coefficients[i, "Pr(>|z|)"],
                                                                 hazard_ratio = hr_tertile[i], hr_lower = hr_confint_lower_tertile[i],
                                                                 hr_upper = hr_confint_upper_tertile[i], AIC = aic_tertile,
                                                                 PH_p_value = ph_p_values_tertile[i], stringsAsFactors = FALSE))
    }
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
    added_km[[length(added_km) + 1]] <- if (km_steps) km_table(surv_curves_tertile, column, "tertile") else km_grid_table(surv_curves_tertile, column, "tertile", day_grid)

    # Pairwise comparisons of the tertiles (Tukey contrasts, single-step P)
    pairwise[[column]] <- cbind(column = column, tukey_single_step(cox_model_tertile), stringsAsFactors = FALSE)
    fits[[column]] <- list(binary = surv_curves_binary, tertile = surv_curves_tertile)
  }
  binary_summary_df <- binary_summary_df[order(binary_summary_df$p_value), ]
  tertile_summary_df <- tertile_summary_df[order(tertile_summary_df$p_value), ]

  # Tertile cut-offs and group sizes
  added_cutoffs <- list()
  for (column in columns_of_interest) {
    filtered_data[[column]] <- as.numeric(filtered_data[[column]])
    lower_cutoff <- quantile(filtered_data[[column]], probs = 1/3, na.rm = TRUE)
    upper_cutoff <- quantile(filtered_data[[column]], probs = 2/3, na.rm = TRUE)
    filtered_data$tertile_group <- cut(filtered_data[[column]], breaks = c(-Inf, lower_cutoff, upper_cutoff, Inf),
                                       labels = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
    tg <- table(filtered_data$tertile_group)
    added_cutoffs[[length(added_cutoffs) + 1]] <- data.frame(
      column = column, n_nonmissing = sum(!is.na(filtered_data[[column]])),
      lower_cutoff = unname(lower_cutoff), upper_cutoff = unname(upper_cutoff),
      n_bottom = unname(tg["Bottom Tertile"]), n_middle = unname(tg["Middle Tertile"]), n_top = unname(tg["Top Tertile"]),
      min = min(filtered_data[[column]], na.rm = TRUE), max = max(filtered_data[[column]], na.rm = TRUE),
      stringsAsFactors = FALSE)
  }

  # Tertile difference analysis at specific time points (summary() without extend)
  all_results <- list()
  for (column in columns_of_interest) {
    filtered_data[[column]] <- as.numeric(filtered_data[[column]])
    lower_cutoff <- quantile(filtered_data[[column]], probs = 1/3, na.rm = TRUE)
    upper_cutoff <- quantile(filtered_data[[column]], probs = 2/3, na.rm = TRUE)
    filtered_data$tertile_group <- cut(filtered_data[[column]], breaks = c(-Inf, lower_cutoff, upper_cutoff, Inf),
                                       labels = c("Low", "Middle", "High"))
    data_high_low <- filtered_data[filtered_data$tertile_group %in% c("Low", "High"),]
    surv_obj <- Surv(data_high_low$adjusted_days_to_death, data_high_low$adjusted_event_indicator)
    fit <- survfit(surv_obj ~ tertile_group, data = data_high_low)
    results <- data.frame(timepoint = TIME_POINTS, high_survival = NA, low_survival = NA, difference = NA,
                          high_n_risk = NA, low_n_risk = NA)
    for (i in seq_along(TIME_POINTS)) {
      surv_summary <- summary(fit, times = TIME_POINTS[i])
      high_idx <- which(surv_summary$strata == "tertile_group=High")
      low_idx <- which(surv_summary$strata == "tertile_group=Low")
      if (length(high_idx) > 0 && length(low_idx) > 0) {
        results$high_survival[i] <- surv_summary$surv[high_idx]
        results$low_survival[i] <- surv_summary$surv[low_idx]
        results$difference[i] <- surv_summary$surv[high_idx] - surv_summary$surv[low_idx]
        results$high_n_risk[i] <- surv_summary$n.risk[high_idx]
        results$low_n_risk[i] <- surv_summary$n.risk[low_idx]
      }
    }
    results$high_survival <- round(results$high_survival * 100, 1)
    results$low_survival <- round(results$low_survival * 100, 1)
    results$difference <- round(results$difference * 100, 1)
    all_results[[column]] <- results
  }
  results_df <- do.call(rbind, lapply(names(all_results), function(name) {
    df <- all_results[[name]]
    data.frame(column = name, timepoint = df$timepoint, high_survival = df$high_survival, low_survival = df$low_survival,
               difference = df$difference, high_n_risk = df$high_n_risk, low_n_risk = df$low_n_risk, stringsAsFactors = FALSE)
  }))

  # dOS curves and fixed time points, tertile and binary
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
      fx <- summary(f_hi, times = TIME_POINTS, extend = TRUE); fy <- summary(f_lo, times = TIME_POINTS, extend = TRUE)
      dos_fixed[[length(dos_fixed) + 1]] <- data.frame(metric = column, scheme = scheme, comparison = paste("High minus", low_label),
        timepoint = TIME_POINTS, high_survival = round(fx$surv * 100, 1), low_survival = round(fy$surv * 100, 1),
        difference = round((fx$surv - fy$surv) * 100, 1), high_n_risk = fx$n.risk, low_n_risk = fy$n.risk, stringsAsFactors = FALSE)
    }
  }
  dos_daily <- do.call(rbind, dos_rows)
  dos_fixed <- do.call(rbind, dos_fixed)
  chk <- merge(results_df, dos_fixed[dos_fixed$scheme == "tertile", ], by.x = c("column", "timepoint"), by.y = c("metric", "timepoint"))
  stopifnot(nrow(chk) == nrow(results_df), isTRUE(all.equal(chk$difference.x, chk$difference.y)))

  # Model statistics per KM panel and the Supplementary Fig. 2 row
  info <- do.call(rbind, added_model_info)
  info$metric_display <- ifelse(info$column == "Tumorpercent", "Tumor_percent", ifelse(info$column == "Lymphoidpercent", "Lymphoid_percent", info$column))
  info$stars_logrank <- stars(info$logrank_P)
  info$stars_cox_wald <- stars(info$cox_wald_P_high_vs_low)
  info$plot_title_as_run <- paste(info$column, ifelse(info$scheme == "binary", "Survival Analysis (Binary Grouping) - p-value:",
                                                      "Survival Analysis (Tertile Grouping) - p-value:"), signif(info$logrank_P, 3))
  tt <- tertile_summary_df[tertile_summary_df$group_comparison == "tertile_groupTop Tertile", ]
  s2 <- info[info$column %in% fig_columns, ]
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

  list(binary_summary = binary_summary_df, tertile_summary = tertile_summary_df, info = info,
       km = do.call(rbind, added_km), pairwise = do.call(rbind, pairwise), cutoffs = do.call(rbind, added_cutoffs),
       tertile_differences = results_df, dos_curve = dos_daily, dos_fixed = dos_fixed, s2 = s2, fits = fits,
       values = filtered_data[, columns_of_interest, drop = FALSE],
       n_rows = nrow(filtered_data), n_events = sum(filtered_data$adjusted_event_indicator))
}

# Supplementary Table 4 layout (as in 02_tcga_comparators.R)
supp_table4_layout <- function(res) {
  info <- res$info
  display_name <- function(x) ifelse(x == "Tumorpercent", "Tumor_percent", ifelse(x == "Lymphoidpercent", "Lymphoid_percent", x))
  b <- res$binary_summary; b$Group <- "Binary"
  b$contrast <- "Top third vs bottom two-thirds"; b$P_type <- "log-rank, 2 groups"
  tt <- res$tertile_summary[res$tertile_summary$group_comparison == "tertile_groupTop Tertile", ]
  tt$Group <- "Tertile"; tt$contrast <- "Top vs bottom tertile (3-level Cox model)"; tt$P_type <- "Cox Wald, top vs bottom"
  keep <- c("column", "Group", "contrast", "P_type", "p_value", "hazard_ratio", "hr_lower", "hr_upper", "AIC", "PH_p_value")
  st4 <- rbind(b[, keep], tt[, keep])
  st4$Metric <- display_name(st4$column)
  st4 <- st4[order(st4$Group, st4$Metric), ]
  m <- info[match(paste(st4$column, tolower(st4$Group)), paste(info$column, info$scheme)), ]
  data.frame(Metric = st4$Metric, Group = st4$Group, contrast = st4$contrast,
             P_value = st4$p_value, P_type = st4$P_type, HR = st4$hazard_ratio, HR_lower = st4$hr_lower, HR_upper = st4$hr_upper,
             AIC = st4$AIC, PH_P_value = st4$PH_p_value, n_model = m$n_model, events_model = m$events_model,
             n_high = m$n_high, n_low = m$n_low, events_high = m$events_high, events_low = m$events_low,
             logrank_P = m$logrank_P, cox_wald_P_high_vs_low = m$cox_wald_P_high_vs_low,
             P_3dp = round(st4$p_value, 3), HR_3dp = round(st4$hazard_ratio, 3), HR_lower_3dp = round(st4$hr_lower, 3),
             HR_upper_3dp = round(st4$hr_upper, 3), AIC_3dp = round(st4$AIC, 3), PH_P_3dp = round(st4$PH_p_value, 3),
             stringsAsFactors = FALSE)
}

# ---- figures (base graphics) -------------------------------------------------------------------------------
step_xy <- function(t, v, tend) {
  k <- length(t)
  list(x = c(t[1], rep(t[-1], each = 2), tend), y = c(rep(v[-k], each = 2), v[k], v[k]))
}
carry_forward <- function(v) { for (i in seq_along(v)[-1]) if (is.na(v[i])) v[i] <- v[i - 1]; v }

draw_km <- function(fit, cols, legend_labels, main, star, note) {
  plot(NA, xlim = c(0, 48), ylim = c(0, 1.02), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", xaxs = "i", yaxs = "i")
  abline(h = seq(0, 1, 0.25), v = c(16, 32), col = "grey92", lwd = 0.6)
  axis(1, at = c(0, 16, 32, 48), cex.axis = 0.75, mgp = c(2, 0.4, 0)); axis(2, at = seq(0, 1, 0.25), labels = sprintf("%.2f", seq(0, 1, 0.25)), las = 1, cex.axis = 0.7, mgp = c(2, 0.5, 0))
  title(main = main, cex.main = 0.85, line = 0.4)
  nst <- length(fit$strata); idx_end <- cumsum(fit$strata); idx_start <- c(1, head(idx_end, -1) + 1)
  for (k in seq_len(nst)) {   # CI ribbons first
    ii <- idx_start[k]:idx_end[k]
    t <- c(0, fit$time[ii]) / DAYS_PER_MONTH; tend <- max(t)
    lo <- carry_forward(c(1, fit$lower[ii])); up <- carry_forward(c(1, fit$upper[ii]))
    a <- step_xy(t, up, tend); b <- step_xy(t, lo, tend)
    polygon(c(a$x, rev(b$x)), c(a$y, rev(b$y)), col = adjustcolor(cols[k], 0.22), border = NA)
  }
  for (k in seq_len(nst)) {
    ii <- idx_start[k]:idx_end[k]
    t <- c(0, fit$time[ii]) / DAYS_PER_MONTH; s <- c(1, fit$surv[ii])
    a <- step_xy(t, s, max(t)); lines(a$x, a$y, col = cols[k], lwd = 1.6)
    cen <- fit$n.censor[ii] > 0
    if (any(cen)) points(fit$time[ii][cen] / DAYS_PER_MONTH, fit$surv[ii][cen], pch = 3, cex = 0.35, col = cols[k])
    below <- which(fit$surv[ii] <= 0.5)   # median survival line ("hv"), when reached
    if (length(below)) { med <- fit$time[ii][below[1]] / DAYS_PER_MONTH
      segments(0, 0.5, med, 0.5, lty = 2, col = "grey40", lwd = 0.7); segments(med, 0, med, 0.5, lty = 2, col = "grey40", lwd = 0.7) }
  }
  text(46, 0.95, star, col = "#6495ED", cex = 1.25, adj = c(1, 0.5), font = 2)
  ly <- 0.05 + 0.075 * (rev(seq_along(legend_labels)) - 1)
  text(46.5, ly, legend_labels, col = cols, cex = 0.7, adj = c(1, 0))
  mtext(note, side = 1, line = 1.35, cex = 0.5, col = "grey25")
}

draw_dos <- function(vals, main, ylab) {
  x <- c(0, TIME_POINTS) / DAYS_PER_MONTH; y <- c(0, vals)
  yl <- range(c(0, y, 10), na.rm = TRUE); yl <- c(min(0, yl[1]), max(yl[2], 10) * 1.15)
  plot(x, y, type = "o", pch = 21, bg = "#6E8FD6", col = "#4F6FBF", xlim = c(0, 48), ylim = yl, xaxt = "n", las = 1,
       xlab = "", ylab = "", cex = 1, cex.axis = 0.65, mgp = c(2, 0.45, 0), bty = "l")
  axis(1, at = c(0, 16, 32, 48), cex.axis = 0.65, mgp = c(2, 0.35, 0))
  abline(h = 0, col = "grey70", lwd = 0.6)
  title(main = main, cex.main = 0.7, line = 0.2, font.main = 1)
  text(x[-1], y[-1], sprintf("%.1f", vals), pos = 3, cex = 0.5, offset = 0.35)
  if (nzchar(ylab)) mtext(ylab, side = 2, line = 1.7, cex = 0.45)
}

plot_suppfig2_row <- function(res, page_title, columns = TL_COLUMNS, radius_labels = NULL) {
  if (length(columns) > 3) {   # generic mode: one page per three metrics
    for (k in seq(1, length(columns), by = 3)) plot_suppfig2_row(res, page_title, columns[k:min(k + 2, length(columns))])
    return(invisible())
  }
  if (is.null(radius_labels)) radius_labels <- ifelse(is.na(radius_of(columns)), columns, as.character(radius_of(columns)))
  layout(matrix(c(1, 2, 3, 7, 1, 2, 3, 8, 1, 2, 3, 9, 4, 5, 6, 10, 4, 5, 6, 11, 4, 5, 6, 12), nrow = 6, byrow = TRUE),
         widths = c(1, 1, 1, 0.78))
  op <- par(oma = c(1.2, 1.5, 2.6, 0.3), mar = c(2.6, 2.6, 1.6, 0.6)); on.exit(par(op))
  for (scheme in c("tertile", "binary")) {
    for (j in 1:3) {
      if (j > length(columns)) { plot.new(); next }
      col <- columns[j]; inf <- res$info[res$info$column == col & res$info$scheme == scheme, ]
      p_star <- if (scheme == "tertile") inf$cox_wald_P_high_vs_low else inf$logrank_P
      note <- if (scheme == "tertile")
        sprintf("n %d, %d deaths; %d/%d/%d; Cox P %s; log-rank P %s",
                inf$n_model, inf$events_model, inf$n_low, inf$n_mid, inf$n_high, signif(inf$cox_wald_P_high_vs_low, 2), signif(inf$logrank_P, 2))
      else sprintf("n %d, %d deaths; top third %d; log-rank P %s", inf$n_model, inf$events_model, inf$n_high, signif(inf$logrank_P, 2))
      draw_km(res$fits[[col]][[scheme]], if (scheme == "tertile") PAL_TERTILE else PAL_BINARY,
              if (scheme == "tertile") c("Low", "Mid", "High") else c("Low", "High"),
              bquote("Tum:Lym (" * G[fx * .(radius_labels[j])] * ")"), stars(p_star), note)
      if (j == 1) mtext("Survival %", side = 2, line = 1.7, cex = 0.6)
    }
  }
  mtext("Time (months); follow-up censored at 1,460 days", side = 1, outer = TRUE, line = 0.1, cex = 0.6, adj = 0.37)
  for (scheme in c("tertile", "binary")) for (j in 1:3) {
    if (j > length(columns)) { plot.new(); next }
    v <- res$dos_fixed[res$dos_fixed$metric == columns[j] & res$dos_fixed$scheme == scheme, "difference"]
    draw_dos(v, bquote(G[fx * .(radius_labels[j])] * "-T:L, " * .(if (scheme == "tertile") "High minus Low" else "High minus (Mid+Low)")),
             if (j == 2) "% difference" else "")
  }
  mtext(page_title, side = 3, outer = TRUE, line = 0.9, cex = 0.8, font = 2)
}

plot_histograms <- function(res, page_title, columns = TL_COLUMNS) {
  radius_labels <- ifelse(is.na(radius_of(columns)), columns, as.character(radius_of(columns)))
  op <- par(mfrow = c(1, min(3, length(columns))), oma = c(0, 0, 2.5, 0), mar = c(4, 4, 2, 1)); on.exit(par(op))
  for (j in seq_along(columns)) {
    x <- res$values[[columns[j]]]; x <- x[!is.na(x)]
    q <- quantile(x, c(1/3, 2/3)); g <- cut(x, c(-Inf, q, Inf), labels = FALSE)
    br <- seq(min(x), max(x), length.out = 31)
    cnt <- sapply(1:3, function(k) if (any(g == k)) hist(x[g == k], breaks = br, plot = FALSE)$counts else rep(0, length(br) - 1))
    plot(NA, xlim = range(br), ylim = c(0, max(rowSums(cnt)) * 1.08), xlab = bquote(G[fx * .(radius_labels[j])] * " Tumor:Lymphoid"),
         ylab = "Count", main = sprintf("n = %d; cut-offs %.3g / %.3g", length(x), q[1], q[2]), cex.main = 0.85, las = 1)
    base <- rep(0, length(br) - 1)
    for (k in 1:3) { rect(br[-length(br)], base, br[-1], base + cnt[, k], col = adjustcolor(PAL_TERTILE[k], 0.75), border = "black", lwd = 0.4); base <- base + cnt[, k] }
    abline(v = q, lty = 2, lwd = 1.2)
  }
  mtext(page_title, side = 3, outer = TRUE, line = 0.8, cex = 0.85, font = 2)
}

# ---- cohort runs ------------------------------------------------------------------------------------------
read_cohort <- function(spec) {
  d <- read.csv(spec$file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(spec$subset)) {
    stopifnot(spec$subset$col %in% names(d))
    d <- d[!is.na(d[[spec$subset$col]]) & d[[spec$subset$col]] == spec$subset$value, , drop = FALSE]
  }
  n_in <- nrow(d); n_dup <- 0L
  if (!is.null(spec$dedup_by)) {   # keep the first row of each key (identifiers are never printed)
    stopifnot(spec$dedup_by %in% names(d))
    n_dup <- sum(duplicated(d[[spec$dedup_by]])); d <- d[!duplicated(d[[spec$dedup_by]]), , drop = FALSE]
  }
  if (!is.null(spec$require_all)) d <- d[stats::complete.cases(d[, spec$require_all, drop = FALSE]), , drop = FALSE]
  out <- data.frame(OS.time = suppressWarnings(as.numeric(d[[spec$time]])) * spec$time_scale,
                    OS = suppressWarnings(as.numeric(d[[spec$event]])))
  for (k in seq_along(spec$metrics)) out[[names(spec$metrics)[k]]] <- suppressWarnings(as.numeric(d[[spec$metrics[[k]]]]))
  attr(out, "n_in") <- n_in; attr(out, "n_dup_removed") <- n_dup
  out
}

run_spec <- function(spec) {
  d <- read_cohort(spec)
  fd <- prepare_os(d)
  res <- comparator_core(fd, names(spec$metrics), names(spec$metrics), day_grid = MONTH_GRID_DAYS, km_steps = FALSE)
  res$spec <- spec; res$n_input <- nrow(d); res$n_dup_removed <- attr(d, "n_dup_removed")
  res$missing_time_set_to_7000 <- sum(is.na(d$OS.time)); res$events_full_followup <- sum(fd$event_indicator)
  cat(sprintf("  %-24s records %d (deaths %d over full follow-up; %d within 1,460 d); missing times set to 7000: %d\n",
              spec$run, res$n_rows, res$events_full_followup, res$n_events, res$missing_time_set_to_7000))
  res
}

tag <- function(df, spec) {
  if (is.null(df) || !nrow(df)) return(df)
  cbind(run = spec$run, cohort = spec$cohort, population = spec$population, df, stringsAsFactors = FALSE)
}

radius_of <- function(col) suppressWarnings(as.integer(ifelse(grepl("\\d", col), sub("^\\D*(\\d+).*$", "\\1", col), NA)))

write_outputs <- function(results, out_dir, make_plots) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  bind <- function(f) do.call(rbind, lapply(results, function(r) tag(f(r), r$spec)))
  rows <- bind(function(r) r$s2)
  rows$radius_um <- radius_of(rows$column)
  rows$P_star <- ifelse(rows$scheme == "tertile", rows$cox_wald_P_high_vs_low, rows$logrank_P)
  rows$P_star_test <- ifelse(rows$scheme == "tertile", "Cox Wald, top vs bottom tertile (3-level model)", "log-rank, top third vs bottom two-thirds")
  rows$star <- stars(rows$P_star)
  rows <- rows[order(match(rows$run, sapply(results, function(r) r$spec$run)), rows$scheme == "binary", rows$radius_um), ]
  lead <- c("run", "cohort", "population", "column", "radius_um", "scheme", "n_model", "events_model", "cutoff_lower", "cutoff_upper",
            "n_low", "n_mid", "n_high", "events_low", "events_mid", "events_high", "HR", "HR_lower", "HR_upper",
            "cox_wald_P_high_vs_low", "logrank_P", "logrank_df", "P_star", "P_star_test", "star",
            "delta_os_pct_day300", "delta_os_pct_day700", "delta_os_pct_day1000", "delta_os_pct_day1300")
  rows <- rows[, c(lead, setdiff(names(rows), c(lead, "metric_display", "stars_logrank", "stars_cox_wald", "plot_title_as_run")))]
  write.csv(rows, file.path(out_dir, "suppfig2_umich_rows.csv"), row.names = FALSE)

  write.csv(bind(function(r) supp_table4_layout(r)), file.path(out_dir, "model_table.csv"), row.names = FALSE)
  write.csv(bind(function(r) r$info), file.path(out_dir, "km_panel_statistics.csv"), row.names = FALSE)
  write.csv(bind(function(r) r$cutoffs), file.path(out_dir, "tertile_cutoffs.csv"), row.names = FALSE)
  write.csv(bind(function(r) r$tertile_differences), file.path(out_dir, "tertile_survival_differences.csv"), row.names = FALSE)
  write.csv(bind(function(r) r$dos_fixed), file.path(out_dir, "delta_os_timepoints.csv"), row.names = FALSE)
  write.csv(bind(function(r) r$dos_curve), file.path(out_dir, "delta_os_monthly.csv"), row.names = FALSE)
  write.csv(bind(function(r) r$km), file.path(out_dir, "km_curves_monthly.csv"), row.names = FALSE)
  write.csv(bind(function(r) r$pairwise), file.path(out_dir, "pairwise_tukey.csv"), row.names = FALSE)
  write.csv(bind(function(r) r$binary_summary), file.path(out_dir, "binary_biomarker_survival_summary.csv"), row.names = FALSE)
  write.csv(bind(function(r) r$tertile_summary), file.path(out_dir, "tertile_biomarker_survival_summary.csv"), row.names = FALSE)

  if (make_plots) for (r in results) {
    pdf(file.path(out_dir, paste0("suppfig2_", r$spec$run, ".pdf")), width = 11, height = 8.5)
    plot_suppfig2_row(r, r$spec$title, names(r$spec$metrics))
    plot_histograms(r, paste0(r$spec$title, ": Gfx distributions and tertile cut-offs"), names(r$spec$metrics))
    dev.off()
  }
  rows
}

manifest_rows <- function(results, files) {
  md5 <- function(f) if (!is.null(f) && file.exists(f)) unname(tools::md5sum(f)) else NA_character_
  rbind(
    do.call(rbind, lapply(names(files), function(k) data.frame(item = c(paste0(k, "_file"), paste0(k, "_md5")),
                                                                  value = c(basename(files[[k]]), md5(files[[k]])), stringsAsFactors = FALSE))),
    do.call(rbind, lapply(results, function(r) data.frame(
      item = paste0(r$spec$run, c("_population", "_records_analysed", "_deaths_full_followup", "_deaths_within_1460d", "_rows_removed_by_dedup")),
      value = c(r$spec$population, r$n_rows, r$events_full_followup, r$n_events, r$n_dup_removed), stringsAsFactors = FALSE))),
    data.frame(item = c("follow_up_cap_days", "missing_time_set_to_days", "dOS_time_points_days", "DEDUP_PATIENTS", "run_time"),
               value = c(CAP_DAYS, MISSING_TIME_DAYS, paste(TIME_POINTS, collapse = "/"), "FALSE (UMICH runs)",
                         format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), stringsAsFactors = FALSE))
}

umich_specs <- function(umich1_file, umich2_file) {
  u1 <- list(file = umich1_file, cohort = "UMICH1", time = "stime", time_scale = DAYS_PER_MONTH, event = "deathstatus",
             metrics = list(Gfx10_Tum_vs_Lym = "Gfx_10um", Gfx20_Tum_vs_Lym = "Gfx_20um", Gfx40_Tum_vs_Lym = "Gfx_40um"))
  u2 <- list(file = umich2_file, cohort = "UMICH2", time = "os_time", time_scale = 1, event = "os_event",
             metrics = list(Gfx10_Tum_vs_Lym = "gfx10", Gfx20_Tum_vs_Lym = "gfx20", Gfx40_Tum_vs_Lym = "gfx40"))
  list(
    modifyList(u1, list(run = "UMICH1_all447", population = "primary: all UMICH1 patients",
                        title = "UMICH1 (training cohort), all 447 patients: Supplementary Fig. 2b-d")),
    modifyList(u2, list(run = "UMICH2_pooled295", population = "primary: all UMICH2 records of the pooled 1,139-record table",
                        title = "UMICH2 (testing cohort), 295 pooled records (20/40 um: 268 with values): Supplementary Fig. 2h-j")),
    modifyList(u2, list(run = "UMICH2_linked268", population = "sensitivity: UMICH2 records with all three radii",
                        require_all = c("gfx10", "gfx20", "gfx40"),
                        title = "UMICH2 sensitivity: the 268 records with Gfx at 10, 20 and 40 um")),
    modifyList(u1, list(run = "UMICH1_oralcavity223", population = "sensitivity: oral-cavity subset",
                        subset = list(col = "site", value = "oral cavity"),
                        title = "UMICH1 sensitivity: oral-cavity patients only (n = 223)")))
}

# ---- main ---------------------------------------------------------------------------------------------------
main <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
  stage_dir <- if (length(script_path) == 1) normalizePath(file.path(dirname(script_path), "..")) else getwd()
  args <- commandArgs(trailingOnly = TRUE)
  opt <- function(name, default = NULL) { i <- which(args == name); if (!length(i)) return(default)
    if (i[1] == length(args)) stop("missing value for ", name); args[i[1] + 1] }
  flag <- function(name) name %in% args
  make_plots <- !flag("--no-plots")
  out_dir <- opt("--out", file.path(stage_dir, "outputs", "07_suppfig2_umich"))

  if (!is.null(opt("--input"))) {   # generic mode: any cohort table
    mets <- strsplit(opt("--metrics", stop("--metrics is required with --input")), ",", fixed = TRUE)[[1]]
    sub_arg <- opt("--subset")
    spec <- list(file = opt("--input"), cohort = opt("--label", "cohort"), run = opt("--label", "cohort"),
                 population = "generic run", title = paste("Generic comparator run:", opt("--label", "cohort")),
                 time = opt("--time", stop("--time is required")), event = opt("--event", stop("--event is required")),
                 time_scale = as.numeric(opt("--time-scale", "1")), dedup_by = opt("--dedup-by"),
                 subset = if (is.null(sub_arg)) NULL else list(col = sub("=.*$", "", sub_arg), value = sub("^[^=]*=", "", sub_arg)),
                 metrics = setNames(as.list(mets), mets))
    res <- list(run_spec(spec))
    write_outputs(res, out_dir, make_plots)
    write.csv(manifest_rows(res, list(input = spec$file)), file.path(out_dir, "run_manifest.csv"), row.names = FALSE)
    writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
    cat("Wrote generic-mode outputs to", out_dir, "\n"); return(invisible())
  }

  umich1_file <- opt("--umich1", Sys.getenv("GFX_UMICH1", file.path(stage_dir, "data", "umich1_recurrence_deid.csv")))
  umich2_file <- opt("--umich2", Sys.getenv("GFX_UMICH2", file.path(stage_dir, "data", "umich2_radius_deid.csv")))
  for (f in c(umich1_file, umich2_file)) if (!file.exists(f))
    stop("Input not found: ", basename(f), ". The UMICH patient-level tables are not distributed with this repository; ",
         "pass --umich1 / --umich2.")

  cat("\n== UMICH runs ==\n")
  specs <- umich_specs(umich1_file, umich2_file)
  results <- lapply(specs, run_spec)
  rows <- write_outputs(results, out_dir, make_plots)
  write.csv(manifest_rows(results, list(umich1 = umich1_file, umich2 = umich2_file)), file.path(out_dir, "run_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))

  cat("\n== Supplementary Fig. 2 UMICH rows: stars (tertile: Cox Wald top vs bottom; binary: log-rank) ==\n")
  show <- rows[, c("run", "scheme", "radius_um", "n_model", "events_model", "HR", "P_star", "star")]
  show$HR <- signif(show$HR, 3); show$P_star <- signif(show$P_star, 3); op <- options(width = 200); on.exit(options(op))
  print(show, row.names = FALSE)
  cat("\nWrote outputs to", out_dir, "\n")
}

if (sys.nframe() == 0L) main()
