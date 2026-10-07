#!/usr/bin/env Rscript
# SPDX-License-Identifier: MIT
#
# 04_figures_fig2b_suppfig3.R
#
# Figure 2b  : Kaplan-Meier overall survival by Gfx10 tertile in UMICH1, TCGA and
#              UMICH2, displayed through 48 months; HR / CI / P are computed from
#              the plotted records.
# Supp Fig 3 : a, pooled within-cohort-standardized Gfx10 tertiles (1,139 records),
#              KM displayed through 60 months; b, thirteen dichotomization cut-points
#              (20th-80th percentile, strict ">" rule); c, restricted cubic (natural)
#              spline relative to the pooled median, with the 2-df non-linearity LRT.
#
# Dependencies: R (tested with 4.4.3; base graphics), survival (tested with 3.8.3)
#               and splines (part of base R).
#
# Usage
#   Rscript 04_figures_fig2b_suppfig3.R [--pooled=PATH] [--tcga=PATH] [--outdir=DIR]
#                                       [--setting=both|tcga_all_rows|pooled_consistent]
#                                       [--conf_type=log-log|log|plain] [--spline_range=p01_p99|full]
#   --outdir        default ../figures relative to this script
#   --setting       Figure 2b cohort populations (see below); default both
#   --conf_type     confidence-interval type of the Kaplan-Meier bands; default log-log
#   --spline_range  x range of Supp Fig 3c: 1st-99th percentile of pooled gfx10_z
#                   (p01_p99, default) or the full range (full)
#
# Inputs
#   --pooled  1,139-record pooled table with columns cohort, gfx10, gfx10_z, os_time
#             (days), os_event (the input of 01_pooled_primary_1139.R). Its UMICH rows
#             are patient-level clinical data, and the table is not distributed with
#             this repository. Default: ../data/pooled_1139_deid.csv relative to this
#             script, or the environment variable GFX_POOLED_1139.
#   --tcga    public TCGA analysis table; default ../data/tcga_analysis_table.csv
#             relative to this script. 540 rows (a TCGA case can have more than one
#             row); columns case, OS, OS.time (days), Gfx10_Tum_vs_Lym.
#
# Cohort populations for Figure 2b (--setting)
#   tcga_all_rows     : UMICH1 447 and UMICH2 295 records from the pooled table;
#                       TCGA = the 539 rows of the public TCGA table that have OS.time,
#                       OS and Gfx10_Tum_vs_Lym (full follow-up).
#   pooled_consistent : 447 / 397 / 295 records exactly as in the 1,139-record pooled
#                       analysis (TCGA = one record per case, case-mean Gfx10, pooled
#                       os_time / os_event).
#   both              : both settings, one PDF each.
#   Tertiles are cohort-specific R type-7 quantiles of Gfx10: low <= lower cut,
#   high > upper cut, middle otherwise. HR = univariable Cox (Efron ties), high vs low
#   with the middle group excluded; P = two-sided Wald. The 3-group log-rank P is also
#   tabulated. Models use full recorded follow-up; only the display is truncated.
#
# Outputs (aggregate only) in --outdir:
#   Fig2b_<setting>.pdf, SuppFig3_pooled1139.pdf,
#   Fig2b_cohort_stats.csv (one row per setting and cohort),
#   SuppFig3_cutpoint_sweep.csv, SuppFig3_summary.csv

suppressPackageStartupMessages(library(survival))
library(splines)

## ---- arguments -------------------------------------------------------------
all_args   <- commandArgs(trailingOnly = FALSE)
file_arg   <- sub("^--file=", "", all_args[grep("^--file=", all_args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
default_pooled <- Sys.getenv("GFX_POOLED_1139",
  file.path(script_dir, "..", "data", "pooled_1139_deid.csv"))
pooled_path <- get_arg("pooled", default_pooled)
tcga_path   <- get_arg("tcga", file.path(script_dir, "..", "data", "tcga_analysis_table.csv"))
out_dir     <- get_arg("outdir", file.path(script_dir, "..", "figures"))
setting_arg <- get_arg("setting", "both")
conf_type   <- get_arg("conf_type", "log-log")
# Supp Fig 3c x-range: 1st-99th percentile of pooled gfx10_z (p01_p99) or the full
# range (full).
spline_range <- get_arg("spline_range", "p01_p99")
stopifnot(spline_range %in% c("p01_p99", "full"))
settings <- if (setting_arg == "both") c("tcga_all_rows", "pooled_consistent") else setting_arg
stopifnot(all(settings %in% c("tcga_all_rows", "pooled_consistent")))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DAYS_PER_MONTH <- 365.25 / 12
COHORTS <- c("UMICH1", "TCGA", "UMICH2")
COHORT_TITLE <- c(UMICH1 = "UMICH1 (training)", TCGA = "TCGA (validation)",
                  UMICH2 = "UMICH2 (testing)")
COL <- c(Low = "#DE3163", Mid = "#4DAF4A", High = "#6495ED")

## ---- data ------------------------------------------------------------------
if (!file.exists(pooled_path)) stop("pooled table not found: ", pooled_path,
  "\nThe patient-level table is not distributed with this repository (the required columns are listed in the header of this script).",
  "\nDe-identified patient-level UMICH data are available from the corresponding author under a data use",
  " agreement (see the manuscript's Data availability statement).")
pooled <- read.csv(pooled_path, stringsAsFactors = FALSE)
stopifnot(all(c("cohort", "gfx10", "gfx10_z", "os_time", "os_event") %in% names(pooled)),
          nrow(pooled) == 1139,
          identical(as.vector(table(pooled$cohort)[COHORTS]), c(447L, 397L, 295L)),
          !anyNA(pooled[, c("gfx10", "gfx10_z", "os_time", "os_event")]))
# gfx10_z must be gfx10 standardized within cohort with the sample SD (denominator
# n - 1), as checked in 01_pooled_primary_1139.R.
z_recomputed <- ave(pooled$gfx10, pooled$cohort, FUN = function(x) (x - mean(x)) / sd(x))
stopifnot(max(abs(z_recomputed - pooled$gfx10_z)) < 1e-12)

tcga <- read.csv(tcga_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(tcga) == 540, all(c("case", "OS", "OS.time", "Gfx10_Tum_vs_Lym") %in% names(tcga)))
tcga$time  <- suppressWarnings(as.numeric(tcga$OS.time))
tcga$event <- suppressWarnings(as.numeric(tcga$OS))
tcga$score <- suppressWarnings(as.numeric(tcga$Gfx10_Tum_vs_Lym))
tcga_curve <- tcga[!is.na(tcga$time) & !is.na(tcga$event) & !is.na(tcga$score),
                   c("time", "event", "score")]

cohort_data <- function(setting) {
  from_pooled <- function(coh) {
    d <- pooled[pooled$cohort == coh, ]
    data.frame(time = d$os_time, event = d$os_event, score = d$gfx10)
  }
  tc <- if (setting == "tcga_all_rows") tcga_curve else from_pooled("TCGA")
  list(UMICH1 = from_pooled("UMICH1"), TCGA = tc, UMICH2 = from_pooled("UMICH2"))
}

## ---- helpers ---------------------------------------------------------------
tertiles <- function(x) {
  q <- quantile(x, c(1/3, 2/3), type = 7, names = FALSE)
  list(group = factor(ifelse(x <= q[1], "Low", ifelse(x > q[2], "High", "Mid")),
                      levels = c("Low", "Mid", "High")), cuts = q)
}
fmt_p <- function(p) if (p < 1e-4) formatC(p, format = "e", digits = 2) else formatC(p, format = "g", digits = 2)
fill_forward <- function(v, init) {
  for (i in seq_along(v)) if (is.na(v[i])) v[i] <- if (i == 1) init else v[i - 1]
  v
}
step_xy <- function(t, y, tmax, y0 = 1) {
  k <- length(t)
  if (k == 0) return(list(x = c(0, tmax), y = c(y0, y0)))
  list(x = c(0, as.vector(rbind(t, t)), tmax),
       y = c(y0, as.vector(rbind(c(y0, y[-k]), y)), y[k]))
}
# Draw one KM curve (with pointwise CI band and censor marks) up to horizon_days.
draw_km <- function(time, event, col, horizon_days) {
  sf <- survfit(Surv(time, event) ~ 1, conf.type = conf_type)
  keep <- sf$time <= horizon_days
  t <- sf$time[keep]; s <- sf$surv[keep]
  lo <- fill_forward(sf$lower[keep], 1); up <- fill_forward(sf$upper[keep], 1)
  tmax <- min(horizon_days, max(time))
  u <- step_xy(t, up, tmax); l <- step_xy(t, lo, tmax); m <- step_xy(t, s, tmax)
  polygon(c(u$x, rev(l$x)) / DAYS_PER_MONTH, c(u$y, rev(l$y)),
          col = adjustcolor(col, alpha.f = 0.18), border = NA)
  lines(m$x / DAYS_PER_MONTH, m$y, col = col, lwd = 1.6)
  cens <- sf$n.censor[keep] > 0
  points(t[cens] / DAYS_PER_MONTH, s[cens], pch = 3, cex = 0.45, col = col)
}
km_panel <- function(d, group, horizon_months, ticks, main, annot) {
  plot(NA, xlim = c(0, horizon_months), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       xaxt = "n", las = 1, xlab = "Months", ylab = "Overall survival", main = main,
       cex.main = 0.95, font.main = 1)
  axis(1, at = ticks)
  abline(v = ticks, h = c(0.25, 0.5, 0.75), col = "grey92", lwd = 0.6)
  for (g in c("Low", "Mid", "High")) {
    idx <- group == g
    if (any(idx)) draw_km(d$time[idx], d$event[idx], COL[g], horizon_months * DAYS_PER_MONTH)
  }
  box()
  n <- table(group)
  legend("bottomleft", bty = "n", cex = 0.75, lwd = 1.6, col = COL[c("High", "Mid", "Low")],
         legend = sprintf("%s (n = %d)", c("High", "Mid", "Low"), as.integer(n[c("High", "Mid", "Low")])))
  mtext(annot, side = 3, line = -2.2, adj = 0.97, cex = 0.62)
}
cohort_stats <- function(d) {
  tg <- tertiles(d$score); d$group <- tg$group
  hl <- d[d$group != "Mid", ]; hl$high <- as.integer(hl$group == "High")
  fit <- coxph(Surv(time, event) ~ high, data = hl, ties = "efron")
  s <- summary(fit)
  lr3 <- survdiff(Surv(time, event) ~ group, data = d)
  lrhl <- survdiff(Surv(time, event) ~ high, data = hl)
  km_at <- function(g) summary(survfit(Surv(time, event) ~ 1, data = d[d$group == g, ]),
                               times = 1460, extend = TRUE)$surv
  list(group = d$group,
       row = data.frame(n = nrow(d), events = sum(d$event),
         n_low = sum(d$group == "Low"), n_mid = sum(d$group == "Mid"), n_high = sum(d$group == "High"),
         lower_cut = tg$cuts[1], upper_cut = tg$cuts[2],
         model_n = fit$n, model_events = fit$nevent,
         HR = s$conf.int[1, "exp(coef)"], lower = s$conf.int[1, "lower .95"],
         upper = s$conf.int[1, "upper .95"], wald_p = s$coefficients[1, "Pr(>|z|)"],
         logrank_3group_p = pchisq(lr3$chisq, length(lr3$n) - 1, lower.tail = FALSE),
         logrank_high_low_p = pchisq(lrhl$chisq, 1, lower.tail = FALSE),
         km_1460d_low = km_at("Low"), km_1460d_high = km_at("High")))
}

## ---- Figure 2b -------------------------------------------------------------
stat_rows <- list()
for (setting in settings) {
  cd <- cohort_data(setting)
  pdf(file.path(out_dir, sprintf("Fig2b_%s.pdf", setting)), width = 11, height = 3.9,
      useDingbats = FALSE)
  par(mfrow = c(1, 3), mar = c(4, 4.2, 2.6, 0.8), mgp = c(2.3, 0.7, 0), cex.axis = 0.85)
  for (coh in COHORTS) {
    st <- cohort_stats(cd[[coh]]); r <- st$row
    stat_rows[[length(stat_rows) + 1]] <- data.frame(setting = setting, cohort = coh, r,
      row.names = NULL)
    annot <- sprintf("High vs low HR %.2f (%.2f-%.2f)\nWald P = %s; log-rank (3 groups) P = %s",
                     r$HR, r$lower, r$upper, fmt_p(r$wald_p), fmt_p(r$logrank_3group_p))
    km_panel(cd[[coh]], st$group, 48, c(0, 16, 32, 48),
             sprintf("%s: %d records, %d deaths", COHORT_TITLE[[coh]], r$n, r$events), annot)
  }
  mtext(sprintf("Fig. 2b, setting = %s (curves to 48 months; Cox on full follow-up; CI %s)",
                setting, conf_type), side = 1, outer = TRUE, line = -1.0, cex = 0.6, col = "grey30")
  dev.off()
}
fig2b <- do.call(rbind, stat_rows)
write.csv(fig2b, file.path(out_dir, "Fig2b_cohort_stats.csv"), row.names = FALSE)

## ---- Supplementary Figure 3 (pooled 1,139) ---------------------------------
d <- pooled
d$cohort <- factor(d$cohort, levels = COHORTS)
# a: pooled standardized-score tertiles
q <- quantile(d$gfx10_z, c(1/3, 2/3), type = 7, names = FALSE)
d$tertile <- factor(ifelse(d$gfx10_z <= q[1], "Low", ifelse(d$gfx10_z > q[2], "High", "Mid")),
                    levels = c("Low", "Mid", "High"))
hl <- d[d$tertile != "Mid", ]; hl$high <- as.integer(hl$tertile == "High")
ft <- coxph(Surv(os_time, os_event) ~ high + strata(cohort), data = hl, ties = "efron")
st <- summary(ft)
tert <- data.frame(n_low = sum(d$tertile == "Low"), n_mid = sum(d$tertile == "Mid"),
                   n_high = sum(d$tertile == "High"), model_n = ft$n, model_events = ft$nevent,
                   HR = st$conf.int[1, 1], lower = st$conf.int[1, 3], upper = st$conf.int[1, 4],
                   p = st$coefficients[1, 5])
# b: thirteen strict cut-points
cuts <- do.call(rbind, lapply(seq(20, 80, 5), function(pct) {
  thr <- quantile(d$gfx10_z, pct / 100, type = 7, names = FALSE)
  d$high <- as.integer(d$gfx10_z > thr)
  s <- summary(coxph(Surv(os_time, os_event) ~ high + strata(cohort), data = d, ties = "efron"))
  data.frame(percentile = pct, threshold_z = thr, n_high = sum(d$high), HR = s$conf.int[1, 1],
             lower = s$conf.int[1, 3], upper = s$conf.int[1, 4], p = s$coefficients[1, 5])
}))
# c: natural cubic spline, knots at the 5th/35th/65th/95th percentiles
k <- quantile(d$gfx10_z, c(.05, .35, .65, .95), type = 7, names = FALSE)
f_lin <- coxph(Surv(os_time, os_event) ~ gfx10_z + strata(cohort), data = d, ties = "efron")
f_spl <- coxph(Surv(os_time, os_event) ~ ns(gfx10_z, knots = k[2:3], Boundary.knots = k[c(1, 4)]) +
                 strata(cohort), data = d, ties = "efron")
lrt <- 2 * (f_spl$loglik[2] - f_lin$loglik[2]); lrt_df <- length(coef(f_spl)) - length(coef(f_lin))
spline_p <- pchisq(lrt, lrt_df, lower.tail = FALSE)
med <- median(d$gfx10_z)
xr <- if (spline_range == "full") range(d$gfx10_z) else quantile(d$gfx10_z, c(.01, .99), type = 7, names = FALSE)
xg <- seq(xr[1], xr[2], length.out = 300)
basis <- function(x) ns(x, knots = k[2:3], Boundary.knots = k[c(1, 4)])
C <- sweep(basis(xg), 2, as.vector(basis(med)))
lp <- as.vector(C %*% coef(f_spl)); se <- sqrt(rowSums((C %*% vcov(f_spl)) * C))
z <- qnorm(0.975)

pdf(file.path(out_dir, "SuppFig3_pooled1139.pdf"), width = 12, height = 4, useDingbats = FALSE)
par(mfrow = c(1, 3), mar = c(4, 4.2, 2.6, 0.8), mgp = c(2.3, 0.7, 0), cex.axis = 0.85)
km_panel(data.frame(time = d$os_time, event = d$os_event), d$tertile, 60, seq(0, 60, 12),
         sprintf("a  All three cohorts pooled: %d records, %d deaths", nrow(d), sum(d$os_event)),
         sprintf("High vs low HR %.2f (%.2f-%.2f), P = %s\n(cohort-stratified Cox, full follow-up)",
                 tert$HR, tert$lower, tert$upper, fmt_p(tert$p)))
plot(cuts$percentile, cuts$HR, log = "y", ylim = c(0.3, 1.1), xlim = c(18, 82), pch = 19,
     las = 1, xlab = "Dichotomization cut-point (percentile)", ylab = "Hazard ratio, high vs low",
     main = "b  Cut-point robustness", cex.main = 0.95, font.main = 1)
abline(v = c(25, 100 / 3, 50, 200 / 3, 75), lty = 3, col = "grey40")
abline(h = 1, lty = 2)
arrows(cuts$percentile, cuts$lower, cuts$percentile, cuts$upper, angle = 90, code = 3,
       length = 0.025)
points(cuts$percentile, cuts$HR, pch = 19)
axis(3, at = c(25, 100 / 3, 50, 200 / 3, 75), labels = c("Q1", "T1", "Median", "T2", "Q3"),
     tick = FALSE, line = -0.9, cex.axis = 0.7)
mtext(sprintf("HR %.2f-%.2f; max P = %s", min(cuts$HR), max(cuts$HR), fmt_p(max(cuts$p))),
      side = 1, line = -1.5, adj = 0.95, cex = 0.62)
plot(xg, exp(lp), type = "n", log = "y", ylim = c(0.3, 2.4), las = 1,
     xlab = "Gfx10 (SD from the within-cohort mean)", ylab = "Hazard ratio vs median",
     main = sprintf("c  Dose-response (natural spline), %s", if (spline_range == "full") "full Gfx10 range" else "1st-99th percentile"),
     cex.main = 0.95, font.main = 1)
polygon(c(xg, rev(xg)), c(exp(lp - z * se), rev(exp(lp + z * se))),
        col = adjustcolor("grey50", alpha.f = 0.3), border = NA)
lines(xg, exp(lp), lwd = 1.6)
abline(h = 1, lty = 2); abline(v = c(k[2], med, k[3]), lty = 3, col = "grey40")
axis(3, at = med, labels = "median", tick = FALSE, line = -0.9, cex.axis = 0.7)
mtext(sprintf("x range: %s; dotted lines: 35th/65th-percentile knots and median", if (spline_range == "full") "full" else "1st-99th percentile"),
      side = 1, line = -2.5, adj = 0.95, cex = 0.55, col = "grey30")
mtext(sprintf("non-linearity P = %.2f (LRT, %d df)", spline_p, lrt_df), side = 1, line = -1.5,
      adj = 0.95, cex = 0.62)
dev.off()

summ <- rbind(
  data.frame(item = "tertile_n_low_mid_high", value = paste(tert$n_low, tert$n_mid, tert$n_high, sep = "/")),
  data.frame(item = "tertile_model_n_events", value = paste(tert$model_n, tert$model_events, sep = "/")),
  data.frame(item = "tertile_HR_CI", value = sprintf("%.10f (%.10f-%.10f)", tert$HR, tert$lower, tert$upper)),
  data.frame(item = "tertile_p", value = format(tert$p, digits = 10)),
  data.frame(item = "sweep_HR_range", value = sprintf("%.10f-%.10f", min(cuts$HR), max(cuts$HR))),
  data.frame(item = "sweep_max_p", value = format(max(cuts$p), digits = 10)),
  data.frame(item = "spline_lrt_df_p", value = sprintf("%.10f / %d / %.10f", lrt, lrt_df, spline_p)),
  data.frame(item = "spline_knots_z", value = paste(sprintf("%.6f", k), collapse = ";")))
write.csv(cuts, file.path(out_dir, "SuppFig3_cutpoint_sweep.csv"), row.names = FALSE)
write.csv(summ, file.path(out_dir, "SuppFig3_summary.csv"), row.names = FALSE)

## ---- console summary -------------------------------------------------------
print(fig2b[, c("setting", "cohort", "n", "events", "n_low", "n_mid", "n_high", "HR", "lower",
                "upper", "wald_p", "logrank_3group_p")], digits = 6)
