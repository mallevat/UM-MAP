#!/usr/bin/env Rscript
# SPDX-License-Identifier: MIT
#
# 05_tcga_genomic_suppfig6.R
#
# Supplementary Fig. 6b-d, lower panels: TCGA-HNSC overall survival by tertile of tumour
# mutation burden (TMB), predicted neoantigen load and PD-L1 (CD274) expression, with a
# 3-group log-rank P.
#
# Method (same data preparation as 02_tcga_comparators.R):
#   - all rows of the input table are analysed (no de-duplication);
#   - a missing OS.time is set to 7000 days, so that row is censored at the cap;
#   - follow-up is capped at 1,460 days (deaths after day 1,460 are censored);
#   - for each metric, rows with a non-missing value; tertiles from R type-7 quantiles:
#     Low <= q(1/3) < Mid <= q(2/3) < High;
#   - P = 3-group log-rank test (survdiff, 2 df).
#
# Usage:
#   Rscript survival/05_tcga_genomic_suppfig6.R [INPUT_CSV] [OUT_DIR] [--no-plots]
#     INPUT_CSV   input table; default data/tcga_analysis_table.csv
#     OUT_DIR     output folder; default outputs/05_tcga_genomic_suppfig6
#                 (both defaults are relative to the folder above this script's folder)
#     --no-plots  skip the PDF plots
#
# Input: CSV with columns case (patient identifier), OS (1 death, 0 censored), OS.time (days),
#   TMB_mut_per_Mb, neoantigen_count and CD274. Other columns are ignored.
#
# Outputs (in OUT_DIR):
#   suppfig6_survival_logrank.csv  per metric: rows, cases, events, tertile cut-offs, group
#                                  sizes and events, log-rank chi-square, df and P
#   suppfig6_km_curves.csv         Kaplan-Meier step functions by tertile
#   suppfig6_km_<metric>.pdf       Kaplan-Meier plots (unless --no-plots)
#   session_info.txt               package versions
#
# Dependencies: R with survival; survminer and ggplot2 for the plots.
#   Tested with R 4.4.3, survival 3.8-3, survminer 0.5.2 and ggplot2 4.0.3.

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
stage_dir <- if (length(script_path) == 1) normalizePath(file.path(dirname(script_path), "..")) else getwd()
args <- commandArgs(trailingOnly = TRUE)
flags <- args[grepl("^--", args)]
pos <- args[!grepl("^--", args)]
IN_FILE <- if (length(pos) >= 1) pos[1] else file.path(stage_dir, "data", "tcga_analysis_table.csv")
OUT_DIR <- if (length(pos) >= 2) pos[2] else file.path(stage_dir, "outputs", "05_tcga_genomic_suppfig6")
MAKE_PLOTS <- !("--no-plots" %in% flags)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

library(survival)
if (MAKE_PLOTS) {
  library(survminer)
  library(ggplot2)
}

CAP_DAYS <- 1460
metrics <- c(TMB_mut_per_Mb = "Tumor mutation burden (mutations per Mb)",
             neoantigen_count = "Predicted neoantigen load",
             CD274 = "PD-L1 (CD274) expression")

raw <- read.csv(IN_FILE, check.names = FALSE, stringsAsFactors = FALSE)
cat("Input:", basename(IN_FILE), "-", nrow(raw), "rows,", length(unique(raw$case)), "cases\n")

# ---- survival data preparation (as in 02_tcga_comparators.R) -------------------
prepare <- function(d, rowset = "all_rows", missing_time = "set7000", cap = CAP_DAYS) {
  if (rowset == "first_row_per_case") d <- d[!duplicated(d$case), ]
  d$OS.time <- as.numeric(as.character(d$OS.time))
  if (missing_time == "drop") d <- d[!is.na(d$OS.time), ] else d$OS.time[is.na(d$OS.time)] <- 7000
  d$event_indicator <- ifelse(d$OS == 1, 1, ifelse(d$OS == 0, 0, NA))
  d <- d[!is.na(d$event_indicator), ]
  d$adjusted_event_indicator <- ifelse(d$OS.time > cap, 0, d$event_indicator)
  d$adjusted_days_to_death <- pmin(d$OS.time, cap)
  d
}
tertiles <- function(x, method = "quantile_cut") {
  if (method == "quantile_cut") {
    q <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE)
    factor(cut(x, breaks = c(-Inf, q, Inf), labels = c("Low", "Mid", "High")), levels = c("Low", "Mid", "High"))
  } else {  # equal-size groups by rank (ties broken by row order), as dplyr::ntile(x, 3) would give
    r <- rank(x, ties.method = "first")
    factor(c("Low", "Mid", "High")[ceiling(3 * r / length(x))], levels = c("Low", "Mid", "High"))
  }
}

# ---- 3-group log-rank test per metric ------------------------------------------------
fd <- prepare(raw)
prim <- list(); km <- list()
for (m in names(metrics)) {
  x <- as.numeric(fd[[m]])
  d <- fd[!is.na(x) & is.finite(x), ]
  d$value <- as.numeric(d[[m]])
  d$tertile <- tertiles(d$value)
  lr <- survdiff(Surv(adjusted_days_to_death, adjusted_event_indicator) ~ tertile, data = d)
  p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  q <- quantile(d$value, probs = c(1/3, 2/3))
  n_g <- table(d$tertile); e_g <- tapply(d$adjusted_event_indicator, d$tertile, sum)
  prim[[m]] <- data.frame(metric = m, label = metrics[[m]], n_rows = nrow(d), n_cases = length(unique(d$case)),
    events = sum(d$adjusted_event_indicator), cutoff_low = unname(q[1]), cutoff_high = unname(q[2]),
    n_low = unname(n_g["Low"]), n_mid = unname(n_g["Mid"]), n_high = unname(n_g["High"]),
    events_low = unname(e_g["Low"]), events_mid = unname(e_g["Mid"]), events_high = unname(e_g["High"]),
    logrank_chisq = lr$chisq, logrank_df = length(lr$n) - 1, logrank_P = p, stringsAsFactors = FALSE)
  fit <- survfit(Surv(adjusted_days_to_death, adjusted_event_indicator) ~ tertile, data = d)
  s <- summary(fit, censored = TRUE)
  km[[m]] <- data.frame(metric = m, group = sub("^tertile=", "", as.character(s$strata)), time_days = s$time,
    n_risk = s$n.risk, n_event = s$n.event, n_censor = s$n.censor, surv = s$surv,
    lower_95 = s$lower, upper_95 = s$upper, stringsAsFactors = FALSE)
  if (MAKE_PLOTS) {
    g <- ggsurvplot(fit, data = d, conf.int = TRUE, risk.table = TRUE, ggtheme = theme_bw(),
                    palette = c('#DE3163', '#64ed6b', '#6495ED'), xlab = "Days", ylab = "Survival probability",
                    legend.labs = c("Low", "Mid", "High"),
                    title = paste0(metrics[[m]], ": log-rank P = ", signif(p, 3)))
    pdf(file.path(OUT_DIR, paste0("suppfig6_km_", m, ".pdf")), width = 6, height = 7)
    print(g)
    dev.off()
  }
}
prim <- do.call(rbind, prim)
write.csv(prim, file.path(OUT_DIR, "suppfig6_survival_logrank.csv"), row.names = FALSE)
write.csv(do.call(rbind, km), file.path(OUT_DIR, "suppfig6_km_curves.csv"), row.names = FALSE)
print(prim[, c("metric", "n_rows", "events", "logrank_P")], digits = 4)

# Host-specific lines are left out of session_info.txt: the BLAS/LAPACK library
# paths are reduced to file names, and the time-zone lines are dropped.
si <- capture.output(sessionInfo())
si <- sub("^(BLAS|LAPACK):(\\s+)\\S*/", "\\1:\\2", si)
si <- si[!grepl("^(time zone|tzcode source):", si)]
si <- si[!(si == "" & c(FALSE, head(si, -1) == ""))]
writeLines(si, file.path(OUT_DIR, "session_info.txt"))
