#!/usr/bin/env Rscript
# SPDX-License-Identifier: MIT
# =============================================================================
# 03_recurrence_umich1.R -- recurrence or persistent disease, UMICH1
#                           (Supplementary Table 5)
# =============================================================================
# Computes the UMICH1 analysis of recurrence or persistent disease
# (Supplementary Table 5): in the complete cases, one multivariable Cox model
# per Gfx measure (10, 20, 30, 40 um and the complete G-cross curve), each
# adjusted for ACE score, the joint anatomic-site/HPV category, AJCC
# 7th-edition stage and smoking status; Schoenfeld tests (cox.zph,
# transform = 'km'); the stage-stratified Gfx10 sensitivity model and two
# other bounded Gfx10 sensitivity models; and, for comparison, the same models
# for overall survival in the same patients.
#
# Usage:  Rscript 03_recurrence_umich1.R [input_csv] [output_dir]
#   input_csv   default ../data/umich1_recurrence_deid.csv (relative to this
#               file)
#   output_dir  default ../outputs (relative to this file)
#
# Input: the UMICH1 table, one row per patient (447 rows). It contains
# patient-level clinical data and is not distributed with this repository.
# Columns used:
#   idnum               record code (must be unique; used only for that check)
#   Gfx_10um, Gfx_20um, Gfx_30um, Gfx_40um, Gfx_completecurve
#                       Tumor:Lymphoid G-cross AUC at 10/20/30/40 um and the
#                       complete-curve summary (hazard ratios are per unit)
#   stime, deathstatus  overall-survival time (months) and death (0/1)
#   rtime, recurstatus  time (months) and event (0/1) for recurrence or
#                       persistent disease
#   persist             1 = persistent disease, 0 = not
#   ace_overall_score   Adult Comorbidity Evaluation: none (reference), mild,
#                       moderate, severe
#   site                joint anatomic-site/HPV category: 3.5 = HPV-positive
#                       oropharynx (reference), hypopharynx, larynx,
#                       oral cavity, oropharynx (HPV-negative or unknown)
#   stage               AJCC 7th edition, 0-4 (stage 0 is grouped with I)
#   smoker              0 = never (reference), 1 = current, 2 = former
#
# Endpoint: the pair (rtime, recurstatus) counts persistent disease as an
# event at ~1 day (0.0329 months) and censors deaths without recorded
# recurrence at the recorded recurrence follow-up. It is "recurrence or
# persistent disease", not recurrence-free survival. stime and rtime are in
# months (days / 30.4375).
#
# Output: <output_dir>/recurrence_umich1_results.csv, aggregate results only:
# one long table with a 'block' column (audit_counts, audit_eligibility,
# coefficients, ph_checks_km_transform, fit_checks). The coefficient rows of
# the main recurrence models carry the Supplementary Table 5 column and row
# labels (supp_table5_column, supp_table5_row). Summaries and sessionInfo()
# are printed to the console.
#
# Dependencies: R (tested with 4.4.3) with the recommended package survival
# (tested with 3.8.3).
# =============================================================================
suppressPackageStartupMessages(library(survival))
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grep('^--file=',args)])
here <- if (length(self)==1) dirname(normalizePath(self)) else getwd()
targs <- commandArgs(trailingOnly=TRUE)
input <- if (length(targs)>=1) targs[1] else file.path(here,'..','data','umich1_recurrence_deid.csv')
out   <- if (length(targs)>=2) targs[2] else file.path(here,'..','outputs')
if (!file.exists(input)) stop('input table not found: ',input,
  '\nThe patient-level table is not distributed with this repository (the required columns are listed in the header of this script).',
  '\nDe-identified patient-level UMICH data are available from the corresponding author under a data use',
  ' agreement (see the manuscript\'s Data availability statement).')
dir.create(out, recursive=TRUE, showWarnings=FALSE)

# ---- data, covariate coding and eligibility ----
source <- read.csv(input, check.names=FALSE, stringsAsFactors=FALSE,
                   na.strings=c("", "NA", "NaN"))
stopifnot(nrow(source)==447, !anyNA(source$idnum), !anyDuplicated(source$idnum))
gfx <- c("Gfx_10um","Gfx_20um","Gfx_30um","Gfx_40um","Gfx_completecurve")
fields <- c(gfx,"stime","deathstatus","rtime","recurstatus","persist",
            "ace_overall_score","site","stage","smoker")
d <- source[,fields]
rm(source)
stopifnot(all(d$deathstatus %in% 0:1), all(d$recurstatus %in% 0:1),
          all(d$persist %in% 0:1), all(d$stime>0), all(d$rtime>0),
          all(d$rtime <= d$stime + 1e-7), all(d$stage %in% 0:4),
          all(na.omit(d$smoker) %in% 0:2))
d$ace <- factor(d$ace_overall_score, levels=c("none","mild","moderate","severe"))
d$joint_site <- factor(d$site,levels=c("3.5","hypopharynx","larynx","oral cavity","oropharynx"))
d$ajcc <- factor(ifelse(d$stage %in% c(0,1),"0_I",as.character(d$stage)),levels=c("0_I","2","3","4"))
d$smoking <- factor(d$smoker,levels=0:2)
complete <- complete.cases(d[,c(gfx,"stime","deathstatus","ace","joint_site","ajcc","smoking")])
# Complete cases: the eligibility rule of the overall-survival models; the
# recurrence times, statuses and persistence flags are complete among them.
stopifnot(sum(complete)==427,
          all(complete.cases(d[complete,c("rtime","recurstatus","persist")])) )
audit <- function(x,label) {
  no_recur_death <- x$deathstatus==1 & x$recurstatus==0
  data.frame(population=label,n=nrow(x),deaths=sum(x$deathstatus),
    stored_recurrence_events=sum(x$recurstatus),persistent_disease=sum(x$persist),
    recurrence_events_without_persistence=sum(x$recurstatus[x$persist==0]),
    deaths_without_recurrence=sum(no_recur_death),
    deaths_without_recurrence_at_recorded_censor=sum(no_recur_death & abs(x$stime-x$rtime)<1e-7),
    deaths_without_recurrence_after_recorded_censor=sum(no_recur_death & x$stime>x$rtime+1e-7),
    persistent_at_one_day=sum(x$persist==1 & x$recurstatus==1 & abs(x$rtime-0.032854209)<1e-8))
}
aud <- rbind(audit(d,"all_source_records"),audit(d[complete,],"same_complete_cases_as_OS"))
# Eligibility audit (counts only; does not change the analysis set).
elig <- data.frame(population=c(rep("all_source_records",8),"same_complete_cases_as_OS"),
  measure=c("records","missing_any_gfx","missing_stime_or_deathstatus","ace_missing_or_unrecognised",
            "site_missing_or_unrecognised","stage_missing","smoking_missing","excluded_incomplete",
            "stage0_recorded"),
  value=c(nrow(d),sum(!complete.cases(d[,gfx])),sum(!complete.cases(d[,c("stime","deathstatus")])),
          sum(is.na(d$ace)),sum(is.na(d$joint_site)),sum(is.na(d$ajcc)),sum(is.na(d$smoking)),
          sum(!complete),sum(d$stage[complete]==0)))
d <- d[complete,]; rownames(d) <- NULL
stopifnot(sum(d$deathstatus)==129, sum(d$recurstatus)==118, sum(d$persist)==36)
models <- list(); rows <- list(); phrows <- list(); checks <- list()
fit_one <- function(label,data,g,time,event,covariates=c("ace","joint_site","ajcc","smoking")) {
  f <- reformulate(c(g,covariates), response=paste0("Surv(",time,", ",event,")"))
  warnings <- character()
  fit <- withCallingHandlers(coxph(f,data=data,ties="efron",na.action=na.fail,
            singular.ok=FALSE,x=TRUE,y=TRUE,model=TRUE),
         warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  stopifnot(all(is.finite(coef(fit))),all(is.finite(vcov(fit))),all(diag(vcov(fit))>0))
  sm <- summary(fit); co <- sm$coefficients; ci <- sm$conf.int
  key <- paste(label,g,sep="__")
  rows[[key]] <<- data.frame(analysis=label,predictor=g,term=rownames(co),n=fit$n,events=fit$nevent,
       coefficient=co[,"coef"],SE=co[,"se(coef)"],HR=co[,"exp(coef)"],
       CI_low=ci[,"lower .95"],CI_high=ci[,"upper .95"],P=co[,"Pr(>|z|)"],row.names=NULL)
  ph <- cox.zph(fit, transform="km", terms=TRUE, global=TRUE)
  phrows[[key]] <<- data.frame(analysis=label,predictor=g,term=rownames(ph$table),
       chisq=ph$table[,"chisq"],df=ph$table[,"df"],P=ph$table[,"p"],row.names=NULL)
  checks[[key]] <<- data.frame(analysis=label,predictor=g,n=fit$n,events=fit$nevent,
       iterations=fit$iter,coefficients=length(coef(fit)),warnings=paste(warnings,collapse=" | "))
  models[[key]] <<- fit
}
for(g in gfx) {
  fit_one("OS",d,g,"stime","deathstatus")
  fit_one("stored_recurrence_or_persistence",d,g,"rtime","recurstatus")
}
# Three bounded, exploratory Gfx10 checks chosen for endpoint/PH concerns,
# not by searching for favourable statistics:
# 1. Remove recorded persistence, which is not a new recurrence.
fit_one("exclude_recorded_persistent_disease",d[d$persist==0,],"Gfx_10um","rtime","recurstatus")
# 2. Count deaths within documented recurrence follow-up; do not impute further
# recurrence-free observation for the two deaths after rtime. This derived
# composite is an endpoint-definition sensitivity analysis.
d$composite_event <- as.integer(d$recurstatus==1 | (d$deathstatus==1 & d$stime<=d$rtime+1e-7))
stopifnot(sum(d$composite_event)==153)
fit_one("add_deaths_within_recorded_recurrence_followup",d,"Gfx_10um","rtime","composite_event")
# 3. The main recurrence fit has a nominal stage PH signal; let stage groups
# have distinct baseline hazards while retaining the other covariates.
fit_one("stored_recurrence_stage_stratified",d,"Gfx_10um","rtime","recurstatus",
        c("ace","joint_site","smoking","strata(ajcc)"))
coefficients <- do.call(rbind,rows)
# ---- combine the tables into one long, aggregate-only table ----

rownames(coefficients) <- NULL
t5_col <- c(Gfx_10um="10 um",Gfx_20um="20 um",Gfx_30um="30 um",Gfx_40um="40 um",Gfx_completecurve="Complete curve")
t5_row <- c(acemild="ACE mild",acemoderate="ACE moderate",acesevere="ACE severe",
  joint_sitehypopharynx="Hypopharynx",joint_sitelarynx="Larynx","joint_siteoral cavity"="Oral cavity",
  joint_siteoropharynx="Oropharynx HPV-negative or unknown",ajcc2="AJCC stage II",ajcc3="AJCC stage III",
  ajcc4="AJCC stage IV",smoking1="Current smoker",smoking2="Former smoker")
main <- coefficients$analysis=="stored_recurrence_or_persistence"
coefficients$supp_table5_column <- ifelse(main,unname(t5_col[coefficients$predictor]),NA)
coefficients$supp_table5_row <- ifelse(main,ifelse(coefficients$term==coefficients$predictor,
  "Gfx score per unit increase",unname(t5_row[coefficients$term])),NA)
stopifnot(!anyNA(coefficients$supp_table5_row[main]), sum(main)==65)
coefficients$block <- "coefficients"
ph_all <- do.call(rbind,phrows); rownames(ph_all) <- NULL; ph_all$block <- "ph_checks_km_transform"
fc <- do.call(rbind,checks); rownames(fc) <- NULL; names(fc)[names(fc)=="coefficients"] <- "n_coefficients"; fc$block <- "fit_checks"
aud_long <- do.call(rbind,lapply(seq_len(nrow(aud)),function(i)
  data.frame(population=aud$population[i],measure=names(aud)[-1],value=as.numeric(unlist(aud[i,-1])),row.names=NULL)))
aud_long$block <- "audit_counts"; elig$block <- "audit_eligibility"
bind_union <- function(l) { cols <- unique(unlist(lapply(l,names)))
 do.call(rbind,lapply(l,function(x) { for(cn in setdiff(cols,names(x))) x[[cn]] <- NA; x[cols] })) }
res <- bind_union(list(aud_long,elig,coefficients,ph_all,fc))
lead <- c("block","population","measure","value","analysis","predictor","term","supp_table5_column","supp_table5_row",
  "n","events","coefficient","SE","HR","CI_low","CI_high","P","chisq","df","iterations","n_coefficients","warnings")
stopifnot(setequal(lead,names(res)))
res <- res[,lead]; rownames(res) <- NULL
outfile <- file.path(out,"recurrence_umich1_results.csv")
write.csv(res,outfile,row.names=FALSE)

print(aud,row.names=FALSE); print(elig[,c("population","measure","value")],row.names=FALSE)
print(coefficients[coefficients$term==coefficients$predictor,c("analysis","predictor","n","events","HR","CI_low","CI_high","P")],digits=8,row.names=FALSE)
print(ph_all[ph_all$analysis=="stored_recurrence_or_persistence" & ph_all$predictor=="Gfx_10um",c("term","chisq","df","P")],digits=8,row.names=FALSE)
cat("Wrote",normalizePath(outfile),"\n")
print(sessionInfo())
