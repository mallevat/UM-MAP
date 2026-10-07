#!/usr/bin/env Rscript
# SPDX-License-Identifier: MIT
# =============================================================================
# 01_pooled_primary_1139.R -- pooled overall-survival analysis
# =============================================================================
# Computes the pooled overall-survival analysis of Gfx10 in the 1,139-record
# UMICH1/TCGA/UMICH2 table, including the statistics shown in Supplementary
# Fig. 3:
#   * cohort-stratified continuous Cox model of gfx10_z (Efron ties)
#   * the same model adjusted for age and HPV status (complete records)
#   * rank-transformed scaled-Schoenfeld tests (cox.zph, transform = 'rank')
#     for both models
#   * median, tertile and quartile high-versus-low contrasts, with group sizes
#     and event counts
#   * 13-threshold cut-point sweep (20th to 80th percentile, 5-point steps)
#   * natural (restricted) cubic spline with knots at the 5th/35th/65th/95th
#     percentiles, tested against the linear model by a 2-df likelihood-ratio test
#   * cohort-specific continuous and tertile models, record/patient/event counts
#     per cohort, and Kaplan-Meier survival at 1,095 and 1,460 days by pooled
#     tertile
#   * the number of administratively assigned censoring times per cohort
#
# Usage:  Rscript 01_pooled_primary_1139.R [input_csv] [output_dir]
#   input_csv   default ../data/pooled_1139_deid.csv (relative to this file)
#   output_dir  default ../outputs (relative to this file)
#
# Input: the pooled table, one row per record (1,139 rows). Its UMICH rows are
# patient-level clinical data, and the table is not distributed with this
# repository. Columns used:
#   cluster_id      patient code; records of the same patient share it
#   cohort          UMICH1, TCGA or UMICH2
#   gfx10           Gfx10, the Tumor:Lymphoid G-cross AUC at 10 um
#   gfx10_z         gfx10 standardized within cohort over all of that cohort's
#                   records: (gfx10 - cohort mean) / cohort sample SD
#                   (denominator n - 1); the script checks this
#   os_time         overall-survival time (days)
#   os_event        1 = death, 0 = censored
#   age             years; NA where not recorded
#   hpv             0/1; NA where not recorded
#   admin_censored  TRUE where the censoring time was assigned administratively
#                   because no last-contact date was recorded (the record is
#                   censored at that time), FALSE otherwise
#
# Output: <output_dir>/pooled_1139_results.csv, aggregate statistics only: one
# long table with a 'block' column (model, cutpoint_sweep, diagnostic,
# cohort_counts, km_pooled_tertile, admin_censoring); 'reported_as' names the
# figure or text item that a row corresponds to. Key results and
# sessionInfo() are printed to the console.
#
# Dependencies: R (tested with 4.4.3) with the recommended packages survival
# (tested with 3.8.3) and splines (part of base R).
# =============================================================================
suppressPackageStartupMessages(library(survival))
library(splines)
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grep('^--file=',args)])
here <- if (length(self)==1) dirname(normalizePath(self)) else getwd()
targs <- commandArgs(trailingOnly=TRUE)
input <- if (length(targs)>=1) targs[1] else file.path(here,'..','data','pooled_1139_deid.csv')
out   <- if (length(targs)>=2) targs[2] else file.path(here,'..','outputs')
if (!file.exists(input)) stop('input table not found: ',input,
  '\nThe patient-level table is not distributed with this repository (the required columns are listed in the header of this script).',
  '\nDe-identified patient-level UMICH data are available from the corresponding author under a data use',
  ' agreement (see the manuscript\'s Data availability statement).')
dir.create(out,recursive=TRUE,showWarnings=FALSE)

name <- 'pooled_1139'
d <- read.csv(input,stringsAsFactors=FALSE)
# Input checks
stopifnot(all(c('cluster_id','cohort','gfx10','gfx10_z','os_time','os_event','age','hpv','admin_censored') %in% names(d)),
          nrow(d)==1139, setequal(unique(d$cohort),c('UMICH1','TCGA','UMICH2')),
          !anyNA(d[,c('cohort','gfx10','gfx10_z','os_time','os_event','admin_censored')]),
          all(d$os_event %in% 0:1), all(d$os_time>0),
          is.logical(d$admin_censored), all(d$os_event[d$admin_censored]==0))

# Input check: gfx10_z must equal gfx10 standardized within cohort (cohort mean
# and sample SD, denominator n - 1) to < 1e-12. The models below use the
# supplied gfx10_z.
z_recomputed <- ave(d$gfx10,d$cohort,FUN=function(x) (x-mean(x))/sd(x))
z_maxdiff <- max(abs(z_recomputed-d$gfx10_z))
stopifnot(z_maxdiff<1e-12)
cat(sprintf('gfx10_z = within-cohort standardized gfx10 (sample SD): max |recomputed - supplied| = %.2g\n',z_maxdiff))

models <- list(); sweeps <- list(); diagnostics <- list(); counts <- list(); kmrows <- list()
extract <- function(fit, term, dataset, model, clustered=FALSE) {
 s <- summary(fit); ci <- s$conf.int[term,]; co <- s$coefficients[term,]
 data.frame(dataset=dataset,model=model,n=fit$n,events=fit$nevent,
  HR=ci['exp(coef)'],lower=ci['lower .95'],upper=ci['upper .95'],p=co['Pr(>|z|)'],
  coefficient=co['coef'],se=co[if(clustered)'robust se' else 'se(coef)'],row.names=NULL)
}
# Group sizes and events for a 0/1/NA grouping (NA = intervening scores
# excluded from a tertile or quartile contrast). Does not enter any model.
groups <- function(high=NULL, event=NULL, lo=NA_real_, hi=NA_real_) {
 if(is.null(high)) return(data.frame(threshold_low=NA_real_,threshold=NA_real_,n_low=NA_integer_,n_mid=NA_integer_,
  n_high=NA_integer_,events_low=NA_integer_,events_mid=NA_integer_,events_high=NA_integer_))
 data.frame(threshold_low=as.numeric(lo),threshold=as.numeric(hi),n_low=sum(high %in% 0),n_mid=sum(is.na(high)),
  n_high=sum(high %in% 1),events_low=sum(event[high %in% 0]),events_mid=sum(event[is.na(high)]),
  events_high=sum(event[high %in% 1]))
}

# ---- models ----
 f <- coxph(Surv(os_time,os_event)~gfx10_z+strata(cohort),data=d,ties='efron',x=TRUE)
 fa <- coxph(Surv(os_time,os_event)~gfx10_z+age+hpv+strata(cohort),data=d,ties='efron',x=TRUE)
 models[[length(models)+1]] <- cbind(extract(f,'gfx10_z',name,'continuous_unadjusted'),groups())
 models[[length(models)+1]] <- cbind(extract(fa,'gfx10_z',name,'continuous_age_hpv'),groups())
 for(coh in unique(d$cohort)) {
  dc <- d[d$cohort==coh,]
  ff <- coxph(Surv(os_time,os_event)~gfx10_z,data=dc,ties='efron')
  models[[length(models)+1]] <- cbind(extract(ff,'gfx10_z',name,paste0(coh,'_continuous')),groups())
  counts[[length(counts)+1]] <- data.frame(dataset=name,cohort=coh,n=nrow(dc),patients=length(unique(dc$cluster_id)),events=sum(dc$os_event))
  q <- quantile(dc$gfx10,c(1/3,2/3),type=7)
  dc$high <- ifelse(dc$gfx10<=q[1],0,ifelse(dc$gfx10>q[2],1,NA))
  ff <- coxph(Surv(os_time,os_event)~high,data=dc,ties='efron')
  models[[length(models)+1]] <- cbind(extract(ff,'high',name,paste0(coh,'_tertile_high_vs_low')),groups(dc$high,dc$os_event))
 }
 for (scheme in c('median','tertile','quartile')) {
  q <- quantile(d$gfx10_z,switch(scheme,median=.5,tertile=c(1/3,2/3),quartile=c(.25,.75)),type=7)
  for(rule in c('gt','ge')) {
   upper <- if(rule=='gt') d$gfx10_z>tail(q,1) else d$gfx10_z>=tail(q,1)
   d$high <- if(scheme=='median') as.integer(upper) else ifelse(d$gfx10_z<=q[1],0,ifelse(upper,1,NA))
   ff <- coxph(Surv(os_time,os_event)~high+strata(cohort),data=d,ties='efron')
   models[[length(models)+1]] <- cbind(extract(ff,'high',name,paste0(scheme,'_high_vs_low_',rule)),groups(d$high,d$os_event,q[1],tail(q,1)))
  }
 }
 for(pct in seq(20,80,5)) {
  threshold <- quantile(d$gfx10_z,pct/100,type=7)
  for(rule in c('strict_greater','greater_or_equal')) {
   d$high <- if(rule=='strict_greater') as.integer(d$gfx10_z>threshold) else as.integer(d$gfx10_z>=threshold)
   ff <- coxph(Surv(os_time,os_event)~high+strata(cohort),data=d,ties='efron')
   r <- extract(ff,'high',name,'cutpoint_unadjusted')
   r$percentile<-pct;r$threshold<-as.numeric(threshold);r$rule<-rule;r$n_high<-sum(d$high)
   r$n_low<-sum(d$high==0);r$events_low<-sum(d$os_event[d$high==0]);r$events_high<-sum(d$os_event[d$high==1])
   sweeps[[length(sweeps)+1]] <- r
  }
 }
 knots <- quantile(d$gfx10_z,c(.05,.35,.65,.95),type=7)
 fs <- coxph(Surv(os_time,os_event)~ns(gfx10_z,knots=knots[2:3],Boundary.knots=knots[c(1,4)])+strata(cohort),data=d,ties='efron',x=TRUE)
 lrt <- 2*(fs$loglik[2]-f$loglik[2]); df <- length(coef(fs))-length(coef(f))
 ph <- cox.zph(f,transform='rank');pha<-cox.zph(fa,transform='rank')
 diagnostics[[length(diagnostics)+1]] <- data.frame(dataset=name,spline_lrt=lrt,spline_df=df,spline_p=pchisq(lrt,df,lower.tail=FALSE),ph_gfx_rank_p=ph$table['gfx10_z','p'],ph_adjusted_gfx_rank_p=pha$table['gfx10_z','p'],ph_adjusted_global_rank_p=pha$table['GLOBAL','p'])
 q <- quantile(d$gfx10_z,c(1/3,2/3),type=7)
 d$tertile <- factor(ifelse(d$gfx10_z<=q[1],'Low',ifelse(d$gfx10_z>q[2],'High','Mid')),levels=c('Low','Mid','High'))
 sf <- survfit(Surv(os_time,os_event)~tertile,data=d)
 for(gr in levels(d$tertile)) {
  dg<-d[d$tertile==gr,];sg<-survfit(Surv(os_time,os_event)~1,data=dg)
  ss<-summary(sg,times=c(1095,1460),extend=FALSE)
  kmrows[[length(kmrows)+1]]<-data.frame(dataset=name,group=gr,n=nrow(dg),events=sum(dg$os_event),survival_1095=ss$surv[1],survival_1460=ss$surv[2])
 }

# Records with an administratively assigned censoring time (admin_censored = TRUE;
# in UMICH2, records with no last-contact date, censored at 1,460 days or, when
# surgery was less than 1,460 days before the database lock, at the lock date).
# Counted per cohort; this enters no model.
adm <- do.call(rbind,lapply(unique(d$cohort),function(coh) {
 a <- d[d$cohort==coh & d$admin_censored,]
 sel <- list(all=rep(TRUE,nrow(a)),censored_at_1460_days=a$os_time==1460,censored_before_1460_days=a$os_time<1460)
 data.frame(dataset=name,cohort=coh,group=names(sel),n=sapply(sel,sum),
  patients=sapply(sel,function(s) length(unique(a$cluster_id[s]))),events=sapply(sel,function(s) sum(a$os_event[s])),row.names=NULL)
}))
adm$block <- 'admin_censoring'; adm$model <- 'administratively_assigned_censoring_times'

# Assemble one long, aggregate-only results table.
bind_union <- function(l) { cols <- unique(unlist(lapply(l,names)))
 do.call(rbind,lapply(l,function(x) { for(cn in setdiff(cols,names(x))) x[[cn]] <- NA; x[cols] })) }
mod <- do.call(rbind,models); mod$block <- 'model'
sw <- do.call(rbind,sweeps); sw$block <- 'cutpoint_sweep'
dg0 <- do.call(rbind,diagnostics)
spl <- data.frame(block='diagnostic',dataset=name,model='spline_vs_linear_LRT',term='ns(gfx10_z)',
 statistic=dg0$spline_lrt,df=dg0$spline_df,p=dg0$spline_p)
phrows <- rbind(
 data.frame(block='diagnostic',dataset=name,model='ph_rank_continuous_unadjusted',term=rownames(ph$table),
  statistic=ph$table[,'chisq'],df=ph$table[,'df'],p=ph$table[,'p'],row.names=NULL),
 data.frame(block='diagnostic',dataset=name,model='ph_rank_continuous_age_hpv',term=rownames(pha$table),
  statistic=pha$table[,'chisq'],df=pha$table[,'df'],p=pha$table[,'p'],row.names=NULL))
stopifnot(isTRUE(all.equal(phrows$p[phrows$model=='ph_rank_continuous_unadjusted' & phrows$term=='gfx10_z'],dg0$ph_gfx_rank_p)),
          isTRUE(all.equal(phrows$p[phrows$model=='ph_rank_continuous_age_hpv' & phrows$term=='gfx10_z'],dg0$ph_adjusted_gfx_rank_p)),
          isTRUE(all.equal(phrows$p[phrows$model=='ph_rank_continuous_age_hpv' & phrows$term=='GLOBAL'],dg0$ph_adjusted_global_rank_p)))
cnt <- do.call(rbind,counts); cnt$block <- 'cohort_counts'; cnt$model <- 'records_patients_events'
km <- do.call(rbind,kmrows); km$block <- 'km_pooled_tertile'; km$model <- 'kaplan_meier'
res <- bind_union(list(mod,sw,spl,phrows,cnt,km,adm))
res$reported_as <- ''
set_rep <- function(i,txt) { res$reported_as[i] <<- txt }
set_rep(res$block=='model' & res$model=='continuous_unadjusted','Results/Methods: primary continuous HR per within-cohort SD (1,139 records, 466 deaths)')
set_rep(res$block=='model' & res$model=='continuous_age_hpv','Results/Statistics: age- and HPV-adjusted HR (867 records, 335 deaths)')
set_rep(res$block=='model' & res$model=='median_high_vs_low_gt','Results: median split HR')
set_rep(res$block=='model' & res$model=='tertile_high_vs_low_gt','Results; Supplementary Fig. 3a: tertile high vs low HR and P')
set_rep(res$block=='model' & res$model=='quartile_high_vs_low_gt','Results: quartile high vs low HR')
set_rep(res$block=='cutpoint_sweep' & res$rule=='strict_greater','Results; Supplementary Fig. 3b: 13 thresholds, high strictly above threshold')
set_rep(res$block=='diagnostic' & res$model=='spline_vs_linear_LRT','Results; Supplementary Fig. 3c: spline non-linearity P (2-df LRT)')
set_rep(res$block=='diagnostic' & res$model=='ph_rank_continuous_unadjusted' & res$term=='gfx10_z','Methods: rank-transformed Schoenfeld P, unadjusted model')
set_rep(res$block=='diagnostic' & res$model=='ph_rank_continuous_age_hpv' & res$term=='gfx10_z','Methods: rank-transformed Schoenfeld P, adjusted model, Gfx10')
set_rep(res$block=='diagnostic' & res$model=='ph_rank_continuous_age_hpv' & res$term=='GLOBAL','Methods: rank-transformed Schoenfeld P, adjusted model, global')
set_rep(res$block=='cohort_counts','Methods: records per cohort (447/397/295)')
set_rep(res$block=='km_pooled_tertile','Supplementary Fig. 3a: pooled tertile group sizes (low/middle/high)')
set_rep(res$block=='admin_censoring' & res$cohort=='UMICH2' & res$group=='all','Methods: 21 administratively assigned UMICH2 censoring times (no last-contact date)')
lead <- c('block','dataset','model','cohort','group','term','rule','percentile','threshold_low','threshold','n','events','patients',
 'n_low','n_mid','n_high','events_low','events_mid','events_high','HR','lower','upper','p','coefficient','se','statistic','df',
 'survival_1095','survival_1460','reported_as')
stopifnot(setequal(lead,names(res)))
res <- res[,lead]; rownames(res) <- NULL
outfile <- file.path(out,'pooled_1139_results.csv')
write.csv(res,outfile,row.names=FALSE)

key <- res$reported_as!='' & res$block %in% c('model','diagnostic')
print(res[key,c('model','term','n','events','n_low','n_mid','n_high','HR','lower','upper','p','statistic','df')],row.names=FALSE,digits=10)
sw_s <- sw[sw$rule=='strict_greater',]
cat(sprintf('Sweep (strict, %d thresholds): HR %.10f-%.10f, max P %.10g\n',nrow(sw_s),min(sw_s$HR),max(sw_s$HR),max(sw_s$p)))
print(km,row.names=FALSE,digits=10)
cat('Records with an administratively assigned censoring time (admin_censored):\n')
print(adm[,c('cohort','group','n','patients','events')],row.names=FALSE)
cat('Wrote',normalizePath(outfile),'\n')
print(sessionInfo())
