#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
#
# 06_tcga_niche_signature_survival.py
#
# TCGA-HNSC overall survival by spatial niche signature score (tumour, lymphoid and interface
# signatures):
#   - continuous Cox model, hazard ratio per SD of the signature score;
#   - median split (High = score > median, Low = score <= median) with a log-rank P;
#   - Kaplan-Meier curves for the median split (Supplementary Fig. 12e).
#
# Method:
#   - Expression: UCSC Xena TCGA.HNSC.sampleMap/HiSeqV2, log2(RSEM normalized count + 1).
#     A signature score is the unweighted mean of these values over the signature's genes.
#   - Genes: for each NicheCompass program selected by SIGNATURE_FILTERS (consensus niche and
#     effect size; for the tumour signature also sample agreement), its 10 genes with the
#     largest decoder weights (Top_10_Genes). Genes absent from the expression matrix are
#     dropped. The lymphoid signature is called "Immune" in the code.
#   - Samples: every RNA-seq sample of a patient with overall-survival data is analysed (primary
#     tumour -01, solid tissue normal -11, metastasis -06).
#   - Survival: GDC clinical data. OS.time = days_to_death, or days_to_last_follow_up if the
#     patient is alive; OS = 1 if vital_status is "Dead". Follow-up is censored at 1,500 days.
#   - Cox: lifelines CoxPHFitter (Efron ties). The score is standardized over the analysed rows
#     (sample SD), so hazard ratios are per SD; P is the Wald test. Log-rank test: lifelines.
#   - Two populations are analysed: all_samples (every scored sample) and primary_tumour_only
#     (sample type -01, one sample per patient).
#
# Usage:
#   python survival/06_tcga_niche_signature_survival.py [--out DIR] [--no-plots]
#       survival analyses from the score table (--scores)
#   python survival/06_tcga_niche_signature_survival.py --expression HiSeqV2.gz [--out DIR] [--no-plots]
#       builds the scores from the Xena matrix first, compares them with the --scores table if
#       that file exists, then runs the survival analyses on the new scores
#   Options (defaults are relative to the folder above this script's folder):
#     --scores      score table; default data/tcga_niche_signature_scores.csv
#     --genes       gene table; default data/tcga_niche_signature_genes.csv
#     --clinical    clinical table; default data/tcga_gdc_clinical_os_snapshot_2025-11-08.csv
#     --programs    program table; default data/niche_program_consensus_top10.csv
#     --expression  Xena expression matrix (runs the scoring step); no default
#     --out         output folder; default outputs/06_tcga_niche_signature_survival
#     --no-plots    skip the PDF figures
#
# Inputs:
#   scores    sample, patient_id, sample_type_code, sample_type, OS.time, OS, tumor_score,
#             lymphoid_score, interface_score (one row per RNA-seq sample)
#   genes     signature, program, consensus_niche, effect_size, n_samples_agreement,
#             rank_in_program, gene, in_expression_matrix (TRUE/FALSE)
#   programs  Program, Consensus_Niche, N_Samples_Agreement, Effect_Size, Top_10_Genes
#             (gene symbols separated by ", ")
#   clinical  submitter_id, OS.time, OS (or days_to_death, days_to_last_follow_up and
#             vital_status, from which OS.time and OS are derived)
#   expression  genes x samples, tab-separated, from
#     https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.HNSC.sampleMap%2FHiSeqV2.gz
#     (version used: 30,432,155 bytes,
#     sha256 ed2b083b327cc1ccb1274b017a9eb1b40ac28fc63b7942b31448e4a6885a596a)
#
# Outputs (in --out):
#   cox_logrank_results.csv     per population and signature: rows, patients, events, HR per SD
#                               (95% CI), Wald P, median cut, group sizes, deaths, log-rank P
#   km_step_functions.csv       Kaplan-Meier step functions of the median-split groups
#   signature_gene_summary.csv  programs and gene counts per signature (when a gene table is
#                               available)
#   km_median_split_all_samples.pdf, km_median_split_primary_tumour_only.pdf  (unless --no-plots)
#   scores_from_expression.csv, signature_genes_from_expression.csv  (scoring step only)
#   run_manifest.txt            arguments, package versions, input sha256, output files
#
# Dependencies: Python 3 with numpy, pandas and lifelines; matplotlib for the figures.
#   Tested with Python 3.11.7, numpy 1.26.4, pandas 2.1.4, lifelines 0.30.0 and matplotlib 3.8.0.

import argparse
import hashlib
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import lifelines

STAGE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = STAGE_DIR / "data"

CENSOR_DAYS = 1500  # follow-up censored at 1,500 days

# Program filters. "Immune" is the lymphoid signature.
SIGNATURE_FILTERS = {
    'Tumor': lambda df: (df['Consensus_Niche'] == 'Malignant') &
                        (df['Effect_Size'] > 0.8) &
                        (df['N_Samples_Agreement'] >= 3),
    'Immune': lambda df: (df['Consensus_Niche'] == 'Lymphoid') &
                         (df['Effect_Size'] > 0),
    'Interface': lambda df: (df['Consensus_Niche'] == 'Combo') &
                            (df['Effect_Size'] > 0.5),
}
# Signature names used in the output tables
SIG_NAME = {'Tumor': 'tumor', 'Immune': 'lymphoid', 'Interface': 'interface'}
SCORE_COL = {k: f"{v}_score" for k, v in SIG_NAME.items()}
SAMPLE_TYPE = {'01': 'Primary Solid Tumor', '06': 'Metastatic', '11': 'Solid Tissue Normal'}


def rel(path):
    """Path relative to the folder above this script's folder, for the manifest (absolute if outside it)."""
    try:
        return str(Path(path).resolve().relative_to(STAGE_DIR))
    except ValueError:
        return str(Path(path).resolve())


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


# =============================================================================
# Scoring step (only with --expression)
# =============================================================================
def score_from_expression(expr_file, clinical_file, programs_file):
    # Xena matrix: genes x samples, log2(norm_count+1).
    expr = pd.read_csv(expr_file, sep='\t', index_col=0)
    if expr.shape[0] < expr.shape[1]:          # orientation guard
        expr = expr.T
    expr.index = [str(g)[:15] for g in expr.index]   # gene names cut to 15 characters
    if not expr.index.is_unique:
        sys.exit("gene names are not unique after the 15-character cut")

    clinical = pd.read_csv(clinical_file)
    if 'submitter_id' in clinical.columns:
        clinical['patient_id'] = clinical['submitter_id'].str[:12]
    if 'OS.time' not in clinical.columns:     # derive OS.time and OS from the GDC fields
        clinical['OS.time'] = clinical.apply(
            lambda x: x['days_to_death'] if pd.notna(x['days_to_death']) else x.get('days_to_last_follow_up'),
            axis=1)
        clinical['OS'] = (clinical['vital_status'] == 'Dead').astype(int)

    # Every RNA-seq sample (column) of a patient with clinical data is kept (primary tumour,
    # solid tissue normal, metastasis), with its full sample barcode.
    samples = pd.DataFrame({'sample': [str(c) for c in expr.columns]})
    samples['patient_id'] = samples['sample'].str[:12]
    overlap = set(samples['patient_id']) & set(clinical['patient_id'])
    keep = samples['patient_id'].isin(overlap).to_numpy()
    expr = expr.loc[:, keep]
    samples = samples.loc[keep].reset_index(drop=True)
    clinical = clinical[clinical['patient_id'].isin(overlap)].copy()
    clinical = clinical[clinical['OS.time'].notna()].copy()

    programs_full = pd.read_csv(programs_file)       # one row per program
    gene_rows, summary_rows = [], []
    for sig_name, filter_func in SIGNATURE_FILTERS.items():
        filtered = programs_full[filter_func(programs_full)]
        all_genes = []
        for _, prog in filtered.iterrows():
            genes10 = prog['Top_10_Genes'].split(', ')
            all_genes.extend(genes10)
            for rank, g in enumerate(genes10, start=1):            # gene table
                gene_rows.append({'signature': SIG_NAME[sig_name], 'program': prog['Program'],
                                  'consensus_niche': prog['Consensus_Niche'],
                                  'effect_size': prog['Effect_Size'],
                                  'n_samples_agreement': prog['N_Samples_Agreement'],
                                  'rank_in_program': rank, 'gene': g,
                                  'in_expression_matrix': g in expr.index})
        unique_genes = list(set(all_genes))
        available = sorted(g for g in unique_genes if g in expr.index)   # sorted: fixed summation order
        score = expr.loc[available].mean(axis=0)                          # unweighted mean over the genes
        samples[SCORE_COL[sig_name]] = score.to_numpy()
        summary_rows.append({'signature': SIG_NAME[sig_name], 'n_programs': len(filtered),
                             'programs': ';'.join(filtered['Program']),
                             'n_gene_listings': len(all_genes), 'n_unique_genes': len(unique_genes),
                             'n_genes_in_matrix_used': len(available),
                             'genes_not_in_matrix': ';'.join(sorted(set(unique_genes) - set(available)))})

    merged = clinical[['patient_id', 'OS.time', 'OS']].merge(samples, on='patient_id')
    merged['sample_type_code'] = merged['sample'].str[13:15]
    merged['sample_type'] = merged['sample_type_code'].map(SAMPLE_TYPE)
    merged = merged[['sample', 'patient_id', 'sample_type_code', 'sample_type', 'OS.time', 'OS',
                     'tumor_score', 'lymphoid_score', 'interface_score']]
    merged = merged.sort_values(['patient_id', 'sample']).reset_index(drop=True)
    return merged, pd.DataFrame(gene_rows), pd.DataFrame(summary_rows)


def write_gene_table(genes, path):
    """Gene table with TRUE/FALSE flags (readable by R and pandas)."""
    g = genes.copy()
    g['in_expression_matrix'] = np.where(g['in_expression_matrix'], 'TRUE', 'FALSE')
    g.to_csv(path, index=False)


def read_gene_table(path):
    g = pd.read_csv(path)
    g['in_expression_matrix'] = g['in_expression_matrix'].astype(str).str.upper().map({'TRUE': True, 'FALSE': False})
    return g


def summarise_genes(genes):
    """Per-signature program and gene counts (a gene listed twice is used once)."""
    rows = []
    for sig in ('tumor', 'lymphoid', 'interface'):
        d = genes[genes['signature'] == sig]
        inm = d['in_expression_matrix'].astype(bool)
        rows.append({'signature': sig, 'n_programs': d['program'].nunique(),
                     'programs': ';'.join(dict.fromkeys(d['program'])),
                     'n_gene_listings': len(d), 'n_unique_genes': d['gene'].nunique(),
                     'n_genes_in_matrix_used': d.loc[inm, 'gene'].nunique(),
                     'genes_not_in_matrix': ';'.join(sorted(set(d.loc[~inm, 'gene'])))})
    return pd.DataFrame(rows)


# =============================================================================
# Survival analyses (applied to one population at a time)
# =============================================================================
def add_censoring(tab):
    tab = tab.copy()
    tab['OS.time_1500'] = tab['OS.time'].apply(lambda x: min(x, CENSOR_DAYS))
    tab['OS_1500'] = tab.apply(lambda row: row['OS'] if row['OS.time'] <= CENSOR_DAYS else 0, axis=1)
    return tab


def analyse(tab, sig_name, population):
    col = SCORE_COL[sig_name]
    merged = tab[['patient_id', 'OS.time_1500', 'OS_1500', col]].rename(columns={col: 'Score'}).dropna()

    # Cox regression (continuous score, standardized)
    cox_df = merged[['OS.time_1500', 'OS_1500', 'Score']].copy()
    cox_df['Score'] = (cox_df['Score'] - cox_df['Score'].mean()) / cox_df['Score'].std()
    cph = CoxPHFitter()
    cph.fit(cox_df, duration_col='OS.time_1500', event_col='OS_1500')
    s = cph.summary.loc['Score']

    # Log-rank (median split)
    median = merged['Score'].median()
    high = merged[merged['Score'] > median]
    low = merged[merged['Score'] <= median]
    lr = logrank_test(high['OS.time_1500'], low['OS.time_1500'], high['OS_1500'], low['OS_1500'])

    return {'population': population, 'signature': SIG_NAME[sig_name], 'name_in_original_code': sig_name,
            'n_rows': len(merged), 'n_patients': merged['patient_id'].nunique(),
            'events_1500d': int(merged['OS_1500'].sum()),
            'hr_per_sd': s['exp(coef)'], 'hr_lower95': s['exp(coef) lower 95%'],
            'hr_upper95': s['exp(coef) upper 95%'], 'cox_wald_p': s['p'],
            'coef': s['coef'], 'se_coef': s['se(coef)'],
            'score_mean': merged['Score'].mean(), 'score_sd': merged['Score'].std(),
            'median_cut': median, 'n_high': len(high), 'n_low': len(low),
            'deaths_high': int(high['OS_1500'].sum()), 'deaths_low': int(low['OS_1500'].sum()),
            'pct_dead_high': 100 * high['OS_1500'].sum() / len(high),
            'pct_dead_low': 100 * low['OS_1500'].sum() / len(low),
            'logrank_chisq': lr.test_statistic, 'logrank_p': lr.p_value}


def km_tables(tab, population):
    """Kaplan-Meier step functions behind the figure (lifelines defaults)."""
    rows = []
    for sig_name in SIGNATURE_FILTERS:
        col = SCORE_COL[sig_name]
        median = tab[col].median()
        for grp, mask in (('High', tab[col] > median), ('Low', tab[col] <= median)):
            kmf = KaplanMeierFitter()
            kmf.fit(tab.loc[mask, 'OS.time_1500'], tab.loc[mask, 'OS_1500'])
            et = kmf.event_table
            ci = kmf.confidence_interval_survival_function_
            sf = kmf.survival_function_
            for t in sf.index:
                rows.append({'population': population, 'signature': SIG_NAME[sig_name], 'group': grp,
                             'time_days': t, 'at_risk': int(et.loc[t, 'at_risk']) if t in et.index else np.nan,
                             'deaths': int(et.loc[t, 'observed']) if t in et.index else np.nan,
                             'censored': int(et.loc[t, 'censored']) if t in et.index else np.nan,
                             'survival': sf.iloc[sf.index.get_loc(t), 0],
                             'ci_lower95': ci.iloc[ci.index.get_loc(t), 0],
                             'ci_upper95': ci.iloc[ci.index.get_loc(t), 1]})
    return pd.DataFrame(rows)


def plot_km(tab, results, out_pdf, title_suffix, p_three_decimals=False):
    """Three-panel Kaplan-Meier figure, median split (days on the x-axis)."""
    import logging
    import matplotlib
    matplotlib.use('Agg')
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)   # no Arial warnings
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.family'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']  # with fallback fonts
    mpl.rcParams['font.size'] = 11
    HIGH_COLOR, LOW_COLOR = '#D6487E', '#5B9BD5'
    titles = {'Tumor': 'Tumor Signature', 'Immune': 'Lymphoid Signature', 'Interface': 'Interface Signature'}
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, sig_name in zip(axes, SIGNATURE_FILTERS):
        col = SCORE_COL[sig_name]
        median = tab[col].median()
        high_mask = tab[col] > median
        low_mask = ~high_mask
        kmf = KaplanMeierFitter()
        kmf.fit(tab.loc[high_mask, 'OS.time_1500'], tab.loc[high_mask, 'OS_1500'],
                label=f'High (n={high_mask.sum()})')
        kmf.plot_survival_function(ax=ax, ci_show=True, color=HIGH_COLOR, linewidth=2.5, alpha=0.9)
        kmf = KaplanMeierFitter()
        kmf.fit(tab.loc[low_mask, 'OS.time_1500'], tab.loc[low_mask, 'OS_1500'],
                label=f'Low (n={low_mask.sum()})')
        kmf.plot_survival_function(ax=ax, ci_show=True, color=LOW_COLOR, linewidth=2.5, alpha=0.9)
        p = results.loc[results['signature'] == SIG_NAME[sig_name], 'logrank_p'].iloc[0]
        p_text = "p < 0.001" if p < 0.001 else (f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}")
        if p_three_decimals and p >= 0.001:   # e.g. 0.050 rather than 0.05
            p_text = f"p = {p:.3f}"
        ax.set_title(titles[sig_name], fontsize=14, fontweight='bold', pad=15)
        ax.set_xlabel('Time (days)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Overall Survival', fontsize=12, fontweight='bold')
        ax.set_xlim(0, 1500)
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
        ax.legend(loc='lower left', frameon=True, fontsize=10, fancybox=False, edgecolor='black', framealpha=1)
        ax.text(0.98, 0.98, p_text, transform=ax.transAxes, fontsize=11, fontweight='bold',
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='black', linewidth=1))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
    plt.tight_layout()
    fig.text(0.5, -0.03, title_suffix, ha='center', fontsize=10)   # population label
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight', format='pdf', transparent=False,
                metadata={'CreationDate': None})
    plt.close(fig)


# =============================================================================
# main
# =============================================================================
def main():
    ap = argparse.ArgumentParser(description='TCGA niche-signature survival (Supplementary Fig. 12e)')
    ap.add_argument('--scores', default=str(DATA_DIR / 'tcga_niche_signature_scores.csv'))
    ap.add_argument('--genes', default=str(DATA_DIR / 'tcga_niche_signature_genes.csv'))
    ap.add_argument('--expression', default=None, help='Xena TCGA.HNSC.sampleMap/HiSeqV2(.gz): run the scoring step')
    ap.add_argument('--clinical', default=str(DATA_DIR / 'tcga_gdc_clinical_os_snapshot_2025-11-08.csv'))
    ap.add_argument('--programs', default=str(DATA_DIR / 'niche_program_consensus_top10.csv'))
    ap.add_argument('--out', default=str(STAGE_DIR / 'outputs' / '06_tcga_niche_signature_survival'))
    ap.add_argument('--no-plots', action='store_true')
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    manifest = [f"script: {Path(__file__).name}", f"run: {datetime.now().isoformat(timespec='seconds')}",
                f"arguments: {' '.join(sys.argv[1:]) or '(defaults)'}", f"output directory: {rel(out)}",
                f"python {platform.python_version()}",
                f"numpy {np.__version__}, pandas {pd.__version__}, lifelines {lifelines.__version__}"]

    gene_summary = None
    if args.expression:
        tab, genes, loop_summary = score_from_expression(args.expression, args.clinical, args.programs)
        gene_summary = summarise_genes(genes)
        if not (gene_summary.set_index('signature')['n_genes_in_matrix_used']
                == loop_summary.set_index('signature')['n_genes_in_matrix_used']).all():
            sys.exit('gene counts from the scoring loop and the gene table disagree')
        tab.to_csv(out / 'scores_from_expression.csv', index=False)
        write_gene_table(genes, out / 'signature_genes_from_expression.csv')
        gene_summary.to_csv(out / 'signature_gene_summary.csv', index=False)
        for f in (args.expression, args.clinical, args.programs):
            manifest.append(f"input sha256 {sha256(f)}  {rel(f)}")
        shipped = Path(args.scores)
        if shipped.exists():   # compare the rebuilt scores with the --scores table
            ref = pd.read_csv(shipped, dtype={'sample_type_code': str}, float_precision='round_trip')
            same_rows = (len(ref) == len(tab) and
                         (ref[['sample', 'patient_id', 'sample_type_code']].to_numpy() ==
                          tab[['sample', 'patient_id', 'sample_type_code']].to_numpy()).all() and
                         np.array_equal(ref[['OS.time', 'OS']].to_numpy(float), tab[['OS.time', 'OS']].to_numpy(float)))
            diff = max(float(np.max(np.abs(ref[c].to_numpy() - tab[c].to_numpy())))
                       for c in ('tumor_score', 'lymphoid_score', 'interface_score')) if same_rows else float('nan')
            msg = f"rebuilt scores vs {shipped.name}: rows identical = {same_rows}; max abs score difference = {diff:.3g}"
            print(msg)
            manifest.append(msg)
    else:
        # float_precision='round_trip' reads the 17-digit scores back exactly
        tab = pd.read_csv(args.scores, dtype={'sample_type_code': str}, float_precision='round_trip')
        manifest.append(f"input sha256 {sha256(args.scores)}  {rel(args.scores)}")
        if Path(args.genes).exists():
            genes = read_gene_table(args.genes)
            manifest.append(f"input sha256 {sha256(args.genes)}  {rel(args.genes)}")
            gene_summary = summarise_genes(genes)
            gene_summary.to_csv(out / 'signature_gene_summary.csv', index=False)

    tab = add_censoring(tab)
    pops = {'all_samples': tab,
            'primary_tumour_only': tab[tab['sample_type_code'] == '01'].copy()}
    if not pops['primary_tumour_only']['patient_id'].is_unique:
        sys.exit('more than one primary-tumour sample for a patient')
    res = pd.DataFrame([analyse(t, s, p) for p, t in pops.items() for s in SIGNATURE_FILTERS])
    res.to_csv(out / 'cox_logrank_results.csv', index=False)
    km = pd.concat([km_tables(t, p) for p, t in pops.items()], ignore_index=True)
    km.to_csv(out / 'km_step_functions.csv', index=False)
    if not args.no_plots:
        a, b = pops['all_samples'], pops['primary_tumour_only']
        types = a['sample_type_code'].value_counts()
        plot_km(a, res[res['population'] == 'all_samples'], out / 'km_median_split_all_samples.pdf',
                f"All samples: {len(a)} RNA-seq samples ({types.get('01', 0)} primary tumour, "
                f"{types.get('11', 0)} solid-tissue normal, {types.get('06', 0)} metastasis) "
                f"from {a['patient_id'].nunique()} patients")
        plot_km(b, res[res['population'] == 'primary_tumour_only'], out / 'km_median_split_primary_tumour_only.pdf',
                f"Primary tumours only: {len(b)} patients", p_three_decimals=True)

    pd.set_option('display.width', 200)
    show = res[['population', 'signature', 'n_rows', 'n_patients', 'events_1500d', 'hr_per_sd', 'hr_lower95',
                'hr_upper95', 'cox_wald_p', 'n_high', 'n_low', 'logrank_p']]
    print(show.to_string(index=False))
    manifest += [f"output {p.name}" for p in sorted(out.iterdir()) if p.is_file() and p.name != 'run_manifest.txt']
    (out / 'run_manifest.txt').write_text('\n'.join(manifest) + '\n')


if __name__ == '__main__':
    main()
