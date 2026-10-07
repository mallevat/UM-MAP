# UM-MAP

Code for cell classification of H&E whole-slide images in QuPath, tumor-lymphoid G-cross scores (Gfx) and the survival analyses of the study "Convergent Spatial Analyses Define a Prognostic Tumor-Immune Interface Niche in Head and Neck Cancer".

## qupath/

`classify_export.groovy` sets QuPath's H&E stain vectors, defines the region to analyze (the whole image, or tissue regions from a pixel classifier), runs watershed cell detection on the hematoxylin channel, classifies every cell with an object classifier and writes one tab-separated table per slide with the columns `Class`, `Centroid X µm` and `Centroid Y µm`. `mkproj.groovy` creates a single-image QuPath project, and `run_slide.sh` runs both scripts on one slide and removes the project afterwards. The defaults are `Combo-12-08-24.json`, the whole image and detection threshold 0.05; the demo below uses `Tumor-Fibroblast-Lymphoid-Myeloid-2.json`, tissue regions from `TissueROI.json` and threshold 0.1.

```bash
export QUPATH=/path/to/QuPath/bin/QuPath
qupath/run_slide.sh <slide.svs> <output_dir>                    # defaults
qupath/run_slide.sh <slide.svs> <output_dir> Tumor-Fibroblast-Lymphoid-Myeloid-2.json TissueROI.json 0.1   # demo settings
```

Arguments are `<slide> <output_dir> [classifier] [tissue_mode] [threshold]`, where `tissue_mode` is `full` or a pixel classifier file name. Classifier files are read from `classifiers/` (or from `$UMMAP_CLASSIFIER_DIR`). The output is `<output_dir>/<slide file name>.tsv`, with every character other than letters and digits replaced by `_` (for example `slide.svs` gives `slide_svs.tsv`). QuPath creates the project in a folder next to the slide, so the slide's folder must be writable. One slide takes about 1 to 6 minutes with 8 CPU cores and 64 GB of memory.

## classifiers/

`Combo-12-08-24.json` (artificial neural network) and `Tumor-Fibroblast-Lymphoid-Myeloid-2.json` (random trees) are QuPath object classifiers. `TissueROI.json` is the QuPath pixel classifier used to detect tissue regions.

## gcross/

`gcross_auc.R` reads a folder of per-cell tables and, for each slide, computes the G-cross function with Kaplan-Meier edge correction (spatstat `Gcross`) within the convex hull of all cells for six cell-type pairs (fibroblast to tumor, fibroblast to lymphoid, tumor to lymphoid, and each class with itself), and the area under each curve from 0 to each radius (default 10, 20 and 40 µm). `Gfx10_Tumor_Lymphoid` is the tumor-to-lymphoid area up to 10 µm. The output has one row per slide with cell counts per class and all areas. Run `Rscript gcross/gcross_auc.R --help` for all options.

```bash
Rscript gcross/gcross_auc.R --input <folder_of_tsv_files> --output gcross_auc.csv
```

## survival/

R scripts (and one Python script) for the survival analyses. Pooled analyses use Gfx10 standardized within each cohort (z-scores). The input tables are not included. The TCGA table for 02, 04 and 05 (public TCGA-HNSC survival and genomic data joined to the study's TCGA Gfx values) is available from the corresponding author; the tables with UMICH patients (the pooled table for 01 and 04, and the UMICH tables for 03 and 07) are available from the corresponding author under a data use agreement and University of Michigan Institutional Review Board approval. The inputs of 06 (not included) are built from public TCGA-HNSC data and the niche-signature gene lists. Column names are listed in each script's header. Scripts 01, 03 and 04 check that their input is the study table (row and cohort counts).

| Script | Analysis |
| --- | --- |
| `01_pooled_primary_1139.R` | Cohort-stratified Cox models of Gfx10 in the pooled dataset: continuous, age- and HPV-adjusted, median, tertile and quartile contrasts, cut-point sweep, restricted cubic spline, proportional-hazards tests |
| `02_tcga_comparators.R` | TCGA survival by Gfx and cell-composition metrics, tertile and top-third groupings (Cox and log-rank) |
| `03_recurrence_umich1.R` | UMICH1 recurrence-or-persistent-disease Cox models |
| `04_figures_fig2b_suppfig3.R` | Kaplan-Meier curves by Gfx10 tertile per cohort and pooled, cut-point and spline plots |
| `05_tcga_genomic_suppfig6.R` | TCGA survival by tertile of tumor mutation burden, neoantigen load and CD274 expression |
| `06_tcga_niche_signature_survival.py` | TCGA survival by niche-signature score (Cox per SD, median-split log-rank) |
| `07_suppfig2_umich_comparators.R` | UMICH1 and UMICH2 survival by Gfx at 10, 20 and 40 µm, tertile and top-third groupings |

```bash
Rscript survival/01_pooled_primary_1139.R <pooled_table.csv> <output_dir>
Rscript survival/02_tcga_comparators.R <tcga_table.csv> <output_dir>
Rscript survival/03_recurrence_umich1.R <umich1_table.csv> <output_dir>
Rscript survival/04_figures_fig2b_suppfig3.R --pooled=<pooled_table.csv> --tcga=<tcga_table.csv> --outdir=<output_dir>
Rscript survival/05_tcga_genomic_suppfig6.R <tcga_table.csv> <output_dir>
python survival/06_tcga_niche_signature_survival.py --scores <scores.csv> --genes <genes.csv> --out <output_dir>
Rscript survival/07_suppfig2_umich_comparators.R --umich1 <umich1_table.csv> --umich2 <umich2_table.csv> --out <output_dir>
```

Script 06 can also rebuild the signature scores from the UCSC Xena TCGA-HNSC expression matrix (`--expression <HiSeqV2.gz> --clinical <clinical.csv> --programs <programs.csv>`). Each script runs in under a few minutes on a desktop computer.

## demo/

`demo/example_input/` holds the per-cell table of the public TCGA-HNSC slide TCGA-HD-A6HZ-01Z-00-DX1 (GDC open access), made with the demo settings above (`Tumor-Fibroblast-Lymphoid-Myeloid-2.json`, `TissueROI.json`, threshold 0.1), and `demo/expected_output/gcross_auc.csv` is the G-cross output for it. The run takes a few seconds:

```bash
Rscript gcross/gcross_auc.R --input demo/example_input --output gcross_auc.csv
```

## Software versions

- QuPath 0.4.3 (qupath/), from <https://github.com/qupath/qupath/releases/tag/v0.4.3>
- G-cross: R 4.4.0, spatstat 3.0-8 (spatstat.geom 3.2-9, spatstat.explore 3.2-7, spatstat.random 3.2-3), pracma 2.4.6; openxlsx 4.2.9 only for `--xlsx`
- Survival (R): R 4.4.3, survival 3.8-3, splines; 02 also needs multcomp 1.4-32, and survminer 0.5.2, ggplot2 4.0.3 and ggpubr 1.0.0 for its plots; 05 needs survminer and ggplot2 for its plots
- Survival (Python, 06): Python 3.11.7, numpy 1.26.4, pandas 2.1.4, lifelines 0.30.0, matplotlib 3.8.0

R packages install from CRAN with `install.packages()` in a few minutes.

## License

MIT (see `LICENSE`).
