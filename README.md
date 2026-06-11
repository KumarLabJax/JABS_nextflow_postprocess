# JABS Nextflow Postprocess

Post-processing pipeline for raw single-mouse outputs from the JABS behavioral tracking pipeline. Takes Nextflow outputs, runs QC screening and data validation, and produces a clean merged dataset ready for downstream analysis.

Both R and Python implementations are maintained in parallel. R scripts are designed for interactive use in RStudio; Python scripts are CLI-first for automation.

---

## Repository Structure

```
JABS_nextflow_postprocess/
├── r/                        R pipeline scripts
│   ├── 2_qc_check.R          compile QC logs, threshold check, missing/dup detection (interactive)
│   ├── 2_qc_check_cli.R      same as above, CLI version for automation
│   ├── 4_combine_batches.R   merge outputs from all batches into unified dataset
│   ├── 5a_outliers.R         outlier detection and QC figures
│   ├── 5b_heatmap.R          phenotype correlation heatmaps
│   └── render_pose.R         utility for rendering pose overlay videos
├── python/                   Python pipeline scripts
│   ├── 2_qc_check.py         compile QC logs, threshold check, missing/dup detection (CLI)
│   └── pose_corner_correction.py  copy pose_v6 files and embed manual corner corrections
├── notebooks/                exploratory analysis (not part of the main pipeline)
│   └── explore_features.py   PCA, clustering, and correlation exploration
├── config/
│   └── QC_params.yaml        default QC thresholds (override via --param flag)
├── pyproject.toml            Python dependencies (managed with uv)
└── renv.lock                 R dependencies (managed with renv)
```

---

## Workflow Overview

The pipeline follows six steps. Steps marked **planned** do not have an implementation yet.

| Step | Description | R | Python |
|------|-------------|---|--------|
| 1 | Visual inspection of videos | — | planned |
| 2+3 | Compile QC logs, check thresholds, flag missing/duplicate entries | `r/2_qc_check.R` | `python/2_qc_check.py` |
| 4 | Combine all outputs from different batches into unified format | `r/4_combine_batches.R` | planned |
| 5 | Check for outliers and plot correlation heatmap | `r/5a_outliers.R`, `r/5b_heatmap.R` | — |
| 6 | View snippets of predicted behaviors | — | planned |

Steps 2 and 3 are combined in a single script since they share the same inputs.

---

## Setup

### R
```r
renv::restore()
```

### Python
```bash
uv sync
```

---

## Step 2+3 — QC Check

**Scripts:** `r/2_qc_check.R` (interactive) · `r/2_qc_check_cli.R` (CLI) · `python/2_qc_check.py` (CLI)

Reads all Nextflow outputs, validates data completeness, screens for duplicates, post-processes gait data, and generates QC figures for fecal boli.

### Inputs

The `--input_dir` should contain one or more batch subdirectories with these files (searched recursively):

| Pattern in filename | Content |
|---------------------|---------|
| `qc_batch_` | QC report CSVs — defines the expected set of videos |
| `gait` | Gait feature CSVs |
| `fecal_boli` | Raw fecal boli CSVs |
| `feature` | JABS behavior feature CSVs |
| `morpho` | Morphometrics CSVs |

### Assumptions

1. Missing gait speed bins mean no gaits were predicted for those bins — not that gait prediction was skipped entirely.
2. All videos to be analyzed appear in at least one `qc_batch_*.csv` file.

### Running (Python CLI)

```bash
python python/2_qc_check.py \
    --input_dir /path/to/NextflowOutput \
    --output_dir /path/to/project_output \
    --param config/QC_params.yaml     # optional, overrides defaults
```

### Running (R CLI)

```bash
Rscript r/2_qc_check_cli.R \
    --input_dir /path/to/NextflowOutput \
    --output_dir /path/to/project_output \
    --param config/QC_params.yaml
```

### QC Parameters (`config/QC_params.yaml`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `expected_length` | 108150 s | Expected video duration |
| `max_tracklet_per_hour` | 6 | Max pose tracklets per hour |
| `max_missing_pose` | 0.5% | Max fraction of missing pose frames |
| `max_missing_segmentation` | 20% | Max fraction of missing segmentation frames |
| `max_missing_keypoint` | 1% | Max fraction of missing keypoint frames |
| `fecal_boli_quantile_plotting` | 0.05 | Lower quantile threshold for fecal boli outlier plots |

### Outputs

```
output_dir/
├── final_nextflow_feature_data/
│   ├── gait_final.csv            wide format; stride counts padded with 0 for missing bins;
│   │                             variances set to NA for bins with <3 strides
│   ├── morphometrics_final.csv
│   ├── JABS_features_final.csv
│   └── fecal_boli_raw.csv        ← requires manual review before proceeding
└── qc/
    ├── nextflow_qc_logs/
    │   ├── qc_all.csv            all videos with pass/fail flags
    │   └── qc_failed.csv         failed videos only
    ├── missing_or_dup_data/
    │   ├── missing_data.csv      in QC report but missing from feature files
    │   ├── videos_not_in_qc_report.csv  in feature files but missing from QC report
    │   └── duplicated_data.xlsx  identical or highly correlated rows
    └── qc_figs/
        └── fecal_boli_qc_figs.pdf
```

### Manual steps before proceeding to Step 4

1. Review `qc/nextflow_qc_logs/qc_failed.csv`. Manually overlay pose on flagged videos, decide which to keep or exclude.
2. Create `videos_to_exclude.txt` in your output directory — one `NetworkFilename` per line.
3. Review `qc/qc_figs/fecal_boli_qc_figs.pdf`. Correct any erroneous fecal boli counts in `final_nextflow_feature_data/fecal_boli_raw.csv`, then save it as `fecal_boli_final.csv` in the same directory.

---

## Step 4 — Combine Batches

**Script:** `r/4_combine_batches.R`

Merges the four final feature files with a `metadata.csv`, ensures all videos are represented across all files, and outputs a single merged dataset.

### Inputs

Expects the outputs of Step 2+3 in `final_nextflow_feature_data/`:
- `gait_final.csv`
- `morphometrics_final.csv`
- `JABS_features_final.csv`
- `fecal_boli_final.csv` (manually corrected from `fecal_boli_raw.csv`)

Also requires:
- `metadata.csv` — must have a `MouseID` column whose values appear as substrings in `NetworkFilename`
- `videos_to_exclude.txt` — created during manual QC in Step 2+3

### Outputs

```
output_dir/
├── merged_nextflow_dataset.csv
└── qc/missing_or_dup_data/
    ├── NetworkFilenames_missing_in_data.csv
    └── mice_missing_in_metadata.csv
```

### Manual steps before proceeding to Step 5

Review the generated QC tables. Mismatches between metadata and data do not stop the merge, but should be resolved or documented before final analysis.

---

## Step 5 — Outliers and Heatmap

**Scripts:** `r/5a_outliers.R` · `r/5b_heatmap.R`

Prepares the final dataset for analysis. `5a_outliers.R` is a template script intended to be customized per project.

### Inputs

- `merged_nextflow_dataset.csv` from Step 4

### Configuration (edit the `Set QC and other Values` section)

- Project name prefix for output files
- Feature substrings to manually exclude (`features.removed.manually`)
- Z-score threshold for outlier detection
- Number of top outlier mice/phenotypes to plot

### Outputs (`5a_outliers.R`)

```
output_dir/
├── YOUR_PROJECT_final_nextflow_dataset.csv
├── features_removed_from_curated_dataset.csv
└── qc/qc_figs/
    ├── scatter_plot_lm_figs.pdf
    └── zscore_boxplots/
        ├── gait_boxplot_outlier_figs.pdf
        ├── JABS_boxplot_outlier_figs.pdf
        ├── morpho_boxplot_outlier_figs.pdf
        ├── highest_z_boxplot_outlier_figs.pdf
        ├── most_freq_outliers_boxplot_outlier_figs.pdf
        └── outlier_videos_summary.csv
```

---

## Pose Corner Correction (utility)

**Script:** `python/pose_corner_correction.py`

Copies `pose_est_v6.h5` files from `failed_corners/` subdirectories and embeds manually corrected corner coordinates from SLEAP annotation files.

```bash
python python/pose_corner_correction.py \
    --input_dir /path/to/NextflowOutput \
    --output_dir /path/to/pose_v6_dir
```

Expects each batch directory inside `--input_dir` to optionally contain:
- `manual_corner_correction.slp` — SLEAP file with corrected corner annotations
- `failed_corners/` — directory containing `*_pose_est_v6.h5` files to correct

---

## Expected Output File Structure

```
/project_output_dir/
├── videos_to_exclude.txt                        (manually created)
├── metadata.csv                                 (manually provided)
├── YOUR_PROJECT_final_nextflow_dataset.csv      (Step 5 output → downstream analysis)
├── features_removed_from_curated_dataset.csv    (Step 5 output)
├── merged_nextflow_dataset.csv                  (Step 4 output)
├── final_nextflow_feature_data/                 (Step 2+3 outputs)
│   ├── gait_final.csv
│   ├── morphometrics_final.csv
│   ├── JABS_features_final.csv
│   ├── fecal_boli_raw.csv
│   └── fecal_boli_final.csv                     (manually corrected from fecal_boli_raw.csv)
└── qc/
    ├── nextflow_qc_logs/
    │   ├── qc_all.csv
    │   └── qc_failed.csv
    ├── missing_or_dup_data/
    │   ├── missing_data.csv
    │   ├── videos_not_in_qc_report.csv
    │   ├── duplicated_data.xlsx
    │   ├── NetworkFilenames_missing_in_data.csv
    │   └── mice_missing_in_metadata.csv
    └── qc_figs/
        ├── fecal_boli_qc_figs.pdf
        ├── scatter_plot_lm_figs.pdf
        └── zscore_boxplots/
            ├── gait_boxplot_outlier_figs.pdf
            ├── JABS_boxplot_outlier_figs.pdf
            ├── morpho_boxplot_outlier_figs.pdf
            ├── highest_z_boxplot_outlier_figs.pdf
            ├── most_freq_outliers_boxplot_outlier_figs.pdf
            └── outlier_videos_summary.csv
```
