# JABS Nextflow Postprocess

Post-processing pipeline for raw single-mouse outputs from the JABS behavioral tracking pipeline. Takes Nextflow outputs, runs QC screening and data validation, and produces a clean merged dataset ready for downstream analysis.

---

## Overview

| Script | What it does |
|--------|-------------|
| `notebooks/1_check_videos.ipynb` | Browse raw video frames on HPC without downloading. Flags likely-empty videos by motion score. |
| `r/2_qc_check.R` · `python/2_qc_check.py` | Read all Nextflow outputs across batches, check each video against QC thresholds (length, pose coverage, tracklet count), flag missing and duplicate entries, post-process gait data into wide format, and plot raw fecal boli curves for manual review. |
| `r/3_combine_batches.R` | Merge the four cleaned feature files (gait, morphometrics, JABS features, fecal boli) with a metadata table. Validates that all videos and mice are accounted for across every file before writing a single merged dataset. |
| `r/4_outliers.R` | Project-specific template for final dataset curation. Removes zero-variance and manually flagged features, runs z-score outlier detection, and generates diagnostic box plots and scatter plots to catch any remaining problematic mice or phenotypes. |
| `r/5_heatmap.R` | Plots pairwise phenotype correlation heatmaps across the final dataset. Under active development. |
| `python/6a_qc_classifiers.py` | For each detected behavior bout in the Nextflow output, samples short video clips at random. Optionally overlays the pose skeleton. Writes one MP4 per clip, organized by behavior label. |
| `python/6b_qc_viewer.py` | Streamlit app that steps through the clips produced by `6a`. Lets you accept, reject, or skip each clip and saves the verdicts to a CSV for downstream use. |
| `python/pose_corner_correction.py` | Utility. When automated arena-corner detection fails for a video, a human re-labels the corners in SLEAP. This script copies the affected pose H5 files and patches them with the corrected corner coordinates. |

---

## Repository Structure

```
JABS_nextflow_postprocess/
├── r/                            R pipeline scripts
│   ├── 2_qc_check.R              compile QC logs, threshold check, missing/dup (interactive)
│   ├── 2_qc_check_cli.R          same as above, CLI version for automation
│   ├── 3_combine_batches.R       merge feature files + metadata into unified dataset
│   ├── 4_outliers.R              outlier detection and QC figures (project-specific template)
│   ├── 5_heatmap.R               phenotype correlation heatmaps
│   ├── render_pose.R             utility — render pose overlay videos
│   └── utils.R                   shared utilities (sourced by other R scripts)
├── python/                       Python pipeline scripts
│   ├── 2_qc_check.py             compile QC logs, threshold check, missing/dup (CLI)
│   ├── 6a_qc_classifiers.py      extract video clips at detected behavior bouts
│   ├── 6b_qc_viewer.py           Streamlit app for reviewing and annotating clips
│   └── pose_corner_correction.py utility — embed manual corner corrections into pose H5 files
├── notebooks/                    interactive tools (run in JupyterLab)
│   ├── 1_check_videos.ipynb      paginated video frame previewer with motion screening
│   └── explore_features.py       PCA, clustering, and correlation exploration
├── src/
│   └── utils.py                  shared Python utilities (video helpers, pose overlay, motion screening)
├── config/
│   └── QC_params.yaml            default QC thresholds (override via --param flag)
├── pyproject.toml                Python dependencies (managed with uv)
└── renv.lock                     R dependencies (managed with renv)
```

Step numbers reflect the recommended order. Utilities (`pose_corner_correction.py`, `utils.R`, `src/utils.py`) are not numbered — they are called by other scripts or used on demand.

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

## Scripts

### 1 — Video Inspection · `notebooks/1_check_videos.ipynb`

Paginated video frame previewer for HPC environments — no GUI or local download required. Samples 5 evenly spaced frames per video and shows 10 videos per page with Next/Previous navigation. Also includes a motion-based screen to flag likely empty videos.

Powered by `src/utils.py`. Open in JupyterLab and run cells interactively.

---

### 2 — QC Check · `r/2_qc_check.R` · `r/2_qc_check_cli.R` · `python/2_qc_check.py`

Reads all Nextflow outputs across batches, validates data completeness, screens for duplicates, post-processes gait data, and generates QC figures for fecal boli. Produces individual cleaned feature files ready for merging.

**Inputs** — `--input_dir` should contain one or more batch subdirectories with these files (searched recursively):

| Pattern in filename | Content |
|---------------------|---------|
| `qc_batch_` | QC report CSVs — defines the expected set of videos |
| `gait` | Gait feature CSVs |
| `fecal_boli` | Raw fecal boli CSVs |
| `feature` | JABS behavior feature CSVs |
| `morpho` | Morphometrics CSVs |

**Assumptions:**
1. Missing gait speed bins mean no gaits were predicted for those bins — not that gait prediction was skipped.
2. All videos to be analyzed appear in at least one `qc_batch_*.csv` file.

**Running (Python):**
```bash
python python/2_qc_check.py \
    --input_dir /path/to/NextflowOutput \
    --output_dir /path/to/project_output \
    --param config/QC_params.yaml     # optional, overrides defaults
```

**Running (R CLI):**
```bash
Rscript r/2_qc_check_cli.R \
    --input_dir /path/to/NextflowOutput \
    --output_dir /path/to/project_output \
    --param config/QC_params.yaml
```

**QC Parameters (`config/QC_params.yaml`):**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `expected_length` | 108150 s | Expected video duration |
| `max_tracklet_per_hour` | 6 | Max pose tracklets per hour |
| `max_missing_pose` | 0.5% | Max fraction of missing pose frames |
| `max_missing_segmentation` | 20% | Max fraction of missing segmentation frames |
| `max_missing_keypoint` | 1% | Max fraction of missing keypoint frames |
| `fecal_boli_quantile_plotting` | 0.05 | Lower quantile threshold for fecal boli outlier plots |

**Outputs:**
```
output_dir/
├── final_nextflow_feature_data/
│   ├── gait_final.csv            wide format; stride counts padded with 0 for missing bins;
│   │                             variances set to NA for bins with <3 strides
│   ├── morphometrics_final.csv
│   ├── JABS_features_final.csv
│   └── fecal_boli_raw.csv        ← requires manual review before merging
└── qc/
    ├── nextflow_qc_logs/
    │   ├── qc_all.csv            all videos with pass/fail flags
    │   └── qc_failed.csv         failed videos only
    ├── missing_or_dup_data/
    │   ├── missing_data.csv      in QC report but missing from feature files
    │   ├── videos_not_in_qc_report.csv
    │   └── duplicated_data.xlsx  identical or highly correlated rows
    └── qc_figs/
        └── fecal_boli_qc_figs.pdf
```

**Manual steps after running:**
1. Review `qc/nextflow_qc_logs/qc_failed.csv`. Overlay pose on flagged videos and decide which to keep or exclude.
2. Create `videos_to_exclude.txt` in your output directory — one `NetworkFilename` per line.
3. Review `qc/qc_figs/fecal_boli_qc_figs.pdf`. Correct erroneous counts in `final_nextflow_feature_data/fecal_boli_raw.csv`, then save as `fecal_boli_final.csv` in the same directory.

---

### 3 — Combine Batches · `r/3_combine_batches.R`

Merges the four cleaned feature files with `metadata.csv`, checks that all videos and mice are represented, and outputs a single merged dataset.

**Inputs** — expects outputs from step 2 in `final_nextflow_feature_data/`:
- `gait_final.csv`, `morphometrics_final.csv`, `JABS_features_final.csv`, `fecal_boli_final.csv`
- `metadata.csv` — must have a `MouseID` column whose values appear as substrings in `NetworkFilename`
- `videos_to_exclude.txt` — created during manual QC in step 2

**Outputs:**
```
output_dir/
├── merged_nextflow_dataset.csv
└── qc/missing_or_dup_data/
    ├── NetworkFilenames_missing_in_data.csv
    └── mice_missing_in_metadata.csv
```

Mismatches between metadata and data do not stop the merge but should be reviewed.

---

### 4 — Outliers · `r/4_outliers.R`

Template script for final dataset preparation — intended to be customized per project. Removes zero-variance phenotypes, generates z-score outlier plots, and produces the final analysis-ready CSV.

**Edit the `Set QC and other Values` section to configure:**
- Project name prefix for output files
- Feature substrings to manually exclude
- Z-score threshold for outlier detection
- Number of top outlier mice/phenotypes to plot

**Inputs:** `merged_nextflow_dataset.csv` from step 3

**Outputs:**
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

### 5 — Heatmap · `r/5_heatmap.R`

Generates phenotype correlation heatmaps. Under active development.

---

### 6a — Behavior Clip Extraction · `python/6a_qc_classifiers.py`

Samples random video clips at detected behavior bouts, with optional pose skeleton overlay. Reads merged behavior CSV tables produced by Nextflow and writes short MP4 clips for review in step 6b.

```bash
# Single behavior, no pose overlay:
python python/6a_qc_classifiers.py \
    --behavior-csv NextflowOutput/batch_aa/merged_behavior_tables/merged_escape_bouts_merged.csv \
    --video-dir NextflowOutput/batch_aa/results/videos/ \
    --output-dir /tmp/qc_clips/ \
    --n-clips 3

# All behaviors in a folder, with pose skeleton:
python python/6a_qc_classifiers.py \
    --behavior-dir NextflowOutput/batch_aa/merged_behavior_tables/ \
    --video-dir NextflowOutput/batch_aa/results/videos/ \
    --output-dir /tmp/qc_clips/ \
    --n-clips 5 --overlay-pose
```

---

### 6b — Behavior Clip Viewer · `python/6b_qc_viewer.py`

Streamlit app for reviewing and annotating the clips produced by `6a_qc_classifiers.py`. Shows one clip at a time with Accept / Reject / Skip buttons; saves verdicts to `annotations.csv`.

```bash
streamlit run python/6b_qc_viewer.py -- \
    --clips-dir /tmp/qc_clips/ \
    --annotations-csv /tmp/annotations.csv
```

---

### Utility — Pose Corner Correction · `python/pose_corner_correction.py`

Copies `pose_est_v6.h5` files and embeds manually corrected corner coordinates from SLEAP annotation files.

```bash
python python/pose_corner_correction.py \
    --input_dir /path/to/NextflowOutput \
    --output_dir /path/to/pose_v6_dir
```

Expects each batch subdirectory to optionally contain:
- `manual_corner_correction.slp` — SLEAP file with corrected corners
- `failed_corners/` — `*_pose_est_v6.h5` files to update

---

## Expected Output File Structure

```
/project_output_dir/
├── videos_to_exclude.txt                        (manually created)
├── metadata.csv                                 (manually provided)
├── YOUR_PROJECT_final_nextflow_dataset.csv      (step 4 output → downstream analysis)
├── features_removed_from_curated_dataset.csv
├── merged_nextflow_dataset.csv                  (step 3 output)
├── final_nextflow_feature_data/                 (step 2 outputs)
│   ├── gait_final.csv
│   ├── morphometrics_final.csv
│   ├── JABS_features_final.csv
│   ├── fecal_boli_raw.csv
│   └── fecal_boli_final.csv                     (manually corrected)
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
