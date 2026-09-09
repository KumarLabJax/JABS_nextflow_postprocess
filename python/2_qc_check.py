#!/usr/bin/env python3
"""
A script to clean NextFlow outputs for final processing
The same as "NextFlow_Output_QC_Postprocess_1.R" but to be run in a script
Developed by Dr. Jake Beierle (don't forget the Dr., it's important)

Documentation:
See comprehensive documentation on the github repository
https://github.com/jacobbeierle/JABS_nextflow_postprocess/tree/main
"""
# %%
import argparse
import os
import sys
from pathlib import Path
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings

warnings.filterwarnings('ignore')
# %%
def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="QC Reporting Pipeline")
    
    parser.add_argument("--input_dir", type=str, required=True,
                        help="Input directory (required)")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Output directory (required)")
    parser.add_argument("--param", type=str, default=None,
                        help="Optional YAML config file to override defaults")
    
    # QC parameters with defaults
    parser.add_argument("--expected_length", type=int, default=60*60*30 + 5*30,
                        help="Expected video length in seconds [default: 108150]")
    parser.add_argument("--max_tracklet_per_hour", type=int, default=6,
                        help="Maximum tracklets per hour [default: 6]")
    parser.add_argument("--max_missing_pose", type=float, default=0.005,
                        help="Maximum fraction of missing pose [default: 0.005]")
    parser.add_argument("--max_missing_segmentation", type=float, default=0.2,
                        help="Maximum fraction of missing segmentation [default: 0.2]")
    parser.add_argument("--max_missing_keypoint", type=float, default=0.01,
                        help="Maximum fraction of missing keypoints [default: 0.01]")
    parser.add_argument("--fecal_boli_quantile_plotting", type=float, default=0.05,
                        help="Quantile for fecal boli plotting [default: 0.05]")
    
    return parser.parse_args()

# %%
def merge_parameters(args):
    """Merge parameter priorities from command line and YAML."""
    params = {
        'expected_length': args.expected_length,
        'max_tracklet_per_hour': args.max_tracklet_per_hour,
        'max_missing_pose': args.max_missing_pose,
        'max_missing_segmentation': args.max_missing_segmentation,
        'max_missing_keypoint': args.max_missing_keypoint,
        'fecal_boli_quantile_plotting': args.fecal_boli_quantile_plotting
    }
    
    # Override with YAML if provided
    if args.param is not None:
        with open(args.param, 'r') as f:
            yaml_vals = yaml.safe_load(f)
            params.update(yaml_vals)
    
    return params

# %%
def create_output_directories(output_dir):
    """Create output directory structure."""
    subdirectories = [
        "final_nextflow_feature_data",
        "qc/nextflow_qc_logs",
        "qc/missing_or_dup_data",
        "qc/qc_figs"
    ]
    
    for subdir in subdirectories:
        dir_path = os.path.join(output_dir, subdir)
        os.makedirs(dir_path, exist_ok=True)

# %%
def read_raw_data(input_dir, pattern):
    """
    Read in multiple CSV files from a directory and harmonize the ID column.
    
    Args:
        input_dir: path to the folder containing CSV files
        pattern: string pattern to match file names
    
    Returns:
        A DataFrame with NetworkFilename as the first column
    """
    # Find all matching files
    files = list(Path(input_dir).rglob(f"*{pattern}*"))
    
    if not files:
        return pd.DataFrame()
    
    # Read and concatenate all files
    dfs = []
    for file in files:
        df = pd.read_csv(file)
        dfs.append(df)
    
    raw_data = pd.concat(dfs, ignore_index=True)
    
    # Use NetworkFilename as the ID column
    if "NetworkFilename" not in raw_data.columns:
        raw_data.rename(columns={raw_data.columns[0]: "NetworkFilename"}, inplace=True)
    
    # Move NetworkFilename to first column
    cols = raw_data.columns.tolist()
    cols.insert(0, cols.pop(cols.index("NetworkFilename")))
    raw_data = raw_data[cols]
    
    # Clean and harmonize NetworkFilename
    raw_data['NetworkFilename'] = (raw_data['NetworkFilename']
                                   .str.replace('_corrected', '', regex=False)
                                   .str.replace('_filtered', '', regex=False)
                                   .str.replace('.avi', '', regex=False)
                                   .str.replace(r'^\.', '', regex=True)
                                   .str.replace(r'^/', '', regex=True))
    
    return raw_data

# %%
def check_missing_and_dup(expected_videos, data_df, corr_thres=0.99):
    """
    Check for missing videos and duplicated rows in a dataset.
    
    Args:
        expected_videos: list of expected NetworkFilename values
        data_df: DataFrame containing NetworkFilename column and data
        corr_thres: threshold for flagging correlated rows
    
    Returns:
        Dictionary with missing_qc, missing_output, and dup_data
    """
    if data_df.empty:
        return {
            'missing_qc': [],
            'missing_output': expected_videos,
            'dup_data': pd.DataFrame()
        }
    
    video_missing_output = list(set(expected_videos) - set(data_df['NetworkFilename']))
    video_missing_in_qc = list(set(data_df['NetworkFilename']) - set(expected_videos))
    
    # Check for duplicate rows
    dup_idx = data_df.iloc[:, 1:].duplicated(keep=False)
    
    # Check for highly correlated rows
    numeric_cols = data_df.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) > 0:
        numeric_data = data_df[numeric_cols].fillna(0)
        if len(numeric_data) > 1:
            # Standardize and compute correlation
            from scipy.stats import zscore
            standardized = numeric_data.apply(zscore, nan_policy='omit')
            corr_mat = standardized.T.corr()
            
            # Find highly correlated pairs
            corr_idx = np.zeros(len(data_df), dtype=bool)
            for i in range(len(corr_mat)):
                for j in range(i+1, len(corr_mat)):
                    if corr_mat.iloc[i, j] > corr_thres:
                        corr_idx[i] = True
                        corr_idx[j] = True
        else:
            corr_idx = np.zeros(len(data_df), dtype=bool)
    else:
        corr_idx = np.zeros(len(data_df), dtype=bool)
    
    # Union of duplicated and correlated rows
    all_idx = dup_idx | corr_idx
    duplicated_rows = data_df[all_idx]
    
    return {
        'missing_qc': video_missing_in_qc,
        'missing_output': video_missing_output,
        'dup_data': duplicated_rows
    }

# %%
def process_qc_logs(input_dir, output_dir, params):
    """Process and publish QC logs."""
    # Read QC files
    qc_files = list(Path(input_dir).rglob("qc_batch_*.csv"))
    
    if not qc_files:
        print("Warning: No QC files found!")
        return []
    
    qc_logs = []
    for file in qc_files:
        df = pd.read_csv(file)
        df['QC_file'] = str(file)
        qc_logs.append(df)
    
    qc_log = pd.concat(qc_logs, ignore_index=True)
    
    # Record why QC failed for each video
    qc_log['passed_duration_QC'] = qc_log['video_duration'] == params['expected_length']
    qc_log['passed_tracklet_QC'] = qc_log['pose_tracklets'] < params['max_tracklet_per_hour'] * params['expected_length'] / 108000
    qc_log['passed_segmentation_QC'] = qc_log['seg_counts'] > (1 - params['max_missing_segmentation']) * params['expected_length']
    qc_log['passed_pose_QC'] = qc_log['pose_counts'] > (1 - params['max_missing_pose']) * params['expected_length']
    qc_log['passed_kp_QC'] = qc_log['missing_keypoint_frames'] < params['max_missing_keypoint'] * params['expected_length']
    
    # Filter failed QC
    passed_cols = [col for col in qc_log.columns if col.startswith('passed_')]
    qc_log_failed = qc_log[~qc_log[passed_cols].all(axis=1)]
    
    # Write QC files
    qc_log.to_csv(os.path.join(output_dir, "qc/nextflow_qc_logs/qc_all.csv"), index=False)
    qc_log_failed.to_csv(os.path.join(output_dir, "qc/nextflow_qc_logs/qc_failed.csv"), index=False)
    
    # Get expected videos
    expected_videos = (qc_log['video_name']
                       .str.replace('_with_fecal_boli', '', regex=False)
                       .str.replace('_filtered', '', regex=False)
                       .str.replace(r'^/', '', regex=True)
                       .unique()
                       .tolist())
    
    return expected_videos

# %%
def plot_fecal_boli_qc(fecal_boli_raw, output_dir, params):
    """Generate fecal boli QC plots."""
    # Pivot data for plotting
    value_cols = [col for col in fecal_boli_raw.columns 
                  if col not in ['NetworkFilename', 'nextflow_version']]
    
    fecal_boli_plot = fecal_boli_raw.melt(
        id_vars=['NetworkFilename', 'nextflow_version'],
        value_vars=value_cols,
        var_name='min',
        value_name='fecal_boli'
    ).dropna(subset=['fecal_boli'])
    
    # Extract minute number
    fecal_boli_plot['min'] = pd.to_numeric(
        fecal_boli_plot['min'].str.extract(r'(\d+)')[0]
    )
    
    # Create PDF with plots
    pdf_path = os.path.join(output_dir, "qc/qc_figs/fecal_boli_qc_figs.pdf")
    
    with PdfPages(pdf_path) as pdf:
        # Plot 1: All mice growth curves
        fig, ax = plt.subplots(figsize=(6, 6))
        for name, group in fecal_boli_plot.groupby('NetworkFilename'):
            ax.plot(group['min'], group['fecal_boli'], alpha=0.5)
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('Fecal Boli Count')
        ax.set_title('Fecal boli growth, all mice')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # Plot 2: Lowest quantile
        max_boli = fecal_boli_plot.groupby('NetworkFilename')['fecal_boli'].max()
        lowest = max_boli.nsmallest(int(len(max_boli) * params['fecal_boli_quantile_plotting']))
        lowest_data = fecal_boli_plot[fecal_boli_plot['NetworkFilename'].isin(lowest.index)]
        
        fig, ax = plt.subplots(figsize=(6, 6))
        for name, group in lowest_data.groupby('NetworkFilename'):
            ax.plot(group['min'], group['fecal_boli'], alpha=0.5)
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('Fecal Boli Count')
        ax.set_title(f'Lowest {params["fecal_boli_quantile_plotting"]*100}% of fecal boli mice')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # Plot 3: Highest quantile
        highest = max_boli.nlargest(int(len(max_boli) * params['fecal_boli_quantile_plotting']))
        highest_data = fecal_boli_plot[fecal_boli_plot['NetworkFilename'].isin(highest.index)]
        
        fig, ax = plt.subplots(figsize=(6, 6))
        for name, group in highest_data.groupby('NetworkFilename'):
            ax.plot(group['min'], group['fecal_boli'], alpha=0.5)
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('Fecal Boli Count')
        ax.set_title(f'Highest {params["fecal_boli_quantile_plotting"]*100}% of fecal boli mice')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # Plot 4: Histogram of final counts
        final_counts = (fecal_boli_plot.sort_values('min', ascending=False)
                        .drop_duplicates('NetworkFilename')['fecal_boli'])
        
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.hist(final_counts, bins=range(int(final_counts.min()), int(final_counts.max())+2))
        ax.set_xlabel('Fecal Boli Count')
        ax.set_ylabel('Count')
        ax.set_title('Fecal boli highest bin, all mice')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

# %%
def process_gait_data(input_dir, output_dir, expected_videos):
    """Process gait data."""
    gait_raw = read_raw_data(input_dir, "gait.csv")
    
    if gait_raw.empty:
        return None, {}
    
    video_level_metrics = ["Distance Traveled", "Body Length", "Speed", 
                          "Speed Variance", "nextflow_version"]
    
    # Remove variance measures from speed bins with fewer than 3 strides
    variance_cols = [col for col in gait_raw.columns if 'Variance' in col 
                     and col not in video_level_metrics]
    
    for col in variance_cols:
        if 'Stride Count' in gait_raw.columns:
            gait_raw.loc[gait_raw['Stride Count'] < 3, col] = np.nan
    
    # Convert to wide format
    id_cols = ['NetworkFilename'] + [col for col in video_level_metrics 
                                     if col in gait_raw.columns]
    
    if 'Speed Bin' in gait_raw.columns:
        value_cols = [col for col in gait_raw.columns 
                     if col not in id_cols + ['Speed Bin']]
        
        gait_wide = gait_raw.pivot_table(
            index=id_cols,
            columns='Speed Bin',
            values=value_cols,
            aggfunc='first'
        ).reset_index()
        
        # Flatten column names
        gait_wide.columns = ['.'.join(map(str, col)).strip('.') if isinstance(col, tuple) 
                            else col for col in gait_wide.columns]
        
        # Fill NA stride counts with 0
        stride_cols = [col for col in gait_wide.columns if 'Stride Count' in col]
        for col in stride_cols:
            gait_wide[col] = gait_wide[col].fillna(0)
    else:
        gait_wide = gait_raw
    
    # Check for missing and duplicated data
    gait_summary = check_missing_and_dup(expected_videos, gait_wide)
    
    # Output to CSV
    gait_wide.to_csv(
        os.path.join(output_dir, "final_nextflow_feature_data/gait_final.csv"),
        index=False
    )
    
    return gait_wide, gait_summary

# %%
def main():
    """Main execution function."""
    args = parse_arguments()
    params = merge_parameters(args)
    
    input_dir = args.input_dir
    output_dir = args.output_dir
    
    # Print configuration
    # Print configuration
    print("=== QC CONFIGURATION ===")
    print(f"Input directory:  {input_dir}")
    print(f"Output directory: {output_dir}")
    print("QC Parameters:")
    for key, value in params.items():
        print(f"  {key}: {value}")
    
    # Create output directories
    create_output_directories(output_dir)
    
    # Process QC logs
    expected_videos = process_qc_logs(input_dir, output_dir, params)
    
    if not expected_videos:
        print("Warning: No expected videos found from QC logs!")
        expected_videos = []
    
    # Process fecal boli data
    print("\nProcessing fecal boli data...")
    fecal_boli_raw = read_raw_data(input_dir, "fecal_boli.csv")
    
    if not fecal_boli_raw.empty:
        fecal_boli_summary = check_missing_and_dup(expected_videos, fecal_boli_raw)
        fecal_boli_raw.to_csv(
            os.path.join(output_dir, "final_nextflow_feature_data/fecal_boli_raw.csv"),
            index=False
        )
        plot_fecal_boli_qc(fecal_boli_raw, output_dir, params)
    else:
        fecal_boli_summary = {'missing_qc': [], 'missing_output': expected_videos, 'dup_data': pd.DataFrame()}
    
    # Process gait data
    print("Processing gait data...")
    gait_wide, gait_summary = process_gait_data(input_dir, output_dir, expected_videos)
    
    # Process JABS features
    print("Processing JABS features...")
    jabs_features = read_raw_data(input_dir, "features.csv")
    
    if not jabs_features.empty:
        jabs_summary = check_missing_and_dup(expected_videos, jabs_features)
        jabs_features.to_csv(
            os.path.join(output_dir, "final_nextflow_feature_data/JABS_features_final.csv"),
            index=False
        )
    else:
        jabs_summary = {'missing_qc': [], 'missing_output': expected_videos, 'dup_data': pd.DataFrame()}
    
    # Process morphometrics
    print("Processing morphometrics...")
    morpho_raw = read_raw_data(input_dir, "morphometrics.csv")
    
    if not morpho_raw.empty:
        morpho_summary = check_missing_and_dup(expected_videos, morpho_raw)
        morpho_raw.to_csv(
            os.path.join(output_dir, "final_nextflow_feature_data/morphometrics_final.csv"),
            index=False
        )
    else:
        morpho_summary = {'missing_qc': [], 'missing_output': expected_videos, 'dup_data': pd.DataFrame()}
    
    # Report and output warnings
    all_missing_data = {
        'fecal_boli': fecal_boli_summary['missing_output'],
        'gait': gait_summary['missing_output'],
        'JABS_features': jabs_summary['missing_output'],
        'Morphometrics': morpho_summary['missing_output']
    }
    
    videos_not_in_qc_report = {
        'fecal_boli': fecal_boli_summary['missing_qc'],
        'gait': gait_summary['missing_qc'],
        'JABS_features': jabs_summary['missing_qc'],
        'morphometrics': morpho_summary['missing_qc']
    }
    
    all_duplicated_data = {
        'fecal_boli': fecal_boli_summary['dup_data'],
        'gait': gait_summary['dup_data'],
        'JABS_features': jabs_summary['dup_data'],
        'morphometrics': morpho_summary['dup_data']
    }
    
    no_missing_output = all(len(v) == 0 for v in all_missing_data.values())
    no_missing_qc = all(len(v) == 0 for v in videos_not_in_qc_report.values())
    no_dups = all(v.empty for v in all_duplicated_data.values())
    
    # Write final reports
    print("\n=== FINAL ERROR REPORT ===")
    
    if no_missing_output and no_missing_qc and no_dups:
        print("No errors to report")
        pd.DataFrame(["No data missing"]).to_csv(
            os.path.join(output_dir, "qc/missing_or_dup_data/missing_data.csv"),
            index=False, header=False
        )
        pd.DataFrame(["No data missing"]).to_csv(
            os.path.join(output_dir, "qc/missing_or_dup_data/videos_not_in_qc_report.csv"),
            index=False, header=False
        )
        pd.DataFrame(["No duplicated data"]).to_csv(
            os.path.join(output_dir, "qc/missing_or_dup_data/duplicated_data.csv"),
            index=False, header=False
        )
    else:
        # Report missing output data
        if not no_missing_output:
            missing_df_list = []
            for output_type, videos in all_missing_data.items():
                if videos:
                    df = pd.DataFrame({'video_path': videos, output_type: True})
                    missing_df_list.append(df)
            
            if missing_df_list:
                missing_df = missing_df_list[0]
                for df in missing_df_list[1:]:
                    missing_df = missing_df.merge(df, on='video_path', how='outer')
                missing_df.to_csv(
                    os.path.join(output_dir, "qc/missing_or_dup_data/missing_data.csv"),
                    index=False
                )
                print("Missing output data for:", ', '.join([k for k, v in all_missing_data.items() if v]))
        else:
            pd.DataFrame(["No data missing"]).to_csv(
                os.path.join(output_dir, "qc/missing_or_dup_data/missing_data.csv"),
                index=False, header=False
            )
        
        # Report missing QC data
        if not no_missing_qc:
            qc_df_list = []
            for output_type, videos in videos_not_in_qc_report.items():
                if videos:
                    df = pd.DataFrame({'video_path': videos, output_type: True})
                    qc_df_list.append(df)
            
            if qc_df_list:
                qc_df = qc_df_list[0]
                for df in qc_df_list[1:]:
                    qc_df = qc_df.merge(df, on='video_path', how='outer')
                qc_df.to_csv(
                    os.path.join(output_dir, "qc/missing_or_dup_data/videos_not_in_qc_report.csv"),
                    index=False
                )
                print("Missing video in QC for:", ', '.join([k for k, v in videos_not_in_qc_report.items() if v]))
        else:
            pd.DataFrame(["No data missing"]).to_csv(
                os.path.join(output_dir, "qc/missing_or_dup_data/videos_not_in_qc_report.csv"),
                index=False, header=False
            )
        
        # Report duplicated data
        if not no_dups:
            with pd.ExcelWriter(
                os.path.join(output_dir, "qc/missing_or_dup_data/duplicated_data.xlsx"),
                engine='openpyxl'
            ) as writer:
                for name, df in all_duplicated_data.items():
                    if not df.empty:
                        df.to_excel(writer, sheet_name=name[:31], index=False)  # Excel sheet name limit
            print("Duplicated data for:", ', '.join([k for k, v in all_duplicated_data.items() if not v.empty]))
        else:
            pd.DataFrame(["No duplicated data"]).to_csv(
                os.path.join(output_dir, "qc/missing_or_dup_data/duplicated_data.csv"),
                index=False, header=False
            )
    
    print("\nProcessing complete!")


if __name__ == "__main__":
    main()