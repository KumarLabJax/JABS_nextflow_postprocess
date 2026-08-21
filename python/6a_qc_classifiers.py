#!/usr/bin/env python3
"""Sample random QC video clips at detected behavior bouts, with optional pose overlay.

Usage examples:

  # Single behavior, no pose:
  python script/qc_classifiers.py \
    --behavior-csv NextflowOutput/batch_aa/merged_behavior_tables/merged_escape_bouts_merged.csv \
    --video-dir NextflowOutput/batch_aa/results/ \
    --output-dir /tmp/qc_clips/ --n-clips 3

  # All behaviors in a folder, with pose skeleton:
  python script/qc_classifiers.py \
    --behavior-dir NextflowOutput/batch_aa/merged_behavior_tables/ \
    --video-dir NextflowOutput/batch_aa/results/ \
    --output-dir /tmp/qc_clips/ --n-clips 5 --overlay-pose

--video-dir must be the 'results' directory containing one namespace/
batch_folder subdirectory per video (e.g. results/videos/batch_1/videos.mp4),
matching the '{namespace} {batch_folder} {video_id}' prefix on video_name.
"""

import argparse
import re
import sys
from pathlib import Path
from urllib.parse import unquote

import cv2
import h5py
import imageio.v2 as imageio
import numpy as np
import pandas as pd
from tqdm import tqdm

# ---------------------------------------------------------------------------
# 12-keypoint JABS mouse skeleton.
# If your model uses a different keypoint ordering, edit KEYPOINT_NAMES and
# SKELETON_EDGES here — these are the only two constants you need to change.
# ---------------------------------------------------------------------------
KEYPOINT_NAMES = [
    "nose",            # 0
    "left_ear",        # 1
    "right_ear",       # 2
    "base_neck",       # 3
    "left_front_paw",  # 4
    "right_front_paw", # 5
    "center_spine",    # 6
    "left_rear_paw",   # 7
    "right_rear_paw",  # 8
    "base_tail",       # 9
    "mid_tail",        # 10
    "tip_tail",        # 11
]
SKELETON_EDGES = [
    (0, 3),                    # nose → base_neck (spine axis)
    (3, 6),                    # base_neck → center_spine
    (6, 4), (6, 5),            # center_spine → front paws
    (6, 9),                    # center_spine → base_tail
    (9, 7), (9, 8),            # base_tail → rear paws
    (9, 10), (10, 11),         # base_tail → mid_tail → tip_tail
]
KEYPOINT_COLOR = (0, 255, 0)   # BGR green dots
EDGE_COLOR = (0, 165, 255)     # BGR orange lines
POSE_CONFIDENCE_THRESHOLD = 0.3


# ---------------------------------------------------------------------------
# CSV parsing
# ---------------------------------------------------------------------------

def parse_behavior_csv(csv_path: Path) -> pd.DataFrame:
    """Read a merged behavior table; return positive-bout rows with behavior_label added."""
    df = pd.read_csv(csv_path, skiprows=2)
    df = df[df["is_behavior"] == 1].copy()
    df["video_name"] = df["video_name"].apply(unquote)
    df["behavior_label"] = _label_from_filename(csv_path.name)
    df["bout_order"] = df.sort_values(["video_name", "start"]).groupby("video_name").cumcount()
    return df


def _label_from_filename(filename: str) -> str:
    m = re.match(r"merged_(.+?)_bouts_merged", filename)
    return m.group(1) if m else Path(filename).stem


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

def _ids_match(video_id: str, stem: str) -> bool:
    """True if stem is video_id (or vice versa) with a '_'/'-'/'.' delimited suffix added/removed."""
    if video_id == stem:
        return True
    shorter, longer = (video_id, stem) if len(video_id) < len(stem) else (stem, video_id)
    if not longer.startswith(shorter):
        return False
    return longer[len(shorter)] in "_-."


def video_name_to_paths(video_name: str, video_dir: Path) -> tuple[Path | None, Path | None]:
    """Return (mp4_path, h5_path) for a decoded video_name string.

    video_name carries a '{namespace} {batch_folder} {video_id}' prefix. Tries an
    exact-filename path first (a single stat call, no directory walk), and only
    falls back to a recursive fuzzy-id search if that miss.
    """
    namespace, batch_folder, video_id = video_name.split(" ", 2)
    search_dir = video_dir / namespace / batch_folder

    mp4 = search_dir / f"{video_id}.mp4"
    if not mp4.exists():
        mp4 = next((p for p in search_dir.rglob("*.mp4") if _ids_match(video_id, p.stem)), None)

    search_root = mp4.parent if mp4 is not None else search_dir
    h5 = search_root / f"{video_id}.h5"
    if not h5.exists():
        h5 = next((p for p in search_root.rglob("*.h5") if _ids_match(video_id, p.stem)), None)

    return mp4, h5


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------

def sample_bouts(df: pd.DataFrame, n_clips: int, seed: int = 42) -> pd.DataFrame:
    """Sample up to n_clips positive bouts per video, reproducibly."""
    rng = np.random.default_rng(seed)
    groups = []
    for _, grp in df.groupby("video_name", sort=False):
        k = min(n_clips, len(grp))
        idx = rng.choice(len(grp), size=k, replace=False)
        groups.append(grp.iloc[sorted(idx)])
    return pd.concat(groups, ignore_index=True) if groups else df.iloc[:0]


# ---------------------------------------------------------------------------
# Pose loading + drawing
# ---------------------------------------------------------------------------

def load_pose_frames(h5_path: Path, frame_start: int, frame_end: int):
    """Return (points, confidence) arrays for frame range [frame_start, frame_end).

    points     shape: (n_frames, 12, 2) float32
    confidence shape: (n_frames, 12)    float32
    """
    with h5py.File(h5_path, "r") as f:
        pts  = f["poseest/points"][frame_start:frame_end, 0, :, :]   # (T,12,2)
        conf = f["poseest/confidence"][frame_start:frame_end, 0, :]   # (T,12)
    return pts.astype(np.float32), conf.astype(np.float32)


def draw_pose(frame: np.ndarray, keypoints_yx: np.ndarray, confidence: np.ndarray) -> np.ndarray:
    """Overlay skeleton edges and keypoint dots onto frame (in-place)."""
    for a, b in SKELETON_EDGES:
        if confidence[a] >= POSE_CONFIDENCE_THRESHOLD and confidence[b] >= POSE_CONFIDENCE_THRESHOLD:
            p1 = (int(keypoints_yx[a, 1]), int(keypoints_yx[a, 0]))
            p2 = (int(keypoints_yx[b, 1]), int(keypoints_yx[b, 0]))
            cv2.line(frame, p1, p2, EDGE_COLOR, 2, cv2.LINE_AA)
    for i in range(len(keypoints_yx)):
        if confidence[i] >= POSE_CONFIDENCE_THRESHOLD:
            cy, cx = int(keypoints_yx[i, 0]), int(keypoints_yx[i, 1])
            cv2.circle(frame, (cx, cy), 4, KEYPOINT_COLOR, -1, cv2.LINE_AA)
    return frame


# ---------------------------------------------------------------------------
# Clip extraction
# ---------------------------------------------------------------------------

def extract_clip(
    mp4_path: Path,
    start_frame: int,
    duration: int,
    padding: int,
    output_path: Path,
    h5_path: Path | None = None,
) -> bool:
    """Write a padded clip [start-padding, start+duration+padding] to output_path.

    If h5_path is provided and exists, overlays the pose skeleton on each frame.
    Returns True on success.
    """
    cap = cv2.VideoCapture(str(mp4_path))
    if not cap.isOpened():
        tqdm.write(f"  [warn] cannot open {mp4_path}")
        return False

    fps   = cap.get(cv2.CAP_PROP_FPS) or 30.0
    total = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))

    clip_start = max(0, start_frame - padding)
    clip_end   = min(total, start_frame + duration + padding)
    n_frames   = clip_end - clip_start

    output_path.parent.mkdir(parents=True, exist_ok=True)
    # libx264 (via imageio-ffmpeg) instead of cv2's mp4v so clips play in Chrome.
    writer = imageio.get_writer(
        str(output_path), fps=fps, codec="libx264", pixelformat="yuv420p"
    )

    pts, conf = None, None
    if h5_path is not None and h5_path.exists():
        try:
            pts, conf = load_pose_frames(h5_path, clip_start, clip_end)
        except Exception as exc:
            tqdm.write(f"  [warn] pose load failed ({h5_path.name}): {exc}")

    cap.set(cv2.CAP_PROP_POS_FRAMES, clip_start)
    for local_i in range(n_frames):
        ok, frame = cap.read()
        if not ok:
            break
        if pts is not None and local_i < len(pts):
            draw_pose(frame, pts[local_i], conf[local_i])
        writer.append_data(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB))

    cap.release()
    writer.close()
    return True


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Sample random QC video clips at detected behavior bouts."
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--behavior-csv",
        type=Path,
        metavar="CSV",
        help="Single merged behavior CSV (e.g. merged_escape_bouts_merged.csv)",
    )
    src.add_argument(
        "--behavior-dir",
        type=Path,
        metavar="DIR",
        help="Directory of merged_*_bouts_merged.csv files",
    )
    parser.add_argument(
        "--video-dir",
        type=Path,
        required=True,
        metavar="DIR",
        help="'results' directory holding namespace/batch_folder subdirectories of "
             "MP4 (and H5 pose) files, matching the '{namespace} {batch_folder} "
             "{video_id}' prefix on video_name",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        metavar="DIR",
        help="Root directory for output clips (one sub-folder per behavior)",
    )
    parser.add_argument(
        "--n-clips",
        type=int,
        default=5,
        metavar="N",
        help="Clips to sample per video (default: 5)",
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=30,
        metavar="FRAMES",
        help="Padding frames before/after each bout (default: 30)",
    )
    parser.add_argument(
        "--overlay-pose",
        action="store_true",
        help="Draw pose skeleton on clips (requires H5 pose file alongside video)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible sampling (default: 42)",
    )
    args = parser.parse_args()

    if args.behavior_csv:
        csv_paths = [args.behavior_csv]
    else:
        csv_paths = sorted(args.behavior_dir.glob("merged_*_bouts_merged.csv"))
        if not csv_paths:
            sys.exit(f"No merged_*_bouts_merged.csv files found in {args.behavior_dir}")

    total_clips = 0
    for csv_path in csv_paths:
        df = parse_behavior_csv(csv_path)
        if df.empty:
            print(f"[skip] {csv_path.name} — no positive bouts")
            continue

        sampled = sample_bouts(df, args.n_clips, seed=args.seed)
        label = df["behavior_label"].iloc[0]
        n_videos = df["video_name"].nunique()

        desc = f"{label} ({n_videos} videos)"
        for _, row in tqdm(list(sampled.iterrows()), desc=desc, unit="clip"):
            mp4_path, h5_path = video_name_to_paths(row["video_name"], args.video_dir)
            if mp4_path is None:
                tqdm.write(f"  [warn] video not found for: {row['video_name']}")
                continue

            namespace, batch_folder, video_id = row["video_name"].split(" ", 2)
            out_name = f"{video_id}_bout{row['bout_order']:03d}_frame{int(row['start'])}.mp4"
            output_path = args.output_dir / label / namespace / batch_folder / out_name

            ok = extract_clip(
                mp4_path=mp4_path,
                start_frame=int(row["start"]),
                duration=int(row["duration"]),
                padding=args.padding,
                output_path=output_path,
                h5_path=h5_path if args.overlay_pose else None,
            )
            if ok:
                total_clips += 1

    print(f"\nDone. {total_clips} clip(s) written to {args.output_dir}")


if __name__ == "__main__":
    main()
