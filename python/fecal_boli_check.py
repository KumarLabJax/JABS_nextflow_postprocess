#!/usr/bin/env python3
"""Render a QC collage of arena-corner + fecal-boli overlays for pose_v6 videos.

For each `_pose_est_v6.h5` file, finds its matching .mp4 video, picks up
to 4 representative frames (biased toward fecal-boli count changes), draws
the arena corners and fecal-boli detections on each, and tiles them into a
single collage PNG.
"""

import argparse
from pathlib import Path
from typing import List, Optional, Tuple

import cv2
import h5py
import numpy as np

H5_SUFFIX = "_pose_est_v6.h5"
LABEL_BAR_HEIGHT = 30

PANEL_LABEL_COLOR = (255, 255, 255)
CORNER_POINT_COLOR = (0, 255, 255)
# Green = corners came from the automatic corner-detection model.
CORNER_LINE_COLOR_AUTO_DETECTED = (0, 220, 0)
# Red = corners did NOT come from the automatic model (cm_per_pixel_source
# != "corner_detection"). This usually means corners that were manually annotated/corrected 
# (e.g. via manual_corner_correction.slp)
CORNER_LINE_COLOR_NOT_AUTO_DETECTED = (0, 0, 255)
BOLI_POINT_COLOR = (0, 140, 255)
NO_CORNERS_COLOR = (0, 0, 0)


def find_video_for_h5(h5_path: Path) -> Optional[Path]:
    """Return the video path corresponding to a pose H5 file, if it exists."""
    video_path = h5_path.with_name(
        h5_path.name.removesuffix(H5_SUFFIX) + ".mp4"
    )
    return video_path if video_path.exists() else None


def load_corners(h5_path: Path) -> Tuple[Optional[np.ndarray], bool]:
    """Find corners coordinate and whether it is automatically detected or not"""
    with h5py.File(h5_path, "r") as f:
        if "static_objects/corners" not in f:
            return None, False
        corners = f["static_objects/corners"][:].astype(np.float32)
        source = f["poseest"].attrs.get("cm_per_pixel_source", b"")
        if isinstance(source, bytes):
            source = source.decode("utf-8", errors="ignore")
        auto_detected = source == "corner_detection"
        return corners, auto_detected


def load_fecal_boli(h5_path: Path) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Extract fecal boli metrics"""
    with h5py.File(h5_path, "r") as f:
        if "dynamic_objects/fecal_boli" not in f:
            return None
        grp = f["dynamic_objects/fecal_boli"]
        counts = grp["counts"][:].reshape(-1).astype(np.int64)
        points = grp["points"][:]
        sample_indices = grp["sample_indices"][:]
        return counts, points, sample_indices


def select_sample_indices(counts: np.ndarray, n_samples: int) -> List[int]:
    """Choose up to 4 frames to plot"""
    prev = np.concatenate(([0], counts[:-1]))
    changes = np.nonzero(counts != prev)[0]

    if len(changes) == 0:
        picks = np.linspace(0, n_samples - 1, min(4, n_samples))
        return sorted(set(int(round(p)) for p in picks))

    if len(changes) <= 4:
        return sorted(changes.tolist())

    picks = np.linspace(0, len(changes) - 1, 4)
    picks = sorted(set(int(round(p)) for p in picks))
    return [int(changes[i]) for i in picks]


def draw_corners(img: np.ndarray, corners: Optional[np.ndarray], auto_detected: bool) -> None:
    """Draw corners on a single frame"""
    if corners is None:
        cv2.putText(
            img, "NO CORNERS", (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.8,
            NO_CORNERS_COLOR, 2, cv2.LINE_AA,
        )
        return

    # Draw individual corner locations.
    for x, y in corners:
        cv2.circle(
            img, (int(x), int(y)), 5,
            CORNER_POINT_COLOR, -1, cv2.LINE_AA,
        )

    # Draw lining
    color = (
        CORNER_LINE_COLOR_AUTO_DETECTED
        if auto_detected
        else CORNER_LINE_COLOR_NOT_AUTO_DETECTED
    )
    hull = cv2.convexHull(corners.astype(np.float32)).reshape(-1, 2)
    cv2.polylines(
        img, [hull.astype(np.int32)], isClosed=True,
        color=color, thickness=2, lineType=cv2.LINE_AA,
    )


def draw_fecal_boli(img: np.ndarray, points_for_sample: np.ndarray, count: int) -> None:
    """Draw fecal boli on a single frame."""
    for j in range(count):
        y, x = points_for_sample[j]
        cv2.circle(img, (int(x), int(y)), 6, BOLI_POINT_COLOR, 2, cv2.LINE_AA)


def render_panel(
    frame_bgr: np.ndarray,
    corners: Optional[np.ndarray],
    auto_detected: bool,
    boli_points: Optional[np.ndarray],
    boli_count: int,
    label_text: str,
) -> np.ndarray:
    """Render a video frame"""
    panel = frame_bgr.copy()
    draw_corners(panel, corners, auto_detected)
    if boli_points is not None:
        draw_fecal_boli(panel, boli_points, boli_count)

    # Add panel label at the bottom of the frame
    h, w = panel.shape[:2]
    labeled = np.zeros((h + LABEL_BAR_HEIGHT, w, 3), dtype=panel.dtype)
    labeled[:h] = panel
    cv2.putText(
        labeled, label_text, (8, h + LABEL_BAR_HEIGHT - 9), cv2.FONT_HERSHEY_SIMPLEX,
        0.55, PANEL_LABEL_COLOR, 1, cv2.LINE_AA,
    )
    return labeled


def blank_panel(size: Tuple[int, int], text: str = "N/A") -> np.ndarray:
    h, w = size
    panel = np.full((h, w, 3), 40, dtype=np.uint8)
    cv2.putText(
        panel, text, (w // 2 - 30, h // 2), cv2.FONT_HERSHEY_SIMPLEX,
        1.0, (150, 150, 150), 2, cv2.LINE_AA,
    )
    return panel


def build_collage(panels: List[np.ndarray]) -> np.ndarray:
    h, w = panels[0].shape[:2]
    while len(panels) < 4:
        panels.append(blank_panel((h, w)))

    top = np.hstack([panels[0], panels[1]])
    bottom = np.hstack([panels[2], panels[3]])
    return np.vstack([top, bottom])


def format_timestamp(frame_idx: int, fps: float) -> str:
    if fps <= 0:
        return "t=?"
    total_seconds = frame_idx / fps
    minutes = int(total_seconds // 60)
    seconds = int(total_seconds % 60)
    return f"{minutes:02d}:{seconds:02d}"


def process_one(h5_path: Path, output_dir: Path) -> bool:
    video_path = find_video_for_h5(h5_path)
    if video_path is None:
        print(f"[skip] no video next to {h5_path}")
        return False

    corners, auto_detected = load_corners(h5_path)
    fecal_boli = load_fecal_boli(h5_path)

    cap = cv2.VideoCapture(str(video_path))
    fps = cap.get(cv2.CAP_PROP_FPS)
    total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))

    if fecal_boli is not None:
        counts, points, sample_indices = fecal_boli
        sample_picks = select_sample_indices(counts, len(counts))
    else:
        sample_picks = None

    panels = []
    if sample_picks is not None:
        for i in sample_picks:
            frame_idx = int(sample_indices[i])
            cap.set(cv2.CAP_PROP_POS_FRAMES, frame_idx)
            ok, frame = cap.read()
            if not ok:
                continue
            count = int(counts[i])
            boli_points = points[i] if count > 0 else None
            label = f"frame {frame_idx} ({format_timestamp(frame_idx, fps)})  boli={count}"
            panels.append(render_panel(frame, corners, auto_detected, boli_points, count, label))
    else:
        picks = np.linspace(0, max(total_frames - 1, 0), 4)
        for frame_idx in picks:
            frame_idx = int(frame_idx)
            cap.set(cv2.CAP_PROP_POS_FRAMES, frame_idx)
            ok, frame = cap.read()
            if not ok:
                continue
            label = f"frame {frame_idx} ({format_timestamp(frame_idx, fps)})  boli=n/a"
            panels.append(render_panel(frame, corners, auto_detected, None, 0, label))

    cap.release()

    if not panels:
        print(f"[skip] could not read any frames for {video_path}")
        return False

    collage = build_collage(panels)
    out_path = output_dir / f"{video_path.stem}_collage.png"
    cv2.imwrite(out_path, collage)
    print(f"[ok] wrote {out_path}")
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--file",
        type=Path,
        help="Single pose H5 file to process.",
    )
    input_group.add_argument(
        "--directory",
        type=Path,
        help="Directory of pose H5 files to process recursively.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory to write collage PNGs into.",
    )
    parser.add_argument(
        "--pattern",
        default=f"*{H5_SUFFIX}",
        help=f"Glob pattern for directory mode (default: *{H5_SUFFIX})",
    )

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.file:
        h5_paths = [args.file]
    else:
        h5_paths = sorted(
            path
            for path in args.directory.rglob(args.pattern)
        )

    if not h5_paths:
        parser.error("No pose H5 files found.")

    ok_count = 0
    for h5_path in h5_paths:
        if process_one(h5_path, args.output_dir):
            ok_count += 1

    print(
        f"\nDone: {ok_count}/{len(h5_paths)} collages written "
        f"to {args.output_dir}"
    )


if __name__ == "__main__":
    main()
