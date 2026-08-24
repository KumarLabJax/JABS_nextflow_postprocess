"""Rescue pose files whose automated arena-corner detection failed.

The JABS pipeline stores arena corner coordinates in each pose H5 file under
``static_objects/corners``.  When automated detection fails, the video lands in
a ``failed_corners/`` subdirectory.  A human then re-labels the four corners in
SLEAP and saves the result as ``manual_corner_correction.slp`` alongside the
batch.

This script takes a single ``manual_corner_correction.slp`` file and, for each
annotated video:

1. Copies the corresponding ``_pose_est_v6.h5`` from ``failed_corners/`` to the
   output directory (decoding URL-encoded ``%20`` spaces in filenames).
2. Overwrites the ``static_objects/corners`` dataset in the copied H5 file with
   the manually labelled corner coordinates.

Usage::

    python python/pose_corner_correction.py \\
        --slp_file /path/to/batch/manual_corner_correction.slp \\
        --output_dir /path/to/pose_v6_dir

By default, the failed pose files are looked up in a ``failed_corners/``
directory next to the ``.slp`` file. Pass ``--failed_pose_dir`` to override.
"""

import argparse
import h5py
import sleap_io as sio
from pathlib import Path
import shutil


def create_pose_v6(slp_correction, failed_pose_dir, results_dir):
    """Copy and correct pose H5 files for one batch.

    Parameters
    ----------
    slp_correction : Path
        SLEAP ``.slp`` file containing one labelled frame per failed video,
        where the single instance holds the four corrected corner coordinates.
    failed_pose_dir : Path
        Directory containing ``*_pose_est_v6.h5`` files that need correction.
    results_dir : str or Path
        Destination directory for the corrected H5 files.
    """
    labels = sio.load_file(slp_correction)

    Path(results_dir).mkdir(parents=True, exist_ok=True)
    for label in labels:
        corners = label.instances[0].numpy()

        file_name = Path(label.video.filename[0]).stem
        pose_file = failed_pose_dir / f"{file_name}_pose_est_v6.h5"

        if pose_file.exists():
            pose_v6_name = f"{file_name.replace('%20', "/")}_pose_est_v6.h5"
            pose_v6_path = Path(results_dir) / pose_v6_name
            pose_v6_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(pose_file, pose_v6_path)

            with h5py.File(pose_v6_path, "r+") as f:
                static_objects = f.require_group("static_objects")

                if "corners" in static_objects:
                    static_objects["corners"][:] = corners
                else:
                    static_objects.create_dataset("corners", data=corners)
            print(f"Wrote {pose_v6_path}")
        else:
            print(f"{pose_file} does not exist")


def main():
    parser = argparse.ArgumentParser(
        description="Copy pose_v6 files and embed manually corrected corner coordinates."
    )
    parser.add_argument(
        "--slp_file", required=True,
        help="Path to a manual_corner_correction.slp file."
    )
    parser.add_argument(
        "--failed_pose_dir",
        help="Directory containing the *_pose_est_v6.h5 files needing correction. "
             "Defaults to a 'failed_corners' directory next to --slp_file."
    )
    parser.add_argument(
        "--output_dir",
        help="Destination directory for corrected pose_v6 files. "
             "Defaults to a 'results' directory next to --slp_file."
    )
    args = parser.parse_args()

    slp_correction = Path(args.slp_file)
    failed_pose_dir = (
        Path(args.failed_pose_dir) if args.failed_pose_dir
        else slp_correction.parent / "failed_corners"
    )

    output_dir = (
        Path(args.output_dir) if args.output_dir
        else slp_correction.parent / "results"
    )

    create_pose_v6(slp_correction, failed_pose_dir, output_dir)


if __name__ == "__main__":
    main()
