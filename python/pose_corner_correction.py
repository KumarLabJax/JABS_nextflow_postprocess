import argparse
import h5py
import sleap_io as sio
from pathlib import Path
import shutil


def create_pose_v6(slp_correction, failed_pose_dir, pose_v6_dir):
    labels = sio.load_file(slp_correction)

    Path(pose_v6_dir).mkdir(parents=True, exist_ok=True)
    for label in labels:
        corners = label.instances[0].numpy()

        file_name = Path(label.video.filename[0]).stem
        pose_file = failed_pose_dir / f"{file_name}_pose_est_v6.h5"

        if pose_file.exists():
            pose_v6_name = f"{file_name.split('%20')[-1]}_pose_est_v6.h5"
            pose_v6_path = Path(pose_v6_dir) / pose_v6_name
            shutil.copy2(pose_file, pose_v6_path)

            with h5py.File(pose_v6_path, "r+") as f:
                static_objects = f.require_group("static_objects")

                if "corners" in static_objects:
                    static_objects["corners"][:] = corners
                else:
                    static_objects.create_dataset("corners", data=corners)


def main():
    parser = argparse.ArgumentParser(
        description="Copy pose_v6 files and embed manually corrected corner coordinates."
    )
    parser.add_argument(
        "--input_dir", required=True,
        help="Nextflow output directory containing batch subdirectories."
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Destination directory for corrected pose_v6 files."
    )
    args = parser.parse_args()

    for batch_dir in Path(args.input_dir).iterdir():
        slp_correction = batch_dir / "manual_corner_correction.slp"
        failed_pose_dir = batch_dir / "failed_corners"
        if slp_correction.exists():
            create_pose_v6(slp_correction, failed_pose_dir, args.output_dir)


if __name__ == "__main__":
    main()
