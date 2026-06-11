qc_failed_file <- read.csv("~/kumar-group/SING-grant/postNextflow/qc/nextflow_qc_logs/qc_failed.csv")
render_dir <- "/home/nguyetu/kumar-group/SING-grant/postNextflow/rendered_poses"

video_files <- file.path(dirname(qc_failed_file$QC_file), "results",
                         gsub("_with_fecal_boli", ".mp4", qc_failed_file$video_name))
video_files <- gsub("_filtered.mp4", ".mp4", video_files)

pose_v6_files <- file.path(dirname(qc_failed_file$QC_file), "results",
                           gsub("_with_fecal_boli.h5", "_pose_est_v6.h5", qc_failed_file$pose_file))
pose_v6_files <- gsub("_filtered.h5", "_pose_est_v6.h5", pose_v6_files)

render_outputs <- file.path(render_dir, gsub("trimmed.mp4", "overlay.avi", basename(video_files)))
commands <- paste("mouse-tracking utils render-pose", video_files, pose_v6_files, render_outputs)

# Ensure both inputs exist
to_render <- file.exists(video_files) & file.exists(pose_v6_files)
cat(commands[to_render], sep = "\n")
