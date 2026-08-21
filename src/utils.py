import cv2
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from pathlib import Path
from ipyfilechooser import FileChooser
from IPython.display import display, clear_output
import ipywidgets as widgets


# ---------------------------------------------------------------------------
# Core video helpers
# ---------------------------------------------------------------------------

def sample_frames(video_path, n_frames=5, start_frame=0, end_frame=None):
    """
    Sample evenly spaced frames from a video file.

    Parameters
    ----------
    video_path : str or Path
    n_frames : int
        Number of frames to sample.
    start_frame : int
        First frame index to sample from (inclusive). Default 0.
    end_frame : int or None
        Last frame index to sample to (inclusive). Default None = last frame.
    """
    cap = cv2.VideoCapture(str(video_path))
    frame_count = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))

    lo = max(0, start_frame)
    hi = (frame_count - 1) if end_frame is None else min(end_frame, frame_count - 1)
    indices = np.linspace(lo, hi, n_frames, dtype=int)

    frames = []
    for idx in indices:
        cap.set(cv2.CAP_PROP_POS_FRAMES, idx)
        ret, frame = cap.read()
        if ret:
            frames.append((idx, cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)))
    cap.release()
    return frames


def show_frames(video_path, n_frames=5, start_frame=0, end_frame=None, axes=None):
    """
    Display sampled frames from a video.

    Parameters
    ----------
    video_path : str or Path
    n_frames : int
    start_frame : int
        First frame index to sample from. Default 0.
    end_frame : int or None
        Last frame index to sample to. Default None = last frame.
    axes : array-like of Axes, optional
        If provided, draw into these axes instead of creating a new figure.
        Must have at least n_frames elements.
    """
    frames = sample_frames(video_path, n_frames, start_frame=start_frame, end_frame=end_frame)
    if not frames:
        print(f"Could not read frames from {video_path}")
        return

    video_path = Path(video_path)
    ncols = 5
    nrows = max(1, -(-n_frames // ncols))  # ceiling division

    if axes is None:
        fig, axes_arr = plt.subplots(nrows, ncols, figsize=(15, 3 * nrows))
        axes_arr = np.array(axes_arr).flatten()
        own_fig = True
    else:
        axes_arr = np.array(axes).flatten()
        own_fig = False

    for ax, (idx, img) in zip(axes_arr, frames):
        ax.imshow(img)
        ax.set_title(f"Frame {idx}")
        ax.axis("off")
    for ax in axes_arr[len(frames):]:
        ax.axis("off")

    if own_fig:
        axes_arr[0].figure.suptitle(f"{video_path.parent.name}/{video_path.name}")
        plt.tight_layout()
        plt.show()


def show_sample_frames_from_directory(directory, keywords=None, n_frames=5,
                                       start_frame=0, end_frame=None,
                                       videos_per_page=10):
    """
    Display sample frames for every video found under *directory*.

    Parameters
    ----------
    directory : str or Path
    keywords : list[str] or None
        Only show videos whose name or parent folder contains any keyword.
    n_frames : int
    start_frame : int
        First frame index to sample from. Default 0.
    end_frame : int or None
        Last frame index to sample to. Default None = last frame.
    videos_per_page : int
        Cap on how many videos to render (first page only).
    """
    folder = Path(directory).expanduser()
    if not folder.is_dir():
        print(f"Directory not found: {folder}")
        return

    all_videos = sorted(folder.rglob("*.mp4")) + sorted(folder.rglob("*.avi"))

    if keywords:
        kws = [k.strip().lower() for k in keywords]
        all_videos = [
            v for v in all_videos
            if any(k in v.name.lower() or k in v.parent.name.lower() for k in kws)
        ]

    if not all_videos:
        print(f"No .mp4 or .avi files found in {folder}")
        return

    for video_path in all_videos[:videos_per_page]:
        show_frames(video_path, n_frames=n_frames, start_frame=start_frame, end_frame=end_frame)


# ---------------------------------------------------------------------------
# Pose helpers  (custom .h5 files from the Kumar-lab pipeline)
# ---------------------------------------------------------------------------

def load_pose(h5_path, conf_threshold=0.1):
    """
    Load a Kumar-lab pose estimation .h5 file.

    Expected datasets
    -----------------
    poseest/points      : (n_frames, n_tracks, n_nodes, 2)  uint16  x/y pixels
    poseest/confidence  : (n_frames, n_tracks, n_nodes)     float32
    poseest/id_mask     : (n_frames, n_tracks)              bool    valid track flag

    Returns
    -------
    trx : np.ndarray, shape (n_frames, n_tracks, n_nodes, 2), float32
        x/y coordinates; NaN where id_mask is False or confidence < conf_threshold.
    skeleton_edges : list of (int, int)
        Always empty — not stored in the h5 format.
    node_names : list[str]
        Generic names ["kp_0", "kp_1", …].
    """
    import h5py

    with h5py.File(str(h5_path), "r") as f:
        pts  = f["poseest/points"][:].astype(np.float32)   # (F, T, N, 2)
        conf = f["poseest/confidence"][:]                   # (F, T, N)

    # Keypoints with zero confidence are missing detections — set to NaN
    missing = conf < conf_threshold                         # (F, T, N)
    trx = pts.copy()
    trx[missing] = np.nan

    n_nodes = trx.shape[2]

    node_names = [f"kp_{i}" for i in range(n_nodes)]
    return trx, [], node_names


def find_pose_file(video_path, pose_dir):
    """
    Search *pose_dir* for a pose .h5 matching *video_path*.

    Matching rule: the video stem must appear as a substring of the h5 stem
    (case-insensitive).  Prefers ``*_pose_est_v6.h5``; falls back to
    ``*_pose_est_v2.h5``; finally any ``*.h5`` match.

    Returns the matched Path or None.
    """
    video_stem = Path(video_path).stem.lower()
    pose_dir = Path(pose_dir)

    candidates = [p for p in pose_dir.rglob("*.h5")
                  if video_stem in p.stem.lower()]

    if not candidates:
        return None

    for suffix in ("_pose_est_v6.h5", "_pose_est_v2.h5"):
        for p in candidates:
            if p.name.endswith(suffix):
                return p

    return candidates[0]


def draw_keypoints(img, keypoints_per_track, skeleton_edges=None,
                   radius=4, alpha=0.85):
    """
    Overlay pose keypoints (and optional skeleton edges) onto *img*.

    Parameters
    ----------
    img : np.ndarray, shape (H, W, 3), uint8, RGB
    keypoints_per_track : np.ndarray, shape (n_tracks, n_nodes, 2)
        x/y coordinates; NaN for missing keypoints.
    skeleton_edges : list of (int, int), optional
    radius : int
        Keypoint dot radius in pixels.
    alpha : float
        Blend factor for the overlay (1.0 = opaque).

    Returns
    -------
    np.ndarray  RGB image with keypoints drawn.
    """
    overlay = img.copy()
    n_tracks = keypoints_per_track.shape[0]

    # One distinct colour per track (BGR for cv2)
    cmap = cm.get_cmap("tab10")
    track_colors = [
        tuple(int(c * 255) for c in cmap(i % 10)[:3])[::-1]
        for i in range(n_tracks)
    ]

    for t, kps in enumerate(keypoints_per_track):
        color = track_colors[t]

        # Draw skeleton edges first (underneath dots)
        if skeleton_edges:
            for src, dst in skeleton_edges:
                p1 = kps[src]
                p2 = kps[dst]
                if np.any(np.isnan(p1)) or np.any(np.isnan(p2)):
                    continue
                cv2.line(overlay,
                         (int(p1[0]), int(p1[1])),
                         (int(p2[0]), int(p2[1])),
                         color, thickness=2, lineType=cv2.LINE_AA)

        # Draw keypoint dots
        for x, y in kps:
            if np.isnan(x) or np.isnan(y):
                continue
            cv2.circle(overlay, (int(x), int(y)), radius, color,
                       thickness=-1, lineType=cv2.LINE_AA)

    return cv2.addWeighted(overlay, alpha, img, 1 - alpha, 0)


def show_frames_with_pose(video_path, pose_path, n_frames=5,
                          start_frame=0, end_frame=None):
    """
    Display sampled frames from *video_path* with SLEAP pose overlay from *pose_path*.

    Falls back to plain frame display if pose data cannot be loaded.

    Parameters
    ----------
    video_path : str or Path
    pose_path : str or Path
    n_frames : int
    start_frame : int
        First frame index to sample from. Default 0.
    end_frame : int or None
        Last frame index to sample to. Default None = last frame.
    """
    video_path = Path(video_path)
    frames = sample_frames(video_path, n_frames, start_frame=start_frame, end_frame=end_frame)
    if not frames:
        print(f"Could not read frames from {video_path}")
        return

    try:
        trx, skeleton_edges, _ = load_pose(pose_path)
        has_pose = True
    except Exception as e:
        print(f"Could not load pose from {pose_path}: {e}")
        has_pose = False

    ncols = 5
    nrows = max(1, -(-n_frames // ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 3 * nrows))
    axes = np.array(axes).flatten()

    for ax, (idx, img) in zip(axes, frames):
        if has_pose and idx < trx.shape[0]:
            kps = trx[idx]  # (n_tracks, n_nodes, 2)
            img = draw_keypoints(img, kps, skeleton_edges)
        ax.imshow(img)
        ax.set_title(f"Frame {idx}")
        ax.axis("off")
    for ax in axes[len(frames):]:
        ax.axis("off")

    fig.suptitle(f"{video_path.parent.name}/{video_path.name}")
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Motion screening helpers
# ---------------------------------------------------------------------------

def motion_score(video_path, n_frames=10):
    """Return mean absolute pixel difference between consecutive sampled frames."""
    cap = cv2.VideoCapture(str(video_path))
    frame_count = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    if frame_count == 0:
        cap.release()
        return 0.0
    indices = np.linspace(0, frame_count - 1, n_frames, dtype=int)

    diffs, prev_gray = [], None
    for idx in indices:
        cap.set(cv2.CAP_PROP_POS_FRAMES, idx)
        ret, frame = cap.read()
        if not ret:
            continue
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY).astype(np.float32)
        if prev_gray is not None:
            diffs.append(np.mean(np.abs(gray - prev_gray)))
        prev_gray = gray
    cap.release()
    return float(np.mean(diffs)) if diffs else 0.0


def is_empty_video(video_path, threshold=3.0):
    """Return (is_empty, score). Videos with score < threshold are flagged empty."""
    score = motion_score(video_path)
    return score < threshold, score


# ---------------------------------------------------------------------------
# Interactive browser — no pose
# ---------------------------------------------------------------------------

def browse_videos(videos_per_page=10):
    """
    Launch an interactive Jupyter widget to browse videos in a chosen folder.
    Folder chooser + keyword filter + frame range controls + paginated previews.
    """
    state = {"videos": [], "page_index": 0}

    chooser = FileChooser()
    chooser.title = "<b>Select Video Directory</b>"
    chooser.show_only_dirs = True

    txt_filter = widgets.Text(
        placeholder="Enter keywords (comma-separated)",
        description="Filter:",
        layout=widgets.Layout(width="400px"),
    )
    btn_apply = widgets.Button(description="Apply Filter")
    filter_row = widgets.HBox([txt_filter, btn_apply])

    int_n_frames = widgets.BoundedIntText(
        value=5, min=1, max=50, description="# frames:", layout=widgets.Layout(width="160px")
    )
    int_start = widgets.IntText(
        value=0, description="Start frame:", layout=widgets.Layout(width="180px")
    )
    int_end = widgets.IntText(
        value=-1, description="End frame:", layout=widgets.Layout(width="180px"),
    )
    frame_help = widgets.Label("(end = -1 means last frame)")
    frame_row = widgets.HBox([int_n_frames, int_start, int_end, frame_help])

    btn_prev = widgets.Button(description="Previous")
    btn_next = widgets.Button(description="Next")
    label_page = widgets.Label(value="Page 0 / 0")
    out = widgets.Output()
    pagination = widgets.HBox([btn_prev, btn_next, label_page])

    display(chooser, filter_row, frame_row, pagination, out)

    def _end_frame():
        v = int_end.value
        return None if v < 0 else v

    def _update():
        videos, page = state["videos"], state["page_index"]
        start, end = page * videos_per_page, (page + 1) * videos_per_page
        total = max(1, -(-len(videos) // videos_per_page))
        label_page.value = f"Page {page + 1} / {total}"
        with out:
            clear_output(wait=True)
            for vp in videos[start:end]:
                show_frames(vp, n_frames=int_n_frames.value,
                            start_frame=int_start.value, end_frame=_end_frame())

    def _load(_=None):
        if not chooser.selected_path:
            return
        folder = Path(chooser.selected_path).expanduser()
        if not folder.is_dir():
            with out:
                clear_output()
                print(f"Directory not found: {folder}")
            return
        all_vids = sorted(folder.rglob("*.mp4")) + sorted(folder.rglob("*.avi"))
        kws = txt_filter.value.strip().lower()
        if kws:
            ks = [k.strip() for k in kws.split(",")]
            all_vids = [v for v in all_vids
                        if any(k in v.name.lower() or k in v.parent.name.lower() for k in ks)]
        state["videos"], state["page_index"] = list(all_vids), 0
        if not state["videos"]:
            with out:
                clear_output()
                print(f"No .mp4 or .avi files found in {folder}")
        else:
            _update()

    def _prev(_):
        if state["page_index"] > 0:
            state["page_index"] -= 1
            _update()

    def _next(_):
        if (state["page_index"] + 1) * videos_per_page < len(state["videos"]):
            state["page_index"] += 1
            _update()

    chooser.register_callback(_load)
    btn_apply.on_click(_load)
    btn_prev.on_click(_prev)
    btn_next.on_click(_next)


# ---------------------------------------------------------------------------
# Interactive browser — with pose overlay
# ---------------------------------------------------------------------------

def browse_videos_with_pose(videos_per_page=10):
    """
    Interactive widget to browse videos with optional JABS pose overlay.

    Two folder choosers:
      - Video directory
      - Pose directory (.slp files; leave unset to skip overlay)

    For each video the browser looks for a .slp file whose stem is a substring
    of (or contains) the video stem.  If found, keypoints are drawn; otherwise
    the plain frame is shown.
    """
    state = {"videos": [], "page_index": 0}

    # --- Video chooser ---
    vid_chooser = FileChooser()
    vid_chooser.title = "<b>Select Video Directory</b>"
    vid_chooser.show_only_dirs = True

    # --- Pose chooser ---
    pose_chooser = FileChooser()
    pose_chooser.title = "<b>Select Pose Directory (.h5 files) — optional</b>"
    pose_chooser.show_only_dirs = True

    chooser_row = widgets.HBox([vid_chooser, pose_chooser])

    # --- Filter ---
    txt_filter = widgets.Text(
        placeholder="Enter keywords (comma-separated)",
        description="Filter:",
        layout=widgets.Layout(width="400px"),
    )
    btn_apply = widgets.Button(description="Apply Filter")
    filter_row = widgets.HBox([txt_filter, btn_apply])

    # --- Frame range ---
    int_n_frames = widgets.BoundedIntText(
        value=5, min=1, max=50, description="# frames:", layout=widgets.Layout(width="160px")
    )
    int_start = widgets.IntText(
        value=0, description="Start frame:", layout=widgets.Layout(width="180px")
    )
    int_end = widgets.IntText(
        value=-1, description="End frame:", layout=widgets.Layout(width="180px")
    )
    frame_help = widgets.Label("(end = -1 means last frame)")
    frame_row = widgets.HBox([int_n_frames, int_start, int_end, frame_help])

    # --- Pagination ---
    btn_prev = widgets.Button(description="Previous")
    btn_next = widgets.Button(description="Next")
    label_page = widgets.Label(value="Page 0 / 0")
    out = widgets.Output()
    pagination = widgets.HBox([btn_prev, btn_next, label_page])

    display(chooser_row, filter_row, frame_row, pagination, out)

    def _end_frame():
        v = int_end.value
        return None if v < 0 else v

    def _render_video(video_path):
        """Render one video — with pose if a matching .slp can be found."""
        pose_dir = pose_chooser.selected_path
        kw = dict(n_frames=int_n_frames.value,
                  start_frame=int_start.value, end_frame=_end_frame())
        if pose_dir:
            slp = find_pose_file(video_path, pose_dir)
            if slp:
                show_frames_with_pose(video_path, slp, **kw)
                return
        show_frames(video_path, **kw)

    def _update():
        videos, page = state["videos"], state["page_index"]
        start, end = page * videos_per_page, (page + 1) * videos_per_page
        total = max(1, -(-len(videos) // videos_per_page))
        label_page.value = f"Page {page + 1} / {total}"
        with out:
            clear_output(wait=True)
            for vp in videos[start:end]:
                _render_video(vp)

    def _load(_=None):
        if not vid_chooser.selected_path:
            return
        folder = Path(vid_chooser.selected_path).expanduser()
        if not folder.is_dir():
            with out:
                clear_output()
                print(f"Directory not found: {folder}")
            return
        all_vids = sorted(folder.rglob("*.mp4")) + sorted(folder.rglob("*.avi"))
        kws = txt_filter.value.strip().lower()
        if kws:
            ks = [k.strip() for k in kws.split(",")]
            all_vids = [v for v in all_vids
                        if any(k in v.name.lower() or k in v.parent.name.lower() for k in ks)]
        state["videos"], state["page_index"] = list(all_vids), 0
        if not state["videos"]:
            with out:
                clear_output()
                print(f"No .mp4 or .avi files found in {folder}")
        else:
            _update()

    def _prev(_):
        if state["page_index"] > 0:
            state["page_index"] -= 1
            _update()

    def _next(_):
        if (state["page_index"] + 1) * videos_per_page < len(state["videos"]):
            state["page_index"] += 1
            _update()

    vid_chooser.register_callback(_load)
    btn_apply.on_click(_load)
    btn_prev.on_click(_prev)
    btn_next.on_click(_next)


# ---------------------------------------------------------------------------
# Interactive browser — motion screening
# ---------------------------------------------------------------------------

def browse_videos_motion_screen(videos_per_page=6, threshold=3.0):
    """
    Browse videos and flag ones with low motion (likely empty OFA arena).

    Parameters
    ----------
    videos_per_page : int
    threshold : float
        Motion score below this value → flagged as potentially empty.
    """
    state = {"videos": [], "scores": {}, "page_index": 0}

    chooser = FileChooser()
    chooser.title = "<b>Select Video Directory</b>"
    chooser.show_only_dirs = True

    int_n_frames = widgets.BoundedIntText(
        value=5, min=1, max=50, description="# frames:", layout=widgets.Layout(width="160px")
    )
    int_start = widgets.IntText(
        value=0, description="Start frame:", layout=widgets.Layout(width="180px")
    )
    int_end = widgets.IntText(
        value=-1, description="End frame:", layout=widgets.Layout(width="180px")
    )
    frame_help = widgets.Label("(end = -1 means last frame)")
    frame_row = widgets.HBox([int_n_frames, int_start, int_end, frame_help])

    btn_prev = widgets.Button(description="Previous")
    btn_next = widgets.Button(description="Next")
    label_page = widgets.Label(value="Page 0 / 0")
    chk_show_all = widgets.Checkbox(value=False, description="Show all videos")
    out = widgets.Output()
    pagination = widgets.HBox([btn_prev, btn_next, chk_show_all, label_page])

    display(chooser, frame_row, pagination, out)

    def _end_frame():
        v = int_end.value
        return None if v < 0 else v

    def _visible():
        if chk_show_all.value:
            return state["videos"]
        return [v for v in state["videos"]
                if state["scores"].get(v, (False, 0))[0]]

    def _update():
        to_show = _visible()
        page = state["page_index"]
        start, end = page * videos_per_page, (page + 1) * videos_per_page
        total = max(1, -(-len(to_show) // videos_per_page)) if to_show else 1
        label_page.value = f"Page {page + 1} / {total}"
        with out:
            clear_output(wait=True)
            if not to_show:
                print("All videos contain movement — no empty OFA detected.")
                return
            for vp in to_show[start:end]:
                empty, score = state["scores"].get(vp, (False, 0.0))
                tag = "  [EMPTY?]" if empty else ""
                print(f"{vp.name}  —  motion score = {score:.2f}{tag}")
                show_frames(vp, n_frames=int_n_frames.value,
                            start_frame=int_start.value, end_frame=_end_frame())

    def _load(_=None):
        if not chooser.selected_path:
            return
        folder = Path(chooser.selected_path).expanduser()
        if not folder.is_dir():
            with out:
                clear_output()
                print(f"Directory not found: {folder}")
            return
        videos = sorted(folder.rglob("*.mp4")) + sorted(folder.rglob("*.avi"))
        if not videos:
            with out:
                clear_output()
                print(f"No .mp4 or .avi files found in {folder}")
            return
        with out:
            clear_output()
            print(f"Scoring {len(videos)} videos for motion…")
        scores = {v: is_empty_video(v, threshold=threshold) for v in videos}
        state["videos"], state["scores"], state["page_index"] = list(videos), scores, 0
        _update()

    def _prev(_):
        if state["page_index"] > 0:
            state["page_index"] -= 1
            _update()

    def _next(_):
        if (state["page_index"] + 1) * videos_per_page < len(_visible()):
            state["page_index"] += 1
            _update()

    chooser.register_callback(_load)
    btn_prev.on_click(_prev)
    btn_next.on_click(_next)
    chk_show_all.observe(lambda _: _update(), names="value")
