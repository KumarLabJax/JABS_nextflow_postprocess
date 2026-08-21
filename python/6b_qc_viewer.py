#!/usr/bin/env python3
"""Streamlit QC viewer for behavior clip snippets produced by qc_classifiers.py.

Launch:
    streamlit run script/qc_viewer.py -- --clips-dir /tmp/qc_clips/
    streamlit run script/qc_viewer.py -- --clips-dir /tmp/qc_clips/ --annotations-csv /tmp/ann.csv
"""

import argparse
import re
from pathlib import Path

import pandas as pd
import streamlit as st


def _parse_args():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--clips-dir", type=Path, required=True)
    parser.add_argument("--annotations-csv", type=Path, default=None)
    args, _ = parser.parse_known_args()
    if args.annotations_csv is None:
        args.annotations_csv = args.clips_dir / "annotations.csv"
    return args


_CLIP_RE = re.compile(r"(.+)_bout(\d+)_frame(\d+)\.mp4$")


def scan_clips(clips_dir: Path) -> list[dict]:
    clips = []
    for mp4 in sorted(clips_dir.rglob("*.mp4")):
        behavior = mp4.relative_to(clips_dir).parts[0]
        m = _CLIP_RE.match(mp4.name)
        if not m:
            continue
        clips.append({
            "path": str(mp4),
            "behavior": behavior,
            "video_id": m.group(1),
            "bout_idx": int(m.group(2)),
            "frame_start": int(m.group(3)),
        })
    return clips


def load_annotations(csv_path: Path) -> dict:
    if not csv_path.exists():
        return {}
    df = pd.read_csv(csv_path, dtype=str).fillna("")
    return {
        row["clip_path"]: {"status": row["status"], "notes": row["notes"]}
        for _, row in df.iterrows()
    }


def save_annotations(annotations: dict, csv_path: Path) -> None:
    rows = [
        {"clip_path": p, "status": a["status"], "notes": a["notes"]}
        for p, a in annotations.items()
    ]
    pd.DataFrame(rows, columns=["clip_path", "status", "notes"]).to_csv(csv_path, index=False)


def record_verdict(clip: dict, status: str, notes: str, csv_path: Path) -> None:
    st.session_state.annotations[clip["path"]] = {"status": status, "notes": notes}
    save_annotations(st.session_state.annotations, csv_path)


def main():
    st.set_page_config(page_title="Behavior QC Viewer", layout="wide")
    args = _parse_args()

    # Initialize session state once
    if "clips" not in st.session_state:
        st.session_state.clips = scan_clips(args.clips_dir)
    if "annotations" not in st.session_state:
        st.session_state.annotations = load_annotations(args.annotations_csv)
    if "idx" not in st.session_state:
        st.session_state.idx = 0

    clips = st.session_state.clips
    if not clips:
        st.error(f"No clips found in {args.clips_dir}")
        st.stop()

    behaviors = sorted({c["behavior"] for c in clips})

    # ------------------------------------------------------------------ Sidebar
    with st.sidebar:
        st.title("Behavior QC")

        selected_behavior = st.selectbox("Behavior filter", ["All"] + behaviors, key="filter")

        # Reset position when filter changes
        if st.session_state.get("_last_filter") != selected_behavior:
            st.session_state.idx = 0
            st.session_state["_last_filter"] = selected_behavior

        filtered = [
            c for c in clips
            if selected_behavior == "All" or c["behavior"] == selected_behavior
        ]
        n_total = len(filtered)
        n_annotated = sum(
            1 for c in filtered if c["path"] in st.session_state.annotations
        )

        st.progress(
            n_annotated / n_total if n_total else 0,
            text=f"{n_annotated} / {n_total} reviewed",
        )

        if st.button("Jump to next unannotated"):
            for i, c in enumerate(filtered):
                if c["path"] not in st.session_state.annotations:
                    st.session_state.idx = i
                    st.rerun()

        st.divider()

        if st.button("💾 Save annotations"):
            save_annotations(st.session_state.annotations, args.annotations_csv)
            st.success(f"Saved to {args.annotations_csv}")

        st.caption(f"`{args.annotations_csv}`")

    # ------------------------------------------------------------------ Guard
    if not filtered:
        st.warning("No clips match the current filter.")
        st.stop()

    st.session_state.idx = max(0, min(st.session_state.idx, n_total - 1))
    clip = filtered[st.session_state.idx]
    ann = st.session_state.annotations.get(clip["path"], {"status": "", "notes": ""})

    # ------------------------------------------------------------------ Header
    st.markdown(
        f"**Clip {st.session_state.idx + 1} of {n_total}** &nbsp;|&nbsp; `{clip['behavior']}`"
    )

    # Show current verdict badge if already annotated
    if ann["status"]:
        color = {"accept": "green", "reject": "red", "skip": "gray"}.get(ann["status"], "gray")
        st.markdown(
            f"<span style='color:{color}'>● Previously marked: <b>{ann['status']}</b></span>",
            unsafe_allow_html=True,
        )

    # ------------------------------------------------------------------ Video
    st.video(clip["path"], loop=True, autoplay=True)

    # ------------------------------------------------------------------ Metadata
    c1, c2, c3, c4 = st.columns(4)
    vid_label = clip["video_id"]
    if len(vid_label) > 24:
        vid_label = vid_label[:21] + "…"
    c1.metric("Behavior", clip["behavior"])
    c2.metric("Video ID", vid_label)
    c3.metric("Bout #", clip["bout_idx"])
    c4.metric("Start frame", clip["frame_start"])

    st.divider()

    # ------------------------------------------------------------------ Notes + verdict
    notes = st.text_input(
        "Notes (optional)",
        value=ann.get("notes", ""),
        key=f"notes_{st.session_state.idx}",
    )

    btn_accept, btn_reject, btn_skip, _, _ = st.columns([1, 1, 1, 2, 2])
    advance = False

    if btn_accept.button("✓ Accept", type="primary", use_container_width=True):
        record_verdict(clip, "accept", notes, args.annotations_csv)
        advance = True
    if btn_reject.button("✗ Reject", type="secondary", use_container_width=True):
        record_verdict(clip, "reject", notes, args.annotations_csv)
        advance = True
    if btn_skip.button("→ Skip", use_container_width=True):
        record_verdict(clip, "skip", notes, args.annotations_csv)
        advance = True

    if advance and st.session_state.idx < n_total - 1:
        st.session_state.idx += 1
        st.rerun()

    # ------------------------------------------------------------------ Navigation
    st.divider()
    nav_prev, nav_counter, nav_next = st.columns([1, 3, 1])

    if nav_prev.button("← Prev", use_container_width=True, disabled=st.session_state.idx == 0):
        st.session_state.idx -= 1
        st.rerun()

    nav_counter.markdown(
        f"<div style='text-align:center;padding-top:8px'>"
        f"{st.session_state.idx + 1} / {n_total}</div>",
        unsafe_allow_html=True,
    )

    if nav_next.button("Next →", use_container_width=True, disabled=st.session_state.idx == n_total - 1):
        st.session_state.idx += 1
        st.rerun()


if __name__ == "__main__":
    main()
