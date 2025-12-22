# -*- coding: utf-8 -*-
"""
Created on Sat Dec 20 14:41:56 2025

@author: inayas1
"""

import shutil
import subprocess
from pathlib import Path
from typing import List, Dict, Any


def process_h264(animal: List[Dict[str, Any]], owr: int = 0, fps: int = 60) -> List[Dict[str, Any]]:
    """
    Convert .h264 camera videos to .mp4 and store in pdir.

    Parameters
    ----------
    animal : list of dict
        Output from get_exp_info() (list of per-animal dicts).
    owr : int, default 0
        Overwrite flag. If 0 and output exists, skip conversion.
    fps : int, default 60
        Input frame rate for ffmpeg (-r).

    Returns
    -------
    animal : list of dict
        Updated with animal[i]["video"]["mp4"][cam] and ["video"]["led"][cam].
    """

    # Check ffmpeg availability
    if shutil.which("ffmpeg") is None:
        raise RuntimeError("ffmpeg not found on PATH. Install ffmpeg or add it to PATH.")

    cams = ["face", "pupil", "paws"]

    for i, a in enumerate(animal):
        video = a.get("video", {})
        h264 = video.get("h264", None)
        if not isinstance(h264, dict):
            continue

        # Ensure output dicts exist
        video.setdefault("mp4", {})
        video.setdefault("led", {})
        a["video"] = video  # write back in case it was missing

        pdir = Path(a["pdir"])
        pdir.mkdir(parents=True, exist_ok=True)

        for cam in cams:
            in_file = h264.get(cam, "")
            if not in_file:
                continue

            in_path = Path(in_file)
            if not in_path.exists():
                print(f"WARNING: Could not find {in_path} for animal index {i} cam {cam}")
                continue

            base = in_path.stem  # fileparts base
            out_mp4_path = pdir / f"{base}.mp4"
            out_csv_path = pdir / f"{base}_intensity.csv"

            # If mp4 exists and overwrite off, skip conversion (like MATLAB)
            if out_mp4_path.exists() and int(owr) == 0:
                print(f"Skipping (already exists): {out_mp4_path.name}")
                # MATLAB: stores full path in the "skip" branch
                video["mp4"][cam] = str(out_mp4_path)
                video["led"][cam] = str(out_csv_path)
                continue

            # Build ffmpeg command (stream copy, no re-encode)
            # cmd = [
            #     "ffmpeg", "-y",
            #     "-r", str(fps),
            #     "-i", str(in_path),
            #     "-c:v", "copy",
            #     str(out_mp4_path),
            # ]

            print(f'Converting {in_path} -> {out_mp4_path} ... (ffmpeg progress hidden)')
            try:
                # capture output so you don't get ffmpeg spam, like your MATLAB note
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError:
                print(f"WARNING: ffmpeg failed for {in_path} (animal index {i}, cam {cam})")
                continue

            # MATLAB stores *just the file name* after converting.
            # Keep the same behavior here:
            video["mp4"][cam] = f"{base}.mp4"
            video["led"][cam] = str(out_csv_path)

    return animal
