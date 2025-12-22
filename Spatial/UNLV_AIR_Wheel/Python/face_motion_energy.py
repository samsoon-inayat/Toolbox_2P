# -*- coding: utf-8 -*-
"""
Created on Sat Dec 20 14:46:54 2025

@author: inayas1
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import cv2
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import json

def select_roi_first_frame(video_path):
    cap = cv2.VideoCapture(str(video_path))
    ok, frame = cap.read()
    cap.release()
    if not ok:
        raise RuntimeError("Could not read first frame.")
    r = cv2.selectROI("Select ROI (ENTER to confirm, ESC to cancel)", frame, showCrosshair=True)
    cv2.destroyAllWindows()
    x, y, w, h = map(int, r)
    if w == 0 or h == 0:
        raise RuntimeError("ROI selection cancelled.")
    return (x, y, w, h)

def motion_energy_roi(
    video_path,
    fps=60.0,
    roi=None,
    downsample=1,
    blur_ksize=5,
    diff_threshold=10,
    use_absdiff=True,
    out_csv=None
):
    """
    Computes motion energy inside ROI by frame-to-frame differences.

    Motion metric options:
      - mean(|frame_t - frame_{t-1}|) after thresholding
      - robust to global brightness changes if ROI is good
    """
    video_path = Path(video_path)
    cap = cv2.VideoCapture(str(video_path))
    if not cap.isOpened():
        raise RuntimeError(f"Could not open video: {video_path}")

    # If fps not provided, try reading from file
    file_fps = cap.get(cv2.CAP_PROP_FPS)
    if file_fps and file_fps > 1:
        fps = float(file_fps)

    # ROI selection if needed
    if roi is None:
        roi = select_roi_first_frame(video_path)
    x, y, w, h = roi

    motion = []
    t_s = []
    frame_idx = 0

    prev_gray = None

    while True:
        ok, frame = cap.read()
        if not ok:
            break

        if downsample > 1 and (frame_idx % downsample != 0):
            frame_idx += 1
            continue

        crop = frame[y:y+h, x:x+w]
        gray = cv2.cvtColor(crop, cv2.COLOR_BGR2GRAY)

        if blur_ksize and blur_ksize > 1:
            gray = cv2.GaussianBlur(gray, (blur_ksize, blur_ksize), 0)

        if prev_gray is None:
            prev_gray = gray
            motion.append(np.nan)
            t_s.append(frame_idx / fps)
            frame_idx += 1
            continue

        if use_absdiff:
            d = cv2.absdiff(gray, prev_gray)
        else:
            d = (gray.astype(np.float32) - prev_gray.astype(np.float32))

        # Threshold to suppress sensor noise
        d = d.astype(np.float32)
        d[d < diff_threshold] = 0

        # Motion energy: mean of remaining diffs
        m = float(np.mean(d))

        motion.append(m)
        t_s.append(frame_idx / fps)

        prev_gray = gray
        frame_idx += 1

    cap.release()

    df = pd.DataFrame({
        "frame": np.arange(len(motion)),
        "time_s": t_s,
        "motion_energy": motion
    })

    if out_csv is None:
        out_csv = video_path.with_suffix("").as_posix() + "_motion_energy.csv"
    out_csv = Path(out_csv)
    df.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")

    return df, roi, fps




def save_or_load_motion_energy(
    animal: List[Dict[str, Any]],
    an_idx: int,
    cam: str,
    df: Optional[pd.DataFrame] = None,
    roi: Optional[Tuple[int, int, int, int]] = None,
    fps: Optional[float] = None,
    base_name: Optional[str] = None,
    owr: int = 0,
    ) -> Tuple[List[Dict[str, Any]], pd.DataFrame, Dict[str, Any]]:
    """
    Save/load motion energy outputs into animal[an_idx]['pdir'] and record names in animal.

    If CSV+JSON exist and owr==0, loads them and returns (animal, df_loaded, meta_loaded).

    Parameters
    ----------
    animal : list of dict
        Your experiment structure.
    an_idx : int
        Which animal entry to operate on.
    cam : str
        'face' / 'pupil' / 'paws' (or any string key).
    df : pandas.DataFrame or None
        Computed motion-energy dataframe (required if saving new).
    roi : (x,y,w,h) or None
        ROI used (required if saving new).
    fps : float or None
        FPS used (required if saving new).
    base_name : str or None
        Base name for files. If None, will try to infer from mp4 filename.
        Example: base_name="face" -> face_motion_energy.csv / face_motion_energy_meta.json
    owr : int
        Overwrite flag. If 0 and outputs exist, load instead of saving.

    Returns
    -------
    animal, df_out, meta_out
    """

    a = animal[an_idx]
    pdir = Path(a["pdir"])
    pdir.mkdir(parents=True, exist_ok=True)

    # Ensure dicts exist to store output references
    a.setdefault("analysis", {})
    a["analysis"].setdefault("motion_energy", {})
    a["analysis"]["motion_energy"].setdefault(cam, {})

    # Infer base_name if not provided: use mp4 path/name if available
    if base_name is None:
        mp4_entry = (
            a.get("video", {})
             .get("mp4", {})
             .get(cam, "")
        )
        if mp4_entry:
            mp4_path = Path(mp4_entry)
            base_name = mp4_path.stem  # "face" or "face_something"
        else:
            base_name = cam

    csv_name = f"{base_name}_motion_energy.csv"
    json_name = f"{base_name}_motion_energy_meta.json"
    csv_path = pdir / csv_name
    json_path = pdir / json_name

    # If already processed and not overwriting: load
    if owr == 0 and csv_path.exists() and json_path.exists():
        df_loaded = pd.read_csv(csv_path)
        with open(json_path, "r") as f:
            meta_loaded = json.load(f)

        # Record in animal
        a["analysis"]["motion_energy"][cam]["csv"] = str(csv_path)
        a["analysis"]["motion_energy"][cam]["meta_json"] = str(json_path)
        a["analysis"]["motion_energy"][cam]["base_name"] = base_name
        a["analysis"]["motion_energy"][cam]["loaded"] = True

        animal[an_idx] = a
        return animal, df_loaded, meta_loaded

    # Otherwise: save new (need df/roi/fps)
    if df is None or roi is None or fps is None:
        raise ValueError("To save new outputs, you must pass df, roi, and fps.")

    # Save CSV
    df.to_csv(csv_path, index=False)

    # Save meta JSON
    meta = {
        "video": a.get("video", {}).get("mp4", {}).get(cam, ""),  # whatever you stored earlier
        "roi": {
            "x": int(roi[0]),
            "y": int(roi[1]),
            "width": int(roi[2]),
            "height": int(roi[3]),
        },
        "fps": float(fps),
    }
    with open(json_path, "w") as f:
        json.dump(meta, f, indent=2)

    # Record in animal
    a["analysis"]["motion_energy"][cam]["csv"] = str(csv_path)
    a["analysis"]["motion_energy"][cam]["meta_json"] = str(json_path)
    a["analysis"]["motion_energy"][cam]["base_name"] = base_name
    a["analysis"]["motion_energy"][cam]["loaded"] = False

    animal[an_idx] = a
    return animal, df, meta

#%%

if __name__ == "__main__":
    # Example usage:
    video = animal[0]['video']['mp4']['face']  # <-- change this
    df, roi, fps = motion_energy_roi(video_path=video, roi=None, out_csv="face_motion_energy.csv")
    print("ROI used:", roi, "FPS:", fps)
    print(df.head())
    
    animal, df_out, meta_out = save_or_load_motion_energy(
        animal=animal,
        an_idx=0,
        cam="face",
        df=df,
        roi=roi,      # (x,y,w,h)
        fps=fps,
        base_name="face",  # or None to infer from mp4 name
        owr=0
        ) 

