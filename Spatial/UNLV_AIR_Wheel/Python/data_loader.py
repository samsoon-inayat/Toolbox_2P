# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict, Any
from process_h264 import process_h264
import pickle


@dataclass
class MouseDirs:
    main_dir: Path
    rdata_dir: Path
    pdata_dir: Path
    adata_dir: Path


def build_dirs(main_dir: str) -> MouseDirs:
    main = Path(main_dir)
    return MouseDirs(
        main_dir=main,
        rdata_dir=main / "RData",
        pdata_dir=main / "PData",
        adata_dir=main / "AData",
    )


def get_exp_info(mD: MouseDirs, animal_list: List[str], date_list: List[str]) -> List[Dict[str, Any]]:
    if len(animal_list) != len(date_list):
        raise ValueError("animal_list and date_list must be the same length.")

    animals: List[Dict[str, Any]] = []

    for animal_id, date_str in zip(animal_list, date_list):
        rdir = mD.rdata_dir / animal_id / date_str
        pdir = mD.pdata_dir / animal_id / date_str
        adir = mD.adata_dir / animal_id / date_str

        # Make output dirs if missing (like MATLAB mkdir)
        pdir.mkdir(parents=True, exist_ok=True)
        adir.mkdir(parents=True, exist_ok=True)

        # Initialize fields (like MATLAB)
        info = {
            "ID": animal_id,
            "date": date_str,
            "rdir": str(rdir),
            "pdir": str(pdir),
            "adir": str(adir),
            "video": {
                "h264": {
                    "face": "",
                    "pupil": "",
                    "paws": "",
                }
            },
            "mat": "",
        }

        # Scan files in rdir (equivalent to dir(fullfile(rdir,'*')))
        if rdir.exists():
            for f in rdir.iterdir():
                if f.is_dir():
                    continue

                name_lower = f.name.lower()

                # ---- H264 CAMERAS ----
                if name_lower.endswith(".h264"):
                    if name_lower.startswith("face"):
                        info["video"]["h264"]["face"] = str(f)
                    elif name_lower.startswith("pupi"):
                        info["video"]["h264"]["pupil"] = str(f)
                    elif name_lower.startswith("video"):
                        info["video"]["h264"]["paws"] = str(f)

                # ---- MAT FILE ----
                elif name_lower.endswith(".mat"):
                    # MATLAB code overwrote this each time; this does the same.
                    # If you want "recording.mat" specifically, see note below.
                    info["mat"] = str(f)
        else:
            # Optional: warn or raise
            # raise FileNotFoundError(f"RData directory not found: {rdir}")
            pass

        animals.append(info)

    return animals

# %%
if __name__ == "__main__":
    main_dir = r"E:\GoogleDrive\InayatSamsoon\UNLV\AIR_Wheel_Methods"
    mD = build_dirs(main_dir)

    # Example lists (same as your MATLAB)
    animal_list = ["NML_GC_01"]
    date_list = ["2025_12_16"]

    animal = get_exp_info(mD, animal_list, date_list)
    print("Done")
    
    animal = process_h264(animal, owr=0)
    data_to_save = {
    "animal": animal,
    "main_dir": main_dir,
    "animal_list": animal_list,
    "date_list": date_list
    }
    # Save the dictionary as one file
    with open("session_data.pkl", "wb") as f:
        pickle.dump(data_to_save, f)

    # with open("animal.pkl", "wb") as f: # Must use 'wb' (write binary)
    #     pickle.dump(animal, f)
    # Optional: show result
    # from pprint import pprint
    # pprint(animal)
