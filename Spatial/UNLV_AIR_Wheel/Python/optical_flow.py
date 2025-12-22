import cv2
import os
import json
import csv
import numpy as np
import pickle
from tqdm import tqdm

def run_comprehensive_motion_analysis(data_list):
    # Fix for OpenCV 4.x/5.x
    try:
        tvl1 = cv2.optflow.DualTVL1OpticalFlow_create()
    except AttributeError:
        print("Error: 'cv2.optflow' not found. Run: pip install opencv-contrib-python")
        return

    targets = ['face', 'paws', 'pupil']
    analysis_queue = []

    # --- PHASE 1: SEQUENTIAL ROI SELECTION ---
    print("\n--- PHASE 1: SELECT ROIs FOR ALL VIDEOS ---")
    for entry in data_list:
        mp4_files = entry['video']['mp4']
        for key in targets:
            original_path = mp4_files.get(key)
            if not original_path: continue

            # Determine path (reduced vs original)
            if key in ['face', 'paws']:
                root, ext = os.path.splitext(original_path)
                video_path = f"{root}_reduced{ext}"
                if not os.path.exists(video_path): video_path = original_path
            else:
                video_path = original_path

            if not os.path.exists(video_path): continue

            cap = cv2.VideoCapture(video_path)
            ret, first_frame = cap.read()
            cap.release()
            if not ret: continue

            # Select ROI
            win_name = f"SELECT ROI: {key} ({os.path.basename(video_path)})"
            roi = cv2.selectROI(win_name, first_frame, fromCenter=False)
            cv2.destroyWindow(win_name)
            
            x, y, w, h = [int(v) for v in roi]
            if w > 0 and h > 0:
                # Save ROI info to JSON immediately
                root, _ = os.path.splitext(video_path)
                roi_json_path = f"{root}_roi.json"
                roi_info = {"x": x, "y": y, "w": w, "h": h, "source": video_path}
                with open(roi_json_path, 'w') as f:
                    json.dump(roi_info, f, indent=4)
                
                # Add to processing queue
                analysis_queue.append({
                    'key': key,
                    'path': video_path,
                    'roi': (x, y, w, h),
                    'output': f"{root}_OF.csv"
                })

    # --- PHASE 2: UNINTERRUPTED PROCESSING ---
    print(f"\n--- PHASE 2: PROCESSING {len(analysis_queue)} ANALYSES ---")
    for task in analysis_queue:
        video_path = task['path']
        x, y, w, h = task['roi']
        
        cap = cv2.VideoCapture(video_path)
        total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
        
        # Read first frame for initialization
        ret, prev_frame = cap.read()
        if not ret: 
            cap.release()
            continue
            
        prev_gray = cv2.cvtColor(prev_frame[y:y+h, x:x+w], cv2.COLOR_BGR2GRAY)
        
        with open(task['output'], mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(['frame_number', 'time_ms', 'avg_u', 'avg_v', 'motion_energy'])

            # Progress bar for the specific video
            pbar = tqdm(total=total_frames, desc=f"Analyzing {task['key']}", unit="fr", leave=True)
            pbar.update(1)
            
            frame_count = 1
            while True:
                ret, frame = cap.read()
                if not ret: break
                
                timestamp = cap.get(cv2.CAP_PROP_POS_MSEC)
                curr_gray = cv2.cvtColor(frame[y:y+h, x:x+w], cv2.COLOR_BGR2GRAY)
                
                # Optical Flow Calculation
                flow = tvl1.calc(prev_gray, curr_gray, None)
                
                # Metrics
                avg_u = np.mean(flow[..., 0])
                avg_v = np.mean(flow[..., 1])
                motion_energy = np.mean(cv2.absdiff(curr_gray, prev_gray))
                
                writer.writerow([frame_count, round(timestamp, 2), avg_u, avg_v, motion_energy])
                
                prev_gray = curr_gray
                frame_count += 1
                pbar.update(1)
            
            pbar.close()
        cap.release()

# Run
with open("session_data.pkl", "rb") as f:
    loaded_data = pickle.load(f)
run_comprehensive_motion_analysis(loaded_data["animal"])
