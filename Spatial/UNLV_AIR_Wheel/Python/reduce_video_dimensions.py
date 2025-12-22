import cv2
import os


def resize_videos_opencv(data_list, scale=0.25):
    targets = ['face', 'paws']
    
    for entry in data_list:
        mp4_files = entry['video']['mp4']
        
        for key in targets:
            input_path = mp4_files.get(key)
            if not input_path or not os.path.exists(input_path):
                continue

            # 1. Setup Input
            cap = cv2.VideoCapture(input_path)
            fps = cap.get(cv2.CAP_PROP_FPS)
            width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
            height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
            
            # 2. Calculate New Dimensions (1/4th)
            new_width = int(width * scale)
            new_height = int(height * scale)
            
            # 3. Setup Output Path and Writer
            file_root, file_ext = os.path.splitext(input_path)
            output_path = f"{file_root}_reduced{file_ext}"
            
            # 'mp4v' is a common codec for MP4 files in OpenCV
            fourcc = cv2.VideoWriter_fourcc(*'mp4v') 
            out = cv2.VideoWriter(output_path, fourcc, fps, (new_width, new_height))

            print(f"Resizing {key}: {width}x{height} -> {new_width}x{new_height}")

            # 4. Process Frame by Frame
            while cap.isOpened():
                ret, frame = cap.read()
                if not ret:
                    break
                
                # Resize using INTER_AREA for best downscaling quality
                resized_frame = cv2.resize(frame, (new_width, new_height), interpolation=cv2.INTER_AREA)
                out.write(resized_frame)

            # 5. Clean up
            cap.release()
            out.release()
            print(f"Saved: {output_path}")

resize_videos_opencv(animal)
