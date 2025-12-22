# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 23:33:30 2025

@author: samso
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle

def plot_animal_motion(animal_list):
    for entry in animal_list:
        animal_id = entry['ID']
        date = entry['date']
        pdir = entry['pdir']
        targets = ['paws']
        
        # Create a figure with 3 subplots (one for each camera view)
        fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True)
        fig.suptitle(f"Motion Analysis: {animal_id} ({date})", fontsize=16)

        for i, key in enumerate(targets):
            # Construct the path to the _reduced_OF.csv file
            original_mp4 = entry['video']['mp4'].get(key)
            print(original_mp4)
            if not original_mp4:
                continue
                
            base_name = os.path.splitext(os.path.basename(original_mp4))[0]
            # Based on your description: video_reduced_OF.csv
            csv_path = os.path.join(pdir, f"{base_name}_reduced_OF.csv")
            
            if not os.path.exists(csv_path):
                # Fallback: check if the file was named without the 'reduced' tag but in the same folder
                csv_path = os.path.join(pdir, f"{base_name}_OF.csv")
            
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                
                # Calculations
                df['speed'] = np.sqrt(df['avg_u']**2 + df['avg_v']**2)
                df['time_sec'] = df['time_ms'] / 1000
                
                ax1 = axes[i]
                color_speed = 'tab:blue'
                ax1.set_ylabel(f'{key.capitalize()}\nSpeed (px/fr)', color=color_speed)
                ax1.plot(df['time_sec'], df['speed'], color=color_speed, alpha=0.8)
                ax1.tick_params(axis='y', labelcolor=color_speed)
                
                # Overlay Motion Energy
                ax2 = ax1.twinx()
                color_energy = 'tab:red'
                ax2.set_ylabel('Energy', color=color_energy)
                ax2.plot(df['time_sec'], df['motion_energy'], color=color_energy, alpha=0.4, linestyle='--')
                ax2.tick_params(axis='y', labelcolor=color_energy)
                
                ax1.grid(True, which='both', linestyle='--', alpha=0.5)
            else:
                axes[i].text(0.5, 0.5, f"CSV not found: {key}", ha='center')

        axes[-1].set_xlabel('Time (seconds)')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        # Save plot to the PData directory
        save_path = os.path.join(pdir, f"{animal_id}_{date}_motion_summary.png")
        plt.savefig(save_path, dpi=300)
        print(f"Saved summary plot to: {save_path}")
        plt.show()

with open("session_data.pkl", "rb") as f:
    loaded_data = pickle.load(f)
# Run the plotting function
plot_animal_motion(loaded_data["animal"])
