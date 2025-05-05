import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl
import os
from PIL import Image
import numpy as np

from PIL import Image, ImageChops
mpl.use("Agg")  # For non-interactive backend

# Input/output directories
base_dir = "/work/cmcc/ag15419/basin_modes"
amp_dir = os.path.join(base_dir, "mode_plots_amp")
pow_dir = os.path.join(base_dir, "mode_plots_pow")
out_dir = os.path.join(base_dir, "mode_plots_all")

#############################

def crop_white_margins(img):
    """Crop white margins from top and bottom of the image."""
    bg = Image.new(img.mode, img.size, (255, 255, 255))
    diff = ImageChops.difference(img, bg)
    bbox = diff.getbbox()
    if bbox:
        # Crop only vertically (keep width)
        left, upper, right, lower = bbox
        return img.crop((0, upper, img.width, lower))
    return img  # return original if no bbox

# Make sure the output directory exists
os.makedirs(out_dir, exist_ok=True)

# Loop over all box indices from 0 to 22
for idx in range(23):
    # Get all files related to this index from the amplitude directory
    all_files = os.listdir(amp_dir)
    files_idx = [f for f in all_files if f"_amp_{idx}_" in f]

    for f in files_idx:
        # Extract the period from the filename
        parts = f.split("_")
        period_with_h = parts[-1]  # e.g., '12h.png'
        period = period_with_h.replace("h.png", "")  # e.g., '12'

        # Construct the expected file paths
        file_ul = os.path.join(amp_dir, f"mode_flag_amp_{idx}_{period}h.png")  # upper left
        file_ur = os.path.join(amp_dir, f"mode_amp_{idx}_{period}h.png")       # upper right
        file_ll = os.path.join(amp_dir, f"mode_ampval_{idx}_{period}h.png")    # lower left
        file_lr = os.path.join(pow_dir, f"mode_powval_{idx}_{period}h.png")    # lower right
        out_file = os.path.join(out_dir, f"modes_all_{idx}_{period}h.png")

        # Check if all required files exist
        if not all(os.path.exists(p) for p in [file_ul, file_ur, file_ll, file_lr]):
            print(f"Some files are missing for IDX={idx}, PERIOD={period}h. Skipping.")
            continue

        # Create a 2x2 subplot figure
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))

        img_ul = crop_white_margins(Image.open(file_ul))
        img_ur = crop_white_margins(Image.open(file_ur))
        img_ll = crop_white_margins(Image.open(file_ll))
        img_lr = crop_white_margins(Image.open(file_lr))

        axs[0, 0].imshow(np.array(img_ul))
        axs[0, 0].axis('off')

        axs[0, 1].imshow(np.array(img_ur))
        axs[0, 1].axis('off')

        axs[1, 0].imshow(np.array(img_ll))
        axs[1, 0].axis('off')

        axs[1, 1].imshow(np.array(img_lr))
        axs[1, 1].axis('off')

        # Add the main title
        #fig.suptitle(f"Mode with period {period} h", fontsize=16)

        # Adjust layout and save the figure
        plt.tight_layout(rect=[0, 0, 1, 0.95])  # Leave space for the title
        #plt.subplots_adjust(wspace=0.05, hspace=0.05)
        plt.savefig(out_file)
        plt.close()

        print(f"Saved: {out_file}")
