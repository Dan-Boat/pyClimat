# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 14:43:34 2023

@author: dboateng
"""

from PIL import Image
import os
import glob

path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/scratch/OCO2_OCO3/OCO2_OCO3_grid_1"
#path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/scratch/OCO2"
# Directory containing PNG files
output_gif = "output.gif"

# Get a list of all PNG files in the directory
png_files = glob.glob(os.path.join(path_to_plots, "*.png"))

# Sort the list of PNG files by filename, if necessary
png_files.sort()

# Create a list to store the image objects
images = []

# Loop through the PNG files and append them to the images list
for png_file in png_files:
    img = Image.open(png_file)
    images.append(img)

# Save the list of images as a GIF
images[0].save(os.path.join(path_to_plots,output_gif), save_all=True, append_images=images[1:], loop=0, duration=100)