# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 13:27:48 2023

@author: dboateng
"""

import imageio
from PIL import Image

def resize_gif(input_path, output_path, factor):
    reader = imageio.get_reader(input_path)
    
    # Get the first frame to get the size
    first_frame = reader.get_next_data()
    width, height = first_frame.shape[1], first_frame.shape[0]

    # Calculate the new size and speed
    new_width = int(width * factor)
    new_height = int(height * factor)

    # Try to get fps from the metadata, use a default value if not available
    try:
        fps = reader.get_meta_data()['fps']
    except KeyError:
        fps = 10  # Set a default value (adjust as needed)

    new_fps = fps / factor

    # Create a writer for the resized GIF
    writer = imageio.get_writer(output_path, fps=new_fps)

    # Resize and write each frame
    for frame in reader:
        resized_frame = Image.fromarray(frame).resize((new_width, new_height), Image.ANTIALIAS)
        writer.append_data(resized_frame)

    writer.close()




if __name__ == "__main__":
    input_gif_path = "MIO_278_climatologies_temp.gif"
    output_gif_path = "resized_output.gif"
    resize_factor = 0.5  # Adjust this factor to decrease size and increase speed

    resize_gif(input_gif_path, output_gif_path, resize_factor)
