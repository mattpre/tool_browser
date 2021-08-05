import imageio
import os

image_folder = '.'
images = [img for img in os.listdir(image_folder) if img.endswith(".jpeg")]

with imageio.get_writer('movie.gif', mode='I') as writer:
    for filename in images[::2]:
        image = imageio.imread(filename)
        writer.append_data(image)
