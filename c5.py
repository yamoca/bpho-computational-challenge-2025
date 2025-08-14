import numpy as np
from PIL import Image

object_image = Image.open("your_image.png")
object_array = np.array(object_image)

height, width, _ = object_array.shape # _ to throw away rgb channels


# empty array for virtual image - same size as original - will need to modify for future challenges
virtual_image_array = np.zeros_like(object_array)

# transform pixel coordinate (f(-x) essentially)
for y in range(height):
    for x in range(width):
        virtual_image_array[y, x] = object_array[y, width - 1 - x]

virtual_image = Image.fromarray(virtual_image_array)
virtual_image.save("virtual_image.png")

