from pyhere import here
from pathlib import Path
import session_info

import tifffile
import numpy as np


img_paths = []

out_path = here('processed-data', '02_image_stitching', 'test_image.tif')

#   Image dimensions, which should match between images
EXPECTED_SHAPE = (100, 100, 3)

#   Array of translations to apply to each image
trans = np.array(
    [[1365, 18], [19, 52], [102, 742], [1236, 797]], dtype = np.float64
)

Path(out_path).parent.mkdir(exist_ok = True)

#   Will be retrieved from spaceranger JSON
SCALEFACTOR = 20

#   Scale to full resolution (values were provided in high resolution)
trans = (trans / SCALEFACTOR).astype(np.uint32)

#   Initialize the combined tiff, for now using floats (we're averaging pixels
#   across several images before fixing the data type)
max0 = EXPECTED_SHAPE[0] + np.max(trans[:, 0])
max1 = EXPECTED_SHAPE[1] + np.max(trans[:, 1])
combined_img = np.zeros((max0, max1, 3), dtype = np.float64)
weights = np.zeros((max0, max1, 3), dtype = np.float64)

#   Define example RGB images just for testing
imgs = [np.zeros(EXPECTED_SHAPE, dtype = np.uint8) for x in range(4)]
imgs[0][:, :, 0] = 255
imgs[1][:, :, 1] = 255
imgs[2][:, :, 2] = 255
imgs[3][:, :, :2] = 127

for i in range(4):
    img = imgs[i]

    #   "Place this image" on the combined image, considering translations
    combined_img[
            trans[i][0]: trans[i][0] + EXPECTED_SHAPE[0],
            trans[i][1]: trans[i][1] + EXPECTED_SHAPE[1]
        ] += img

    #   Keep track of how many times each pixel was added to
    weights[
            trans[i][0]: trans[i][0] + EXPECTED_SHAPE[0],
            trans[i][1]: trans[i][1] + EXPECTED_SHAPE[1]
        ] += 1

#   Rescale image, noting that "untouched" pixels should be divided by 1, not 0
weights[weights == 0] = 1
combined_img = (combined_img / weights).astype(np.uint8)

with tifffile.TiffWriter(out_path) as tiff:
    tiff.write(combined_img)

session_info.show()
