from pyhere import here
from pathlib import Path
import session_info

import numpy as np
import pandas as pd
import json
from loopy.sample import Sample

img_paths = []

sample_info_path = here('processed-data', '02_image_stitching', 'sample_info_clean.csv')
out_dir = here('processed-data', '02_image_stitching', 'combined_Br8325')

spot_diameter_m = 55e-6 # 55-micrometer diameter for Visium spot

#   Image dimensions, which should match between images
EXPECTED_SHAPE = (100, 100, 3)

#   Array of translations to apply to each image
trans = np.array(
    [[1365, 18], [19, 52], [102, 742], [1236, 797]], dtype = np.float64
)

Path(out_dir).mkdir(parents = True, exist_ok = True)

sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info.query('Brain == "Br8325"')

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

#   For now, just read in one JSON file and assume all samples are roughly on
#   the same spatial scale (spots have the same diameter in pixels)
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']

tissue_positions_list = []

################################################################################
#   Read in and translate images and spot coordinates
################################################################################

for i in range(4):
    img = imgs[i]

    #   Define "sample_id"; "tissue_path" for this sample

    #   Read in the tissue positions file to get spatial coordinates. Index by
    #   barcode + sample ID
    tissue_positions = pd.read_csv(
        tissue_path,
        header = None,
        names = ["in_tissue", "row", "col", "y", "x"], # Note the switch of x and y
        index_col = 0
    )
    tissue_positions.index = tissue_positions.index + sample_id

    #   Translate coordinates as specified
    tissue_positions['x'] = tissue_positions['x'] + trans[i][0]
    tissue_positions['y'] = tissue_positions['y'] + trans[i][1]

    tissue_positions_list.append(tissue_positions)

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

################################################################################
#   Use the Samui API to create the importable directory for this combined
#   "sample"
################################################################################

this_sample = Sample(name = Path(out_dir).name, path = out_dir)

tissue_positions = pd.concat(tissue_positions_list)[['x', 'y']]
this_sample.add_coords(
    tissue_positions, name = "coords", mPerPx = m_per_px, size = spot_diameter_m
)

this_sample.write()

session_info.show()
