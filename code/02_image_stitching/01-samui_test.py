from pyhere import here
from pathlib import Path
import session_info

import numpy as np
import pandas as pd
import json
import os
from loopy.sample import Sample
import tifffile

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
out_dir = here('processed-data', '02_image_stitching', 'combined_Br8325')
img_out_path = here('processed-data', '02_image_stitching', 'combined_Br8325.tif')

#   55-micrometer diameter for Visium spot
SPOT_DIAMETER_M = 55e-6

#   Array defining an affine transformation:
#   [x0', y0'] = trans[0] @ [x0, y0, 1]
trans = np.array(
    [
        [[1, 0, 1365], [0, 1, 18]],
        [[1, 0, 19], [0, 1, 52]],
        [[1, 0, 102], [0, 1, 742]],
        [[1, 0, 1236], [0, 1, 797]]
    ],
    dtype = np.float64
)

Path(out_dir).mkdir(parents = True, exist_ok = True)

#   Read in sample info and subset to the samples of interest
sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info = sample_info.loc[
    (sample_info['Brain'] == "Br8325") &
    (sample_info['Sample #'] != "32v_svb"),
    :
]

#   Loop through all samples and rescale translations to be at full resolutions
#   (the initial definition of 'trans' uses high-resolution values). Grab
#   image dimensions so we can compute the dimensions of the combined image
img_shapes = []
for i in range(sample_info.shape[0]):
    #   Grab high-res scale factor and scale translations accordingly
    json_path = os.path.join(
        sample_info['spaceranger_dir'].iloc[i], 'scalefactors_json.json'
    )
    with open(json_path, 'r') as f:
        spaceranger_json = json.load(f)
    
    trans[i, :, 2] /= spaceranger_json['tissue_hires_scalef']

    #   Grab image dimensions without loading into memory
    tif = tifffile.TiffFile(sample_info['raw_image_path'].iloc[i])
    img_shapes.append(tif.pages[0].shape)
    tif.close()

img_shapes = np.array(img_shapes)

#   Translations must be an integer number of pixels
trans[:, :, 2] = np.round(trans[:, :, 2])

#   Initialize the combined tiff, for now using floats (we're averaging pixels
#   across several images before fixing the data type). 16 bits to save memory
max0, max1 = np.max(trans[:, :, 2] + img_shapes[:, :2], axis = 0).astype(int)
combined_img = np.zeros((max0, max1, 3), dtype = np.float16)
weights = np.zeros((max0, max1, 1), dtype = np.float16)

#   Define example RGB images just for testing
# imgs = [np.zeros(EXPECTED_SHAPE, dtype = np.uint8) for x in range(4)]
# imgs[0][:, :, 0] = 255
# imgs[1][:, :, 1] = 255
# imgs[2][:, :, 2] = 255
# imgs[3][:, :, :2] = 127

#   For now, just read in one JSON file and assume all samples are roughly on
#   the same spatial scale (spots have the same diameter in pixels)
json_path = os.path.join(
    sample_info['spaceranger_dir'].iloc[0], 'scalefactors_json.json'
)
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = SPOT_DIAMETER_M / spaceranger_json['spot_diameter_fullres']

tissue_positions_list = []

################################################################################
#   Read in and translate images and spot coordinates
################################################################################

for i in range(sample_info.shape[0]):
    img = tifffile.imread(sample_info['raw_image_path'].iloc[i])

    #   Read in the tissue positions file to get spatial coordinates. Index by
    #   barcode + sample ID, and subset to only spots within tissue
    tissue_path = os.path.join(
        sample_info['spaceranger_dir'].iloc[i], 'tissue_positions_list.csv'
    )
    tissue_positions = pd.read_csv(
            tissue_path,
            header = None,
            # Note the switch of x and y
            names = ["in_tissue", "row", "col", "y", "x"],
            index_col = 0
        ).query('in_tissue == 1')
    tissue_positions.index = tissue_positions.index + '_' + sample_info.index[i]

    #   Apply affine transform of coordinates
    tissue_positions[['x', 'y']] = (
        trans[i] @ 
        np.array(tissue_positions.assign(ones = 1)[['x', 'y', 'ones']]).T
    ).T

    tissue_positions_list.append(tissue_positions)

    #   "Place this image" on the combined image, considering translations
    combined_img[
            int(trans[i, 0, 2]): int(trans[i, 0, 2] + img.shape[0]),
            int(trans[i, 1, 2]): int(trans[i, 1, 2] + img.shape[1]),
            :
        ] += img

    #   Keep track of how many times each pixel was added to
    weights[
            int(trans[i, 0, 2]): int(trans[i, 0, 2] + img.shape[0]),
            int(trans[i, 1, 2]): int(trans[i, 1, 2] + img.shape[1])
        ] += 1

#   Rescale image, noting that "untouched" pixels should be divided by 1, not 0.
#   Note that all channels (RGB) have the same weights
weights[weights == 0] = 1
combined_img = (combined_img / weights).astype(np.uint8)

with tifffile.TiffWriter(img_out_path, bigtiff = True) as tiff:
    tiff.write(combined_img)

################################################################################
#   Use the Samui API to create the importable directory for this combined
#   "sample"
################################################################################

this_sample = Sample(name = Path(out_dir).name, path = out_dir)

tissue_positions = pd.concat(tissue_positions_list)[['x', 'y']].astype(int)
this_sample.add_coords(
    tissue_positions, name = "coords", mPerPx = m_per_px, size = SPOT_DIAMETER_M
)

this_sample.add_image(tiff = img_out_path, channels = 'rgb', scale = m_per_px)

sample_df = pd.DataFrame(
    {
        'sample_id': [
            x.split('_')[-2] + '_' + x.split('_')[-1]
            for x in tissue_positions.index
        ]
    },
    index = tissue_positions.index
)
sample_df['sample_id'] = sample_df['sample_id'].astype('category')

this_sample.add_csv_feature(
    sample_df, name = "Sample Info", coordName = "coords"
)

this_sample.write()

session_info.show()
