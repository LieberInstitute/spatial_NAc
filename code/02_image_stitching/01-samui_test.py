from pyhere import here
from pathlib import Path
import session_info

import numpy as np
import pandas as pd
import json
import os
import sys
from loopy.sample import Sample
import tifffile
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

#   Grab the slide number and 
this_slide, file_suffix = sys.argv[1:]
assert file_suffix in ('initial', 'adjusted')

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
estimate_path = here(
    'processed-data', '02_image_stitching', 'transformation_estimates_{}.csv'
)
out_dir = here('processed-data', '02_image_stitching', 'combined_{}_{}')
img_out_path = here('processed-data', '02_image_stitching', 'combined_{}_{}.tif')

#   55-micrometer diameter for Visium spot
SPOT_DIAMETER_M = 55e-6

#   Read in sample info and subset to slide of interest, with capture areas in
#   order
sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info = sample_info.loc[sample_info['Slide #'] == this_slide, :]
sample_info = sample_info.loc[sample_info.index.sort_values()]

#   Determine sample and run dependent paths (is this to display the initial
#   transformation from ImageJ or the adjusted transformation?)
estimate_path = Path(str(estimate_path).format(this_slide))
is_adjusted = file_suffix == 'adjusted'

out_dir = Path(str(out_dir).format(this_slide, file_suffix))
img_out_path = Path(str(img_out_path).format(this_slide, file_suffix))

################################################################################
#   Define the transformation matrix
################################################################################

#   Check if the adjusted estimates exist (i.e. if we already ran this script to
#   display in Samui, annotated ROIs, then ran '03-adjust_transform.py' to
#   refine the initial transformations from ImageJ). If not, use the initial
#   estimate from ImageJ
if is_adjusted:
    #   Read in the adjusted transformations
    estimates_df = pd.read_csv(estimate_path)
    estimates_df = estimates_df.loc[estimates_df['adjusted']]
else:
    estimates_df = sample_info.rename(
        {
            'initial_transform_x': 'x',
            'initial_transform_y': 'y',
            'initial_transform_theta': 'theta'
        },
        axis = 1
    )

#   x and y switch, so angles invert
estimates_df['theta'] *= -1
theta = np.array(estimates_df['theta'])

#   Array defining an affine transformation:
#   [x0', y0'] = trans[0] @ [x0, y0, 1]
trans = np.array(
    [
        [
            [
                np.cos(estimates_df['theta'].iloc[i]),
                -1 * np.sin(estimates_df['theta'].iloc[i]),
                estimates_df['x'].iloc[i]
            ],
            [
                np.sin(estimates_df['theta'].iloc[i]),
                np.cos(estimates_df['theta'].iloc[i]),
                estimates_df['y'].iloc[i]
            ]
        ]
        for i in range(estimates_df.shape[0])
    ],
    dtype = np.float64
)

if not is_adjusted:
    #   The initial translations are in high resolution. Scale to full
    #   resolution
    for i in range(sample_info.shape[0]):
        #   Grab high-res scale factor and scale translations accordingly
        json_path = os.path.join(
            sample_info['spaceranger_dir'].iloc[i], 'scalefactors_json.json'
        )
        with open(json_path, 'r') as f:
            spaceranger_json = json.load(f)
        
        trans[i, :, 2] /= spaceranger_json['tissue_hires_scalef']

#   Flip x and y to follow Samui conventions. Translations must be an integer
#   number of pixels
trans[:, :, 2] = np.flip(np.round(trans[:, :, 2]), axis = 1)

out_dir.mkdir(parents = True, exist_ok = True)

################################################################################
#   Determine image shapes and initialize the combined image
################################################################################

#   Loop through all samples and grab image dimensions so we can compute
#   the dimensions of the combined image. Doesn't load the images into memory
img_shapes = []
for i in range(sample_info.shape[0]):
    tif = tifffile.TiffFile(sample_info['raw_image_path'].iloc[i])
    img_shapes.append(tif.pages[0].shape)
    tif.close()

img_shapes = np.array(img_shapes)

#   Initialize the combined tiff
max0, max1 = np.max(trans[:, :, 2] + img_shapes[:, :2], axis = 0).astype(int)
combined_img = np.zeros((max0, max1, sample_info.shape[0]), dtype = np.uint8)

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
    img = Image.open(sample_info['raw_image_path'].iloc[i]).convert('L')

    #   Rotate about the top left corner of the image
    this_theta = 180 * theta[i] / np.pi # '.rotate' uses degrees, not radians
    img = np.array(img.rotate(this_theta, expand = False, center = (0, 0)))

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
    tissue_positions[['y', 'x']] = (
        trans[i] @ 
        np.array(tissue_positions.assign(ones = 1)[['y', 'x', 'ones']]).T
    ).T

    tissue_positions_list.append(tissue_positions)

    #   "Place this image" on the combined image, considering translations. Use
    #   separate channels for each image
    combined_img[
            int(trans[i, 0, 2]): int(trans[i, 0, 2] + img.shape[0]),
            int(trans[i, 1, 2]): int(trans[i, 1, 2] + img.shape[1]),
            i : (i + 1)
        ] += img.reshape((img.shape[0], img.shape[1], 1))

with tifffile.TiffWriter(img_out_path, bigtiff = True) as tiff:
    tiff.write(combined_img)

################################################################################
#   Use the Samui API to create the importable directory for this combined
#   "sample"
################################################################################

this_sample = Sample(name = out_dir.name, path = out_dir)

tissue_positions = pd.concat(tissue_positions_list)[['x', 'y']].astype(int)
this_sample.add_coords(
    tissue_positions, name = "coords", mPerPx = m_per_px, size = SPOT_DIAMETER_M
)

this_sample.add_image(
    tiff = img_out_path, channels = [x.split('_')[-1] for x in sample_info.index], scale = m_per_px
)

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
