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

#   Parse command-line arguments, ultimately determining:
#       1. whether to produce the combined image of the initial or adjusted
#          estimates
#       2. the sample IDs (slide + arrays) to include in the combined image
this_slide, arrays, file_suffix = sys.argv[1:]
sample_ids = [f'{this_slide}_{x}' for x in arrays.split('_')]
assert file_suffix in ('initial', 'adjusted')

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
estimate_path = Path(
    here(
        'processed-data', '02_image_stitching',
        f'transformation_estimates_{this_slide}_{arrays}.csv'
    )
)
out_dir = Path(
    here(
        'processed-data', '02_image_stitching',
        f'combined_{this_slide}_{arrays}_{file_suffix}2'
    )
)
img_out_path = Path(
    here(
        'processed-data', '02_image_stitching',
        f'combined_{this_slide}_{arrays}_{file_suffix}.tif'
    )
)
tissue_out_path = Path(
    here(
        'processed-data', '02_image_stitching',
        f'tissue_positions_{this_slide}_{arrays}.csv'
    )
)

#   55-micrometer diameter for Visium spot
SPOT_DIAMETER_M = 55e-6

#   Read in sample info and subset to samples of interest
sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info = sample_info.loc[sample_ids, :]

################################################################################
#   Functions
################################################################################

#   Given the size (shape) of a PIL Image 'img', simulate a rotation
#   counter-clockwise by theta (radians) about the top left of the image (0, 0),
#   and return the translation (shape (2,) numpy array) needed to add to the
#   indices of the numpy array given by np.array(img.rotate([theta in degrees]))
#   to correctly apply the rotation about the top left of the image
def adjust_trans(img_shape, theta):
    #   (Negative) angle the top edge of the image makes with the diagonal,
    #   and the length of that diagonal. We'll need these to calculate
    #   where the bottom-right corner of the image moves as a function of
    #   theta
    theta_r = np.arctan(img_shape[1] / img_shape[0])
    r = np.linalg.norm(img_shape)

    #   Consider the image as a rectangle, and rotate that rectangle about the
    #   top-left corner (which is considered the origin). Track where all 4
    #   corners of the rectangle go, since that determines where the originally
    #   top-left corner ends up relative to the expanded-dimension numpy array's
    #   top-left corner. 'adjust' is simply the vector difference between those
    #   corners
    adjust = np.array(
        [
            -1 * max(
                0,
                img_shape[0] * np.sin(theta),
                img_shape[1] * np.sin(theta - np.pi / 2),
                r * np.sin(theta - theta_r)
            ),
            min(
                0,
                img_shape[0] * np.cos(theta),
                img_shape[1] * np.cos(theta - np.pi / 2),
                r * np.cos(theta - theta_r)
            )
        ]
    )

    return adjust

#   Unit tests of 'adjust_trans'
def is_equal(a, b, tol):
    return np.sum((a - b) @ (a - b)) < tol

tol = 1e-10
assert is_equal(adjust_trans((100, 150), 2 * np.pi), np.zeros(2), tol)
assert is_equal(adjust_trans((87, 13), np.pi), np.array([-13, -87]), tol)
assert is_equal(
    adjust_trans((4, 4), np.pi / 4),
    np.array([-2 * np.sqrt(2), 0]),
    tol
)

#   Given a numpy array of shape (n, 2) representing the shapes of n images
#   (stored as numpy arrays), and 'trans', the shape (n, 2, 3) numpy array
#   used to represent affine transformations of n images in this script,
#   return a shape (n, 2) numpy array representing the shapes of n images
#   after applying .rotate(expand = True) to the PIL Image version of each
#   image
def rotate_shapes(shapes, trans):
    #   For each image shape, form a list of coordinates for the corners of the
    #   rectangle representing the image
    start_points = np.array(
        [
            [(0, 0), (0, shapes[i, 1]), (shapes[i, 0], 0), shapes[i]]
            for i in range(shapes.shape[0])
        ],
        dtype = np.float64
    )

    #   Apply the appropriate rotation to each of those points for each image
    end_points = np.array(
        [
            (trans[i, :, :2] @ start_points[i, :, :].T).T
            for i in range(shapes.shape[0])
        ],
        dtype = np.float64
    )

    #   Return the image dimensions of the rotated image by forming a rectangle
    #   that fits all 4 corners of each image
    end_shape = np.ceil(
        np.max(end_points, axis = 1) - np.min(end_points, axis = 1)
    ).astype(int)
    return end_shape

################################################################################
#   Define the transformation matrix
################################################################################

#   Check if the adjusted estimates exist (i.e. if we already ran this script to
#   display in Samui, annotated ROIs, then ran '03-adjust_transform.py' to
#   refine the initial transformations from ImageJ). If not, use the initial
#   estimate from ImageJ
if file_suffix == 'adjusted':
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
rotated_shapes = rotate_shapes(img_shapes[:, :2], trans)

trans_img = np.array(
    [
        trans[i, :, 2] +
        adjust_trans(np.flip(img_shapes[i, :2]), theta[i])
        for i in range(sample_info.shape[0])
    ],
    dtype = np.float64
)
trans_img = np.round(trans_img).astype(int)

#   Equally translate images and spot coords such that the minimum pixel
#   occupied by any image is at (0, 0) (also potentially saving memory)
trans[:, :, 2] -= np.min(trans_img, axis = 0)
trans_img -= np.min(trans_img, axis = 0)

#   Initialize the combined tiff. Determine the boundaries by computing the
#   maximal coordinates in each dimension of each rotated and translated image
max0, max1 = np.max(rotated_shapes + trans_img, axis = 0)
combined_img = np.zeros(
    (max0 + 1, max1 + 1, sample_info.shape[0]), dtype = np.uint8
)

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
    img_pil = Image.open(sample_info['raw_image_path'].iloc[i]).convert('L')

    #   Rotate about the top left corner of the image
    theta_deg = 180 * theta[i] / np.pi # '.rotate' uses degrees, not radians
    img = np.array(img_pil.rotate(theta_deg, expand = True), dtype = np.uint8)

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
            trans_img[i, 0]: trans_img[i, 0] + img.shape[0],
            trans_img[i, 1]: trans_img[i, 1] + img.shape[1],
            i : (i + 1)
        ] += img.reshape((img.shape[0], img.shape[1], 1))

with tifffile.TiffWriter(img_out_path, bigtiff = True) as tiff:
    tiff.write(combined_img)

#   Also export tissue positions with SpatialExperiment-friendly colnames
tissue_positions_r = tissue_positions.rename(
    {
        'row': 'array_row', 'col': 'array_col', 'y': 'pxl_row_in_fullres',
        'x': 'pxl_row_in_fullres'
    },
    axis = 1
)
tissue_positions_r.index.name = 'key'

tissue_positions_r.to_csv(tissue_out_path, index = True)

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
    tiff = img_out_path,
    channels = [x.split('_')[-1] for x in sample_info.index],
    scale = m_per_px
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
