from pyhere import here
from pathlib import Path
import session_info

import numpy as np
import pandas as pd
import json
import os
import sys

SPOT_DIAMETER_M = 55e-6

this_slide = sys.argv[1]

#   Inputs
roi_json_path = here(
    'processed-data', '02_image_stitching', 'rois_combined_{}.json'
)

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)

#   Outputs
estimate_path = here(
    'processed-data', '02_image_stitching', 'transformation_estimates_{}.csv'
)

pairwise_path = here(
    'processed-data', '02_image_stitching', 'pairwise_errors_{}.csv'
)

#   Capture areas we'll compare ROIs (in pairs) for in order to deduce
#   adjustments to the translations + rotations. The convention here will be to
#   rotate and translate the second element of the pair to match the first
array_pairs = {
    'V11U08-082': [('B1', 'C1'), ('A1', 'D1')],
    'V11U23-406': [('A1', 'B1')],
    'V11U23-404': []
}

################################################################################
#   Function definitions
################################################################################

#   Form a sort of polygon where each ROI is connected to exactly 2 other ROIs.
#   Compute an angle between the edges formed by points (rows) in 'a' vs. those
#   in 'b'. Estimate theta via a weighted average of each estimate for a pair of
#   edges; weights are the geometric mean of the lengths of each pair of edges.
#   Intuitively, this means a pair of long edges gives a more accurate estimate
#   of theta than a pair of short edges.
# 
#   Return theta in radians (the angle such that rotating points in 'b'
#   counter-clockwise by theta would result in points that line up angularly
#   with points in 'a'.
def get_theta(a, b, tol = 1e-8):
    weights = np.zeros(a.shape[0])
    theta = 0

    #   Loop through every pair of points such that a given ROI connects to a
    #   total of 2 edges used to compute theta. Note that the order here is
    #   arbitrary, and results can in theory change (very subtly) if we didn't
    #   traverse the pairs of points row by row (e.g., we could do rows {2, 1},
    #   {4, 2}, {3, 4}, {1, 3} instead of {2, 1}, {3, 2}, {4, 3}, {1, 4}).
    for i in range(a.shape[0]):
        #   Take 2 pairs of points. Determine edges, which are expected to form
        #   a very subtle X shape (with lines nearly parallel)
        a_edge = a[i, :] - a[(i + 1) % a.shape[0], :]
        b_edge = b[i, :] - b[(i + 1) % b.shape[0], :]

        #   Geometric mean of length of both edges. Used to compute theta and
        #   weight the estimate of theta 
        avg_mag = np.sqrt((a_edge @ a_edge) * (b_edge @ b_edge))
        
        weights[i] = avg_mag

        #   Use the definition of a dot product to solve for theta (the angle
        #   between the 2 edges). Note that by convention, theta is positive!
        #   a @ b = |a||b|cos(theta)
        theta_mag = np.arccos(a_edge @ b_edge / avg_mag)
        assert theta_mag >= 0, theta_mag

        #   If necessary, switch the sign of theta to represent the clockwise
        #   rotation necessary to line up with a when applied to b. Here,
        #   we simply check if applying the rotation in this way creates
        #   alignment (a dot product very close to 1), and negate theta
        #   otherwise
        rot = np.array(
            [
                [np.cos(theta_mag), -1 * np.sin(theta_mag)],
                [np.sin(theta_mag), np.cos(theta_mag)]
            ]
        )
        if abs((rot @ a_edge) @ b_edge / avg_mag - 1) < tol:
            theta_mag *= -1

        #   Add one term of the weighted sum
        theta += avg_mag * theta_mag
    
    #   Return the weighted sum of theta estimates by edge length
    return theta / np.sum(weights)

#   Return the root mean squared euclidean distance between ROIs in number of
#   spot diameters. Here 'a' and 'b' are each supposed to contain ROIs such
#   that row 1 in a is the point paired with row 1 of b
def get_rmse(a, b, M_PER_PX, SPOT_DIAMETER_M):
    #   Mean of the squares of the edge lengths between each pair of ROIs, then
    #   take the square root of the result
    rmse = np.sqrt(
        ((b - a)[:, 0] @ (b - a)[:, 0] + (b - a)[:, 1] @ (b - a)[:, 1]) /
        a.shape[0]
    )

    #   Convert from full-resolution pixels to number of spot diameters
    return rmse * M_PER_PX / SPOT_DIAMETER_M

#   Return the average euclidean distance between ROIs in 'a' and 'b', in
#   spot diameters
def get_avg_distance(a, b, M_PER_PX, SPOT_DIAMETER_M):
    err = np.mean(
        #   For each pair, calculate the length between ROI in a and b 
        np.sqrt(
            (b - a)[:, 0] ** 2 + (b - a)[:, 1] ** 2
        )
    )

    #   Convert from full-resolution pixels to number of spot diameters
    return err * M_PER_PX / SPOT_DIAMETER_M

#   Return a row-ordering of b that lines up with the rows of a. The algorithm
#   assumes the same ROI between a and b is closer than any other combination
#   of ROIs between a and b
def arrange_b(a, b):
    indices = np.zeros(a.shape[0], dtype = np.uint16)
    #   For each ROI in a, find the index of the closest ROI in b
    for i in range(a.shape[0]):
        indices[i] = np.argmin(np.sum((b - a[i, :]) * (b - a[i, :]), axis = 1))
    
    return b[indices, :].copy()

################################################################################
#   Process initial transformation estimates to be at full resolution
################################################################################

#   Read in sample info and subset to slide of interest, with capture areas in
#   order
sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info = sample_info.loc[sample_info['Slide #'] == this_slide, :]
sample_info = sample_info.loc[sample_info.index.sort_values()]

#   Adjusts paths for this slide
roi_json_path = Path(str(roi_json_path).format(this_slide))
estimate_path = Path(str(estimate_path).format(this_slide))
pairwise_path = Path(str(pairwise_path).format(this_slide))

#   Initial estimates of (by row): x translation, y translation, theta in
#   radians (counterclockwise relative to canvas)
init_estimate = np.array(
    sample_info.loc[
        :,
        [
            'initial_transform_x', 'initial_transform_y',
            'initial_transform_theta'
        ]
    ],
    dtype = np.float64
)

#   Loop through all samples and rescale translations to be at full resolutions
#   (the initial definition of 'init_estimate' uses high-resolution values).
for i in range(sample_info.shape[0]):
    #   Grab high-res scale factor and scale translations accordingly
    json_path = os.path.join(
        sample_info['spaceranger_dir'].iloc[i], 'scalefactors_json.json'
    )
    with open(json_path, 'r') as f:
        spaceranger_json = json.load(f)
    
    init_estimate[i, :2] /= spaceranger_json['tissue_hires_scalef']

#   Initialize a DataFrame of transformations for each sample. Include a
#   before and after (not adjusted vs. adjusted)
estimate_df = pd.DataFrame(init_estimate, columns = ['x', 'y', 'theta'])
estimate_df['adjusted'] = False
estimate_df['sample_id'] = sample_info.index
estimate_df = pd.concat([estimate_df, estimate_df.assign(adjusted = True)])

################################################################################
#   Read in ROIs from Samui and use to deduce adjustments to initial
#   transformation estimates
################################################################################

with open(roi_json_path, 'r') as f:
    roi_json = json.load(f)

#   Form a DataFrame with x and y coordinates (in full-resolution pixels) and
#   the associated label. x and y follow ImageJ coordinate system
roi_df = pd.DataFrame(
    {
        'x': pd.Series([x['coords'][0] for x in roi_json['rois']]) / roi_json['mPerPx'],
        'y': -1 * pd.Series([x['coords'][1] for x in roi_json['rois']]) / roi_json['mPerPx'],
        'label': [x['label'] for x in roi_json['rois']]
    }
)

#   Initialize a DataFrame of pairwise relationships between capture areas,
#   involving translations and rotations required for the second member of the
#   pair to best line up with the first, and associated error metrics for the
#   resulting alignment (before and after).
pairwise_df = pd.DataFrame(
    array_pairs[this_slide], columns = ['capture_area1', 'capture_area2']
)
other_columns = (
    'x', 'y', 'theta', 'before_err_rmse', 'before_err_avg', 'after_err_rmse',
    'after_err_avg'
)
for col in other_columns:
    pairwise_df[col] = None

#   Loop through pairs of capture areas to compute transformations and error
#   information, adding it to 'pairwise_df'
for pair in array_pairs[this_slide]:
    #   Extract ROIs of each label into separate arrays. Every ROI is paired,
    #   so a and b should be identical in dimension
    label1 = pair[0] + '_artifact_centroid'
    label2 = pair[1] + '_artifact_centroid'
    a = np.array(roi_df.loc[roi_df['label'] == label1, ['x', 'y']])
    b = np.array(roi_df.loc[roi_df['label'] == label2, ['x', 'y']])
    assert a.shape == b.shape

    #   Ensure row i in a refers to the same ROI as row i in b for all i
    b = arrange_b(a, b)

    #   Calculate initial error metrics
    init_err_avg = get_avg_distance(a, b, roi_json['mPerPx'], SPOT_DIAMETER_M)
    init_err_rmse = get_rmse(a, b, roi_json['mPerPx'], SPOT_DIAMETER_M)
    pairwise_df.loc[
            (pairwise_df['capture_area1'] == pair[0]) &
            (pairwise_df['capture_area2'] == pair[1]),
            ['before_err_avg', 'before_err_rmse']
        ] = [init_err_avg, init_err_rmse]

    #   Determine how much to rotate 'b' and apply that rotation. Center
    #   rotation at the top left of the image so that rotation parameters
    #   aren't dependent on the initial translation estimates
    theta = get_theta(a, b)
    rot = np.array(
        [
            [np.cos(theta), -1 * np.sin(theta)],
            [np.sin(theta), np.cos(theta)]
        ]
    )

    #   The translation from the origin to the top left of the capture area
    origin_to_corner = np.array(
        estimate_df.loc[
            (estimate_df['sample_id'] == f"{this_slide}_{pair[1]}") &
            ~estimate_df['adjusted'],
            ['x', 'y']
        ]
    ).reshape((1, 2))

    #   Rotate about the top-left corner of the capture area
    b = (rot @ (b - origin_to_corner).T).T + origin_to_corner

    #   Now translate 'b' to minimize error between ROIs of 'a' and 'b'
    trans = np.mean(a - b, axis = 0)
    b += trans

    #   Add error metrics after applying transformation, and the transformation
    #   itself (translation and theta)
    err_avg = get_avg_distance(a, b, roi_json['mPerPx'], SPOT_DIAMETER_M)
    err_rmse = get_rmse(a, b, roi_json['mPerPx'], SPOT_DIAMETER_M)
    pairwise_df.loc[
            (pairwise_df['capture_area1'] == pair[0]) &
            (pairwise_df['capture_area2'] == pair[1]),
            ['after_err_avg', 'after_err_rmse', 'x', 'y', 'theta']
        ] = [err_avg, err_rmse, trans[0], trans[1], theta]


#   Use information about adjustments to form a final table of (x, y)
#   translations and theta values, that can be used to prepare for reviewing
#   in Samui
for i in range(pairwise_df.shape[0]):
    estimate_df.loc[
        estimate_df['sample_id'] == f"{this_slide}_{pairwise_df['capture_area2'].iloc[i]}" &
        estimate_df['adjusted'],
        ['x', 'y', 'theta']
    ] += pairwise_df.loc[:, ['x', 'y', 'theta']].iloc[i]

#   Write estimates and pairwise errors to disk
estimate_df.to_csv(estimate_path, index = False)
pairwise_df.to_csv(pairwise_path, index = False)

session_info.show()
