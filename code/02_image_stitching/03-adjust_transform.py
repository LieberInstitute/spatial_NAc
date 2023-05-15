from pyhere import here
from pathlib import Path
import session_info

import numpy as np
import pandas as pd
import json
import os

SPOT_DIAMETER_M = 55e-6

roi_json_path = here(
    'processed-data', '02_image_stitching', 'rois_combined_Br8325.json'
)

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)

#   Initial estimates of (by row): x translation, y translation, theta in
#   radians (counterclockwise relative to canvas)
init_estimate = np.array(
    [
        [1365, 18, 0],
        [19, 52, 0],
        [102, 742, 0],
        [1236, 797, 0]
    ]
)

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

################################################################################
#   Process initial transformation estimates to be at full resolution
################################################################################

#   Read in sample info and subset to the samples of interest
sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info = sample_info.loc[
    (sample_info['Brain'] == "Br8325") &
    (sample_info['Sample #'] != "32v_svb"),
    :
]

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

################################################################################
#   Read in ROIs from Samui and use to deduce adjustments to initial
#   transformation estimates
################################################################################

with open(roi_json_path, 'r') as f:
    roi_json = json.load(f)

#   Form a DataFrame with x and y coordinates (in full-resolution pixels) and
#   the associated label
roi_df = pd.DataFrame(
    {
        'x': pd.Series([x['coords'][0] for x in roi_json['rois']]) / roi_json['mPerPx'],
        'y': pd.Series([x['coords'][1] for x in roi_json['rois']]) / roi_json['mPerPx'],
        'label': [x['label'] for x in roi_json['rois']]
    }
)

#   Extract ROIs of each label into separate arrays. Every ROI is paired,
#   so a and b should be identical in dimension
label1 = 'B1_artifact_centroid'
label2 = 'C1_artifact_centroid'
a = np.array(roi_df.loc[roi_df['label'] == label1, ['x', 'y']])
b = np.array(roi_df.loc[roi_df['label'] == label2, ['x', 'y']])
assert a.shape == b.shape

print("Error of initial guess for aligment of B1 and C1:")
init_err = get_avg_distance(a, b, roi_json['mPerPx'], SPOT_DIAMETER_M)
print(f"{round(init_err, 2)} spot diameters")

print("Error after simply translating:")
err = get_avg_distance(a, b - np.mean(b - a, axis = 0), roi_json['mPerPx'], SPOT_DIAMETER_M)
print(f"{round(err, 2)} spot diameters")

#   Determine how much to rotate 'b' and apply that rotation
theta = get_theta(a, b)
rot = np.array(
    [
        [np.cos(theta), -1 * np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ]
)
b = (rot @ b.T).T

#   Now translate 'b' to minimize error between ROIs of 'a' and 'b'
b -= np.mean(b - a, axis = 0)

print("Error after rotating and translating:")
err = get_avg_distance(a, b, roi_json['mPerPx'], SPOT_DIAMETER_M)
print(f"{round(err, 2)} spot diameters")