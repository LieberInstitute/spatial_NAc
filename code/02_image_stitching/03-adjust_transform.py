from pyhere import here
from pathlib import Path
import session_info

import numpy as np
import pandas as pd
import json
import os
from loopy.sample import Sample
import tifffile
from PIL import Image

SPOT_DIAMETER_M = 55e-6

roi_json_path = here(
    'processed-data', '02_image_stitching', 'rois_combined_Br8325.json'
)

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

#   Form a sort of polygon where each ROI is connected to exactly 2 other ROIs.
#   Compute an angle between the edges formed by points with label1 vs. those of
#   label2. Estimate theta via a weighted average of each estimate for a pair of
#   edges; weights are the geometric mean of the lengths of each pair of edges.
#   Intuitively, this means a pair of long edges gives a more accurate estimate
#   of theta than a pair of short edges.
# 
#   Return (translation, theta) where theta is in radians (the angle such
#   that rotating points of label2 counter-clockwise by theta would result in
#   points that line up angularly with points of label1). 'translation' is the
#   vector to be added to points of label2 to minimize the resulting average
#   distance between pairs of ROIs.
def get_transform(roi_df, label1, label2):
    #   Extract ROIs of each label into separate arrays. Every ROI is paired,
    #   so a and b should be identical in dimension
    a = np.array(roi_df.loc[roi_df['label'] == label1, ['x', 'y']])
    b = np.array(roi_df.loc[roi_df['label'] == label2, ['x', 'y']])
    assert a.shape == b.shape

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
        #   between the 2 edges):
        #   a @ b = |a||b|cos(theta)
        theta += avg_mag * np.arccos(a_edge @ b_edge / avg_mag)
    
    #   Return (translation, angle). Here the angle (theta) is the weighted sum
    #   of theta estimates by edge length
    return (np.mean(a - b, axis = 0), theta / np.sum(weights))

#   Return the root mean squared euclidean distance between ROIs in number of
#   spot diameters. Here 'a' and 'b' are each supposed to contain ROIs such
#   that row 1 in a is the point paired with row 1 of b
def get_error(a, b, M_PER_PX, SPOT_DIAMETER_M):
    #   Sum of the squares of the edge lengths between each pair of ROIs, then
    #   take the square root of the result
    rmse = np.sqrt(
        (b - a)[:, 0] @ (b - a)[:, 0] + (b - a)[:, 1] @ (b - a)[:, 1]
    )

    #   Convert from full-resolution pixels to number of spot diameters
    return rmse * M_PER_PX / SPOT_DIAMETER_M

print("Error of initial guess for aligment of B1 and C1:")
err = get_error(
    np.array(roi_df.loc[roi_df['label'] == 'B1_artifact_centroid', ['x', 'y']]),
    np.array(roi_df.loc[roi_df['label'] == 'C1_artifact_centroid', ['x', 'y']]),
    roi_json['mPerPx'], SPOT_DIAMETER_M
)
print(f"{round(err, 2)} spot diemeters")