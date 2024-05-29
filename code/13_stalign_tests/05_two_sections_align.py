import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import pandas as pd
import torch
from STalign import STalign
from pathlib import Path
from pyhere import here
from PIL import Image
import re
import os
import json
import pickle

#   Set a non-interactive backend for writing figures to a file
import matplotlib
matplotlib.use('agg')

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
coldata_path = here(
    'processed-data', '13_stalign_tests', 'raw_coldata.csv'
)
out_dir = Path(here('processed-data', '13_stalign_tests', 'two_sections'))
plot_dir = Path(here('plots', '13_stalign_tests', 'two_sections'))
capture_areas = ['V12D07-078_B1', 'V12D07-078_D1']

plot_dir.mkdir(parents = False, exist_ok = True)

sample_info = pd.read_csv(sample_info_path, index_col = 0)

################################################################################
#   Read in landmarks as done in the tutorial
################################################################################

pointsIlist = np.load(
        out_dir / f'src_image_points.npy', allow_pickle=True
    ).tolist()
pointsJlist = np.load(
        out_dir / f'target_image_points.npy', allow_pickle=True
    ).tolist()

pointsI = []
pointsJ = []
for i in pointsIlist.keys():
    for j in range(len(pointsIlist[i])):
        pointsI.append([pointsIlist[i][j][1], pointsIlist[i][j][0]])

for i in pointsJlist.keys():
    for j in range(len(pointsJlist[i])):
        pointsJ.append([pointsJlist[i][j][1], pointsJlist[i][j][0]])

pointsI = np.array(pointsI)
pointsJ = np.array(pointsJ)

################################################################################
#   Read in tissue positions for source capture area and scale, if applicable
################################################################################

spaceranger_json_list = []
tissue_positions_list = []
for i in range(2):
    pattern = re.compile(r'^tissue_positions(_list|)\.csv$')
    this_dir = sample_info.loc[capture_areas[i], 'spaceranger_dir']
    tissue_path = [
            Path(os.path.join(this_dir, x)) for x in os.listdir(this_dir)
            if pattern.match(x)
        ][0]
    if '_list' in tissue_path.name: 
        tissue_positions = pd.read_csv(
            tissue_path,
            header = None,
            # Note the switch of x and y
            names = ["in_tissue", "row", "col", "y", "x"],
            index_col = 0
        )
    else:
        tissue_positions = pd.read_csv(
            tissue_path,
            skiprows = 1,
            # Note the switch of x and y
            names = ["in_tissue", "row", "col", "y", "x"],
            index_col = 0
        )
    #
    sf_json_path = os.path.join(this_dir, 'scalefactors_json.json')
    with open(sf_json_path, 'r') as f:
        spaceranger_json_list.append(json.load(f))
    #
    #   Scale spatial coordinates to high or lowres
    tissue_positions['x'] *= spaceranger_json_list[-1]['tissue_hires_scalef']
    tissue_positions['y'] *= spaceranger_json_list[-1]['tissue_hires_scalef']
    #
    tissue_positions.index += '_' + capture_areas[i]
    tissue_positions_list.append(tissue_positions)

################################################################################
#   Compute affine transformation and apply to source landmarks and spatial
#   coords
################################################################################

# set device for building tensors
if torch.cuda.is_available():
    torch.set_default_device('cuda:0')
else:
    torch.set_default_device('cpu')

# compute initial affine transformation from points
L,T = STalign.L_T_from_points(pointsI,pointsJ)

#   Doesn't work because in general this is an affine transform, not an affine
#   transform without scaling

#   Determine theta by using arccos, but negate it when sin(theta) is negative,
#   to account for the convention that arccos returns positive angles
theta_cos = np.arccos(L[0, 0]) * (2 * (L[1, 0] > 0) - 1)


A = STalign.to_A(torch.tensor(L),torch.tensor(T))