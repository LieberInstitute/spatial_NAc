import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import pandas as pd
import torch
import plotly
import requests
from STalign import STalign
from pathlib import Path
from pyhere import here
from PIL import Image
import re
import os
import json

#   Set a non-interactive backend for writing figures to a file
import matplotlib
matplotlib.use('agg')

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
out_dir = Path(here('processed-data', '13_stalign_tests'))
plot_dir = Path(here('plots', '13_stalign_tests'))
capture_areas = ['V12D07-078_B1', 'V12D07-078_D1']

out_dir.mkdir(parents = False, exist_ok = True)
plot_dir.mkdir(parents = False, exist_ok = True)

################################################################################
#   Read in source and target images and normalize
################################################################################

sample_info = pd.read_csv(sample_info_path, index_col = 0)
img_paths = [
    str(here(x, 'tissue_hires_image.png'))
    for x in sample_info.loc[capture_areas, 'spaceranger_dir']
]
img_arrs = [
    STalign.normalize(np.array(Image.open(x), dtype = np.float64) / 256)
    for x in img_paths
]

I = np.moveaxis(img_arrs[0], 2, 0)
J = np.moveaxis(img_arrs[1], 2, 0)
XI = np.array(range(I.shape[2])) * 1.
YI = np.array(range(I.shape[1])) * 1.
XJ = np.array(range(J.shape[2])) * 1.
YJ = np.array(range(J.shape[1])) * 1.
extentJ = STalign.extent_from_x((YJ,XJ))
extentI = STalign.extent_from_x((YI,XI))

#   Save images in format expected by 'point_annotator.py'
np.savez(out_dir / 'src_image', x = XI, y = YI, I = I)
np.savez(out_dir / 'target_image', x = XJ, y = YJ, I = J)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Then run annotate_landmarks.sh interactively
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

################################################################################
#   Read in landmarks as done in the tutorial
################################################################################
 
pointsIlist = np.load(out_dir / 'src_image_points.npy', allow_pickle=True).tolist()
pointsJlist = np.load(out_dir / 'target_image_points.npy', allow_pickle=True).tolist()

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
#   Plot capture areas with landmarks overlayed
################################################################################

fig,ax = plt.subplots(1,2)
ax[0].imshow((I.transpose(1,2,0).squeeze()), extent=extentI)
ax[1].imshow((J.transpose(1,2,0).squeeze()), extent=extentJ)

trans_offset_0 = mtransforms.offset_copy(
    ax[0].transData, fig=fig, x=0.05, y=-0.05, units='inches'
)
trans_offset_1 = mtransforms.offset_copy(
    ax[1].transData, fig=fig, x=0.05, y=-0.05, units='inches'
)

ax[0].scatter(pointsI[:,1],pointsI[:,0], c='red', s=10)
ax[1].scatter(pointsJ[:,1],pointsJ[:,0], c='red', s=10)

for i in pointsIlist.keys():
    for j in range(len(pointsIlist[i])):
        ax[0].text(pointsIlist[i][j][0], pointsIlist[i][j][1],f'{i}{j}', c='red', transform=trans_offset_0, fontsize= 8)

for i in pointsJlist.keys():
    for j in range(len(pointsJlist[i])):
        ax[1].text(pointsJlist[i][j][0], pointsJlist[i][j][1],f'{i}{j}', c='red', transform=trans_offset_1, fontsize= 8)

ax[0].set_title(capture_areas[0], fontsize=15)
ax[1].set_title(capture_areas[1], fontsize=15)

plt.savefig(str(plot_dir / '_'.join(capture_areas)) + '_with_landmarks.png')
plt.clf()

################################################################################
#   Read in tissue positions for source capture area and scale to highres
################################################################################

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
        spaceranger_json = json.load(f)
    #
    #   Scale spatial coordinates to highres (like the images we're using here)
    tissue_positions['x'] *= spaceranger_json['tissue_hires_scalef']
    tissue_positions['y'] *= spaceranger_json['tissue_hires_scalef']
    #
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
A = STalign.to_A(torch.tensor(L),torch.tensor(T))

#   Apply affine transform to source spatial coordinates
affine = np.matmul(
    np.array(A.cpu()),
    np.array(
        [
            tissue_positions_list[0]['y'],
            tissue_positions_list[0]['x'],
            np.ones(len(tissue_positions_list[0]['x']))
        ]
    )
)

xIaffine = affine[1,:]
yIaffine = affine[0,:]

#   Apply affine transform to source landmarks
ypointsI = pointsI[:,0]
xpointsI = pointsI[:,1]
affine = np.matmul(np.array(A.cpu()),np.array([ypointsI, xpointsI, np.ones(len(ypointsI))]))

xpointsIaffine = affine[1,:]
ypointsIaffine = affine[0,:]
pointsIaffine = np.column_stack((ypointsIaffine,xpointsIaffine))

################################################################################
#   Plot each image with its spatial coords and landmarks (aligned, for the
#   source capture area)
################################################################################

fig, ax = plt.subplots(1, 2)

#   Images
ax[0].imshow((I).transpose(1,2,0),extent=extentI)
ax[1].imshow((J).transpose(1,2,0),extent=extentJ)

#   Spatial coords
ax[0].scatter(xIaffine,yIaffine,s=1,alpha=0.1, label = 'source spatial coords')
ax[1].scatter(
    tissue_positions_list[1]['x'],
    tissue_positions_list[1]['y'],
    s=1, alpha=0.2, label = 'target spatial coords'
)

#   Landmarks and visual settings
for this_ax in ax:
    #   Landmarks
    this_ax.scatter(pointsIaffine[:,1],pointsIaffine[:,0],c="blue", label='source landmarks aligned', s=100, marker = "1")
    this_ax.scatter(pointsJ[:,1],pointsJ[:,0], c='red', label='target landmarks', s=100, marker = "2")
    #
    this_ax.legend()
    this_ax.set_aspect('equal')

#   Add titles
ax[0].set_title(f'Source ({capture_areas[0]})', fontsize=15)
ax[1].set_title(f'Target ({capture_areas[1]})', fontsize=15)

plt.savefig(str(plot_dir / '_'.join(capture_areas)) + '_aligned.png')
plt.clf()
