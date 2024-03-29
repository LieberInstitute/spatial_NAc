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