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

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
out_dir = Path(here('processed-data', '13_stalign_tests'))
capture_areas = ['V12D07-078_B1', 'V12D07-078_D1']

out_dir.mkdir(parents = False, exist_ok = True)

#   Read in source and target images and normalize
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

#   Save images in format expected by 'point_annotator.py'
np.savez(
    out_dir / 'src_image',
    x = np.array(range(I.shape[2])) * 1.,
    y = np.array(range(I.shape[1])) * 1.,
    I = I
)
np.savez(
    out_dir / 'target_image',
    x = np.array(range(J.shape[2])) * 1.,
    y = np.array(range(J.shape[1])) * 1.,
    I = J
)

#   Then run annotate_landmarks.sh interactively