from pyhere import here
from pathlib import Path
import session_info

import numpy as np
import pandas as pd
import json
import os
import sys
import glob

tissue_paths = glob.glob(
    str(
        here(
            'processed-data', '02_image_stitching','tissue_positions*.csv'
        )
    )
)

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)

#   55-micrometer diameter for Visium spot
SPOT_DIAMETER_M = 55e-6

coords = pd.concat([pd.read_csv(x, index_col = 'key') for x in tissue_paths])
coords['sample_id'] = ['_'.join(x[-2:]) for x in coords.index.str.split('_')]
# coords.rename({'pxl_row_in_fullres.1': 'pxl_col_in_fullres'}, axis = 1, inplace = True)

sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info.index.name = 'sample_id'

#   Subset to 1 donor (TODO: not by slide; make more general)
sample_ids = sample_info.index[sample_info['Slide #'] == 'V12D07-074'].tolist()
this_sample_info = sample_info.loc[sample_ids, :].copy()

#   Add in some info from spaceranger
spaceranger_df_list = []
for spaceranger_dir in this_sample_info.loc[:, 'spaceranger_dir']:
    #   Grab high-res scale factor
    json_path = os.path.join(spaceranger_dir, 'scalefactors_json.json')
    with open(json_path, 'r') as f:
        spaceranger_df = pd.DataFrame({k:[v] for k,v in json.load(f).items()})
    
    spaceranger_df_list.append(spaceranger_df)

spaceranger_info = pd.concat(spaceranger_df_list)
spaceranger_info.index = this_sample_info.index

this_sample_info = pd.merge(
    this_sample_info, spaceranger_info, how = 'left', on = 'sample_id'
)

these_coords = coords.loc[coords['sample_id'].isin(sample_ids), :]

#   Goal: adjust these_coords.loc[:, ['array_row', 'array_col']] to most
#   appropriately integrate the capture areas, given their overlaps after
#   transforming/ stitching

#   One option: write a function to calculate best array_row (array_col)
#   given pxl_row_in_fullres (pxl_col_in_fullres) and relevant spaceranger info,
#   almost treating the merged image as one big capture area. Then apply the
#   function to overwrite 'array_row' and 'array_col'. If there are 2+ spots per
#   value, double the density of the spot array until there aren't?
