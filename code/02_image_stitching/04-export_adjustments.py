from pyhere import here
import session_info

import pandas as pd
import json
import os
import glob

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)

estimate_paths = glob.glob(
    str(
        here(
            'processed-data', '02_image_stitching',
            'transformation_estimates_*.csv'
        )
    )
)

out_path = here(
    'processed-data', '02_image_stitching', 'adjusted_transforms_highres.csv'
)

#   Read in sample info and subset to relevant columns
sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info = sample_info.loc[
    :, ['Brain', 'Slide #', 'Array #', 'spaceranger_dir']
]
sample_info.index.name = 'sample_id'

#   Read in and concatenate all adjusted estimates that exist so far
estimates = pd.concat(
    [pd.read_csv(x, index_col = 'sample_id') for x in estimate_paths]
)

estimates = estimates.loc[
    estimates['adjusted'], estimates.columns != 'adjusted'
]

#   Add in some sample info (particularly we'll need the spaceranger paths to
#   scale to high resolution)
estimates = pd.merge(estimates, sample_info, how = 'left', on = 'sample_id')

#   Scale translations to high resolution
for sample_id in estimates.index:
    #   Grab high-res scale factor
    json_path = os.path.join(
        estimates.loc[sample_id, 'spaceranger_dir'], 'scalefactors_json.json'
    )
    with open(json_path, 'r') as f:
        spaceranger_json = json.load(f)
    
    #   Scale to high resolution
    estimates.loc[sample_id, ['x', 'y']] *= spaceranger_json['tissue_hires_scalef']

(
    estimates
        .drop('spaceranger_dir', axis = 1)
        .to_csv(out_path)
)

session_info.show()
