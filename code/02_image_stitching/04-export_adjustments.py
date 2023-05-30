from pyhere import here
import session_info

import pandas as pd
import json
import os

slides = ['V11U08-082', 'V11U23-406']

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)

estimate_path = here(
    'processed-data', '02_image_stitching', 'transformation_estimates_{}.csv'
)

out_path = here(
    'processed-data', '02_image_stitching', 'adjusted_transforms_highres.csv'
)

#   Read in sample info and subset to what we need
sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info = sample_info.loc[
    sample_info['Slide #'].isin(slides),
    ['Brain', 'Slide #', 'Array #', 'spaceranger_dir']
]
sample_info.index.name = 'sample_id'

#   Read in and concatenate adjusted estimates for the selected slides
estimates = pd.concat(
    [
        pd.read_csv(str(estimate_path).format(slide)) for slide in slides
    ]
)
estimates.set_index('sample_id', inplace = True)

estimates = estimates.loc[
    estimates['adjusted'], estimates.columns != 'adjusted'
]

estimates = pd.merge(estimates, sample_info, how = 'left', on = 'sample_id')

#   Scale translations to high resolution
for sample_id in estimates.index:
    #   Grab high-res scale factor
    json_path = os.path.join(
        estimates.loc[sample_id, 'spaceranger_dir'], 'scalefactors_json.json'
    )
    with open(json_path, 'r') as f:
        spaceranger_json = json.load(f)
    
    estimates.loc[sample_id, ['x', 'y']] *= spaceranger_json['tissue_hires_scalef']

(
    estimates
        .drop('spaceranger_dir', axis = 1)
        .to_csv(out_path)
)

session_info.show()
