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

unadjusted_donors = ['Br8667']

#   Read in sample info and only grab samples in analysis
sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info = sample_info.loc[
    ~sample_info['In analysis'].isna() & sample_info['In analysis'], :
]
sample_info.index.name = 'sample_id'

#   Read in and concatenate all adjusted estimates that exist so far
estimates = pd.concat(
    [pd.read_csv(x, index_col = 'sample_id') for x in estimate_paths]
)

estimates = estimates.loc[
    estimates['adjusted'],
    estimates.columns != 'adjusted'
]

#   Add in some sample info (particularly we'll need the spaceranger paths to
#   scale to high resolution)
estimates = pd.merge(estimates, sample_info, how = 'right', on = 'sample_id')

#   Use initial transforms where adjusted don't exist
for colname in ['x', 'y', 'theta']:
    #   Adjusted transformation estimates should be missing only for unadjusted
    #   donors
    assert all(
        estimates['Brain'].isin(unadjusted_donors) == estimates[colname].isna()
    )
    
    estimates.loc[estimates[colname].isna(), colname] = estimates.loc[
        estimates[colname].isna(), 'initial_transform_' + colname
    ]

#   Scale translations to high resolution
for sample_id in estimates.index:
    #   Grab high-res scale factor
    json_path = os.path.join(
        estimates.loc[sample_id, 'spaceranger_dir'], 'scalefactors_json.json'
    )
    with open(json_path, 'r') as f:
        spaceranger_json = json.load(f)
    #
    #   Scale to high resolution
    estimates.loc[sample_id, ['x', 'y']] *= spaceranger_json['tissue_hires_scalef']

#   Subset to important columns and export
(
    estimates
        .loc[:, ['x', 'y', 'theta', 'Brain', 'Slide #', 'Array #']]
        .to_csv(out_path)
)

session_info.show()
