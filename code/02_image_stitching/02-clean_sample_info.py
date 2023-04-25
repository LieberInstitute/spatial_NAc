from pyhere import here
from pathlib import Path
import session_info
import os

import pandas as pd

out_path = here('processed-data', '02_image_stitching', 'sample_info_clean.csv')

sample_info_paths = [
    here('raw-data', 'sample_info', 'Visium-samples-for-seq-1v_svb-16v_svb.xlsx'),
    here('raw-data', 'sample_info', 'Visium_NAc_Round3_072822_SVB-complete_info.xlsx'),
    here('raw-data', 'sample_info', 'Visium_NAc_Round4_081022_SVB-complete_info_final.xlsx')
]

sr_info_path = here('code', '01_spaceranger', 'spaceranger_parameters.txt')

#   Read in the individual sample sheets and concatenate
sample_info_list = [
        (
            pd.read_excel(x)
                .drop(columns = "NOTE", errors = "ignore")
                .dropna()
        )
            for x in sample_info_paths
    ]
sample_info = pd.concat(sample_info_list)

#   Make brain number consistent
sample_info['Brain'] = (
    sample_info['Brain']
        .astype(str)
        #   Add 'Br' if it's missing
        .replace(to_replace = r'^([0-9])', value = r'Br\1', regex = True)
        #   Remove decimals
        .replace(to_replace = '\.0$', value = '', regex = True)
)

#   TODO: add round as a column, and whatever information is necessary to in 
#   general find raw image paths and spaceranger outputs

sr_info = pd.read_csv(sr_info_path, header = 0, sep = '\t')

