from pyhere import here
from pathlib import Path
import session_info
import os
import glob

import pandas as pd

out_path = here('processed-data', '02_image_stitching', 'sample_info_clean.csv')

sample_info_paths = [
    here('raw-data', 'sample_info', 'Visium-samples-for-seq-1v_svb-16v_svb.xlsx'),
    here('raw-data', 'sample_info', 'Visium_NAc_Round3_072822_SVB-complete_info.xlsx'),
    here('raw-data', 'sample_info', 'Visium_NAc_Round4_081022_SVB-complete_info_final.xlsx')
]

sr_info_path = here('code', '01_spaceranger', 'spaceranger_parameters.txt')

################################################################################
#   Functions
################################################################################

#   Given the sample_info DataFrame and numeric index i, return the path to the
#   full-resolution/ raw image, or None if it can't be found
def get_img_path(sample_info, i):
    img_dir = here('raw-data', 'images', 'CS2')
    
    a = glob.glob(
        os.path.join(
            img_dir,
            'round*',
            '*_' + sample_info['Slide #'].iloc[i] + '*' +
                sample_info['Brain'].iloc[i] + '_*' +
                sample_info['Array #'].iloc[i] + '.tif'
        )
    )
    
    if len(a) == 1:
        return a[0]
    else:
        return None

################################################################################
#   Gather sample info
################################################################################

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

sample_info.index = sample_info['Brain'] + '_' + sample_info['Array #']

################################################################################
#   Determine the path to the full-resolution/raw images and spaceranger output
#   directories
################################################################################

#   Grab the image paths for samples run through spaceranger
sr_info = pd.read_csv(sr_info_path, header = 0, sep = '\t', index_col = 0)
sample_info['raw_image_path'] = sr_info.iloc[:, 2]

#   Try a glob-based method to find image paths for all samples
b = [get_img_path(sample_info, i) for i in range(sample_info.shape[0])]

#   We'll prefer the glob-based method. Just verify it agrees with the
#   spaceranger info where both are defined
non_null_indices = ~sample_info['temp'].isna() & \
    ~sample_info['raw_image_path'].isna()
assert all(sample_info['temp'][non_null_indices] == sample_info['raw_image_path'][non_null_indices])

#   Verify the spaceranger info is a subset of the glob-based paths (actually,
#   there's one path that is known to be bad)
sample_info['raw_image_path'][~sample_info['raw_image_path'].isna() & sample_info['temp'].isna()]

#   We'll just use the paths retrieved with the glob method
sample_info['raw_image_path'] = sample_info['temp']

sample_info['spaceranger_dir'] = [
    Path(x)
    for x in here(
        'processed-data', '01_spaceranger', sample_info.index, 'outs', 'spatial'
    )
]

sample_info.to_csv(out_path)

session_info.show()
