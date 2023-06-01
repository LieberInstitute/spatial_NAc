from pyhere import here
from pathlib import Path
import session_info
import os
import glob
import re
import json

import pandas as pd
import numpy as np

out_path = here('processed-data', '02_image_stitching', 'sample_info_clean.csv')

sample_info_paths = [
    here('raw-data', 'sample_info', 'Visium-samples-for-seq-1v_svb-16v_svb.xlsx'),
    here('raw-data', 'sample_info', 'Visium_NAc_Round3_072822_SVB-complete_info.xlsx'),
    here('raw-data', 'sample_info', 'Visium_NAc_Round4_081022_SVB-complete_info_final.xlsx'),
    here('raw-data', 'sample_info', 'SPage_20230225.xlsx')
]

xml_map_path = here('raw-data', 'sample_key_spatial_NAc.csv')

#   Which "round" a given slide was part of: this info is necessary to uniquely
#   determine where the raw image is for a given sample. In some cases, there
#   isn't enough information to determine round, so it's left NA
rounds = {
    'V11U08-082': 1,
    'V11U08-083': 1,
    'V11U23-404': 2,
    'V11U23-406': 2,
    'V11U23-403': 3,
    'V11U23-405': 3,
    'V11D01-385': np.NaN,
    'V12D07-074': 5,
    'V11D01-386': np.NaN,
    'V12N28-334': np.NaN
}

################################################################################
#   Functions
################################################################################

#   Given the sample_info DataFrame and numeric index i, return the path to the
#   full-resolution/ raw image, or None if it can't be found
def get_img_path(sample_info, i):
    img_dir = here('raw-data', 'images', 'CS2')
    
    #   For samples with no known round, return None
    if sample_info['Round'].isna().iloc[i]:
        return None
    
    a = glob.glob(
        os.path.join(
            img_dir,
            'round' + str(int(sample_info['Round'].iloc[i])),
            '*' + sample_info['Brain'].iloc[i] + '_*' +
                sample_info['Array #'].iloc[i] + '.tif'
        )
    )

    assert len(a) == 1, "Did not find exactly one raw image"
    return a[0]

#   Get a numpy array representing transformation matrices from ImageJ (of shape
#   (n, 2, 3) for n capture areas), return theta (of shape (n,)), the angles
#   determined by the rotation-matrix portion of 'mat'
def theta_from_mat(mat):
    #   Determine theta by using arccos. Here the real angle might by theta_cos
    #   or -1 * theta_cos
    theta_cos = np.arccos(mat[:, 0, 0])

    #   Use 'theta_cos', but negate it when sin(theta) is negative, to account
    #   for the convention that arccos returns positive angles
    return theta_cos * (2 * (mat[:, 1, 0] > 0) - 1)

################################################################################
#   Gather sample info
################################################################################

#   Read in the individual sample sheets and concatenate
sample_info_list = [
        (
            pd.read_excel(x)
                .loc[:, ['Tissue', 'Brain', 'Slide #', 'Array #']]
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
        #   Remove prefix of "Hs_" if it's present
        .replace('^Hs_', '', regex = True)
        #   Remove decimals
        .replace(to_replace = '\.0$', value = '', regex = True)
)

#   Make tissue (brain region) consistent
sample_info['Tissue'] = sample_info['Tissue'].str.lower()
sample_info.loc[sample_info['Tissue'] == 'nac', 'Tissue'] = 'NAc'
sample_info.loc[sample_info['Tissue'] == 'dacc', 'Tissue'] = 'dACC'

sample_info.index = sample_info['Slide #'] + '_' + sample_info['Array #']

#   Add "round" where that information can be uniquely deduced
sample_info['Round'] = [rounds[x] for x in sample_info['Slide #']]

################################################################################
#   Determine the path to the full-resolution/raw images and spaceranger output
#   directories
################################################################################

#   Try a glob-based method to find image paths for all samples
sample_info['raw_image_path'] = [
    get_img_path(sample_info, i) for i in range(sample_info.shape[0])
]

sample_info['spaceranger_dir'] = [
    Path(x)
    for x in here(
        'processed-data', '01_spaceranger_reorg', sample_info.index, 'outs',
        'spatial'
    )
]

################################################################################
#   Add the initial transformation estimates for image stitching from ImageJ
################################################################################

xml_map = pd.read_csv(xml_map_path, index_col = 'Slide')
xml_map.index.name = 'sample_id'
xml_map['Slide #'] = [x.split('_')[0] for x in xml_map.index]

transform_df_list = []
for slide in xml_map['Slide #'].unique():
    assert slide in sample_info['Slide #'].unique(), f"{slide} not found in 'sample_info'"

    #   Open the ImageJ XML output for this slide. Assumes every array in a slide
    #   is in the same XML file associated with each slide
    imagej_xml_path = Path(
        here() /
        xml_map.loc[xml_map['Slide #'] == slide, 'XML file name'].iloc[0]
    )
    with open(imagej_xml_path) as f:
        imagej_xml = f.read()
    
    #   Clean the file; make sure new lines only separate XML elements
    imagej_xml = re.sub('\n', '', imagej_xml)
    imagej_xml = re.sub('>', '>\n', imagej_xml)

    array_nums = [
        re.sub('[_/]', '', x)
        for x in re.findall(rf'_[ABCD]1/', imagej_xml)
    ]
    
    #   Grab the transformation matrices and import as numpy array with the
    #   structure used in 01-samui_test.py
    matrices = re.findall(r'matrix\(.*\)".*file_path=', imagej_xml)
    trans_mat = [re.sub(r'matrix\((.*)\).*', '\\1', x).split(',') for x in matrices]
    trans_mat = np.transpose(
        np.array(
            [[int(float(y)) for y in x] for x in trans_mat], dtype = np.float64
        )
            .reshape((-1, 3, 2)),
        axes = [0, 2, 1]
    )
    assert trans_mat.shape[1:3] == (2, 3), "Improperly read ImageJ XML outputs"

    #   Grab sample info for this slide, ordered how the array numbers
    #   appear in the ImageJ output
    this_sample_info = sample_info.loc[
        [slide + '_' + array_nums[i] for i in range(len(array_nums))]
    ]

    #   Adjust translations to represent pixels in full resolution
    for i in range(len(array_nums)):
        #   Grab high-res scale factor and scale translations accordingly
        json_path = os.path.join(
            this_sample_info['spaceranger_dir'].iloc[i],
            'scalefactors_json.json'
        )
        with open(json_path, 'r') as f:
            spaceranger_json = json.load(f)
        
        trans_mat[i, :, 2] /= spaceranger_json['tissue_hires_scalef']

    #   Compile transformation information for this slide in a DataFrame
    transform_df_list.append(
            pd.DataFrame(
            {
                'initial_transform_x': trans_mat[:, 0, 2],
                'initial_transform_y': trans_mat[:, 1, 2],
                'initial_transform_theta': theta_from_mat(trans_mat),
                'Array #': array_nums,
                'Slide #': slide
            }
        )
    )

#   Combine across slides, then merge into sample_info
transform_df = pd.concat(transform_df_list)
index = sample_info.index
sample_info = pd.merge(
    sample_info, transform_df, how = 'left', on = ['Slide #', 'Array #']
)
sample_info.index = index

sample_info.to_csv(out_path)

session_info.show()
