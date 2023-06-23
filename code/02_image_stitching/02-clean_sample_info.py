from pyhere import here
from pathlib import Path
import session_info
import os
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

################################################################################
#   Functions
################################################################################

#   Given the contents of the spaceranger '_invocation' file as a single string
#   'file_string', return the path to the image used as input to spaceranger
def parse_sr_invocation(file_string):
    #   Find the line containing the image path, and confirm there is only
    #   one image 
    match = re.findall(r'tissue_image_paths.*=.*\[".*"\]', file_string)
    assert len(match) == 1

    #   Grab just the path and return
    path = Path(re.sub(r'(tissue_image_paths|[ "\[\]=,])', '', match[0]))
    assert path.exists()
    return path

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

################################################################################
#   Merge in info about where initial transforms are
################################################################################

xml_map = pd.read_csv(xml_map_path, index_col = 'Slide')
xml_map.index.name = 'sample_id'
xml_map = xml_map.loc[:, ['Brain', 'XML file name']].copy()
xml_map['Slide #'] = [x.split('_')[0] for x in xml_map.index]
xml_map['Array #'] = [x.split('_')[1] for x in xml_map.index]

#   Merge sample_info and xml_map by index, keeping the union of their columns.
#   There should be a more elegant way to do this...
sample_info = pd.merge(
    sample_info, xml_map, how = "outer", left_index = True, right_index = True
)
sample_info['Brain'] = sample_info['Brain_x'].combine_first(sample_info['Brain_y'])
sample_info['Slide #'] = sample_info['Slide #_x'].combine_first(sample_info['Slide #_y'])
sample_info['Array #'] = sample_info['Array #_x'].combine_first(sample_info['Array #_y'])
sample_info.drop(
    ['Brain_x', 'Slide #_x', 'Brain_y', 'Slide #_y', 'Array #_x', 'Array #_y'],
    axis = 1, inplace = True
)

################################################################################
#   Determine the path to the full-resolution/raw images and spaceranger output
#   directories
################################################################################

sample_info['spaceranger_dir'] = [
    Path(x)
    for x in here(
        'processed-data', '01_spaceranger_reorg', sample_info.index, 'outs',
        'spatial'
    )
]

#   Grab the full-resolution ("raw") image paths by extracting which images were
#   used as input to spaceranger
raw_image_paths = []
for spaceranger_dir in sample_info['spaceranger_dir']:
    invocation_path = spaceranger_dir.parent.parent / '_invocation'
    
    #   Spaceranger has not been run for all samples
    if invocation_path.exists():
        with open(invocation_path, 'r') as f:
            invocation_str = f.read()
        
        raw_image_paths.append(parse_sr_invocation(invocation_str))
    else:
        raw_image_paths.append(None)

sample_info['raw_image_path'] = raw_image_paths

################################################################################
#   Add the initial transformation estimates for image stitching from ImageJ
################################################################################

transform_df_list = []
for slide in xml_map['Slide #'].unique():
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
