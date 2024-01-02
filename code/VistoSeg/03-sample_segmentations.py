from pyhere import here
import pandas as pd
import numpy as np
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

vs_inputs_path = here('processed-data', 'VistoSeg', 'VistoSeg_inputs.csv')
n_clusters = 5
subregion_ratio = 10
comp_ratio = 8
num_pix_buffer = 10

vs_inputs = pd.read_csv(vs_inputs_path)

for sample_id in vs_inputs.loc[:, 'sample_id']:
    #   Grab the path to the raw image used as input to VNS
    raw_path = vs_inputs.loc[
        vs_inputs.loc[:, 'sample_id'] == sample_id, 'raw_image_path'
    ][0]

    #   Subset and compress segmentations for each cluster, appending the
    #   result to 'img_list'
    img_list = []
    for i in range(n_clusters):
        #   Read in the image for one cluster and convert to grayscale
        img = np.array(
            Image
                .open(raw_path.replace('.tif', f'_cluster{i}.tif'))
                .convert('L')
        )
        
        #   Subset image to a subregion in the center of the image, whose size
        #   is determined by 'subregion_ratio'
        bottom = (subregion_ratio / 2 - 0.5) / subregion_ratio
        top = (subregion_ratio / 2 + 0.5) / subregion_ratio
        img = img[
            int(bottom * img.shape[0]): int(top * img.shape[0]),
            int(bottom * img.shape[1]): int(top * img.shape[1])
        ]

        #   Compress image according to 'compression_ratio'
        img = np.array(
            Image
                .fromarray(img)
                .resize(
                    (int(img.shape[0] / comp_ratio), int(img.shape[1] / comp_ratio))
                )
        )
        img_list.append(img)
    
    #   New images for each cluster should all be the same size
    assert all([x.shape == img_list[0].shape for x in img_list])

    #   Initialize the final image, which is a horizontal concatenation of the
    #   subsampled images
    total_buffer_size = num_pix_buffer * (len(img_list) - 1)
    final_img = np.zeros(
        (
            img_list[0].shape * len(img_list) + total_buffer_size,
            img_list[1].shape * len(img_list) + total_buffer_size
        ),
        dtype = np.uint8
    )

    #   Place all images in a horizontal series on the final image
    for i, img in enumerate(img_list):
        start = (img.shape[0] + num_pix_buffer) * i
        final_img[start: start + img.shape[0], :] = img
    
    #   Save results in the same directory
    (
        Image
            .fromarray(final_img)
            .save(raw_path.replace('.tif', '_all_clusters_subsampled.tif'))
    )
