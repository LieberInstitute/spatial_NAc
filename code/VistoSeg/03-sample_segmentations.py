from pyhere import here
import pandas as pd
import numpy as np
from PIL import Image
from pathlib import Path

Image.MAX_IMAGE_PIXELS = None

vs_inputs_path = Path(here('processed-data', 'VistoSeg', 'VistoSeg_inputs.csv'))
out_dir = Path(here('processed-data', 'VistoSeg'))

n_clusters = 5
subregion_size = 300
comp_ratio = 1

num_pix_buffer = int(subregion_size / comp_ratio / 10)
half_sub_size = int(subregion_size / 2)
assert subregion_size % 2 == 0

vs_inputs = pd.read_csv(vs_inputs_path)

#   Initialize the final combined image, which stacks all clusters vertically
#   and all samples horizontally so the nuclei cluster can be easily determined
#   for every sample in a single image. Note "n_clusters + 1", which is done
#   because the raw image is also included once vertically
final_img = np.zeros(
    (
        subregion_size * (n_clusters + 1) + num_pix_buffer * n_clusters,
        subregion_size * vs_inputs.shape[0] +
            num_pix_buffer * (vs_inputs.shape[0] - 1)
    ),
    dtype = np.uint8
)

#   Loop through all capture areas
for sample_index, sample_id in enumerate(vs_inputs.loc[:, 'sample_id']):
    print(f"Processing sample {sample_id} ({sample_index} of {vs_inputs.shape[0]})...")

    #   Grab the path to the raw image used as input to VNS
    raw_path = vs_inputs.loc[:, 'raw_image_path'].iloc[sample_index]

    #   Subset and compress segmentations for each cluster, placing the
    #   result on 'final_img'
    for cluster_index in range(n_clusters + 1):
        #   Actually, the 0th cluster will represent the raw image
        if cluster_index == 0:
            img_path = raw_path
        else:
            img_path = raw_path.replace('.tif', f'_cluster{cluster_index}.tif')
        
        #   Read in the image for one cluster and convert to grayscale
        img = np.array(Image.open(img_path).convert('L'))
        
        #   Subset image to a square subregion in the center of the image, whose
        #   length is 'subregion_size'
        center = (int(img.shape[0] / 2), int(img.shape[1] / 2))
        img = img[
            center[0] - half_sub_size: center[0] + half_sub_size,
            center[1] - half_sub_size: center[1] + half_sub_size
        ]
        
        #   Compress image according to 'compression_ratio'
        img = np.array(
            Image
                .fromarray(img)
                .resize(
                    (int(img.shape[0] / comp_ratio), int(img.shape[1] / comp_ratio))
                )
        )
        
        #   Place the image for this sample and cluster at the appropriate
        #   position on the final combined image
        start = (
            (img.shape[0] + num_pix_buffer) * cluster_index,
            (img.shape[1] + num_pix_buffer) * sample_index
        )
        final_img[
            start[0]: start[0] + img.shape[0],
            start[1]: start[1] + img.shape[1]
        ] = img

#   Save final image
(
    Image
        .fromarray(final_img)
        .save(out_dir / 'all_VNS_outputs_subsampled.tif')
)
