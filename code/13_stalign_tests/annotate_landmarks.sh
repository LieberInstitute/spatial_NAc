resolution=lowres

module load stalign/1.0.1
python $(which point_annotator.py) \
    ../../processed-data/13_stalign_tests/src_image_${resolution}.npz \
    ../../processed-data/13_stalign_tests/target_image_${resolution}.npz
