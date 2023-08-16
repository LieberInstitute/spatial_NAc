The scripts in this directory take initial transformation estimates, calculated in ImageJ, and refine them using Samui, ultimately producing spatial coordinates and images that precisely align/ merge multiple capture areas into a single sample per donor. This is in preparation for creating a `SpatialExperiment` object, done in `05_harmony_BayesSpace`.

Scripts are run in the following order:

- `02-clean_sample_info.*`: Gather and clean existing sample-info sheets into a single CSV. Read in initial transformation estimates from ImageJ
- `01-samui_test.*`: Create the Samui-compatible directory for each donor (set one `donor` and `mode="initial"` in `01-samui_test.sh`)
- `03-adjust_transform.*`: After annotating landmark centroids in Samui, a CSV is exported and used here to determine the optimal adjustment to make to initial transformations
- `01-samui_test.*`: Create the Samui-compatible directory for the adjusted-transformed samples as a visual check of quality. Export transformed spot coordinates, scale factors, and the high-res merged image. Set `mode="adjusted"` in `01-samui_test.sh`
- `04-export_adjustments.*`: Export adjusted transformations in high resolution, for use in Visium Stitcher downstream
- `06-spatial_coords.*`: Compute the `array_row` and `array_col` pieces of the adjusted spot coordinates and export as a CSV (with pixel coordinates as well)
