* `01-r_to_python.*`: Convert `SpatialExperiment` and snRNA-seq `SingleCellExperiment` to `AnnData`s. Save a copy of the modified SCE as an R object as well. This is a processing step prior to running `cell2location`. Finally, run `getMeanRatio2` to rank genes as markers, and save the resulting object. Convert R objects to `AnnData` objects for use in python, and rank all genes as markers.
* `02-find_markers.*`: find 20 "mean ratio" markers for each cell type and produce many exploratory plots to assess the quality of these markers
* `03-c2l_prepare_anndata.*`: finish proper conversion begun in `01-r_to_python.*`, and otherwise prepare the `AnnData` objects for running `cell2location`
* `04-c2l_registration.*`: run `cell2location` and export results to CSV
* `shared_functions.R`: partially deprecated code-- `spot_plot` has been replaced with `spatialNAcUtils::spot_plot`, but the other functions help form plots used in `02-find_markers.*`
