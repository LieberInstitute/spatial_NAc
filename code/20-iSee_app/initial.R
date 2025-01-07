initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "pp-GLM-PCA", XAxis = 1L,
                                          YAxis = 2L, FacetRowByColData = "sort", FacetColumnByColData = "sort",
                                          ColorByColumnData = "Barcode", ColorByFeatureNameAssay = "logcounts",
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "sort",
                                          SizeByColumnData = "sum", TooltipColumnData = character(0),
                                          FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "None",
                                          ColorByDefaultColor = "#000000", ColorByFeatureName = "KCNIP4",
                                          ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
                                          ColorBySampleName = "1_AAACCCAAGATGAAGG-1", ColorBySampleSource = "---",
                                          ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
                                          SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
                                          VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
                                          ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
                                          Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
                                          CustomLabelsText = "1_AAACCCAAGATGAAGG-1", FontSize = 1,
                                          LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
                                          LabelCenters = FALSE, LabelCentersBy = "sort", LabelCentersColor = "#000000",
                                          VersionInfo = list(iSEE = structure(list(c(2L, 17L, 3L)), class = c("package_version",
                                                                                                              "numeric_version"))), PanelId = c(ReducedDimensionPlot = 1L),
                                          PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE,
                                          RowSelectionSource = "---", ColumnSelectionSource = "---",
                                          DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                          RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                          SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "None",
                                      XAxisColumnData = "Barcode", XAxisFeatureName = "KCNIP4",
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
                                      YAxisFeatureName = "KCNIP4", YAxisFeatureSource = "---",
                                      YAxisFeatureDynamicSource = FALSE, FacetRowByColData = "sort",
                                      FacetColumnByColData = "sort", ColorByColumnData = "Barcode",
                                      ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
                                      ShapeByColumnData = "sort", SizeByColumnData = "sum", TooltipColumnData = character(0),
                                      FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "None",
                                      ColorByDefaultColor = "#000000", ColorByFeatureName = "KCNIP4",
                                      ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
                                      ColorBySampleName = "1_AAACCCAAGATGAAGG-1", ColorBySampleSource = "---",
                                      ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
                                      SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
                                      VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
                                      ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
                                      Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
                                      CustomLabelsText = "1_AAACCCAAGATGAAGG-1", FontSize = 1,
                                      LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
                                      LabelCenters = FALSE, LabelCentersBy = "sort", LabelCentersColor = "#000000",
                                      VersionInfo = list(iSEE = structure(list(c(2L, 17L, 3L)), class = c("package_version",
                                                                                                          "numeric_version"))), PanelId = c(FeatureAssayPlot = 1L),
                                      PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE,
                                      RowSelectionSource = "---", ColumnSelectionSource = "---",
                                      DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                      RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                      SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE,
                                        CustomRowsText = "KCNIP4", ClusterRows = FALSE, ClusterRowsDistance = "spearman",
                                        ClusterRowsMethod = "ward.D2", DataBoxOpen = FALSE, VisualChoices = "Annotations",
                                        ColumnData = character(0), RowData = character(0), CustomBounds = FALSE,
                                        LowerBound = NA_real_, UpperBound = NA_real_, AssayCenterRows = FALSE,
                                        AssayScaleRows = FALSE, DivergentColormap = "purple < black < yellow",
                                        ShowDimNames = "Rows", LegendPosition = "Bottom", LegendDirection = "Horizontal",
                                        VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10,
                                        ShowColumnSelection = TRUE, OrderColumnSelection = TRUE,
                                        VersionInfo = list(iSEE = structure(list(c(2L, 17L, 3L)), class = c("package_version",
                                                                                                            "numeric_version"))), PanelId = c(ComplexHeatmapPlot = 1L),
                                        PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE,
                                        RowSelectionSource = "---", ColumnSelectionSource = "---",
                                        RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                        RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                        SelectionHistory = list())
