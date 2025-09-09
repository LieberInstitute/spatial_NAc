initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "tSNE_HARMONY", XAxis = 1L, 
                                          YAxis = 2L, FacetRowByColData = "Sample", FacetColumnByColData = "Sample", 
                                          ColorByColumnData = "CellType.Final", ColorByFeatureNameAssay = "logcounts", 
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Sample", 
                                          SizeByColumnData = "sum", TooltipColumnData = character(0), 
                                          FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data", 
                                          ColorByDefaultColor = "#000000", ColorByFeatureName = "KCNIP4", 
                                          ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
                                          ColorBySampleName = "1_AAACCCAAGACCAACG-1", ColorBySampleSource = "---", 
                                          ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
                                          SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
                                          VisualBoxOpen = TRUE, VisualChoices = "Color", ContourAdd = FALSE, 
                                          ContourColor = "#0000FF", FixAspectRatio = FALSE, ViolinAdd = TRUE, 
                                          PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                          CustomLabels = FALSE, CustomLabelsText = "1_AAACCCAAGACCAACG-1", 
                                          FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                          HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "Sample", 
                                          LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                            c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                            ))), PanelId = c(ReducedDimensionPlot = 1L), PanelHeight = 500L, 
                                          PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---", 
                                          ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, 
                                          ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE, 
                                          ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data", 
                                      XAxisColumnData = "CellType.Final", XAxisFeatureName = "KCNIP4", 
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE, 
                                      YAxisFeatureName = "DRD1", YAxisFeatureSource = "---", YAxisFeatureDynamicSource = FALSE, 
                                      FacetRowByColData = "Sample", FacetColumnByColData = "Sample", 
                                      ColorByColumnData = "CellType.Final", ColorByFeatureNameAssay = "logcounts", 
                                      ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Sample", 
                                      SizeByColumnData = "sum", TooltipColumnData = character(0), 
                                      FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data", 
                                      ColorByDefaultColor = "#000000", ColorByFeatureName = "KCNIP4", 
                                      ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
                                      ColorBySampleName = "1_AAACCCAAGACCAACG-1", ColorBySampleSource = "---", 
                                      ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
                                      SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
                                      VisualBoxOpen = TRUE, VisualChoices = "Color", ContourAdd = FALSE, 
                                      ContourColor = "#0000FF", FixAspectRatio = FALSE, ViolinAdd = TRUE, 
                                      PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                      CustomLabels = FALSE, CustomLabelsText = "1_AAACCCAAGACCAACG-1", 
                                      FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                      HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "Sample", 
                                      LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                        c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                        ))), PanelId = c(FeatureAssayPlot = 1L), PanelHeight = 500L, 
                                      PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---", 
                                      ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, 
                                      ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE, 
                                      ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE, 
                                        CustomRowsText = "GJA1\nDRD1\nDRD2\nCLDN5\nFOXJ1\nSLC17A7\nIL1RAPL2\nVIP\nGLP1R\nCHAT\nSST\nKCNC2\nARHGAP15\nST18\nPDGFRA", 
                                        ClusterRows = FALSE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2", 
                                        DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = "CellType.Final", 
                                        RowData = character(0), CustomBounds = TRUE, LowerBound = -5L, 
                                        UpperBound = 5L, AssayCenterRows = TRUE, AssayScaleRows = TRUE, 
                                        DivergentColormap = "blue < white < red", ShowDimNames = "Rows", 
                                        LegendPosition = "Bottom", LegendDirection = "Horizontal", 
                                        VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10, 
                                        ShowColumnSelection = TRUE, OrderColumnSelection = TRUE, 
                                        VersionInfo = list(iSEE = structure(list(c(2L, 18L, 0L)), class = c("package_version", 
                                                                                                            "numeric_version"))), PanelId = c(ComplexHeatmapPlot = 1L), 
                                        PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE, 
                                        RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                        RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                        RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                        SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "MIR1302-2HG", Search = "", SearchColumns = c("", 
                                                                                                           "", "", "", "", "", "", "", "", "", ""), HiddenColumns = character(0), 
                                  VersionInfo = list(iSEE = structure(list(c(2L, 18L, 0L)), class = c("package_version", 
                                                                                                      "numeric_version"))), PanelId = 1L, PanelHeight = 500L, PanelWidth = 6L, 
                                  SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                  DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                  RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                  SelectionHistory = list())

################################################################################
# Settings for Column data table 1
################################################################################

initial[["ColumnDataTable1"]] <- new("ColumnDataTable", Selected = "1_AAACCCAAGACCAACG-1", Search = "", 
                                     SearchColumns = c("", "", "", "", "", "", "", "", "", "", 
                                                       "", "", "", "", "", "", "", "", "", "", ""), HiddenColumns = character(0), 
                                     VersionInfo = list(iSEE = structure(list(c(2L, 18L, 0L)), class = c("package_version", 
                                                                                                         "numeric_version"))), PanelId = 1L, PanelHeight = 500L, PanelWidth = 6L, 
                                     SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                     DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                     RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                     SelectionHistory = list())
