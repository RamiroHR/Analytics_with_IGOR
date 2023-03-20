# IGOR procedures

## ACL_WindowDesktops
Allows to create multiple desktops within a single IGOR run for better data and plots organization. Include a Help file.

## Base
version: _RRv1
Basic set of functions for data analysis and visualisation.
Updated XYZ3Matrix function to remove I_offset.

## IGOR_Functions_RR
version: v5
Includes diverse functions created by me, such as:\

Math
- [x] Compute ellipse area and perimeter.
- [x] Equivalent series capa.
- [x] LC circuit resonant frequency.

For AFM data analysis:
- [x] ExtractLayer: Extract all Layers from the 3D matrix
- [x] xtract1Layer: Extract a selected Layer from the 3D matrix
- [x] xtractZoom: Extract a selected area from a 2D-plot
- [x] AreaJJ_v3: Measure JJ Area from a AFM picture.
- [x] AreaJJ_v4: Measure JJ Area (AFM pic) within designed zone (interactive)
- [x] MeanValues
- [x] PlotLinesProfiles
- [x] AddRetrapingToMap
- [x] AddRetrapingToMap_Max
- [x] AlignPeak
- [x] GetV_when_Iths
- [x] HideRetraping
- [x] HideSwitching
- [x] AlignAtMaxV (wV,Voffset)
- [x] CutRangeVJ_v0
- [x] CutRangeVJ_v1

## JJSpectro
version: v4
From script started by JG. Continue adding functionalities.
- [x] Menu JJ spectro
- [x] SeriesR
- [x] PlotIVs
- [x] GetLRC
- [x] FitPeakLowCoupling
- [x] FitPeakHighCoupling
- [x] ToBeSolvedResPeak
- [x] PeakHighCoupling, ..2, ..3
- [x] ToBeIntegratedqp
- [x] Iqp(V,alpha,Delta)
- [x] qpexcitation(w,V)
- [x] RemoveHalfPoints(w1,w2)


## Procedure_noise_analysis
version: v3
Functions for computing spectral noise density and related signals.
- [x] Function PlotPSDs: compute PSD of many waves.

## ykDL750
version: _1.19v2
Procedure for controlling sources and measurement equipments

