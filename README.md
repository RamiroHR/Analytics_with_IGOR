# IGOR procedures

## ACL_WindowDesktops
Allows to create multiple desktops within a single IGOR run for better data and plots organization. Include a Help file.

## Base
version: _RRv1
Basic set of functions for data analysis and visualisation.
Updated XYZ3Matrix function to remove I_offset.

## Functions_RR
version: v7
Includes diverse functions created by me, such as:  

Math
- Compute ellipse area and perimeter.
- Equivalent series capa.
- LC circuit resonant frequency.

For AFM data analysis:
- ExtractLayer: Extract all Layers from the 3D matrix
- xtract1Layer: Extract a selected Layer from the 3D matrix
- xtractZoom: Extract a selected area from a 2D-plot
- AreaJJ_v3: Measure JJ Area from a AFM picture.
- AreaJJ_v4: Measure JJ Area (AFM pic) within designed zone (interactive)
- MeanValues
- PlotLinesProfiles
- AddRetrapingToMap
- AddRetrapingToMap_Max
- AlignPeak
- GetV_when_Iths
- HideRetraping
- HideSwitching
- AlignAtMaxV (wV,Voffset)
- CutRangeVJ_v0
- CutRangeVJ_v1
- PlotTracesFromRange #plot IV map traces selected by an index range [Ni,Nf] 
- AppendTracesFromRange  # similar to previous but append to an existing plot


## JJSpectro
version: v10_JG
From script started by JG (<v10). I continue adding functionalities (>v10).
- Menu JJ spectro
- InitializeConstants
- FitRN
- GetGoodGainMismatch
- UserCursorAdjust_ContButtonProc
- CorrectTilt
- SeriesR
- PlotIVs
- GetLRC
- FitPeakLowCoupling
- FitPeakOffLoop
- ToBeIntegratedq
- qpexcitation
- Iqp
- FitPeakHighCouplingQP
- RemoveHalfPoints
- TrueMod
- IcLagr
- IcLagrsameL
- IcLagrsameLm
- L1
- L2
- L3
- ToBeSolvedPlas
- Plasma
- GetVmax
- GetImax
- GetIsw
- GetIswm
- SuperIsw
- PlasmaRF
- GetIp
- GetVp
- FindPeaks
- FindPeaksMax
- GetVname
- CopyMatrixf
- NiceIV
- Legend0Pi
- MakeNicePlotFromYkdl
- VoltToHz
- NiceMap
- PlotAll2DWaves
- AddFreqAxis
- CurrentToFlux
- TopAxisFlux
- UserCursorAdjust
- SuperXYZ
- GetRvsT
- GetRvsTemp
- GetRvsTemp1
- AmpereToHz
- AddFreqAxisCurrent


## Procedure_noise_analysis
version: v3
Functions for computing spectral noise density and related signals.
- Function PlotPSDs: compute PSD of many waves.

## ykDL750
version: _1.19_v3
Procedure for controlling sources and measurement equipments

