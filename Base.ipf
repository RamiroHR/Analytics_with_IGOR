 //  PHI0: Base Utilities for Data/Graph Loading and Manipulation
//
//
//
#pragma rtGlobals = 1
#pragma version=2.03	// 2017.08.22 
							// Added EliminateOutliers: removes points captured during switching in IV curves

#pragma IgorVersion = 6.3

#include <KBColorizeTraces>
#include <Multi-peak fitting 2.0>
#include <All IP Procedures> 
#include <Image Saver>


Constant kElectronCharge = 1.602e-19
Constant kPhi0 = 3.291e-16
Constant khbar = 1.055e-34
 
StrConstant kPackagePath = "root:Packages:GQ:Base"

// BEGIN Menu
//
//
//

Menu "Base"
	"Copy Waves/1", /Q,CopyWaves()
	"Paste Waves/2",/Q,PasteWaves() 
	"-"
	"Symmetrize Graph/3",/Q, SymmetrizeGraph()
	"Center Graph/4",/Q,CenterGraph() 
	"Offset Y-axis",OffsetY()
	"Offset X-axis",OffsetX()
	"Combine Traces", CombineTraces()
	"-"
	"Same Z Scale", SameZScale("")
	"Get Background", GetBackground()
	"-"
	"Add Title Textbox",/Q, AddTitle()
	"Correct Bias Resistance", CorrectBiasResistor()
	"Clear Graph", ClearGraph("")
	"Shift Wave", WaveShift("")
	"-"
	"Load Matrix Panel",LoadMatrixPanel()
	"Load Matrix Partial Panel",LoadMatrixPanelPartial()
	"Load File Panel", LoadFilePanel()
	"Load File Panel BR", LoadFilePanel_BR()
	"Data Load Panel", DataLoadPanel()
	"Subtractor Panel", Subtractor()
End

Menu "GraphMarquee"
	"-"
	"Flatten Image Vert", FlattenImage()
	"Flatten Image Horz", FlattenImage(rank=1)
End


//
//
//
// END Menu


// BEGIN Function
//
//
//

Function Initialize_Base(reinit)
	Variable reinit
	String savDF= GetDataFolder(1) // Save current DF for restore.
	
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S GQ // Make sure this exists.
	
	if(!DataFolderExists("Base")  || reinit) // Already created?
		NewDataFolder/O/S Base // Our stuff goes in here.
		Make/O/N=(256,256) m, dm
//		Make/O/B/U/N=0 datafoldersel
		Make/O/T/N=1 datafolderlist
		Variable/G firstpoint=0, lastpoint, istart,istartBR,nbtoskip, iend, Vstartslow, Vfinalslow, xscale1 = 1, yscale1 = 1, yscale = 1, zscale = 1,npnts, xstep, xstart, iref, ifile
		Variable/G Vcenter, Vperiod, Delta, harmonic, nptsslow, channel_no, DeltaEjj
		Variable/G addtograph, correctr, Xcol=1, Ycol=2, skipLines=0
		Variable/G rbias = 2000
		Variable/G slowaxis = 0
		String/G filePathStr = "X:\ManipP08\BR3\Timedomain", filePrefix = "BR", xunits = "V", yunits = "A", zunits = "A", activeMatrix
		String/G xscale1str = "1", yscale1str = "1", xscalestr = "1", yscalestr = "1", zscalestr = "1"
	endif
	
	NewDataFolder/O root:Data
	SetDataFolder root:Packages:GQ:Base
	Make/O/T/N=1 datafolderlist
	
	SetDataFolder savDF // Restore current DF.
End

Function Load2Matrix(filePrefix, filePathStr, istart, iend,[extend])
	// Load Data Series into a matrix
	// filePathStr : data file path
	// filePrefix : prefix string for filenames, i.e. JPS3 for JPS3nnnn
	// istart : starting file index
	// iend : ending file index

	String filePrefix, filePathStr
	Variable istart, iend, extend

	String fileName, wName
	NewPath/O filePath filePathStr 
	sprintf wName, "%s%04d", filePrefix, istart
	
	if(ParamIsDefault(extend))
		extend = 0
		KillWaves/Z dm
	endif
	
	Variable i
	for(i=istart;i<=iend;i+=1)
		sprintf fileName, "%s%04d.DAT", filePrefix, i
		LoadWave/J/M/Q/N=$"m"/K=1/V={"\t","$",0,1}/P=filePath fileName
		
		Concatenate/NP=0 {m0}, dm
		
		If (mod(i-1,10) == 0)
			if(mod(i-1, 100) == 0)
				printf "%3g\r", i
			else
				printf "%3g ", i
			endif
		EndIf
	EndFor
	
	printf "\nDone. %s loaded, %d points total, %d columns per file.\r", wName, DimSize(dm,0), DimSize(dm,1)
	Duplicate/O dm $wName
	KillPath filePath
	return 0 // Signifies success.
End

Function SameZScale(graphName)
	String graphName
	
	String imagelist = ImageNameList(graphName,";")
	Variable i = 0, numimages = ItemsInList(imagelist)
	String rec = WMGetRECREATIONFromInfo(ImageInfo(graphName, StringFromList(0, imagelist),0))
	rec = RemoveEnding(rec)
	Print "ModifyImage " + StringFromList(0,imagelist) + " ctab = " + StringByKey("ctab",rec,"=",";",0)
	rec = ReplaceString(";",rec,",")
	For(i=1;i<numimages;i+=1)
		Execute "ModifyImage " + StringFromList(i,imagelist) + " " + rec
	EndFor

	return 0	
End

Function ChangeUnits(w,unit,unitstr)
	WAVE w
	String unit
	String unitstr
	
	//String wName = NameOfWave(w)
	Variable offset, delta
		
	strswitch(unit)	// numeric switch
		case "x":		// execute if case matches expression
			offset = DimOffset(w,0)
			delta = DimDelta(w,0)
			SetScale/P x,offset,delta,unitstr, w
			break				// exit from switch
		case "y":		// execute if case matches expression
			offset = DimOffset(w,1)
			delta = DimDelta(w,1)
			SetScale/P x,offset,delta,unitstr, w
			break				// exit from switch
		case "d":	  	// execute if case matches expression
			SetScale/P d,0,0,unitstr, w
			break				// exit from switch
	endswitch
End

Function LoadToMatrix(m, filePrefix, filePathStr, istart, iend,[quiet,skipLines,extension])
	// Load Data Series into a matrix
	// m : Wave reference
	// filePathStr : data file path
	// filePrefix : prefix string for filenames, i.e. JPS3 for JPS3nnnn
	// istart : starting file index
	// iend : ending file index
	// quiet : default equal to 0 shows progress
	// skipLines : lines to skip when loading data file
	// Modif HP 4/11/2013 to allow 5 digits in index
	WAVE m
	String filePrefix, filePathStr,extension
	Variable istart, iend
	Variable quiet, skipLines

	String fileName, wName,loadStr
	
	if(ParamIsDefault(extension))
		extension = ".dat"
	endif
	
	loadStr= filePrefix + extension
	NewPath/Z/Q/O filePath filePathStr 
	sprintf wName, loadStr, istart
	
	if(V_flag != 0)
		DoAlert/T="LoadToMatrix Error" 0,"Path not found!"
		return -1
	endif
	
	Variable refNum = 0
	Open/Z/R/P=filePath refNum as wName
	if(V_flag != 0)
		DoAlert/T="LoadToMatrix Error" 0,"File not found!"
		return -1
	endif
	Close refNum

	sprintf wName, filePrefix, istart
	if(!quiet)
		Printf "LoadToMatrix:\tloading %s to "+filePrefix+" in %s:\r\t", wName, iend, filePathStr
	endif
	
	if(ParamIsDefault(skipLines))
		skipLines = 0
	endif
		
	Variable i
	for(i=istart;i<=iend;i+=1)
		sprintf fileName, loadStr, i
		LoadWave/O/L={0, skipLines, 0, 0, 0 }/J/M/Q/N=$"m"/K=1/V={"\t","$",0,1}/P=filePath fileName

		//LoadWave/O/J/M/Q/N=$"m"/K=1/V={"\t","$",0,1}/P=filePath fileName
		
		if(i == istart)
			WAVE m0 = m0
			If(DimSize(m0,1) != DimSize(m,1))
				Redimension/N=(DimSize(m0,0),DimSize(m0,1)) m
				m = m0
			Endif
		else
			Concatenate/NP=0 {m0}, m
		endif
		
		if(!quiet)
			If (mod(i,10) == 0)
				printf "%5u ", i
				if(mod(i, 100) == 0)
					Printf "\r\t"
				endif
			EndIf
		endif
	EndFor
	
	sprintf wName, filePrefix, istart
	if(!quiet)
		printf "\r  Done. %s(%d,%d) loaded.\r", wName, DimSize(m,0), DimSize(m,1)
	endif
	
	KillPath filePath
	return 0 // Signifies success.
End


Function CenterGraph()
	If(strlen(CsrInfo(A)))
		WAVE wx = CsrXwaveRef(A)
		WAVE wy = CsrWaveRef(A)
	
		Variable dx = hcsr(A), dy = vcsr(A)
		wx -= dx
		wy -= dy
		Cursor/K A
	EndIf
End

Function OffsetY()
	If(strlen(CsrInfo(A)))
		WAVE w = CsrWaveRef(A)
		Variable dy = vcsr(A)
		w -= dy
	EndIf
End

Function OffsetX()
	If(strlen(CsrInfo(A)))
		WAVE w = CsrXwaveRef(A)
		Variable dx = hcsr(A)
		w -= dx
	EndIf
End

Function CombineTraces()
	String wysumname
	wysumname = NameOfWave(WaveRefIndexed("",0,1))+"_SUM"
	Concatenate/O/NP TraceNameList("",";",1+4), $wysumname
//	Print TraceNameList("",";",1+4)
	Duplicate/O $wysumname, wxsum
	wxsum = 0
	Variable i = 0, j = 0
	For(i = 0; WaveExists(WaveRefIndexed("",i,2)); i+=1)
		WAVE wx = WaveRefIndexed("",i,2)
		wxsum[j,j+numpnts(wx)-1] = wx[p-j]
		j += numpnts(wx)
	EndFor
	
	Duplicate/O wxsum, $(NameOfWave(WaveRefIndexed("",0,2))+"_SUM")
	DoWindow/K $wysumname
	Display/K=1/N=wysumname $(NameOfWave(WaveRefIndexed("",0,1))+"_SUM") vs $(NameOfWave(WaveRefIndexed("",0,2))+"_SUM")
End

Function RemoveAllTraces([graphName, matchStr])
	String graphName, matchStr
	
	if(ParamIsDefault(graphName))
		graphName = ""
	endif
	
	if(ParamIsDefault(matchStr))
		matchStr = "*"
	endif
	
	String traces = ListMatch(TraceNameList(graphName,";",7),matchStr)
	Variable i, n = itemsinlist(traces)
	For(i=0;i<n;i+=1)
		RemoveFromGraph/Z/W=$graphName $(StringFromList(i,traces))
	EndFor
	
	return 0
End

Function Unwrap2D(m,[direction,nocenter])
	WAVE m
	Variable direction, nocenter
	
	Duplicate/FREE m, m0
	
	if(ParamIsDefault(direction))
		MatrixTranspose m0
	endif			
	
	Variable nx = DimSize(m0,0), ny = DimSize(m0,1)
	Variable i
	Make/FREE/N=(nx) tempw
	
	For(i=0;i<ny;i+=1)
		tempw = m0[p][i]
		Unwrap 360, tempw
		m0[][i] = tempw[p]		
	EndFor
	
	if(ParamIsDefault(direction))
		MatrixTranspose m0
	endif
	
	
	if(ParamIsDefault(nocenter))
		Variable center = 360*round(mean(m0)/360)
		m0 -= center
	endif
	
	m = m0
End


Function GraphAllInDataFolder(dfStr, windowName, [hidelegend, wtypes, copy])
	String dfStr
	String windowName
	Variable hidelegend
	Variable wtypes // 1 for 1D waves, and 2 for 2D waves, and 3 for both
	Variable copy // copies plotted waves to root
	
	String savDF = GetDataFolder(1)
	SetDataFolder dfStr
	
	String wList, wName, prepend = ""
	Variable i, n
	
	if(ParamIsDefault(wtypes))
		wtypes = 3
	endif
	
	if(ParamIsDefault(hidelegend))
		hidelegend = 0
	endif

	if(ParamIsDefault(copy))
		copy = 0
	endif
	
	if(copy == 1)
		prepend = "root:"
	endif

	wList = WaveList("*",";","DIMS:1")
	n =  ItemsInList(wList)
		
	if(n && (wtypes & 1))
		// check if subwindow
		if(!StringMatch(windowName, "*#*"))
			DoWindow/K $(windowName)
			Display/N=$(windowName)
		endif
		
		For(i=0;i<n;i+=1)
			wName = StringFromList(i,wList)
			
			if(copy)
				Duplicate/O $wName $(prepend+wName)
				wName = prepend + wName
			endif
		
			AppendToGraph/W=$(windowName) $(wName)
		EndFor
			
		if(!hidelegend)
			if(n>0)
				Legend/W=$(windowName)/C/N=DataLoadPanelLegend
			else
				Legend/W=$(windowName)/K/N=DataLoadPanelLegend
			endif
		endif
	endif

	wList = WaveList("*",";","DIMS:2")
	n =  ItemsInList(wList)

	if(n && (wtypes & 2))
		For(i=0;i<n;i+=1)
			wName = StringFromList(i,wList)
			
			DoWindow/K $(wName+"_WIN");DoUpdate
			Display/N=$(wName+"_WIN")
	
			if(copy)
				Duplicate/O $(wName) $(prepend+wName)
			endif	
	
			AppendImage/W=$(wName+"_WIN") $(prepend+wName)
			
			if(!hidelegend)
				TextBox/C/N=WaveNameBox/A=RT wName
			endif

		EndFor
	endif
	
	SetDataFolder savDF
	
	return 0
End

Function FindPhase(w, f,[nsearch])
	WAVE w
	Variable f, nsearch
	
	if(ParamIsDefault(nsearch))
		nsearch = 128
	endif
	
//	Duplicate/FREE w, cosw, fitw, selw
	Duplicate/O w, cosw, fitw, selw
	Make/FREE/N=(nsearch) phiw, fitness
//	SetScale/I x,leftx(w),rightx(w), phiw, fitness
	phiw = -pi+2*pi*p/nsearch	

	Variable i=0, sv, wmax, wmin
	sv = sqrt(variance(w))
	wmax = wavemax(w)
	wmin = wavemin(w)
	selw = (abs(w) > sv && (w > (wmax-sv) || w < (wmin+sv))) ? 1 : 0
	For(i=0;i<nsearch;i+=1)
		cosw = sin(f*x+phiw[i]+pi/2)
		FastOp fitw = (1/sv)*w*cosw
//		FastOp fitw = (1/sv)*w*fitw
		FastOp fitw = fitw*fitw
		FastOp fitw = selw*fitw
//		fitness[i] = log(abs(sum(fitw)))
//		WaveStats/Q fitw
//		fitness[i] = wavemax(fitw)
		fitness[i] = sum(fitw)
	Endfor
	
	WaveStats/Q fitness
	Variable phi = phiw[V_minrowloc]
	phiw = phi-pi/nsearch+2*pi/(nsearch)^2*p
		
	For(i=0;i<nsearch;i+=1)
		cosw = sin(f*x+phiw[i]+pi/2)
		FastOp fitw = (1/sv)*w*cosw
//		FastOp fitw = (1/sv)*w*fitw
		FastOp fitw = fitw*fitw
		FastOp fitw = selw*fitw
//		fitness[i] = log(abs(sum(fitw)))
//		fitness[i] = wavemax(fitw)
		fitness[i] = sum(fitw)
	Endfor

	WaveStats/Q fitness
	phi = phiw[V_minrowloc]
	Duplicate/O fitness, fitnessw

	
	Return mod(phi+pi/2,pi)
End

Function CorrectMCoeff(dataName,peakpos,method)
	String dataName
	WAVE peakpos
	Variable method
	
	WAVE cw = $(dataName + "_coeffs")
	Duplicate/O cw, $(dataName + "_cV")
	WAVE cw = $(dataName + "_cV")
	Print cw[0][1], cw[0][2], cw[0][3], cw[0][10], cw[0][11], cw[0][12]
//	Make/N=(DimSize(cw,0))/FREE tw
//	tw = 1/2*(cw[][peakpos[0]] +cw[p][peakpos[1]])
//	Variable Vc = 
	cw[][peakpos[0]+1,peakpos[0]+2] +=cw[p][q+(peakpos[1]-peakpos[0])]
	cw[][peakpos[0]] -= cw[p][peakpos[1]]
	cw[][peakpos[0],peakpos[0]+2] /= 2
	//cw[][peakpos[0]+2] /= x
	Print cw[0][1], cw[0][2], cw[0][3]
	CorrectCoefficients(cw,peakpos[0]+1,peakpos[0],method=method)
	CorrectCoefficients(cw,peakpos[0]+2,peakpos[0],method=method)
		
	Make/O/N=(DimSize(cw,0)) $(dataName+"_width")
	WAVE widths = $(dataName+"_width")
	CopyScales cw, widths
	SetScale d,0,0,"V",widths
	widths = cw[p][peakpos[0]+1]
		
	Make/O/N=(DimSize(cw,0)) $(dataName+"_mcurr")
	WAVE mcurr = $(dataName+"_mcurr")
	CopyScales cw, mcurr
	SetScale d,0,0,"A",mcurr
	mcurr = 2*cw[p][peakpos[0]+2]/cw[p][peakpos[0]+1]/pi
	
	DoWindow/F $(dataName+"Fit")
	if(V_flag == 0)
		Display/N=$(dataName+"Fit") widths as (dataName+" Fit Width and Current")
		AppendToGraph/L=currL mcurr
		Execute/Q "FitStyleCoeff()"
	endif	

	Duplicate/O cw, $(dataName + "_cHz")
	WAVE cw = $(dataName + "_cHz")

	cw[][peakpos[0]+1] *= (kElectronCharge/(2*pi*khbar))
	cw[][peakpos[0]+2] /= (2*kElectronCharge*x)
	SetScale d,0,0,"Hz",cw
	Make/O/N=(DimSize(cw,0)) $(dataName+"_gamma"), $(dataName+"_trans")
	WAVE decay = $(dataName+"_gamma")
	WAVE trans = $(dataName+"_trans")
	CopyScales cw, decay, trans
	decay = cw[p][peakpos[0]+1]
	trans = cw[p][peakpos[0]+2]
	
	DoWindow/F $(dataName+"FitRate")
	if(V_flag == 0)
		Display/N=$(dataName+"FitRate") decay as (dataName+" Fit Rates")
		AppendToGraph trans
		Execute/Q "FitStyleRates()"
	endif
	
	return 0
End
	
Proc FitStyleRates() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z width=283.465,height={Aspect,1}
	ModifyGraph/Z mode=3
	ModifyGraph/Z marker[1]=8
	ModifyGraph/Z rgb[0]=(0,52224,26368),rgb[1]=(0,0,0)
	ModifyGraph/Z mrkThick=1
	ModifyGraph/Z log=1
	Label/Z left "Rates"
	Label/Z bottom "V\\BJJ\\M (\\U)"
	SetAxis/Z/A/N=1 left
	Legend/C/N=text0/A=RC/J "\\s(#0) Decay Rate\r\\s(#1) Transition Rate"
EndMacro

Proc FitStyleCoeff() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z width=283.465,height={Aspect,1}
	ModifyGraph/Z mode=3
	ModifyGraph/Z marker[1]=8
	ModifyGraph/Z rgb[0]=(0,52224,26368),rgb[1]=(0,0,0)
	ModifyGraph/Z mrkThick=1
	ModifyGraph/Z lblPos(left)=46,lblPos(currL)=46
	ModifyGraph/Z freePos(currL)={0,bottom}
	ModifyGraph/Z axisEnab(left)={0,0.48}
	ModifyGraph/Z axisEnab(currL)={0.52,1}
	Label/Z left "Width (\\U)"
	Label/Z bottom "V\\BJJ\\M (\\U)"
	Label/Z currL "Max Current (\\U)"
	SetAxis/Z/A/N=1/E=1 left
	SetAxis/Z/A/N=1/E=1 currL
	SetAxis/Z/A/N=1 bottom
EndMacro


Function CorrectCoefficients(coefw, c1, c2,[method])
	// Correct coefw by multiplying column c1 by derivative of column c2
	// method: specifies method to smooth derivative
	// 	method = 1 does no smoothing, method = 0 is linear fit
	//	method = n does boxcar smoothing with window n
	WAVE coefw
	Variable c1, c2, method
	
	if(ParamIsDefault(method))
		method = 0
	endif	
	
	Make/N=(DimSize(coefw,0))/FREE wc2 = coefw[p][c2]
	//Make/O/N=(DimSize(coefw,0)) wc2 = coefw[p][c2]
	
	CopyScales coefw wc2
	
	Differentiate wc2
	//Smooth/B 3, wc2
	if (method == 0)
		Variable slope = mean(wc2)
		coefw[][c1] /= abs(slope)
//		Print slope
	else
		Smooth/E=3 method, wc2
		coefw[][c1] /= abs(wc2[p])
//		Print wc2[0]
	endif
	
	return 0
End

Function CorrectMap(dataName, [offsets, phase, phaseIndex, nx, ny, neighbors, nanfactor, vbiasScale, rbias,singlesided])
	// CorrectMap:	interpolate flux maps onto regular grid 
	// 
	// Required parameter:
	//	dataName: string containing base name of data file, such as "JT66666" for "JT66666_21"
	//
	// Optional parameters: 
	//	Ex: CorrectMap("JT66666",offsets=offsetwave,nx=256,rbias = 2000)
	//	offsets: Name of offset wave, defaults to "offset" located in root. the x wave scaling should correspond to the applied bias voltages 
	//	phase: phase correction for sinusoidal flux biasing. defaults to automatic phase finder, which is not always reliable.
	//	phaseIndex: column index corresponding to trace to use for phase correction, defaults to last column.
	//	nx, ny: number of points to interpolate in the x and y direction.  They default to the smaller of 1024 and the respective dimensions of the wave
	//		You should set ny larger than DimSize(dataWave, 1) for maps below the plasma frequency.  I use ny = 256.  For maps above the plasma frequency, you can
	//		leave to default or increase 1.5x for extra smoothing.
	//	neighbors: neighborhood for interpolation. Can leave to default of 3 in most cases.
	//	nanfactor: IMPORTANT: in order not to have interpolation across the plasma frequency, pick a number such that nanfactor*(dimension delta)
	//		is smaller than the gap in data at the plasma frequency.  For standard nx, ny, this is usually nanfactor=8.
	//	vbiasScale: defaults to 1108, the spectrometry JJ voltage bias scaling 
	//	rbias: defaults to 2000, the specJJ bias resistance
	//	singlesided: to do a single sided interpolation of the data, set to 1. defaults to double-sided interpolation. i tried single-sided, but the results
	//		were worse
	//
	// Usage:
	//	imagine you have loaded data "JT66666_21" which is a map above the plasma frequency.  you have also loaded an itx with the offset
	//	information.  imagine column 33 has a strong signal (use the cursors to identify the column number). then:
	//
	//	CorrectMap("JT66666",phaseIndex=33)
	//	
	//	Ex. Output:
	//	Phase is x.xxxx
	//	nx = 1024
	//	ny = 134
	//	
	//	this will automatically find the phase and do the interpolation.  it will output the default values it has determined for the parameters.
	//	it will also open a window, "Phase Tuner" which allows you to manually find
	//	the phase.  Use the sliders to precisely determine the phase. now rerun CorrectMap:
	//
	//	CorrectMap("JT66666",phase=y.yyyyy,nanfactor=8)
	//
	//	Do "Display; AppendImage JT66666m" to see the results.  if you want more smoothing, try CorrectMap("JT66666",phase=y.yyyyy,ny=200).
	//	
	//	If JT66666_21 were a map below the frequency, I would use from ny=128 to ny=256.  Try it out to find a good value.  If there are lots of
	//	missing data, it will take a longer time to run.  Be patient.
	String dataName
	Wave offsets
	Variable phase, phaseIndex, nx, ny, neighbors, nanfactor, vbiasScale, rbias, singlesided

	if(ItemsInList(dataName) == 0)
		Print "ERROR: input wave nonexistent"
		return -1
	endif
	
	// to add multiple file incorporation later

	If(!WaveExists($(dataName + "_21")))
		Print "ERROR: input wave nonexistent"
		return -1
	endif

	WAVE wz = $(dataName + "_21")
	
	Variable n0 = DimSize(wz,0)
	Variable n1 = DimSize(wz,1)
	
	If(!WaveExists($(dataName+"_11")))
		Print dataname+"_11 Does not exist!  Creating..."
		WAVE/Z worig = root:Data:$(dataName):$(dataName)
		
		if(!WaveExists(worig))
			WAVE/Z worig = $(dataName)
		endif
		
		if(!WaveExists(worig))
			Print "ERROR: input wave nonexistent"
			return -1
		endif
		
		Duplicate/O wz $(dataName+"_11")
		WAVE wx = $(dataName+"_11")
		Variable maxpnts = max(n0,n1)
		//Print maxpnts
		wx = worig[p+q*maxpnts][0]
	else
		WAVE wx = $(dataName + "_11")
	endif

	Duplicate/FREE wz, wy
	
	If(ParamIsDefault(offsets))
		if(WaveExists(offset))
			WAVE offsets = offset
			WAVE bias = bias
		//	SetScale/I x,bias[0],bias[numpnts(bias)-1],"V", offsets, bias
			SetScale/I x,wavemax(bias),wavemin(bias),"V", offsets, bias
		//	Interpolate2/Y=offset_L/T=1 bias, offsets
		//	WAVE offsets = offset_L
		else
			Make/O/FREE/N=(n1) offsets = 0
		endif
	EndIf
	
	Smooth/M=(NaN) 4, offset
	
	If(ParamIsDefault(phaseIndex))
		phaseIndex = n1-1
		Print "phaseIndex = " + num2str(phaseIndex)
	Endif
	
	If(ParamIsDefault(nx))
		nx = min(round(n0),1024)
		Print "nx = " + num2str(nx)
	Endif
	
	If(ParamIsDefault(ny))
		ny = min(round(n1),1024)
		Print "ny = " + num2str(ny)
	Endif
	
	if(ParamIsDefault(neighbors))
		neighbors = 3
	endif
	
	if(ParamIsDefault(nanfactor))
		nanfactor = round(0.1*ny)
	endif
	
	if(ParamIsDefault(vbiasScale))
		vbiasScale = 1108   // NOTE: used to be 1069, then 1117
	endif
	
	if(ParamIsDefault(rbias))
		rbias = 2000		// Probably 10% smaller
	endif

	wy = y/vbiasScale-offsets(y)-rbias*wz
	
	SetScale d,0,0,"V", wx, wy
	
	KillWaves/Z W_coef
	CurveFit/Q/W=0 sin, wx[][0]
	WAVE coef = W_coef

	Make/FREE/N=(n0) w = wz[p][phaseIndex]
	CopyScales wz, w

	If(ParamIsDefault(phase))
		phase = FindPhase(w, coef[2])
		Print "Phase is", phase
	Endif
	
	TunePhase(w,coef[2],phase)
	Duplicate/FREE wx wxnew
	wxnew = coef[0]+coef[1]*sin(coef[2]*x+phase)
	Duplicate/FREE wy wynew
	Duplicate/FREE wz wznew
	
	if(!ParamIsDefault(singlesided))
		if(singlesided == 1)
			Redimension/N=(n0/2,n1) wxnew, wynew, wznew
			nx = nx/2
		endif
	endif
	
	String wName = dataName+"m"
	//XYZ2Matrix(wx,wy,wz,0,0,nx,0,0,ny,neighbors,wName)
	XYZ3Matrix(wxnew,wynew,wznew,wName,nx=nx,ny=ny,nanfactor=nanfactor,neighbors=neighbors)
//	Duplicate/O wy, vjj
//	Duplicate/O wx, vbob
	return 0
End

Function TunePhase(w,freq,phase)
	WAVE w
	Variable freq, phase
	
//	String dfSav = GetDataFolder(1)
//	SetDataFolder root:Packages:GQ:Base
//	NewDataFolder/O/S TunePhase
	Variable/G phi = phase, frequency = freq
	Duplicate/O w, sinw, datasin
	SetFormula sinw "sin(x*frequency+phi)"
	
	if(strlen(WinList("TunePhaseWin*",";","")) == 0)
		Display/N=$"TunePhaseWin" datasin vs sinw as "Phase Tuner"
	
		ModifyGraph margin(top)=100
	
		Slider TunePhase_sliderQ vert=0,pos={60,9},size={410,49},limits={-3.14159,3.14159,0},ticks=32,variable=phi,proc=TunePhaseSliderProc
		Slider TunePhase_sliderZ pos={65,56},size={410,49},vert=0,variable=phi;DelayUpdate
		Slider TunePhase_sliderZ limits={phi-pi/128,phi+pi/128,0},ticks=16
		SetVariable phi_display bodyWidth=60,value=phi,format="%.4f",noedit=0;DelayUpdate
		SetVariable phi_display limits={-inf,inf,0}, pos={237,114}
	endif
	
	return 0
End

//Function SetPhiVarProc(sva) : SetVariableControl
//	STRUCT WMSetVariableAction &sva
//
//	switch( sva.eventCode )
//		case 1: // mouse up
//		case 2: // Enter key
//		case 3: // Live update
//			Variable dval = sva.dval
//			String sval = sva.sval
//			Print "Phase is", dval
//			break
//		case -1: // control being killed
//			break
//	endswitch
//
//	return 0
//End

Function TunePhaseSliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa

	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				String slidername = sa.ctrlName
				slidername[16,17] = "Z"
				 
				Slider $(slidername) limits={curval-pi/128,curval+pi/128,0};DelayUpdate
			endif
			break
	endswitch

	return 0
End

Function ClearGraph(graphName)
	String graphName
	
	RemoveAllTraces(graphName = graphName)
	
	String annlist = AnnotationList(graphName)
	Variable i, n = ItemsInList(annlist)
	For(i=0;i<n;i+=1)
		TextBox/K/W=$graphName/N=$(StringFromList(i,annlist))
	EndFor
	
	return 0
End

Function WaveShift(graphName)
	String graphName	// Name of graph or "" for top graph
	Variable numSmoothingPasses	// Parameter to Smooth operation, e.g., 15
	
	// Get list of all traces in graph
	String traceList = TraceNameList(graphName, ";", 5), traceName
	Variable numTraces = ItemsInList(traceList)
	Variable traceIndex
	String info, offsetStr, mulStr
	Variable xMult, yMult, xOffset, yOffset
	
	for(traceIndex=0; traceIndex<numTraces; traceIndex+=1)
		traceName = StringFromList(traceIndex, traceList)
		
		xMult = 0
		yMult = 0
		xOffset = 0
		yOffset = 0
		
		// info will be "" if no offsets or something like "offset(x)={10,20}"
		info = TraceInfo(graphName, traceName, 0)

		offsetStr = StringByKey("offset(x)", info, "=")	// e.g., "{10,20}"
		mulStr = StringByKey("muloffset(x)",info,"=")
		
		if (strlen(offsetStr) > 0)
			sscanf offsetStr, "{%g,%g}", xOffset, yOffset
		endif
		
		if (strlen(mulStr) > 0)
			sscanf mulStr, "{%g,%g}", xMult, yMult
		endif
		
		ModifyGraph offset($traceName) = {0, 0}
		ModifyGraph muloffset($traceName) = {0, 0}
		
		WAVE wy = TraceNameToWaveRef(graphName,traceName)
		WAVE/Z wx = XWaveRefFromTrace(graphName,traceName)
		
		if (yOffset != 0)
			wy += yOffset
			Note wy, "yOffset = " + num2str(yOffset)
		endif

		if (yMult != 0)
			wy *= yMult
			Note wy, "yMult = " + num2str(yMult)
		endif
	
		if(WaveExists(wx))
			if (xOffset != 0)
				wx += xOffset
				Note wx, "xOffset = " + num2str(xOffset)
			endif

			if (xMult != 0) 
				wx *= xMult 
				Note wx, "xMult = " + num2str(xMult)
			endif
		else
			if (xOffset != 0)
				SetScale/P x,leftx(wy)+xOffset,deltax(wy),waveunits(wy,0), wy
				Note wy, "xOffset = " + num2str(xOffset)
			endif

			if (xMult != 0) 
				SetScale/P x,leftx(wy)+xOffset,deltax(wy)*xOffset,waveunits(wy,0), wy
				Note wy, "xMult = " + num2str(xMult)
			endif
		endif
	Endfor
end



//Function TraceToWave(graphName)
//	String graphName
//	
//	String traceName = StringByKey("TNAME",CsrInfo(A,graphName))
//	
//	Variable instance = 0
//	
//	if(StringMatch(traceName,"*#*")) 
//		// Trace Instance
//		String baseTraceName = StringFromList(0,traceName,"#")
//		instance = str2num(StringFromList(1,traceName,"#"))
//	else
//		baseTraceName = traceName
//	endif
//
//	WAVE wy = TraceNameToWaveRef(graphName,baseTraceName)
//	WAVE/Z wx = XWaveRefFromTrace(graphName,baseTraceName)
//		
//	Duplicate/FREE wy, wytrace
//	
//	Print tracename
//	
//	String shiftstr = StringByKey("offset(x)",TraceInfo(graphName,traceName,instance),"=",";")
//	String multstr = StringByKey("muloffset(x)", TraceInfo(graphName,traceName,instance),"=",";")
//	
//	Variable xmult, ymult, xshift, yshift
//	sscanf shiftstr, "{%f,%f}", xshift, yshift
//	sscanf multstr, "{%f,%f}", xmult, ymult
//	
//	Print xshift, yshift, xmult, ymult
//	
//	wytrace *= ymult
//	wytrace += yshift
//	
//	if(WaveExists(wx))
//		Duplicate/FREE wx, wxtrace
//		wxtrace *= xmult
//		wxtrace += xshift
//	else
//		SetScale/P x,leftx(wytrace)+xshift,deltax(wytrace)*xmult,waveunits(wytrace,0), wytrace
//	endif
//	
//	Duplicate/O wytrace, $(baseTraceName + "_" + num2str(instance))
//
//	if(WaveExists(wx))
//		Duplicate/O wxtrace, $(NameOfWave(wx) + "_" + num2str(instance))
//	endif
//	
//
//End
//



Function SymmetrizeGraph([wName, res, show])
	String wName
	Variable res, show

	if(ParamIsDefault(wName))
		wName = WinName(0,1)
	endif
	
	if(ParamIsDefault(res))
		res = 1024
	endif
	
	WAVE wx = WaveRefIndexed(wName,0,2)
	WAVE wy = WaveRefIndexed(wName,0,1)

	Variable xmin, xmax, ymin, ymax

	GetAxis/W=$(wName)/Q bottom
	xmin = V_min
	xmax = V_max
	GetAxis/W=$(wName)/Q left
	ymin = V_min
	ymax = V_max
	Symmetrize(wx,wy,res,xmin,xmax,ymin,ymax)
	
	If(show)
		Duplicate/O wx, wxs
		Duplicate/O wy, wys
		wxs *= -1
		wys *= -1
	
		Duplicate/O wxs, $(NameOfWave(wx)+"_SYM")
		Duplicate/O wys, $(NameOfWave(wy)+"_SYM")
	
	
		RemoveFromGraph/W=$wName/Z $(NameOfWave(wy)+"_SYM")
		AppendToGraph/W=$wName $(NameOfWave(wy)+"_SYM") vs $(NameOfWave(wy)+"_SYM")
		ModifyGraph/W=$wName rgb($(NameOfWave(wy)+"_SYM"))=(0,0,0)
	
		WAVE wxdif = $(NameOfWave(wx)+"_DIF"), wydif = $(NameOfWave(wy)+"_DIF")
		//	If (WaveExists(wxdif))
		//	RemoveFromGraph/W=$wName/Z $NameOfWave(wydif)
		//	AppendToGraph/W=$wName wydif vs wxdif; DoUpdate
		//	ModifyGraph/W=$wName rgb($NameOfWave(wydif))=(0,0,32768)
		//     EndIf		
	EndIf
End

Function Symmetrize(wx, wy, res,xmin,xmax,ymin,ymax)
       WAVE wx, wy
       Variable res, xmin, xmax, ymin, ymax
   
       Make/FREE/O/N=(res,res) w2d = 0

       SetScale/I x,xmin,xmax,WaveUnits(wx,1), w2d
       SetScale/I y,ymin,ymax,WaveUnits(wy,1), w2d

       Variable i, r,c,npnts = numpnts(wx)
       For(i=0;i<npnts;i+=1)		
		r = floor((wx[i]-xmin)/(xmax-xmin)*res)
		c = floor((wy[i]-ymin)/(ymax-ymin)*res)

		if(r >= 0 && r < res && c >= 0 && c < res)
			w2d[r][c] = 1
		endif
	EndFor
       
       Duplicate/FREE/O w2d, w2dt
       ImageRotate/O/H w2dt
       ImageRotate/O/V w2dt
       
       Variable dx, dy
      
	MatrixOp/O corrw = correlate(w2d,w2dt,4)
	MatrixFilter/O/N=4 gauss corrw
	WaveStats/Q corrw
	
	dx = (V_maxRowLoc > res/2 ? V_maxRowLoc - res : V_maxRowLoc)*DimDelta(w2d,0)/2-(xmax+xmin)/2
	dy = (V_maxColLoc > res/2 ? V_maxColLoc - res : V_maxColLoc)*DimDelta(w2d,1)/2-(ymax+ymin)/2

   	wx += dx
	wy += dy
	
	KillWaves/Z corrw
	
	return 0
 End
 

Function RemoveOffsets(m,n1,ncol, scaling, slowaxis, fname, pathStr)
	// Removes offsets in amplifiers using Offset data in file "fname" at path "pathStr"
	// m: matrix containing data from which offsets are removed
	// n1: index of first data file for offset removal
	// ncol: column number in offset data file
	// scaling: scaling of data
	// slowaxis: whether slow axis is Y-axis
	WAVE m
	Variable n1, ncol, scaling, slowaxis
	String fname, pathStr
	
	NewPath pathName, pathStr
	Variable fref
	Open /R/P=pathName fref as fname
	
	if (fref == 0)
		Print "ERROR: File not found."
		return -1
	endif
	
	String buffer
	Variable i = 0, len, tmp
	
	do
		FReadLine fref, buffer
		len = strlen(buffer)
		if (len == 0)
			Print "ERROR: Offset data not found."
			KillPath pathName
			Close fref
			return -1
		endif
		sscanf buffer, "%d", tmp
		if (tmp == (n1+1))
			break
		endif
	while (1)
	
	Variable n = DimSize(m,0)
	
	if(slowaxis) 
		n = DimSize(m,1)
	endif
	
	for(i=0;i<n;i+=1)
		FReadLine fref, buffer
		len = strlen(buffer)
		if (len == 0)
			break						// No more lines to be read
		endif
		if (slowaxis)
			m[][i] -= str2num(StringFromList(ncol-1,buffer,"\t"))/scaling
		else
			m[i][] -= str2num(StringFromList(ncol-1,buffer,"\t"))/scaling
		endif
		//Print str2num(StringFromList(ncol-1,buffer,"\t"))/scaling
		FReadLine fref, buffer
		FReadLine fref, buffer
		FReadLine fref, buffer
	endfor
		
		
	KillPath pathName
	Close fref

	Return 0
End

Function ReadOffsets(fname,pathStr, ncol)
	String fname, pathStr
	Variable ncol
	Make/O/N=(10000,2) offsets
	NewPath pathName, pathStr
	Variable fref
	Open /R/P=pathName fref as fname
	
	if (fref == 0)
		return -1
	endif
	
	String buffer
	Variable i = 0, len, tmp
	do
		FReadLine fref, buffer
		len = strlen(buffer)
		if (len == 0)
			break						// No more lines to be read
		endif
		sscanf buffer, "%d", tmp
		offsets[i][0] = tmp
		FReadLine fref, buffer
		//Print buffer
		tmp = str2num(StringFromList(ncol-1,buffer,"\t"))
		//Print tmp
		//break
		offsets[i][1] = tmp
		FReadLine fref, buffer
		FReadLine fref, buffer
		i += 1
	while (1)
	
	Redimension/N=(i,2) offsets
	KillPath pathName
	Close fref
	return 0
End

Function SubtractData(y1,x1,y2,x2,n,wname, filePath)
	// calculates difference between two XY waves
	WAVE y1, x1, y2, x2
	Variable n
	String wname, filePath
	
	Interpolate2/T=1/N=(n)/E=1/Y=y3a x1,y1
	Interpolate2/T=1/N=(n)/E=1/Y=y3b x2,y2
	
	Duplicate/O y3b y3
	
	y3 = y3b(x) - y3a(x)
	
	Duplicate/O y3, $(wname)
	
	SaveForMARFit(y3,wname, filePath)
	KillWaves y3a, y3b, y3
End

Function SaveForMARFit(w,wname,filePath)
	WAVE w
	String wname, filePath
	Duplicate w, wtmp
	wtmp = x
	NewPath pathname, filePath
	Variable refNum
	Open/P=pathname refNum as (wname+".dat")			// open file for write
	wfprintf refNum, "%.5e\t%.5e\r\n", wtmp, w	// print 6 values each
	Close refNum
	KillPath pathname
	Duplicate/O wtmp, $(wname+"x")
	KillWaves wtmp
End

Function FlattenImageOld([rank])
	Variable rank

	String winrec = WinRecreation("",0)
	WAVE m = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
	
	if(ParamIsDefault(rank) || (rank == 0))
		GetMarquee/K left
	
		Duplicate/O m, $(NameOfWave(m) + "_F0")
		WAVE mf = $(NameOfWave(m) + "_F0")
		
		Variable q1 = (V_bottom - DimOffset(m,1))/DimDelta(m,1)
		Variable dq = floor(abs((V_top - V_bottom)/DimDelta(m,1)))
	
		Make/FREE/N=(DimSize(m,0),dq) subm
	
		subm = m[p][q+q1]
	
		MatrixOp/O subm = sumRows(subm)/numCols(subm)
	
		mf -= subm[p]
	else
		GetMarquee/K bottom
	
		Duplicate/O m, $(NameOfWave(m) + "_F1")
		WAVE mf = $(NameOfWave(m) + "_F1")
		
		Variable p1 = (V_left - DimOffset(m,0))/DimDelta(m,0)
		Variable dp = floor(abs((V_right - V_left)/DimDelta(m,0)))
	
		Make/FREE/N=(dp,DimSize(m,1)) subm
	
		subm = m[p+p1][q]
	
		MatrixOp/O subm = sumCols(subm)/numRows(subm)
		MatrixTranspose subm
		mf -= subm[q]
	endif
	
	Execute/Q winrec
	
	AppendImage mf
	RemoveImage $(NameOfWave(m))
	
	return 0
end

Function GetBackground()
	WAVE m = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))

	GetMarquee/K left, bottom 
	Variable p1 = round(max((V_left - DimOffset(m,0))/DimDelta(m,0),0))
	Variable p2 = round(min((V_right - DimOffset(m,0))/DimDelta(m,0),DimSize(m,0)-1))
	
	Variable q1 = round(max((V_bottom - DimOffset(m,1))/DimDelta(m,1),0))
	Variable q2 = round(min((V_top - DimOffset(m,1))/DimDelta(m,1),DimSize(m,1)-1))
	
	Variable dp = abs(p2-p1)
	Variable dq = abs(q2-q1)
	
	Variable tmp
	if(q2 < q1)
		tmp = q2
		q2 = q1
		q1 = tmp
	endif
	
	if(p2 < p1)
		tmp = p2
		p2 = p1
		p1 = tmp
	endif
	
	Make/FREE/N=(dp,dq) subm
	//Make/O/N=(dp,dq) subm
	subm = m[p+p1][q+q1]

	Variable/G V_background = mean(subm)
	return V_background
End

Function HasNan(w)
	WAVE w
	
	Variable i, n = numpnts(w)
	For(i=0;i<n;i+=1)
		if(numtype(w[i])==2)
			return 1
		endif
	EndFor
	
	return 0
End

Function MeanNan(w)
	WAVE w
	
	Duplicate/FREE w, wtmp
	
	Redimension/N=(numpnts(w)) wtmp
	
	WaveTransform zapNans, wtmp
	
	return mean(wtmp)
End

Function FlattenImage([rank, zerobkg])
	Variable rank, zerobkg
	
	String winrec = WinRecreation("",0)
	WAVE m = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))

	Variable firsttime = 1, swapXY
	
	swapXY = stringmatch(winrec,"*swapXY=1*")		
	
	if(swapXY)
		MatrixTranspose m
	endif
	
	if(stringmatch(NameOfWave(m),"*_F"))
		firsttime = 0
		WAVE mf = m
	elseif(ItemsInList(WinList(NameOfWave(m) + "_FW",";","WIN:1")))
		firsttime = 0
		Duplicate/O m, $(NameOfWave(m) + "_F")
		WAVE mf = $(NameOfWave(m) + "_F")
	else
		Duplicate/O m, $(NameOfWave(m) + "_F")
		WAVE mf = $(NameOfWave(m) + "_F")
	endif

	GetMarquee/K left, bottom 
	
	Variable p1 = 0, p2 = DimSize(m,0)
	Variable q1 = 0, q2 = DimSize(m,1)

	// Bug fixed in version 1.20: +1 was missing in definitions of dp and dq
	Variable dp = abs(p2-p1)+1
	Variable dq = abs(q2-q1)+1
	
	if(V_flag)
		p1 = round(max((V_left - DimOffset(m,0))/DimDelta(m,0),0))
		p2 = round(min((V_right - DimOffset(m,0))/DimDelta(m,0),DimSize(m,0)-1))
	
		q1 = round(max((V_bottom - DimOffset(m,1))/DimDelta(m,1),0))
		q2 = round(min((V_top - DimOffset(m,1))/DimDelta(m,1),DimSize(m,1)-1))
			
		dp = abs(p2-p1)+1
		dq = abs(q2-q1)+1

		Variable tmp
		if(q2 < q1)
			tmp = q2
			q2 = q1
			q1 = tmp
		endif
	
		if(p2 < p1)
			tmp = p2
			p2 = p1
			p1 = tmp
		endif
	endif
//	
//	

//	if(ParamIsDefault(meanbkg))
//		NVAR background = V_background
//		if(NVAR_Exists(background))
//			bkg = background
//		else
//			bkg = 0
//		endif
//	endif
	

	Make/FREE/N=(dp,dq) subm
	Duplicate/FREE mf, bkg
	//Make/O/N=(dp,dq) subm
	subm = m[p+p1][q+q1]
	
	Variable meanbkg
	
	if(HasNan(subm))
		MatrixFilter/N=(3) NanZapMedian subm
	endif
	
	meanbkg = mean(subm)
	
	if(ParamIsDefault(rank) || (rank == 0)) // FlattenImage Vert
		if(dp == DimSize(m,0))
			bkg = meanbkg
		elseif(dq == DimSize(m,1))
			//	MatrixOp/O/FREE bkg = sumCols(subm)/numRows(subm)
			//Redimension/N=(DimSize(m,1)) bkg 
			//bkg = (subm[dp-1][p]+subm[0][p])/2 + meanbkg
			bkg = meanbkg
		else
			Redimension/N=(dq) bkg
			bkg = (subm[dp-1][p]+subm[0][p])/2
			meanbkg = mean(bkg)
			bkg = meanbkg
		endif
		
		MatrixOp/O/FREE subm = sumRows(subm)/numCols(subm)
		
		if(ParamIsDefault(zerobkg))
			mf[p1,p2][] -= subm[p-p1]
		else
			mf[p1,p2][] -= subm[p-p1] - bkg[q]
		endif
		
	else  							//FlattenImage Horz
		if(dq == DimSize(m,1))
			bkg = meanbkg
		elseif(dp == DimSize(m,0))
			//MatrixOp/O/FREE bkg = sumRows(subm)/numCols(subm)
			//Redimension/N=(DimSize(m,0)) bkg
			//bkg = (subm[p][0]+subm[p][dq-1])/2
			bkg = meanbkg
		else
			Redimension/N=(dp) bkg
			bkg = (subm[p][0]+subm[p][dq-1])/2
			meanbkg = mean(bkg)
			bkg = meanbkg
		endif

		MatrixOp/O/FREE subm = sumCols(subm)/numRows(subm)
		MatrixTranspose subm

		if(ParamIsDefault(zerobkg))
			mf[][q1,q2] -= subm[q-q1]
		else
			mf[][q1,q2] -= subm[q-q1] - bkg[p]
		endif
	endif
	
	if(swapXY)
		MatrixTranspose mf
		MatrixTranspose m
	endif
	
	if(firsttime)
		Execute/Q winrec
		ReplaceWave image=$(NameOfWave(m)), $(NameOfWave(m)+"_F")
		String windowName = winName(0,1)
		String newWindowName = NameOfWave(m) + "_FW"
		RenameWindow $(windowName), $(newWindowName)
	endif
	
	KillWaves/Z W_WaveList
	
	return 0
end	

Function FlattenImagePeak([rank, bkgorder])
	Variable rank, bkgorder
	
	if(ParamIsDefault(bkgorder))
		bkgorder = 2
	endif
	
	String winrec = WinRecreation("",0)
	WAVE m = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))

	Variable firsttime = 1
	
	if(stringmatch(NameOfWave(m),"*_FP"))
		firsttime = 0
		WAVE mf = m
	else
		Duplicate/O m, $(NameOfWave(m) + "_FP")
		WAVE mf = $(NameOfWave(m) + "_FP")
	endif

	GetMarquee/K left, bottom 
	Variable p1 = round(max((V_left - DimOffset(m,0))/DimDelta(m,0),0))
	Variable p2 = round(min((V_right - DimOffset(m,0))/DimDelta(m,0),DimSize(m,0)-1))
	
	Variable q1 = round(max((V_bottom - DimOffset(m,1))/DimDelta(m,1),0))
	Variable q2 = round(min((V_top - DimOffset(m,1))/DimDelta(m,1),DimSize(m,1)-1))

	Variable dp = abs(p2-p1)
	Variable dq = abs(q2-q1)
	
	Variable tmp
	if(q2 < q1)
		tmp = q2
		q2 = q1
		q1 = tmp
	endif
	
	if(p2 < p1)
		tmp = p2
		p2 = p1
		p1 = tmp
	endif

	Make/O/N=(dp,dq) subm
	//Duplicate/O subm, bkg
	Make/FREE/N=(4) lorcoeff
	subm = m[p+p1][q+q1]
	//SetScale/P x,p1*DimDelta(m,0)+DimOffset(m,0),DimDelta(m,0),WaveUnits(m,0),subm
	//SetScale/P y,q1*DimDelta(m,1)+DimOffset(m,1),DimDelta(m,1),WaveUnits(m,1),subm

	if(ParamIsDefault(rank) || (rank == 0)) // FlattenImage Vert
		MatrixOp/O/FREE subm = sumRows(subm)/numCols(subm)
		Redimension/N=(dp) subm
		CurveFit/Q lor, kwCWave=lorcoeff, subm
		Duplicate/O subm, bkg
		if(bkgorder==0)
			bkg = lorcoeff[1]/((x-lorcoeff[2])^2+lorcoeff[3])
			
		endif
	//	if(bkgorder==1)
	//		FuncFit lorpoly1, cwaveName, waveName  [ flag parameters ]
	//	elseif(
		mf[p1,p2][] -= bkg[p-p1]
	else  							//FlattenImage Horz
		MatrixOp/O subm = sumCols(subm)/numRows(subm)
		MatrixTranspose subm
		//Print DimSize(subm,0), DimSize(subm,1), dp, dq
		
		Redimension/N=(dq) subm
		Duplicate/O subm,bkg
		CurveFit/Q lor, kwCWave=lorcoeff, subm
		if(bkgorder==0)
			bkg = lorcoeff[1]/((p-lorcoeff[2])^2+lorcoeff[3])
		elseif(bkgorder==1)
			Duplicate/FREE lorcoeff, fitcoeff
			Redimension/N=5 fitcoeff
			fitcoeff[0] = lorcoeff[1]
			fitcoeff[1] = lorcoeff[3]
			fitcoeff[2] = lorcoeff[2]
			fitcoeff[3] = lorcoeff[0]
			fitcoeff[4] = 0
			FuncFit/Q lorpoly1, fitcoeff, subm
			bkg = fitcoeff[0]/((p-fitcoeff[2])^2+fitcoeff[1])
		else
			Duplicate/FREE lorcoeff, fitcoeff
			Redimension/N=6 fitcoeff
			fitcoeff[0] = lorcoeff[1]
			fitcoeff[1] = lorcoeff[3]
			fitcoeff[2] = lorcoeff[2]
			fitcoeff[3] = lorcoeff[0]
			fitcoeff[4,5] = 0
			FuncFit/Q lorpoly2, fitcoeff, subm
			bkg = fitcoeff[0]/((p-fitcoeff[2])^2+fitcoeff[1])
		endif


		mf[][q1,q2] -= bkg[q-q1]
	endif
	
//	Print lorcoeff
//	firsttime = 0
	if(firsttime)
		Execute/Q winrec
		AppendImage mf
		RemoveImage $(NameOfWave(m))
	endif
	
	return 0
end	

Function/WAVE RunLengths(w,[threshold])
	WAVE w
	Variable threshold
	
	Variable i, n = numpnts(w), counts = 0, countIndex = 0, dt = deltax(w)

	if(ParamIsDefault(threshold))
		threshold = FindBimodalThreshold(w)
	endif

	Make/FREE/N=(n) wruns = NaN

	For(i=0;i < n;i+=1)
		if(w[i] > threshold)
			counts += 1
		elseif(counts > 0)
			wruns[countIndex] = counts*dt
			countIndex += 1
			counts = 0
		endif
	EndFor
	
	Redimension/N=(countIndex+1) wruns
	SetScale d,0,0,WaveUnits(w,0),wruns
	Duplicate/O wruns, $(NameOfWave(w)+"_runs")
		
	WAVE wrunsref = $(NameOfWave(w)+"_runs")

	return wrunsref
End


Function FindBimodalThreshold(w)
	Wave w

	Variable n = numpnts(w)
	
	Variable middle = (wavemax(w)+wavemin(w))/2
	Variable counts = 0
	
	Duplicate/FREE w, wlower, wupper
	
	wlower = (wlower[p] < middle ? wlower[p] : NaN)
	wupper = (wupper[p] >= middle ? wupper[p] : NaN)
	
	WaveTransform zapNans wlower
	WaveTransform zapNans wupper
	
	//Print (mean(wupper)+mean(wlower))/2
	
	return (mean(wupper)+mean(wlower))/2
End

Function MatrixMap(m, f)
	// Maps a function f which operates on column waves onto a matrix
	// the input is a 2D wave m ( np x nq ) and a function f : 1D (np) --> 1D (npOut)
	// the output is a 2D wave mOut ( npout x nq )
	// the size of f(mcol), npOut, must be the same for each column mcol of m
	WAVE m
	FUNCREF protofunction f

	Variable np = DimSize(m,0), nq = DimSize(m,1)
	String funcName = StringByKey("NAME", FuncRefInfo(f))
	Variable i
	
	Make/FREE/N=(np) mcol = m[p][0]
	CopyScales m, mcol
	
	WAVE w = f(mcol)
	Variable npOut = DimSize(w,0)
	
	Redimension/N=(npOut,nq) w
	
	Duplicate/O w, $(NameOfWave(m)+"_"+funcName)
	WAVE mOut = $(NameOfWave(m)+"_"+funcName)

	KillWaves $NameOfWave(w)

	For(i=1;i<nq;i+=1)
		mcol = m[p][i]
		WAVE w = f(mcol)
		mOut[][i] = w[p]
		KillWaves $NameOfWave(w)
	EndFor

	WaveClear w

	SetScale/P y,DimOffset(m,1),DimDelta(m,1),WaveUnits(m,1),mOut
	
	return 0
End

Function RunLengthHistOfMatrix(m,nhist)
	WAVE m
	Variable nhist
	// Matrix m produced by matrixmap(min,RunLengths)
	
	Variable np = DimSize(m,0), nq = DimSize(m,1)
	
	Variable i
	
	Make/FREE/N=(np) w
	Make/FREE/N=(nhist,nq) wout
	Make/FREE/N=(nhist) whist
	
	SetScale/P x,DimOffset(m,0),DimDelta(m,0),WaveUnits(m,0), w
	SetScale/P y,DimOffset(m,1),DimDelta(m,1),WaveUnits(m,1), wout
	SetScale d,0,0,WaveUnits(m,0),whist,wout
	SetScale/I x,0,wavemax(m),WaveUnits(m,0),whist,wout
	
	For(i=1;i<nq;i+=1)
		w = m[p][i]
		
		Histogram/B=2/C w, whist
		
		wout[][i] = whist[p]
	EndFor
	
	Duplicate/O wout, $(CleanupName(NameOfWave(m)+"_hist",0))
	
	return 0
End

Function/WAVE protofunction(w,[optparam])
	WAVE w
	Variable optparam
	
	Duplicate/O w, $(NameOfWave(w)+"_2")
	
	WAVE wout = $(NameOfWave(w)+"_2")
	
	wout *= 2
	
	return wout
End

Function ShiftIVs(iw,vw,istart,iend,nb,vp1,vp2,voffset,ip1,ip2,ioffset)
	WAVE iw, vw
	Variable istart, iend, nb, vp1, vp2, voffset, ip1,ip2, ioffset
	
	Make/O/N=(nb) wx, wy
	
	Variable i, j = 0, offset
	For(i=istart;i<iend;i+=1)
		wx = vw[p+j]
		wy = iw[p+j]
		
		if(vp2 > 0)
			// to correct for voltage offsets
			offset = mean(wx,vp1,vp2)
			wx -= offset
			// to change offset value for voltage away from zero:
			wx += voffset
		endif
	
		if(ip2 > 0)
			// to correct for current offsets:
			offset = mean(wy, ip1, ip2)
			wy -= offset
			wy += ioffset
		endif
	
		vw[j,j+nb] = wx[p-j]
		iw[j,j+nb] = wy[p-j]
		j += nb
	EndFor
End

Function RemoveZeros(iw,vw,nb)
	WAVE iw, vw
	Variable nb
	
	Variable i, np = numpnts(iw)
	For(i=0;i<np;i+=1)
		If(iw[i] == 0)
			iw[i] = iw[i-nb]
			vw[i] = vw[i-nb]
		Endif
	EndFor
End

//#if Exists("MPF2_DoMPFit")

//Structure MPFitInfoStruct
//EndStructure
//#endif
// Function designed to take as input the name of a Multi-peak Fit 2
// set folder found in root:Packages:MultiPeakFit2:, along with a
// list of input data waves. This function loops through the list using
// the coefficient waves found in the data folder as the initial guesses
// for each of the data waves found in ywavelist, wxavelist.
Function MultiPeakBatchFit(datafolderName, matrix,[istart,iend,show,fitpnts,transpose])
	String datafolderName
	WAVE matrix
	Variable istart, iend, show, fitpnts, transpose
	
	if(ParamIsDefault(istart))
		istart = 0
	endif
	
	if(ParamIsDefault(iend))
		iend = DimSize(matrix,1)-1
	endif
	
	if(ParamIsDefault(show))
		show = 0
	endif
	
	if(!ParamIsDefault(transpose))
		MatrixTranspose matrix
	endif
	
	if(istart >= iend)
		Print "Range error!"
		return -1
	endif
	
	//Duplicate/FREE/R=[][istart,iend] data, matrix
	
	Variable nfits = iend-istart+1
	
	// The structure required by MPF2_DoMPFit()
	STRUCT MPFitInfoStruct MPStruct

	Variable i, err
	Variable numcoeffs = 0

	String saveDF = GetDataFolder(1)
	SetDataFolder dataFolderName
		// Look up the function type list saved by Multi-peak Fit
		SVAR SavedFunctionTypes
		MPStruct.ListOfFunctions = SavedFunctionTypes

		Wave/Z cw = 'Baseline Coefs'
		if (WaveExists(cw))
			MPStruct.ListOfCWaveNames = NameOfWave(cw)
			numcoeffs += numpnts(cw)
		endif
		MPStruct.ListOfCWaveNames += ";"
	
		Variable	npeaks = 0
		do
			Wave/Z cw = $("Peak "+num2istr(npeaks)+" Coefs")
			if (!WaveExists(cw))
				break
			endif
			MPStruct.ListOfCWaveNames += NameOfWave(cw)+";"
			numcoeffs += numpnts(cw)
			npeaks += 1
		while (1)
		
		NVAR XPointRangeBegin, XPointRangeEnd, XPointRangeReversed
		MPStruct.XPointRangeBegin = XPointRangeBegin
		MPStruct.XPointRangeEnd = XPointRangeEnd
	SetDataFolder saveDF
	MPStruct.NPeaks = npeaks

	if(ParamIsDefault(fitpnts))
		//		fitpnts = XPointRangeEnd-XPointRangeBegin-+1
		fitpnts = 256
	endif


	// fill in more members of the structure
	//MPStruct.XPointRangeBegin = 0
	// this function always fits the entire data set
	MPStruct.FitCurvePoints = fitpnts
	MPStruct.fitOptions = 0	// you might want to change this to 4
	
	// Create the list of coefficient waves and hold strings
	// This function won't actually use holds
	MPStruct.ListOfCWaveNames = ""
	MPStruct.ListOfHoldStrings = ""
	//MPStruct.FuncListString = "/D=fitc/R=resc"
	
	Variable nr = DimSize(matrix,0), nc = DimSize(matrix,1)
	Make/O/N=(nfits,numcoeffs) $(NameOfWave(matrix)+"_coeffs")
	WAVE coeffs = $(NameOfWave(matrix)+"_coeffs")

	Make/O/N=(fitpnts,nfits) $(NameOfWave(matrix)+"_fit")
	Make/O/N=(fitpnts,nfits) $(NameOfWave(matrix)+"_res")
	WAVE fitw = $(NameOfWave(matrix)+"_fit")
	WAVE resw = $(NameOfWave(matrix)+"_res")
	
	Variable lx = DimOffset(matrix,0)+DimDelta(matrix,0)*XPointRangeBegin
	Variable rx = DimOffset(matrix,0)+DimDelta(matrix,0)*XPointRangeEnd

	Variable lxo = DimOffset(matrix,0)
	Variable rxo = DimOffset(matrix,0)+DimDelta(matrix,0)*nr
	
	Variable ly = DimOffset(matrix,1)+DimDelta(matrix,1)*istart
	Variable ry = DimOffset(matrix,1)+DimDelta(matrix,1)*iend

	SetScale/I x,ly,ry,WaveUnits(matrix,1), coeffs
	SetScale/I x,lx,rx,WaveUnits(matrix,0), fitw
	SetScale/I x,lxo,rxo,WaveUnits(matrix,0), resw
	SetScale/I y,ly,ry,WaveUnits(matrix,1), fitw, resw
	SetScale d,0,0,WaveUnits(matrix,-1), fitw, resw

	// loop through each of the data sets, calling MPF2_DoMPFit
	// on each one. After each call, the results are checked and acted
	// on appropriately.
	String BatchDFName = "BatchPeakFit"
	DoWindow/K MultiPeakBatchFitWindow
	KillDataFolder/Z $BatchDFName
	DuplicateDataFolder $dataFolderName, $BatchDFName

	SetDataFolder BatchDFName	

	Make/D/O/N=(nr) tw = matrix[p][istart]
	CopyScales matrix, tw
	//	Make/D/O/N=(fitpnts) fit_tw, Res_tw
	//SetScale/I x,lx,rx,WaveUnits(matrix,0), fit_tw, Res_tw
	//	SetScale d,0,0,WaveUnits(matrix,-1), fit_tw, Res_tw
	Make/D/O/N=(1) fit_tw, Res_tw


	if (show)
		Display /N=MultiPeakBatchFitWindow tw,fit_tw as "MultiPeakBatchFit"
		AppendToGraph/L=Res_Left res_tw
		Execute/Q "MultiPeakBatchFitWindowStyle()"
	// 	DoAlert 0,"OK?"
	endif
	
	for (i = 0; i < nfits; i += 1)
			// Look up the function type list saved by Multi-peak Fit
//			SVAR SavedFunctionTypes
//			MPStruct.ListOfFunctions = SavedFunctionTypes

			Wave/Z cw = 'Baseline Coefs'
			if (WaveExists(cw))
				MPStruct.ListOfCWaveNames = NameOfWave(cw)
			endif
			MPStruct.ListOfCWaveNames += ";"
		
			npeaks = 0
			do
				Wave/Z cw = $("Peak "+num2istr(npeaks)+" Coefs")
				if (!WaveExists(cw))
					break
				endif
				MPStruct.ListOfCWaveNames += NameOfWave(cw)+";"
				npeaks += 1
			while (1)
			
			//KillWaves $(WaveList("fit_*",",",""))
			//KillWaves fit_tw, Res_tw
			//Print i, WaveList("fit_*",";",""); DelayUpdate
		
		MPStruct.NPeaks = npeaks
		//Print MPStruct.FitCurvePoints
		// Most of the contents of the structure don't need to change.
		// Naturally, we have to put in a wave reference to the
		// data set being fit this time through the loop
		tw = matrix[p][i+istart]
		Wave MPStruct.yWave = tw

		// the X wave for the current data set won't exist if the data
		// is waveform data.
		Wave/Z MPStruct.xWave = $""

		// just in case the new wave has a different number of points
		//MPStruct.XPointRangeEnd = numpnts(MPStruct.yWave)-1
		
		err = MPF2_DoMPFitLocal(MPStruct, saveDF+BatchDFName+":")//,doUpdates=1,doAutoDest=1,doAutoResid=1)

		if (err)
			// error return from MPF2_DoMPFit generally indicates
			// a programmer error
			DoAlert 0, "Error calling MPF2_DoMPFit: "+num2str(err)
			SetDataFolder saveDF
			return err
		endif
		if (MPStruct.fitQuitReason == 2)
			// if the user aborts a fit, we assume that the whole process
			// should be aborted
			DoAlert 0, "User aborted batch fit"
			SetDataFolder saveDF
			return -1
		endif
		if (MPStruct.fitError)
			// Something went wrong with the current fit, and it
			// failed to converge. We note this fact to the user via
			// an alert, then move on to the next one. You may wish
			// to do something more sophisticated than this.
			
			// Avoid a long line by concatenating the message
			// in small pieces
			String alertMsg = "Error doing fit to "
			alertMsg += "Fit "+num2str(i)+": "
			alertMsg += num2str(MPStruct.fitError)
			alertMsg += "; continuing with next fit."
			DoAlert 0, alertMsg
		endif
		
		Variable j = 0
		// Fill in coeffs matrix

			Wave/Z cw = 'Baseline Coefs'
			if (WaveExists(cw))
				coeffs[i][0,numpnts(cw)-1] = cw[q]
				j += numpnts(cw)
			endif
		
			npeaks = 0
			do
				Wave/Z cw = $("Peak "+num2istr(npeaks)+" Coefs")
				if (!WaveExists(cw))
					break
				endif
			//	coeffs[i][1+npeaks*3,1+(npeaks+1)*3] = cw[q-(1+npeaks*3)]
				coeffs[i][j,j+numpnts(cw)-1] = cw[q-j]
				npeaks += 1
				j += numpnts(cw)
			while (1)
			
			// Get Fit curve
			//WAVE fittw = fit_tw
			//WAVE restw = Res_tw

			if(XPointRangeReversed)
				Reverse fit_tw, Res_tw
			endif
			
			fitw[][i] = fit_tw[p]
			resw[][i] = Res_tw[p]
	endfor
	
	SetDataFolder saveDF


	if(show)
		Variable/G fitindex = i-1
		Variable/G i0 = istart
		String/G dataname = NameOfWave(matrix)
		Duplicate/O tw, dataw
		Duplicate/O resw, datares
		Duplicate/O fitw, datafit
		Setformula dataw, NameOfWave(matrix)+"[p][fitindex+"+num2str(istart)+"]"
		Setformula datares, NameOfWave(matrix)+"_res[p][fitindex]"
		Setformula datafit, NameOfWave(matrix)+"_fit[p][fitindex]"
		AppendToGraph dataw, datafit
		RemoveFromGraph tw, fit_tw
		AppendToGraph/L=Res_Left datares
		RemoveFromGraph res_tw
		Execute/P/Q "MultiPeakBatchFitWindowStyle()"	
//		ModifyGraph mode(dataw)=2
//		ModifyGraph marker(dataw)=8
//		ModifyGraph lSize(dataw)=1.5
//		ModifyGraph rgb(dataw)=(0,0,0)
//		ModifyGraph lblPos(left)=52,lblPos(Res_Left)=52
//		ModifyGraph freePos(Res_Left)=0
//		ModifyGraph axisEnab(left)={0,0.75}
//		ModifyGraph axisEnab(Res_Left)={0.8,1}
		SetVariable setfitindex pos={11,322},size={88,16},value=fitindex,format="%4u";DelayUpdate
		SetVariable setfitindex limits={0,(nfits-1),1}
		Button SpawnFitButton title="Spawn",pos={8,295},proc=SpawnFitButtonProc
	endif
	
	if(!ParamIsDefault(transpose))
		MatrixTranspose matrix
	endif

	
//	KillWaves/Z tw, fit_tw, res_tw
	return 0
end

Function SpawnFitButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			NVAR i = fitindex, i0 = i0
			SVAR dname = dataname
			WAVE tw = $(dname)
			WAVE fitw = $(dname+"_fit")
			WAVE resw = $(dname+"_res")
			Display /W=(34.8,42.2,430.2,249.8) tw[][i+i0], fitw[][i] as ("Fit " + dname + "[][" + num2str(i+i0) + "]");DoUpdate
			AppendToGraph/L=Res_Left resw[][i]	
			Execute/Q/P "MultiPeakBatchFitWindowStyle()"	
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Proc MultiPeakBatchFitWindowStyle() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode[0]=2
	ModifyGraph/Z marker[0]=8
	ModifyGraph/Z lSize[0]=1.5
	ModifyGraph/Z rgb[0]=(0,0,0)
	ModifyGraph/Z lblPos(left)=52,lblPos(Res_Left)=52
	ModifyGraph/Z freePos(Res_Left)=0
	ModifyGraph/Z axisEnab(left)={0,0.75}
	ModifyGraph/Z axisEnab(Res_Left)={0.8,1}
EndMacro

// Function to do a multi-peak fit from client code.
// Input: a structure with information about the fit, plus the name of a data folder containing coefficient waves,
// one for each peak, plus a coefficient wave for the baseline (unless the baseline function is "none").
// The contents of DataFolderName is a full path to the data folder, with final ":"
// Returns the function list string used in FuncFit sum-of-functions list.
//
// LOCAL: changed to make use of MPstruct.FitCurvePoints
Function MPF2_DoMPFitLocal(MPstruct, DataFolderName [, doUpdates, doAutoDest, doAutoResid])
	STRUCT MPFitInfoStruct &MPstruct
	String DataFolderName
	Variable doUpdates, doAutoDest, doAutoResid
	
	if (ParamIsDefault(doUpdates))
		doUpdates = 1
	endif
	
	if (ParamIsDefault(doAutoDest))
		doAutoDest = 1
	endif
	
	if (ParamIsDefault(doAutoResid))
		doAutoResid = 1
	endif
	
	Variable npeaks = MPstruct.NPeaks
	Wave yw = MPstruct.yWave
	Wave/Z xw = MPstruct.xWave
	
	if (ItemsInList(MPstruct.ListOfFunctions) != npeaks+1)		// +1 for the baseline function
		return MPF2_Err_BadNumberOfFunctions
	endif
	
	if (ItemsInList(MPstruct.ListOfCWaveNames) != npeaks+1)
		return MPF2_Err_BadNumberOfCWaves
	endif
	
	MPstruct.FuncListString = ""
	String holdString
	
	String BL_TypeName = StringFromList(0, MPstruct.ListOfFunctions)
	Variable doBaseLine = CmpStr(BL_TypeName, "None") != 0
	if (doBaseLine)
		String BL_FuncName
		Variable nBLParams
		
		FUNCREF MPF2_FuncInfoTemplate blinfo = $(BL_typename + BL_INFO_SUFFIX)
		BL_FuncName = blinfo(BLFuncInfo_BaselineFName)
		nBLParams = ItemsInList(blinfo(BLFuncInfo_ParamNames))
		if (nBLParams == 0)
			return MPF2_Err_NoSuchBLType
		endif
		
		STRUCT MPF2_BLFitStruct BLStruct
		if (WaveExists(xw))
			BLStruct.xStart = xw[MPstruct.XPointRangeBegin]
			BLStruct.xEnd = xw[MPstruct.XPointRangeEnd]
		else
			BLStruct.xStart = pnt2x(yw, MPstruct.XPointRangeBegin)
			BLStruct.xEnd = pnt2x(yw, MPstruct.XPointRangeEnd)
		endif
		String blcoefwname = StringFromList(0, MPstruct.ListOfCWaveNames)
		Wave/Z blcoefwave = $(DataFolderName+PossiblyQuoteName(blcoefwname))
		if (!WaveExists(blcoefwave))
			return MPF2_Err_BLCoefWaveNotFound
		endif
		MPstruct.FuncListString += "{"+BL_FuncName+", "+GetWavesDataFolder(blcoefwave,2)
Duplicate/O blcoefwave, blepswave
blepswave = 1e-6
MPstruct.FuncListString += ", EPSW="+GetWavesDataFolder(blepswave,2)
		holdString = StringFromList(0, MPstruct.ListOfHoldStrings)
		if (strlen(holdString) > 0)
			MPstruct.FuncListString += ", HOLD=\""+holdString+"\""
		endif
		MPstruct.FuncListString += ", STRC=BLStruct}"
	endif
	
	Variable i
	for (i = 0; i < nPeaks; i += 1)
		String PeakTypeName = StringFromList(i+1, MPstruct.ListOfFunctions)
		
		FUNCREF MPF2_FuncInfoTemplate infoFunc=$(PeakTypeName+PEAK_INFO_SUFFIX)
		String PeakFuncName = 	infoFunc(PeakFuncInfo_PeakFName)
		if (strlen(PeakFuncName) == 0)
			return MPF2_Err_NoSuchPeakType
		endif

		String pwname = StringFromList(i+1, MPstruct.ListOfCWaveNames)
		pwname = PossiblyQuoteName(pwname)
		pwname = DataFolderName + pwname
		Wave/Z coefw = $pwname
		if (!WaveExists(coefw))
			return MPF2_Err_PeakCoefWaveNotFound
		endif
		
		MPstruct.FuncListString += "{"+PeakFuncName+","+pwname
Duplicate/O coefw, $(NameOfWave(coefw)+"eps")
Wave epsw = $(NameOfWave(coefw)+"eps")
epsw = 1e-6
MPstruct.FuncListString += ", EPSW="+GetWavesDataFolder(epsw,2)
		holdString = StringFromList(i+1, MPstruct.ListOfHoldStrings)			// i+1 to account for the fact that the first hold string goes with the baseline
		if (strlen(holdString) > 0)
			MPstruct.FuncListString += ", HOLD=\""+holdString+"\""
		endif
		MPstruct.FuncListString += "}"
	endfor

	Variable V_FitQuitReason = 0
	Variable V_FitMaxIters=500
	Variable V_FitOptions=MPStruct.fitOptions

//print MPstruct.FuncListString
	MPstruct.fitErrorMsg = ""
	
//	MPF2_BackupCoefWaves(MPstruct.ListOfCWaveNames, DataFolderName)

	Variable errorCode=0
	DebuggerOptions
	Variable doDebugOnError = V_debugOnError
	DebuggerOptions debugOnError=0
	Variable fs = MPstruct.fitCurvePoints
	try
		FuncFit/Q=1/L=(fs)/N=(doUpdates==0?1:0)/M=2 {string=MPstruct.FuncListString} yw[MPstruct.XPointRangeBegin, MPstruct.XPointRangeEnd] /X=xw[MPstruct.XPointRangeBegin, MPstruct.XPointRangeEnd]/W=MPstruct.weightWave[MPstruct.XPointRangeBegin, MPstruct.XPointRangeEnd]/I=1/M=MPstruct.maskWave[MPstruct.XPointRangeBegin, MPstruct.XPointRangeEnd] /AD=(doAutoDest)/AR=(doAutoResid)/A=1/NWOK;AbortOnRTE
	catch
		MPstruct.fitErrorMsg = GetRTErrMessage()
		Variable semiPos = strsearch(MPstruct.fitErrorMsg, ";", 0)
		if (semiPos >= 0)
			MPstruct.fitErrorMsg = MPstruct.fitErrorMsg[semiPos+1, inf]
		endif
		errorCode = GetRTError(1)
	endtry
	DebuggerOptions debugOnError=doDebugOnError
	
	MPstruct.dateTimeOfFit = DateTime
	MPstruct.fitPnts = V_npnts
	MPstruct.chisq = V_chisq
	MPstruct.fitError = errorCode
	MPstruct.fitQuitReason = V_FitQuitReason
		
	return MPF2_Err_NoError
end

//#endif


Function Fit2PeaksInMatrix(matrix, p1, p2)
	// Fit two peaks in matrix
	Wave matrix
	Variable p1, p2
	
	Variable i, nr = DimSize(matrix, 0), nc = DimSize(matrix,1)

	Make/O/N=(nr,7) coeffs
	Make/O/N=(nc) tw
	SetScale/P x,DimOffset(matrix,0),DimDelta(matrix,0), coeffs
	SetScale/P x,DimOffset(matrix,1),DimDelta(matrix,1), tw
	
	Make/O/N=7 W_coef={0,32e-6,64e-6,4e-6,6e-6,20e-9,11e-9}
	For(i = 0; i < nr; i+=1)
		tw = matrix[i][p]
		FuncFit /Q/M=2/W=0  twolor, W_coef, tw[p1,p2]
		coeffs[i][] = W_coef[q]
	EndFor
	
	//Display coeffs
	//KillWaves tw
	Rename coeffs, $(Nameofwave(matrix) + "_coeffs")
End		

Function XYZtoMatrix2(mname,wx,wy,wz,nb,nx,ny)
	// for example call XYZtoMatrix("JPS31190",nw,vw,iw,500,1689-1190,512)
	// nb = number of points per curve
	// nx = total number of curves along x axis
	// ny = number of points in interpolation of y wave
	String mname
	WAVE wx, wy, wz
	Variable nb, nx, ny
	
	Make/O/N=(nx,ny) matrix
	WaveStats/Q wx
	SetScale/I x,V_min,V_max,WaveUnits(wx,-1), matrix
	WaveStats/Q wy
	SetScale/I y,V_min,V_max,WaveUnits(wy,-1), matrix
	SetScale d,0,0,WaveUnits(wz,-1), matrix
	
	Make/O/N=(nb) twy, twz
	Make/O/N=(ny) dwz
	SetScale/I x,V_min,V_max,WaveUnits(wy,-1),dwz
	
	Variable i, j
	j = 0
	For(i=0;i<nx;i+=1)
		twy = wy[p+j]
		twz = wz[p+j]

		Interpolate2/I=3/J=2/T=1/Y=dwz twy, twz
		
		matrix[i][] = dwz[q]
		j += nb
	EndFor
	
	//NewImage/F matrix
	
	KillWaves twy, twz, dwz
//	if (WaveExists($mname)) 
//		KillWaves $mname
//	endif
	Duplicate/O matrix $mname
End

Function XYZtoContourZ(wx,wy,wz, nx, ny)
	WAVE wx, wy, wz
	Variable nx, ny
	
	Display; AppendXYZContour wz vs {wx, wy}; DelayUpdate
	ModifyContour wz autoLevels={*,*,0}
	
	ModifyContour wz nullValue=NaN
	
	Make/O/N=(nx,ny) matrix

	WaveStats/Q wx
	SetScale/I x,V_min,V_max,WaveUnits(wx,-1), matrix
	WaveStats/Q wy
	SetScale/I y,V_min,V_max,WaveUnits(wy,-1), matrix
	matrix = ContourZ("","mv",0,x,y)
End

Function XYZ2Matrix(wx,wy,wz,xmin,xmax,nx,ymin,ymax,ny,neighbors,mname)	
	WAVE wx, wy, wz
	Variable xmin, xmax, nx, ymin, ymax, ny
	Variable neighbors
	String mname
	
	if( xmin == xmax)
		xmin = wavemin(wx)
		xmax = wavemax(wx)
	endif
	
	if( ymin == ymax)
		ymin = wavemin(wy)
		ymax = wavemax(wy)
	endif

	if(nx <= 0 || ny <= 0)
		return -1
	endif
		
	
	Make/O/N=(nx,ny) mhits = 0, m = 0
	
	Variable i, r, c	
	Variable npnts = min(numpnts(wx), min(numpnts(wy),numpnts(wz)))	
	
	For(i=0;i<npnts;i+=1)		
		r = floor((wx[i]-xmin)/(xmax-xmin)*nx)
		c = floor((wy[i]-ymin)/(ymax-ymin)*ny)
		
		if(r >= 0 && r < nx && c >= 0 && c < ny)
			mhits[r][c] += 1
			m[r][c] += wz[i]
		endif
	EndFor

	m = (mhits[p][q] > 0) ? m[p][q]/mhits[p][q] : NaN

	if(neighbors > 0)
		MatrixFilter/N=(neighbors) NanZapMedian m
	EndiF
	
	SetScale d,0,0,WaveUnits(wz,-1),m	
	SetScale/I x,xmin,xmax,WaveUnits(wx,-1),m
	SetScale/I y,ymin,ymax,WaveUnits(wy,-1),m	
	
	Duplicate/O m, $mname	
	KillWaves mhits,m
End

Function EliminateOutliers(wx,wy,[dx,dy])
	WAVE wx, wy
	Variable dx, dy
	
	If(ParamisDefault(dx))
		dx = 4e-6
	Endif
	
	If(ParamisDefault(dy))
		dy = 5e-9
	Endif
	
	Variable i, n = numpnts(wx)
	
	Duplicate/FREE wx, wxl
	Duplicate/FREE wy, wyl
	
	For(i=1;i<n-2;i+=1)
		if((abs(wxl[i]-wxl[i+1]) + abs(wxl[i]-wxl[i-1]) > dx) || (abs(wyl[i]-wyl[i+1]) + abs(wyl[i]-wyl[i-1]) > dy))
			wx[i] = NaN
			wy[i] = NaN
		endif
	Endfor
End


Function XYZ3Matrix(wx,wy,wz,mname,[xmin,xmax,nx,ymin,ymax,ny,neighbors,nanfactor]) //,thresh])	
	// XYZ3Matrix takes three matrices wx wy and wz and creates a third matrix mname (size nx x ny) whose
	// rows are taken from values of wx
	// columns are taken from values of wy
	// and z-values (color) from wz
	//
	// This function is handy for making IV maps as a function of a control parameter (flux, voltage, etc)
	//
	// Parameters in brackets [] are optional:
	// xmin, xmax: minimum and maximum values for the rows of mname.  any data in wx outside this range
	//             is discarded
	// ymin, ymax: same as above for wy
	// nx: number of rows of matrix mname
	// ny: number of columns of matrix mname
	// neighbors: number of neighboring elements in mname to include when averaging over nans (recommended 3)
	// nanfactor: how many neighboring nan elements to include in nanzapping (ask Caglar)
	// (thresh: if there are less than thresh points in a matrix element, this element is ignored (NaN)
	//				this is useful to eliminate isolated points, such as when switching.  
	//				NOTE: this functionality was replaced with EliminateOutliers function.)
	
	WAVE wx, wy, wz
	String mname
	Variable xmin, xmax, nx, ymin, ymax, ny
	Variable neighbors, nanfactor  //, thresh
	
	if(ParamIsDefault(xmin) || ParamIsDefault(xmax))
		xmin = wavemin(wx)
		xmax = wavemax(wx)
		//Print xmin, xmax
	endif
	
	if(ParamIsDefault(ymin) || ParamIsDefault(ymax))
		ymin = wavemin(wy)
		ymax = wavemax(wy)
	endif

	if(ParamIsDefault(nx))
		Duplicate/FREE wx, dwx
		Differentiate wx/D=dwx 
		dwx = abs(dwx)
		nx = abs(xmax-xmin)/wavemin(dwx)
		nx = min(nx,1024)
		Print "nx =", nx
	endif
		
	if(ParamIsDefault(ny))
		Duplicate/FREE wy, dwy
		Differentiate wy/D=dwy
		dwy = abs(dwy)
		ny = abs(ymax-ymin)/wavemin(dwy)
		ny = min(ny,1024)
		Print "ny =", ny
	endif

	Make/O/N=(nx,ny)/FREE mhits = 0, m = 0
	
	SetScale d,0,0,WaveUnits(wz,-1), m, mhits
	SetScale/I x,xmin,xmax,WaveUnits(wx,-1),m, mhits
	SetScale/I y,ymin,ymax,WaveUnits(wy,-1),m	, mhits
	
	Variable i, r, c
	Variable npnts = min(numpnts(wx), min(numpnts(wy),numpnts(wz)))	
	
//#if (IgorVersion() >= 6.3)
//	ImageFromXYZ {wx,wy,wz}, m, mhits
//	m /= mhits
//#else
	For(i=0;i<npnts;i+=1)		
		r =  (wx[i]-xmin)/(xmax-xmin)*(nx-1)
		c = (wy[i]-ymin)/(ymax-ymin)*(ny-1)
	
		if(r >= 0 && r < nx && c >= 0 && c < ny)
			mhits[r][c] += 1
			m[r][c] += wz[i]
		endif
	EndFor
	
m = (mhits[p][q] > 2) ? m[p][q]/mhits[p][q] : NaN
//#endif

//	m /= mhits
	
//	if(ParamIsDefault(nanoption))
//		nanoption = 1
//	endif

	if(!ParamIsDefault(nanfactor))
		if(ParamIsDefault(neighbors))
			neighbors = 3
		endif
		
		neighbors = (mod(neighbors,2) == 1) ? neighbors : neighbors+1
		
		Duplicate/FREE mhits, mhitimage
		Redimension/B/U mhitimage
		mhitimage = (mhits[p][q] > 0) ? 255 : 0
		ImageTransform/N={(2*nanfactor+2),(2*nanfactor+2)} padimage mhitimage
		WAVE/B/U tempM = M_PaddedImage
		ImageTransform/IOFF={nanfactor+1,nanfactor+1,0} offsetImage tempM
		WAVE/B/U tempM = M_offsetImage
		Make/FREE/B/U/N=(nanfactor,nanfactor) sel = 255
		//		sel = (((p-(nanfactor-1)/2)^2+(q-(nanfactor-1)/2)^2) <= ((nanfactor-1)/2)^2) ? sel[p][q] : 0
		// For circular neighborhood		
		ImageMorphology /S=sel closing tempM
		WAVE/B/U tempM = M_ImageMorph
		ImageTransform/IOFF={-nanfactor-1,-nanfactor-1,0} offsetImage tempM
		WAVE/B/U mhitsclosed = M_offsetImage
		Redimension/N=(DimSize(mhitimage,0),DimSize(mhitimage,1)) mhitsclosed
		MatrixFilter/N=(neighbors)/R=mhitsclosed NanZapMedian m
		
		KillWaves/Z M_ImageMorph, M_offsetImage, M_PaddedImage		
	else
		if(!ParamIsDefault(neighbors))
			MatrixFilter/N=(neighbors)/P=1 NanZapMedian m
		endif
	endif
	
	Duplicate/O m, $mname	
End

Function MakeImage(matrix,xmin,xmax)
	WAVE matrix
	Variable xmin, xmax
	
	SetScale/I x,xmin,xmax,matrix
	
	String wname = NameofWave(matrix)
	
	//DoWindow/K 
	Display/N=wname
	AppendImage /W=wname matrix 
	
	ColorScale /W=wname/A=RT widthPct=5,heightPct=20
	
	Variable logcolor = 0
	
	If (logcolor == 1)
		ModifyImage matrix log=1
	EndIf
End
	
Function AddTitle()
	String wavenames
	
	wavenames = RemoveEnding(TraceNameList("","\r",0))
	wavenames += RemoveEnding(ImageNameList("","\r"))
	TextBox/C/N=$("TB"+WinName(0,1))/A=RB wavenames
End

Function SetGainIV()
	
	WAVE vw = WaveRefIndexed("",0,2)
	// Voltage gain
	vw /= 1000
	
	WAVE iw = WaveRefIndexed("",0,1)
	iw /= 1.33e6
	SetScale d,0,0,"A",iw
End

Function DownSample(w,n)
	WAVE w // wave to downsample
	Variable n // number of points
	
	Variable m = numpnts(w), nlen = floor(m/n)
	
	if(n > m)
		return -1
	endif
	
	Make/FREE/N=(n) wout
	
	wout = mean(w,p*nlen,(p+1)*nlen)
	
	Duplicate/O wout, $(NameOfWave(w)+"ds")
	
	return 0
End

Function AddChargingLine(m, R)
	WAVE m
	Variable R
	
	Duplicate/O m, mr
	Variable dy = DimDelta(m,1)
	
	Variable i
	mr = ((m[p][q]-m[p][q-1])/dy-1/R) > 0 ? m[p][q] : mr[p][q-1]-1/R*dy

	Duplicate/O mr, $(NameOfWave(m)+"_R")
End


Function FitPeaksInMatrix(matrix, [p1, p2, rank])
	// Fits peaks in matrix
	Wave matrix
	Variable p1, p2, rank
	
	if(!ParamIsDefault(rank))
		MatrixTranspose matrix
	endif

	Variable i, nr = DimSize(matrix, 0), nc = DimSize(matrix,1)


	if(ParamIsDefault(p1))
		p1 = 0
	endif
	
	if(ParamIsDefault(p2))
		p2 = nr
	endif	
	
	Make/O/N=(nr,3)/FREE coeffs 
	coeffs[0][] = {0,100,15}
	Make/O/N=(nc)/FREE tw
	SetScale/P x,DimOffset(matrix,0),DimDelta(matrix,0), coeffs
	SetScale/P x,DimOffset(matrix,1),DimDelta(matrix,1), tw
	
	Make/O/N=3 W_coef = coeffs[0][p]
	For(i = 0; i < nr; i+=1)
		tw = matrix[i][p]
		CurveFit/Q/M=2/W=0/N/H="100" exp, kwCWave = W_coef, tw[p1,p2]
		coeffs[i][] = W_coef[q]
	EndFor
	
	//Display coeffs
	//KillWaves tw
	Duplicate/O coeffs, $(Nameofwave(matrix) + "_coeffs")
	
	If(!ParamIsDefault(rank))
		MatrixTranspose matrix
	endif
End		

Function XYtoImage(wx,wy,[mOut, nx, ny])
	WAVE wx, wy, mOut
	Variable nx, ny
	
	if(ParamIsDefault(mOut))
		WAVE mOut = $(NameOfWave(wy)+"_2D")
		
		if(WaveExists(mOut))
			mOut = 0
		else
			Variable wxmax = wavemax(wx)
			Variable wxmin = wavemin(wx)
			Variable wymax = wavemax(wy)
			Variable wymin = wavemin(wy)
		
			if(ParamIsDefault(nx))
				nx = (wxmax-wxmin)/(3.49*sqrt(variance(wx))*numpnts(wx)^(-1/3))
			endif
		
			if(ParamIsDefault(ny))
				ny = (wymax-wymin)/(3.49*sqrt(variance(wy))*numpnts(wy)^(-1/3))
			endif

			Make/O/N=(nx,ny) $(NameOfWave(wy)+"_2D") = 0
			WAVE mOut = $(NameOfWave(wy)+"_2D")
			SetScale/I x,wxmin,wxmax,WaveUnits(wx,-1),mOut
			SetScale/I y,wymin,wymax,WaveUnits(wy,-1),mOut
		endif
	
	endif
	
	Duplicate/FREE wx cnts
	cnts = 1
	
	Duplicate/FREE mOut, mCnts
	
	ImageFromXYZ {wx,wy,cnts}, mOut, mCnts
	
	return 0
End

Function FitPeakSeries(m)
	WAVE m
	
	Duplicate/O m, mt
	Redimension/N=(DimSize(m,0)) mt
	
	Make/D/O fitcoefs = {0,0,2e-6,-180e-3,8e-3}
	
	Duplicate/O fitcoefs fitc
	Redimension/N=(numpnts(fitcoefs),DimSize(m,1)) fitc
	SetScale/P y,DimOffset(m,1),DimDelta(m,1),"V",fitc
	SetScale d,0,0,"V",fitc
	Duplicate/O m, mfit
	
	Variable i, n = DimSize(m,1)
	For(i=0;i<n;i+=1)
		mt = m[p][i]
		FuncFit/Q/W=0 Gausslin, fitcoefs, mt/D
		fitc[][i] = fitcoefs[p]
		fitcoefs[0] = 0
		fitcoefs[1] = 0
		mfit[][i] = Gausslin(fitcoefs,x)
	EndFor
	
	fitc[4][] = abs(fitc[4][q])
	return 0	
End

Function CorrectPhase(wx, wy,wName)
	//
	// Takes lock-in X and Y data and corrects any error in phase setting
	//
	WAVE wx, wy
	String wName
	
	if (numpnts(wx) != numpnts(wy))
		Print "ERROR: X and Y waves unequal size."
		return -1
	endif
	
	Duplicate/O wy, wc, wt
	
	wt = atan(wy/wx)
	
	Make/N=1000/O wt_Hist;DelayUpdate
	Histogram/P/C/B=1 wt,wt_Hist
	CurveFit/Q/M=2/W=0 lor, wt_Hist/D
	
	WAVE W_coef = W_coef
	Variable phi = W_coef[2]
	printf "Angle is %.1f°\r", phi/pi*180
	
	wc = wx*cos(phi)+wy*sin(phi)
	
	Duplicate/O wc, $wName
	
	KillWaves wc, wt
	return 0
End

Function SubtractMin(m,c,ns)
	// Subtracts minimum value along rows or columns of a matrix
	// c is direction (0 rows; 1 columns)
	// ns is smoothing factor
	WAVE m
	Variable c
	Variable ns

	Variable i, mint
	
	if(c==0)
		Make/O/N=(DimSize(m,0)) tmp
		
		For(i=0;i<DimSize(m,1);i+=1)
			tmp = m[p][i]
			if(ns > 5)
				Smooth/S=2 25, tmp
			endif
			mint = WaveMin(tmp)
			m[][i] -= mint
		EndFor
	else
		Make/O/N=(DimSize(m,1)) tmp
		
		For(i=0;i<DimSize(m,0);i+=1)
			tmp = m[i][p]
			if(ns > 5)
				Smooth/S=2 25, tmp
			endif
			mint = WaveMin(tmp)
			m[i][] -= mint
		EndFor
	endif
	
	return 0
End


Function CorrectBiasResistor([voltage, current, resistance])
	WAVE voltage, current
	Variable resistance
	
	if(ParamIsDefault(resistance))
		NVAR/Z res = $(kPackagePath + ":rbias")
		if(NVAR_exists(res))
			resistance = res
		else
			Variable/G $(kPackagePath + ":rbias")
			NVAR res = $(kPackagePath + ":rbias")
			resistance = 50
		endif
		Prompt resistance, "Enter resistance:"
		DoPrompt "Correct Bias Resistor", resistance
		if(V_flag)
			return -1
		endif
		res = resistance
	endif
	
	if(ParamIsDefault(current) || ParamIsDefault(voltage))
		String wName = WinName(0,1)
		WAVE/Z current  = WaveRefIndexed(wName,0,1)
		WAVE/Z voltage = WaveRefIndexed(wName,0,2)
	endif	
		
	if(WaveExists(voltage) && WaveExists(current))
		voltage -= resistance*current
		return 0
	elseif(WaveExists(current))
		Duplicate/FREE current, voltage, newcurrent
		voltage = x
		voltage -= resistance*current
		Interpolate2/T=1/I=3/Y=newcurrent voltage, current
		current = newcurrent
	else
		Print "Error!"
		return -1
	endif
End

Function CDiGamma(xvalue,x0,y0,g1,g2)
	Variable xvalue, x0, y0, g1, g2
	
	if(xvalue > x0)
		return y0 + (1-y0)*((xvalue-x0)/(1-x0))^g2
	else
		return y0*(xvalue/x0)^g1
	endif
End

Function Initialize_GammaFy()
	DoWindow GammaFy
	if(V_Flag == 1)
		return 1
	else
		WAVE/Z gammplot = gammplot
		NVAR/Z x0 = x0
		if(!WaveExists(gammplot) || !NVAR_Exists(x0))
			Variable/G x0 = 0.25, y0=0.75, g1 = 0.5, g2 = 1.5
			Make/O/N=(1024) gammplot
			SetScale x,0,1,gammplot
			Execute/P/Q	"gammplot := CDiGamma(x,x0,y0,g1,g2)"
		endif
		return 0
	endif
End

Window GammaFy() : Graph
	PauseUpdate; Silent 1		// building window...
	
	if(Initialize_GammaFy())
		DoWindow/F GammaFy
		return -1
	endif
	
	Display /W=(80.25,394.25,597.75,648.5) gammplot
	ModifyGraph margin(left)=283
	ModifyGraph axOffset(left)=-0.571429
	SetAxis left 0,1
	ShowTools/A
	Slider g1,pos={189,33},size={50,285},limits={0,50,0},variable= g1,ticks= 5
	Slider g2,pos={267,33},size={50,285},limits={0,10,0},variable= g2,ticks= 5
	Slider x0,pos={15,39},size={56,279},limits={0,1,0},variable= x0,ticks= 5
	Slider y0,pos={106,39},size={56,279},limits={0,1,0},variable= y0,ticks= 5
	SetVariable x0setvar,pos={6,9},size={63,16}
	SetVariable x0setvar,limits={-inf,inf,0},value= x0
	SetVariable y0setvar,pos={93,9},size={63,16}
	SetVariable y0setvar,limits={-inf,inf,0},value= y0
	SetVariable g1setvar,pos={179,9},size={63,16}
	SetVariable g1setvar,limits={-inf,inf,0},value= g1
	SetVariable g2setvar,pos={262,8},size={63,16}
	SetVariable g2setvar,limits={-inf,inf,0},value= g2
	SetDrawLayer UserFront
EndMacro


Function CopyWaves([wName])
	String wName
	
	if(ParamIsDefault(wName))
		wName = WinName(0,1)
	endif
	
	if(!DataFolderExists("root:Packages:GQ:DataUtils"))
		NewDataFolder/O root:Packages
		NewDataFolder/O root:Packages:GQ
		NewDataFolder/O root:Packages:GQ:DataUtils
	endif
	
	String/G root:Packages:GQ:DataUtils:wName = wName
End

Function PasteWaves([wName])
	String wName
	
	if(ParamIsDefault(wName))
		wName = WinName(0,1)
	endif
	
	if(!DataFolderExists("root:Packages:GQ:DataUtils"))
		Print "PasteWaves ERROR: no waves to copy"
		return -1
	endif
	
	SVAR oldwName = root:Packages:GQ:DataUtils:wName
	
	if(strlen(WinList(oldwName,"","")) == 0)
		Print "PasteWaves ERROR: source graph not open"
		return -1
	endif
	
	if(strlen(CsrInfo(A,oldwName)) > 0)
		WAVE wy = CsrWaveRef(A,oldwName)
		
		if(WaveExists(CsrXWaveRef(A,oldwName)))
			WAVE wx = 	CsrXWaveRef(A,oldwName)
			AppendToGraph/W=$wName wy vs wx
		else
			AppendToGraph/W=$wName wy
		endif
	else
		String traceList = TraceNameList(oldwName,";",5), traceName
		Variable i, n = ItemsInList(traceList,";")
		For(i=0;i<n;i+=1)
			traceName = StringFromList(i,traceList)
			WAVE wy = TraceNameToWaveRef(oldwName, traceName)
		
			if(WaveExists(XWaveRefFromTrace(oldwName, traceName)))
				WAVE wx = 	XWaveRefFromTrace(oldwName, traceName)
				AppendToGraph/W=$wName wy vs wx
			else
				AppendToGraph/W=$wName wy
			endif
		EndFor
	endif
	return 0
End

//Function/S WaveToString(w)
//	WAVE w
//	String wavestr = NameOfWave(w) + " = {" + num2str(w[0])
//	Variable i
//	For(i=1;i<numpnts(w);i+=1)
//		wavestr += "," + num2str(w[i])
//	EndFor
//	wavestr += "}"
////	TextBox/C/N=text0/F=0/A=MC taustr
//	
//	return wavestr
//End

Function/S WaveToString(w)
	WAVE w
	String wavestr = ""
	Variable i, n = numpnts(w)
	For(i=0;i<n;i+=1)
		wavestr += ";" + num2str(w[i])
	EndFor
//	TextBox/C/N=text0/F=0/A=MC taustr
	
	return wavestr
End


// BEGIN FitFunc
//
// 
//

Function PlasmaRes(w,x) : FitFunc
	Wave w
	Variable x
	
	Variable phi0 = 3.291E-16 // in volts/angular freq
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = Lp*C*I0*x
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = Lp
	//CurveFitDialog/ w[1] = C
	//CurveFitDialog/ w[2] = I0
	//CurveFitDialog/ w[3] = xperiod
	//CurveFitDialog/ w[4] = xoffset
	
	return phi0/sqrt(w[0]*w[1])*sqrt(1+abs(w[0]*w[2]/phi0*cos(pi*(x-w[4])/w[3])))
End

Function twolor(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0+A1/(((x-x1)/w1)^2+1)+A2/(((x-x2)/w2)^2+1)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 7
	//CurveFitDialog/ w[0] = y0
	//CurveFitDialog/ w[1] = x1
	//CurveFitDialog/ w[2] = x2
	//CurveFitDialog/ w[3] = w1
	//CurveFitDialog/ w[4] = w2
	//CurveFitDialog/ w[5] = A1
	//CurveFitDialog/ w[6] = A2

	return w[0]+w[5]/(((x-w[1])/w[3])^2+1)+w[6]/(((x-w[2])/w[4])^2+1)
End

Function Gausslin(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0+x*y1+A+x0+width
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = y0
	//CurveFitDialog/ w[1] = y1
	//CurveFitDialog/ w[2] = A
	//CurveFitDialog/ w[3] = x0
	//CurveFitDialog/ w[4] = width

	return w[0]+x*w[1]+w[2]*exp(-((x-w[3])/w[4])^2)
End

//
//
//
// END FitFunc




// BEGIN Panel Controls
//
//
//

Function ListBoxProc(lba) : ListBoxControl
	STRUCT WMListboxAction &lba

	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave
//	Print lba.eventCode
	String datadf = "root:Data"
	if(!DataFolderExists(datadf))
		listWave = ""
		return -1
	endif

	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 1: // mouse down
			Variable numdfs = CountObjects(datadf,4)
			Redimension/N=(numdfs) listWave
			listWave = GetIndexedObjName(datadf,4,p)
			Sort/R listWave, listWave
			//Print listWave, selWave
			break
		case 3: // double click
			break
		case 4: // cell selection
		case 5: // cell selection plus shift key
			ClearGraph("DataLoadPanel#G0")
			if(row < numpnts(listWave))		
				GraphAllInDataFolder(datadf+":"+listWave[row],"DataLoadPanel#G0",wtypes=1)
			endif
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
		case 13: // checkbox clicked (Igor 6.2 or later)
			break
	endswitch

	return 0
End


Function DataLoaderGraphData(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	DFREF saveDFR = GetDataFolderDFR()

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			ControlInfo DataList
			
			Variable row = V_Value
			WAVE/T listWave = $(S_DataFolder+S_Value)
			
			String dfstr = "root:Data:"+listWave[row]
			String winrecpath = dfstr+":winrec"
			
			
			if(exists(winrecpath)==2)
				Execute/P/Q ("Execute/Q " + winrecpath)
			else
				GraphAllInDataFolder(dfstr, listwave[row], copy=V_Value)
			endif

// TODO Copycheck copy waves to root!
//			ControlInfo CopyCheck
			
//			SetDataFolder dfstr
//			Duplicate 
			
			break
		case -1: // control being killed
			break
	endswitch
	
	SetDataFolder saveDFR

	return 0
End

Function SubtractButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	String dfSav = GetDataFolder(1)
	SetDataFolder root:Packages:GQ:Base
	switch( ba.eventCode )
		case 2: // mouse up
			NVAR iref, ifile
			SVAR filePrefix, filePathStr
			String wnameref
			sprintf wnameref, "root:%s%04d", filePrefix, iref
			String wnamerefY=wnameref+"_2"
			String wnamerefX=wnameref+"_1"
			String wnamefile
			sprintf wnamefile, "root:%s%04d", filePrefix, ifile
			String wnamefileY=wnamefile+"_2"
			String wnamefileX=wnamefile+"_1"
			
			String filename
			sprintf filename, "%s%04dSUB", filePrefix, ifile
			
			WAVE wxb = $wnamerefX
			SubtractData($wnamerefY,$wnamerefX,$wnamefileY,$wnamefileX,1024,filename,filePathStr)
			break
		case -1: // control being killed
			break
	endswitch
	SetDataFolder dfSav
	return 0
End

Function LoadFile()
String dfSav = GetDataFolder(1)
	SetDataFolder root:Packages:GQ:Base
	Variable err = 0
	NVAR istart, xscale1, yscale1,addtograph,xcol,ycol, rbias, correctr
	SVAR filePrefix, filePathStr, xunits, yunits
	String wname
	sprintf wname, "%s%05d", filePrefix, istart
	String fname
	sprintf fname, "%s\\%s%05d.dat", filePathStr, filePrefix, istart
			
	GetFileFolderInfo/Q/Z=1 fname
	if(V_flag != 0)
		sprintf fname, "%s\\%s%05d.itx", filePathStr, filePrefix, istart
		GetFileFolderInfo/Q/Z=1 fname
		if(V_flag == 0)
			// do itx load
			LoadWave/T/O fname
			DoUpdate
			if (correctr == 1)
				WAVE wy = WaveRefIndexed("",0,1)
				WAVE wx = WaveRefIndexed("",0,2)
				wx -= wy*rbias
			endif
		endif
	endif

	LoadWave/N=tempM/M/J/D/O/K=1 fname

	WAVE m = tempM0
	Make/FREE/O/N=(dimSize(m,0)) wx, wy
	wx = m[p][xcol-1]
	wy = m[p][ycol-1]
						
	SetScale d 0,0,xunits,wx
	wx /= xscale1
	SetScale d 0,0,yunits,wy
	wy /= yscale1
			
	if (correctr == 1)
		wx -= wy*rbias
	endif
			
	NewDataFolder/O/S root:Data:$(wname)
	Duplicate /O wx, $(wname+"_"+num2str(xcol))
	Duplicate /O wy, $(wname+"_"+num2str(ycol))
	WAVE wx = $(wname+"_"+num2str(xcol))
	WAVE wy = $(wname+"_"+num2str(ycol))
			
	if (addtograph == 1 && (ItemsInList(WinList("*",";","WIN:1")) != 0))
		AppendToGraph wy vs wx
	else
		Display/N=$wname wy vs wx
	endIf
			
	String/G winrec = WinRecreation("",0)
	SetDataFolder dfSav
End

Function LoadFileButton(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	String dfSav = GetDataFolder(1)
	SetDataFolder root:Packages:GQ:Base
	Variable err = 0
	switch( ba.eventCode )
		case 2: // mouse up
			LoadFile()
		case -1: // control being killed
			break
	endswitch

	SetDataFolder dfSav
	return err
End

Function LoadButtonProc2(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	String dfSav = GetDataFolder(1)
	SetDataFolder root:Packages:GQ:Base
	switch( ba.eventCode )
		case 2: // mouse up
			NVAR firstpoint, lastpoint, istart, iend, npnts, Vstartslow, Vfinalslow, yscale, zscale,xcol, ycol,nptsslow,xstart,xstep
			NVAR slowaxis
			SVAR filePrefix, filePathStr, xunits, yunits, zunits, activeMatrix, extension
			xstart = Vstartslow
			xstep = (Vfinalslow-Vstartslow)/(nptsslow-1)

			String wname
			sprintf wname, "%s%04d_%d%d", filePrefix, istart, ycol, xcol
			
			Make/O dm
			LoadToMatrix(dm, filePrefix, filePathStr, istart, iend,extension=extension)
			
			Variable i = 0
			For(i=(nptsslow-1);i>-1;i-=1)
				DeletePoints/M=0 i*npnts,firstpoint,dm
				DeletePoints/M=0 (i*npnts+lastpoint+1),(npnts-lastpoint-1),dm
			EndFor
			Print (npnts-lastpoint-1)
			Make/O/N=(DimSize(dm,0)) yw, zw
			Print DimSize(dm,0)/nptsslow
			zw = dm[p][ycol-1]
			yw = dm[p][xcol-1]
			
			zw /= zscale
			yw /= yscale
			
			Redimension/N=((lastpoint-firstpoint+1),numpnts(zw)/(lastpoint-firstpoint+1)) zw
			SetScale/P y,xstart,xstep,xunits,zw
			SetScale/I x,yw[0],yw[lastpoint-firstpoint],yunits,zw
			
			if (slowaxis)
				MatrixTranspose zw
			endif
			
			SetScale d,0,0,zunits,zw

			Duplicate/O zw, m, root:$wname
			
			sprintf activeMatrix, "root:%s", wname
	
			break
		case -1: // control being killed
			break
	endswitch
	SetDataFolder dfSav

	return 0
End




Function SpawnButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	String dfSav = GetDataFolder(1)
	SetDataFolder root:Packages:GQ:Base

	switch( ba.eventCode )
		case 2: // mouse up
			SVAR activeMatrix
			WAVE m = $(activeMatrix)
			Display; AppendImage m
			TextBox/A=LT/N=title0 NameOfWave(m)
			//ModifyGraph width={Aspect,0.33}
			//ColorScale/N=colorscale0/A=RT heightPct=10
			// click code here
			break
		case -1: // control being killed
			break
	endswitch
	SetDataFolder dfSav

	return 0
End

Function SetScaling(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	String dfSav = GetDataFolder(1)
	SetDataFolder root:Packages:GQ:Base

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
			String scale_var = sva.Vname
			String cmd
			Variable nsv = strlen(scale_var)-4
			scale_var = scale_var[0,nsv]
			sprintf cmd, "%s = %s; %s = num2str(%s)", scale_var, sva.sval, sva.Vname, scale_var
			Execute cmd
			//	sva.sval = num2str(yscale)
			break
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			break
		case -1: // control being killed
			break
	endswitch

	SetDataFolder dfSav
	return 0
End

//
//
//
// END Panel Controls

// BEGIN Panel 
//
//
//

Window LoadFilePanel() : Panel
	PauseUpdate; Silent 1		// building window...
	String savDF= GetDataFolder(1)
	Initialize_Base(0)
	SetDataFolder root:Packages:GQ:Base
	NewPanel /W=(1005,98,1582,267) as "Load file"
	SetVariable filePath,pos={28,14},size={524,16},title="File Path"
	SetVariable filePath,limits={-inf,inf,0},value=filePathStr
	SetVariable filePrefix,pos={28,36},size={104,16},title="Data Series"
	SetVariable filePrefix,value=filePrefix
	Button LoadMatrix,pos={462,87},size={84,35},proc=LoadFileButton,title="Load"
	SetVariable istart,pos={141,35},size={112,16},title="File Index",format="%04d"
	SetVariable istart,limits={0,inf,1},value=istart
	SetVariable xscale,pos={28,79},size={126,16},proc=SetScaling,title="X-axis scaling"
	SetVariable xscale,limits={-inf,inf,0},value=xscale1str
	SetVariable yscale,pos={28,102},size={130,16},proc=SetScaling,title="Y-axis scaling"
	SetVariable yscale,limits={-inf,inf,0},value=yscale1str
	SetVariable xunits,pos={171,79},size={69,16},title="X Units"
	SetVariable xunits,value=xunits
	SetVariable yunits,pos={171,102},size={69,16},title="Y Units"
	SetVariable yunits,value=yunits
	CheckBox Add_to_active_graph,pos={446,70},size={111,14},title="Add to active graph"
	CheckBox Add_to_active_graph,variable=addtograph
	SetVariable Xcol,pos={256,79},size={112,16},title="X-column",format="%d"
	SetVariable Xcol,limits={1,inf,1},value=xcol
	SetVariable Ycol,pos={255,101},size={112,16},title="Y-column",format="%d"
	SetVariable Ycol,limits={1,inf,1},value=ycol
	CheckBox UseBiasRes,pos={167,140},size={122,14},title="Correct Bias Resistor?"
	CheckBox UseBiasRes,variable=correctr
	SetVariable BiasRes,pos={34,138},size={119,16},title="Bias Resistor"
	SetVariable BiasRes,limits={0,inf,0},value=rbias
	SetDataFolder savDF
EndMacro

//Window DataLoadPanel() : Panel
//	PauseUpdate; Silent 1		// building window...
//	Initialize_Base(0)
//	NewPanel /W=(739,138,1216,450) as "Data Loader"
//	SetDrawLayer UserBack
//	DrawText 150,33,"Data Preview"
//	DrawText 20,34,"root:Data:"
//	ListBox DataList,pos={16,41},size={113,218},proc=ListBoxProc
//	ListBox DataList,listWave=root:Packages:GQ:Base:datafolderlist,mode= 1,selRow= 0
//	Button DataLoaderB,pos={15,270},size={116,30},proc=DataLoaderGraphData,title="Graph Data"
//	Display/W=(150,41,458,299)/HOST=#
//	RenameWindow #,G0
//	SetActiveSubwindow ##
//EndMacro

Window DataLoadPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	Initialize_Base(0)
	NewPanel /W=(1084,137,1561,473) as "Data Loader"
	ShowTools/A
	SetDrawLayer UserBack
	DrawText 150,33,"Data Preview"
	DrawText 20,34,"root:Data:"
	ListBox DataList,pos={16,41},size={113,218},proc=ListBoxProc
	ListBox DataList,listWave=root:Packages:GQ:Base:datafolderlist,mode= 1,selRow= 0
	Button DataLoaderB,pos={15,270},size={116,30},proc=DataLoaderGraphData,title="Graph Data"
	CheckBox CopyCheck,pos={28,309},size={94,14},title="Copy to \"root:\"?",value= 0
	Display/W=(150,41,458,326)/HOST=# 
	RenameWindow #,G0
	SetActiveSubwindow ##
EndMacro

Window LoadMatrixPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1207,112,1746,733) as "Load Matrix"
	SetVariable filePath,pos={50,14},size={338,16},title="File Path"
	SetVariable filePath,limits={-inf,inf,0},value= root:Packages:GQ:Base:filePathStr
	SetVariable filePrefix,pos={21,61},size={136,16},title="Data Identifier"
	SetVariable filePrefix,help={"Format specifier for filenames, see help for sprintf for explanation"}
	SetVariable filePrefix,value= root:Packages:GQ:Base:filePrefix
	Button LoadMatrix,pos={30,527},size={133,33},proc=LoadButtonProc,title="Load"
	SetVariable istart,pos={35,119},size={112,16},title="Initial Index"
	SetVariable istart,format="%05d"
	SetVariable istart,limits={0,99999,1},value= root:Packages:GQ:Base:istart
	SetVariable iend,pos={38,142},size={112,16},title="Final Index",format="%05d"
	SetVariable iend,limits={0,99999,1},value= root:Packages:GQ:Base:iend
	SetVariable Vstartslow,pos={17,213},size={138,16},title="Slow ramp start"
	SetVariable Vstartslow,format="%g"
	SetVariable Vstartslow,limits={-inf,inf,0},value= root:Packages:GQ:Base:Vstartslow
	SetVariable Vfinalslow,pos={14,236},size={139,16},title="Slow Ramp End"
	SetVariable Vfinalslow,format="%g"
	SetVariable Vfinalslow,limits={-inf,inf,0},value= root:Packages:GQ:Base:Vfinalslow
	SetVariable yscale,pos={19,460},size={168,16},proc=SetScaling,title="Y-axis scaling"
	SetVariable yscale,limits={-inf,inf,0},value= root:Packages:GQ:Base:yscalestr
	SetVariable zscale,pos={20,481},size={167,16},proc=SetScaling,title="Z-axis scaling"
	SetVariable zscale,limits={-inf,inf,0},value= root:Packages:GQ:Base:zscalestr
	SetVariable xunits,pos={50,382},size={69,16},title="X Units"
	SetVariable xunits,value= root:Packages:GQ:Base:xunits
	SetVariable yunits,pos={50,404},size={69,16},title="Y Units"
	SetVariable yunits,value= root:Packages:GQ:Base:yunits
	SetVariable zunits,pos={50,425},size={69,16},title="Z Units"
	SetVariable zunits,value= root:Packages:GQ:Base:zunits
	SetVariable npnts,pos={39,166},size={112,16},title="No. Points",format="%d"
	SetVariable npnts,limits={0,inf,1},value= root:Packages:GQ:Base:npnts
	Button spawn_image,pos={31,576},size={130,33},proc=SpawnButtonProc,title="Spawn Image"
	SetVariable npntsslow,pos={19,257},size={112,16},title="No. Points slow"
	SetVariable npntsslow,format="%d"
	SetVariable npntsslow,limits={0,inf,1},value= root:Packages:GQ:Base:nptsslow
	SetVariable Xcol,pos={31,186},size={92,16},title="Fast Column",format="%d"
	SetVariable Xcol,limits={1,inf,1},value= root:Packages:GQ:Base:xcol
	SetVariable Ycol,pos={32,327},size={112,16},title="Z-column",format="%d"
	SetVariable Ycol,limits={1,inf,1},value= root:Packages:GQ:Base:ycol
	CheckBox SlowAxis,pos={30,279},size={82,14},title="Slow axis is X"
	CheckBox SlowAxis,variable= root:Packages:GQ:Base:slowaxis,side= 1
	SetVariable SkipLines,pos={19,99},size={118,16},title="No. Lines Skip"
	SetVariable SkipLines,help={"Number of lines to skip when loading data file"}
	SetVariable SkipLines,limits={0,inf,1},value= root:Packages:GQ:Base:skipLines
	SetVariable extension2,pos={42,79},size={75,16},title="extension"
	SetVariable extension2,value= root:Packages:GQ:Base:extension
	SetWindow kwTopWin,userdata(ACL_desktopNum)=  "17"
	Display/W=(206,50,474,633)/HOST=# 
	AppendImage :Packages:GQ:Base:m
	ModifyImage m ctab= {*,*,Grays,0}
	ModifyImage m ctabAutoscale=1
	ModifyGraph mirror=2
	RenameWindow #,G0
	SetActiveSubwindow ##
EndMacro

Window Subtractor() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(705,247,1272,334) as "Subtract reference from file and save"
	String savDF= GetDataFolder(1)
	Initialize_Base(0)

	SetDataFolder root:Packages:GQ:Base

	SetVariable filePath,pos={28,14},size={524,16},title="File Path"
	SetVariable filePath,limits={-inf,inf,0},value=filePathStr
	SetVariable filePrefix,pos={28,36},size={104,16},title="Data Series"
	SetVariable filePrefix,value=filePrefix
	SetVariable iref,pos={140,35},size={168,16},title="Reference File Index"
	SetVariable iref,format="%04d"
	SetVariable iref,limits={0,inf,1},value=iref
	SetVariable istart1,pos={195,57},size={112,16},title="File Index",format="%04d"
	SetVariable istart1,limits={0,inf,1},value=ifile
	Button subtract,pos={323,35},size={106,40},proc=SubtractButtonProc,title="Subtract"
	SetDataFolder savDF
EndMacro


Window LoadMatrixPanelPartial() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1038,305,1544,957) as "Load Matrix Partial"
	String savDF = GetDataFolder(1)
	Initialize_Base(0)

	SetDataFolder root:Packages:GQ:Base
	
	SetVariable filePath,pos={50,14},size={338,16},title="File Path"
	SetVariable filePath,limits={-inf,inf,0},value=filePathStr
	SetVariable filePrefix,pos={34,64},size={104,16},title="Data Series"
	SetVariable filePrefix,value=filePrefix
	Button LoadMatrix,pos={30,529},size={133,33},proc=LoadButtonProc2,title="Load"
	SetVariable istart,pos={35,91},size={112,16},title="Initial Index",format="%04d"
	SetVariable istart,limits={0,10000,1},value=istart
	SetVariable iend,pos={38,114},size={112,16},title="Final Index",format="%04d"
	SetVariable iend,limits={0,10000,1},value=iend
	SetVariable Vstartslow,pos={17,213},size={138,16},title="Slow ramp start"
	SetVariable Vstartslow,format="%g"
	SetVariable Vstartslow,limits={-inf,inf,0},value=Vstartslow
	SetVariable Vfinalslow,pos={14,236},size={139,16},title="Slow Ramp End"
	SetVariable Vfinalslow,format="%g"
	SetVariable Vfinalslow,limits={-inf,inf,0},value=Vfinalslow
	SetVariable yscale,pos={19,460},size={168,16},proc=SetScaling,title="Y-axis scaling"
	SetVariable yscale,limits={-inf,inf,0},value=yscalestr
	SetVariable zscale,pos={20,481},size={167,16},proc=SetScaling,title="Z-axis scaling"
	SetVariable zscale,limits={-inf,inf,0},value=zscalestr
	SetVariable xunits,pos={50,382},size={69,16},title="X Units"
	SetVariable xunits,value=xunits
	SetVariable yunits,pos={50,404},size={69,16},title="Y Units"
	SetVariable yunits,value=yunits
	SetVariable zunits,pos={50,425},size={69,16},title="Z Units"
	SetVariable zunits,value=zunits
	SetVariable npnts,pos={39,138},size={112,16},title="No. Points",format="%d"
	SetVariable npnts,limits={0,inf,1},value=npnts
	Button spawn_image,pos={31,576},size={130,33},proc=SpawnButtonProc,title="Spawn Image"
	SetVariable npntsslow,pos={19,257},size={112,16},title="No. Points slow"
	SetVariable npntsslow,format="%d"
	SetVariable npntsslow,limits={0,inf,1},value=nptsslow
	SetVariable Xcol,pos={37,176},size={92,16},title="Fast Column",format="%d"
	SetVariable Xcol,limits={1,inf,1},value=xcol
	SetVariable Ycol,pos={32,327},size={112,16},title="Z-column",format="%d"
	SetVariable Ycol,limits={1,inf,1},value=ycol
	CheckBox SlowAxis,pos={30,279},size={82,14},title="Slow axis is X"
	CheckBox SlowAxis,variable=slowaxis,side= 1
	SetVariable npnts1,pos={13,155},size={95,16},title="First point",format="%04d"
	SetVariable npnts1,limits={0,inf,1},value=firstpoint
	SetVariable npnts2,pos={107,154},size={95,16},title="Last point",format="%04d"
	SetVariable npnts2,limits={0,inf,1},value=lastpoint
	Display/W=(206,50,474,633)/HOST=# 
	AppendImage root:Packages:GQ:Base:m
	ModifyImage m ctab= {*,*,Grays,0}
	ModifyImage m ctabAutoscale=1
	ModifyGraph mirror=2
	RenameWindow #,G0
	SetActiveSubwindow ##
	SetDataFolder savDF
EndMacro

//
//
//
// END Panel

// BEGIN GraphStyle
//
//
//

Proc IVStyle() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z width={Aspect,1}
	ModifyGraph/Z mode=2
	ModifyGraph/Z lSize=1.5
	ModifyGraph/Z rgb[0]=(0,26112,39168),rgb[1]=(0,0,0)
	ModifyGraph/Z grid=2
	ModifyGraph/Z gridRGB=(56576,56576,56576)
	Label/Z left "\\Z12\\f01\\U"
	Label/Z bottom "\\Z12\\f01\\U"
	SetAxis/Z/A/N=1/E=2 left
	SetAxis/Z/A/N=1/E=2 bottom
EndMacro

Proc FluxMap() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z width=283.465,height=425.197,gbRGB=(65280,48896,48896)
//	ModifyGraph/Z lSize[10]=1.5
//	ModifyGraph/Z lStyle[10]=8
//	ModifyGraph/Z rgb[3]=(0,0,0),rgb[4]=(0,0,0),rgb[5]=(0,0,0),rgb[6]=(0,0,0),rgb[7]=(0,0,0)
//	ModifyGraph/Z rgb[8]=(0,0,0),rgb[9]=(0,0,0),rgb[10]=(0,60652,60652)
	ModifyGraph/Z mirror=2
	Label/Z left "V\\BJJ\\M (\\U)"
	Label/Z bottom "V\\Bbob\\M (\\U)"
	SetAxis/Z left 0,0.00018
	
	//SetAxis/Z bottom -1.27703399950125,1.38537324738155
End

Proc FitStyle() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode=4
	ModifyGraph/Z marker[1]=8
	ModifyGraph/Z rgb[1]=(0,0,62976)
EndMacro

Proc WaterfallStyle() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z rgb=(0,0,0)
	ModifyGraph/Z negRGB=(0,0,65535)
	ModifyGraph/Z lblMargin(left)=9,lblMargin(right)=52
	ModifyGraph/Z lblLatPos(left)=4,lblLatPos(right)=21
	ModifyGraph/Z lblRot(right)=90
	Label/Z left "I\\BJJ\\M (\\U)"
	Label/Z bottom "Reduced Flux (2\\F'Arial Greek'ðÖ/Ö\\B0\\M)"
	Label/Z right "V\\BJJ\\M (\\U)"
EndMacro

Proc S21withPhase() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z rgb[0]=(0,0,65280),rgb[1]=(0,0,65280)
	ModifyGraph/Z grid(left)=1,grid(bottom)=1
	ModifyGraph/Z minor(bottom)=1
	ModifyGraph/Z gridRGB(left)=(56576,56576,56576),gridRGB(bottom)=(56576,56576,56576)
	ModifyGraph/Z lblPos(left)=48,lblPos(L1)=48
	ModifyGraph/Z freePos(L1)={0,bottom}
	ModifyGraph/Z axisEnab(left)={0.5,1}
	ModifyGraph/Z axisEnab(L1)={0,0.5}
	Label/Z left "S21 (\\U)"
	Label/Z bottom "Frequency (\\U)"
	Label/Z L1 "Phase (\\U)"
EndMacro

Proc IVMap() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z margin(left)=28,margin(bottom)=28,margin(top)=14,margin(right)=14
	ModifyGraph/Z width=113.386,height={Aspect,4}
	ModifyGraph/Z rgb=(1,4,52428)
	ModifyGraph/Z live=1
	ModifyGraph/Z quickdrag=1
	ModifyGraph/Z mirror=2
	ModifyGraph/Z nticks(left)=10,nticks(bottom)=4
	ModifyGraph/Z minor=1
	ModifyGraph/Z fSize=8
	ModifyGraph/Z standoff=0
	ModifyGraph/Z tkLblRot(left)=90
	ModifyGraph/Z btLen=3
	ModifyGraph/Z tlOffset=-2
EndMacro

Proc StandardIV1() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z width=170.079,height=170.079
	ModifyGraph/Z grid=1
	ModifyGraph/Z tick=2
	ModifyGraph/Z mirror=1
	ModifyGraph/Z axOffset(left)=-1.14286
	Label/Z left "Current (\\U)"
	Label/Z bottom "Voltage (\\U)"
EndMacro


Proc PresentationStyle() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z margin(top)=14,margin(right)=14,expand=-0.667,gFont="Arial",gfSize=16
	ModifyGraph/Z height={Aspect,1.618}
	ModifyGraph/Z mirror=2
	Label/Z left "\\Z202eV\\BJ\\M\\u#2"
	Label/Z bottom "\\Z20\\F'Arial Greek'ö"
	SetAxis/Z left 0,0.00018
	SetAxis/Z bottom 0.0249999999999999,2.725
EndMacro

//
//
//
// END GraphStyle





Function LoadButtonFunc()
	SetDataFolder root:Packages:GQ:Base
	
	NVAR istart, iend, npnts, Vstartslow, Vfinalslow, yscale, zscale,xcol, ycol,nptsslow,xstart,xstep,skipLines
	NVAR slowaxis
	SVAR filePrefix, filePathStr, xunits, yunits, zunits, activeMatrix, extension
			
	xstart = Vstartslow
	xstep = (Vfinalslow-Vstartslow)/(nptsslow-1)

	String wname, wname2
			
	Make/FREE/O dm
	LoadToMatrix(dm,filePrefix, filePathStr, istart, iend,skipLines=skipLines,extension=extension)
			
	sprintf wname, filePrefix, istart
	NewDataFolder/O root:Data
	NewDataFolder/O root:Data:$wName
	Duplicate/O dm, root:Data:$(wName):$(wName)
			
	sprintf wname2,fileprefix+"_%d%d", istart, ycol, xcol
			
	Make/FREE/O/N=(DimSize(dm,0)) yw, zw
			
	zw = dm[p][ycol-1]
	yw = dm[p][xcol-1]
			
	zw /= zscale
	yw /= yscale
			
	Redimension/N=(npnts,numpnts(zw)/npnts) zw
	SetScale/P y,xstart,xstep,xunits,zw
	SetScale/I x,yw[0],yw[npnts-1],yunits,zw
	//Print yw[0],yw[numpnts(yw)-2],numpnts(yw)-1
	if (slowaxis)
		MatrixTranspose zw
	endif
			
	SetScale d,0,0,zunits,zw

	Duplicate/O zw, m, root:$(wname2), root:Data:$(wName):$(wname2)
			
	sprintf activeMatrix, "root:%s", wname2
End
	

Function LoadButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	String dfSav = GetDataFolder(1)
	SetDataFolder root:Packages:GQ:Base
	switch( ba.eventCode )
		case 2: // mouse up
			LoadButtonFunc()
			break
		case -1: // control being killed
			break
	endswitch
	SetDataFolder dfSav

	return 0
End



Window LoadFilePanel_BR() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1005,308,1588,477) as "Load file BR"
	SetVariable filePath,pos={28,14},size={524,16},title="File Path"
	SetVariable filePath,limits={-inf,inf,0},value= root:Packages:GQ:Base:filePathStr
	SetVariable filePrefix,pos={28,36},size={104,16},title="Data Series"
	SetVariable filePrefix,value= root:Packages:GQ:Base:filePrefix
	SetVariable istart,pos={141,35},size={112,16},title="File Index",format="%d"
	SetVariable istart,value= root:Packages:GQ:Base:istartBR
	CheckBox Add_to_active_graph,pos={239,103},size={111,14},title="Add to active graph"
	CheckBox Add_to_active_graph,variable= root:Packages:GQ:Base:addtograph
	Button LoadMatrix,pos={323,39},size={84,35},proc=LoadFileButtonBR,title="Load"
	PopupMenu popupx,pos={54,85},size={100,21},title="X wave"
	PopupMenu popupx,mode=1,popvalue="wave0",value= #"root:Packages:GQ:Base:waveNames"
	PopupMenu popupy,pos={52,114},size={103,21},title="Y Wave"
	PopupMenu popupy,mode=2,popvalue="wave1",value= #"root:Packages:GQ:Base:waveNames"
	Button PlotButton,pos={175,82},size={60,56},proc=PlotDataButton,title="Plot"
	SetVariable nbtoskip1,pos={144,55},size={127,16},title="Nb of lines to skip"
	SetVariable nbtoskip1,limits={0,inf,1},value= root:Packages:GQ:Base:nbtoskip
EndMacro



Function PlotDataButton(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	String dfSav = GetDataFolder(1)
	SetDataFolder root:Packages:GQ:Base
	SVAR filePrefix, waveNames
	NVar istartBR, addtograph
	SetDataFolder dfSav

	switch( ba.eventCode )
		case 2: // mouse up
			String wname
			sprintf wname, filePrefix, istartBR
			String fname
			String wavenameX, wavenameY
			ControlInfo popupx
			String labelX=S_value
			sprintf wavenameX,"%s_%s",wname,S_value
			ControlInfo popupy
			String labelY=S_value
			sprintf wavenameY,"%s_%s",wname,S_value
			
			if (addtograph == 1)
				String s = WinName(0,1)
				if (numtype(str2num(s)) != 2)
					AppendToGraph  $(wavenameY) vs $(wavenameX)
				else
					Display/N=$wname  $(wavenameY) vs $(wavenameX)
				endif
			elseif (addtograph == 0)
				Display/N=$wname  $(wavenameY) vs $(wavenameX)
			EndIf
				Label bottom labelX
				Label left labelY
				Legend/C/N=text0/A=MC			
			//Display $(wavenameY) vs $(wavenameX)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function JointHistogram1(w0,w1,hist)
	wave w0,w1,hist
 
	variable bins0=dimsize(hist,0)
	variable bins1=dimsize(hist,1)
	variable n=numpnts(w0)
	variable left0=dimoffset(hist,0)
	variable left1=dimoffset(hist,1)
	variable right0=left0+bins0*dimdelta(hist,0)
	variable right1=left1+bins1*dimdelta(hist,1)
 
	// Scale between 0 and the number of bins to create an index wave.  
	if(ThreadProcessorCount<4) // For older machines, matrixop is faster.  
		matrixop /free idx=round(bins0*(w0-left0)/(right0-left0))+bins0*round(bins1*(w1-left1)/(right1-left1))
	else // For newer machines with many cores, multithreading with make is faster.  
		make /free/n=(n) idx
		multithread idx=round(bins0*(w0-left0)/(right0-left0))+bins0*round(bins1*(w1-left1)/(right1-left1))
	endif
 
	// Compute the histogram and redimension it.  
	histogram /b={0,1,bins0*bins1} idx,hist
	redimension /n=(bins0,bins1) hist // Redimension to 2D.  
	setscale x,left0,right0,hist // Fix the histogram scaling in the x-dimension.  
	setscale y,left1,right1,hist // Fix the histogram scaling in the y-dimension.  
End

Function KillAllGraphs()
	string fulllist = WinList("*", ";","WIN:1")
	string name, cmd
	variable i
 
	for(i=0; i<itemsinlist(fulllist); i +=1)
		name= stringfromlist(i, fulllist)
		sprintf  cmd, "Dowindow/K %s", name
		execute cmd		
	endfor
end

Function LoadFileButtonBR(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	String dfSav = GetDataFolder(1)
	SetDataFolder root:Packages:GQ:Base
	NVAR istartBR, xscale1, yscale1,addtograph,xcol,ycol, rbias, correctr
	NVAR nbtoskip
	SVAR filePrefix, filePathStr, xunits, yunits
	SetDataFolder dfSav

	switch( ba.eventCode )
		case 2: // mouse up
			String wname
			sprintf wname, filePrefix, istartBR
			String fname
			sprintf fname, "%s%s.txt", filePathStr, wname

			LoadWave/A/V={"\t"," ",0,0}/D/W/J/L={nbtoskip-1,nbtoskip,0,0,0}/K=1/O fname
			Variable i, n = V_Flag
			String wavenamein, wavenameout			
			for(i=0;i<n;i+=1) 
			     wavenamein = StringFromList(i,S_waveNames)
			     //WAVE w = $(wavename2)
			     sprintf  wavenameout,"%s_%s",wname,wavenamein
			     KillWaves/Z $wavenameout
			     Rename $wavenamein, $wavenameout
			     Print wavenameout
			 endfor

		SetDataFolder root:Packages:GQ:Base
		String/G waveNames 
		waveNames = S_waveNames
		SetDataFolder dfSav


			break
		case -1: // control being killed
			break
	endswitch

	return 0
End
