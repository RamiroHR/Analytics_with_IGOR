#pragma rtGlobals=1		// Use modern global access method.
#pragma version=1.21

// CHANGES:
// 1.21 added SweepScopeFlux to sweep flux in a better way
// 1.20 added "Initialize_ykDL750(0)" to ykPanel in order to force config directory creation
				
#pragma IgorVersion = 6.3
		
//#pragma IndependentModule = ykDL750

//	ykDL750.ipf
//
//	Package ykDL750 interfaces Yokogawa DL 750 Oscilloscope to Igor Pro through GPIB.
//
//	requires VISA xop and GQ Base.ipf 
//
//	2016.08.12

//#include ":Base"

//
// CHANGES
//


// 2018.05.25
// Added option to control two Yoko 7651 sources using the "Sweep and Scope Twice" Panel: SweepScopePanel_Twice()
//


//	GPIB Address

#if Exists("viOpen") 		// Only include this code if the VISA xop and Base.ipf is available
						// otherwise few functions will work.
						
//#include <Execute Cmd On List>  
// BEGIN Macros 
//
Menu "ykDL750"
	"Yoko DL750",ykPanel()
	"Sweep Scope",SweepScopePanel()
	"Sweep Scope Flux", SweepScopeFluxPanel()
	"DoubleSweep Scope",SweepScopePanel_Twice()
	"Rigol Sweep Scope",SweepScopePanel_R()

End
//
// END Macros

// BEGIN Constants
// 
//StrConstant kGPIBn = "GPIB1"
StrConstant kGPIBn = "GPIB0"

// Yokogawa DL750 with DSP
//StrConstant kAddress = "GPIB1::8::INSTR"
StrConstant kAddress = "GPIB0::11::INSTR"  
Constant kNumDSP = 6

// Yokogawa DL750 without DSP
// StrConstant kAddress = "GPIB1::10::INSTR"
// Constant kNumDSP = 0

// Max number of acquisition channels
Constant kNumAcq = 16

// Max number of DSP channels
// Set to ZERO if there is no DSP option in DL750
// Constant kNumDSP = 6
// DEFINED ABOVE AT GPIB ADDRESS

// Max number of math channels
Constant kNumMath = 8

// Sum of all channels
Constant kNumChan = 30

Constant kChunkSize = 1e5
//
//END Constants

// ReportVISAError (name, session, status)
// See documentation for viStatusDesc function for an example.
static Function ReportVISAError(name, session, status)
	String name					// VISA function name, e.g., viRead or other identifier
	Variable session			// Session ID obtained from viOpen
	Variable status			// Status code from VISA library
	
	String desc
	
	viStatusDesc(session, status, desc)
	Printf "%s error (%x): %s\r", name, status, desc
	Beep
End
 

Function Initialize_ykDL750(reinit,[gpibaddress])
	Variable reinit
	String gpibaddress
	
	if(ParamIsDefault(gpibaddress))
		gpibaddress = kAddress
	endif

	String savDF= GetDataFolder(1) // Save current DF for restore.
	
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S GQ	// Make sure this exists.
	
	if(!DataFolderExists("ykDL750")  || reinit) // Already created?
		NewDataFolder/O/S ykDL750 // Our stuff goes in here.
		Make/O tracew

		String/G address = gpibaddress
		Variable/G fileIndex = 0, showLegend = 1, graphStyle = 1, scaleVar = 1, appending = 0
		Variable/G kNumChan = kNumAcq + kNumDSP + kNumMath
		String/G filePathStr = "U:Joël:Data:20160812 - Hf&Nb", basename = "Hf"
		Make/O/T/N=(kNumChan,3) traceSettings

		Variable i
		For(i=0;i<kNumAcq;i+=1)
			SetDimLabel 0,i,$("CH"+num2str(i+1)),traceSettings
		EndFor

		For(i=0;i<kNumDSP;i+=1)
			SetDimLabel 0,i+kNumAcq,$("DSP"+num2str(i+1)),traceSettings
		EndFor

		For(i=0;i<kNumMath;i+=1)
			SetDimLabel 0,i+kNumAcq+kNumDSP,$("MTH"+num2str(i+1)),traceSettings
		EndFor
		
		SetDimLabel 1,0,Gain,traceSettings
		SetDimLabel 1,1,Units,traceSettings
		SetDimLabel 1,2,ID,traceSettings
		
		traceSettings[][0] = "1"
		traceSettings[][1] = "V"
		traceSettings[][2] = num2str(p+1)

		// Prepare YOKO 7651 1

		Variable/G yoko7651address = 16
		Variable/G yoko7651currentsrc = 1
		Variable/G yoko7651initialvalue = 0
		Variable/G yoko7651finalvalue = 3e-4
		Variable/G yoko7651step = 1e-6
		Variable/G yoko7651numpnts = 301
		Variable/G yoko7651currpnt = NaN
	
		Make/O/N=(yoko7651numpnts) yoko7651loopwave = yoko7651initialvalue + p*yoko7651step
		SetScale d,0,0,"A",yoko7651loopwave
		

		// Prepare YOKO 7651 2

		Variable/G yoko7651address2 = 5
		Variable/G yoko7651currentsrc2 = 1
		Variable/G yoko7651initialvalue2 = 0
		Variable/G yoko7651finalvalue2 = 3e-4
		Variable/G yoko7651step2 = 1e-6

		Make/O/N=(yoko7651numpnts) yoko7651loopwave2 = yoko7651initialvalue2 + p*yoko7651step2
		SetScale d,0,0,"A",yoko7651loopwave2
	
	
		// For Flux Sweeps
		Variable/G yoko7651phizero1 = 0
		Variable/G yoko7651phi2pi1 = 1e-6
		Variable/G yoko7651phi01 = 0
		Variable/G yoko7651phi11 = 1
		
		Variable/G yoko7651phizero2 = 0
		Variable/G yoko7651phi2pi2 = 1e-6
		Variable/G yoko7651phi02 = 0
		Variable/G yoko7651phi12 = 1

		// For Flux compensation
		Variable/G yoko7651gamma12 = 0
		
		// For sinusoidal flux compensation
		Variable/G yoko7651sinamp = 0	
		
				
		// Prepare Rigol 1032Z 1

		Variable/G rigolCH1amp = 1.0
		Variable/G rigolCH2amp = 1.0		
		Variable/G rigolnumpnts = 128
		Variable/G rigolcurrpnt = NaN
		String/G rigoladdress = "USB0::0x1AB1::0x0642::DG1ZA183802718::0::INSTR"
		
		Make/O/N=(rigolnumpnts,2) rigolloopwave
		rigolloopwave[][0] = rigolCH1amp*cos(p/rigolnumpnts*2*pi)
		rigolloopwave[][1] = rigolCH1amp*sin(p/rigolnumpnts*2*pi)
		SetScale d,0,0,"V",rigolloopwave
		
		// Rigol cannot output waves with amplitude less than 2mV
		
		NewDataFolder/O root:Data
	endif
	
	SetDataFolder root:Packages:GQ:ykDL750
	String/G address = gpibaddress
	WAVE tracew = tracew
	if(QueryIDN(address) != 0)
		Beep
		Print address + ": communication error!\r Package ykDL750 may not work."
		return -1
	else
		ykAvailableTraces(address, tracew)
		//Redimension/N=(numpnts(traceTw),3) traceCW
	endif
	
	SetDataFolder savDF // Restore current DF.
End

//Function ReportVISAError(name, session, status)
//	String name					// VISA function name, e.g., viRead or other identifier
//	Variable session			// Session ID obtained from viOpen
//	Variable status			// Status code from VISA library
//	
//	String desc
//	
//	viStatusDesc(session, status, desc)
//	Printf "%s error (%x): %s\r", name, status, desc
//	Beep
//End
//
//Function ClearGraph(graphName)
//	String graphName
//	
//	RemoveAllTraces(graphName = graphName)
//	
//	String annlist = AnnotationList(graphName)
//	Variable i, n = ItemsInList(annlist)
//	For(i=0;i<n;i+=1)
//		TextBox/K/W=$graphName/N=$(StringFromList(i,annlist))
//	EndFor
//	
//	return 0
//End
//
//Function RemoveAllTraces([graphName, matchStr])
//	String graphName, matchStr
//	
//	if(ParamIsDefault(graphName))
//		graphName = ""
//	endif
//	
//	if(ParamIsDefault(matchStr))
//		matchStr = "*"
//	endif
//	
//	String traces = ListMatch(TraceNameList(graphName,";",1),matchStr)
//	Variable i, n = itemsinlist(traces)
//	For(i=0;i<n;i+=1)
//		RemoveFromGraph/Z/W=$graphName $(StringFromList(i,traces))
//	EndFor
//	
//	return 0
//End


Function QueryIDN(resourceName)
	String resourceName
	
	Variable defaultRM, instr
	
	Variable status
	Variable numErrors = 0

	status = viOpenDefaultRM(defaultRM)
	//	Printf "DefaultRM=%d\r", defaultRM
	if (status < 0)
		ReportVISAError("viOpenDefaultRM", defaultRM, status)
		numErrors += 1
	endif

	status = viOpen(defaultRM, resourceName, 0, 0, instr)
	//	Printf "instr=%d\r", instr
	if (status < 0)
		ReportVISAError("viOpen", instr, status)
		numErrors += 1
	endif
	
	String str1 = VISABinQuery(instr, "*IDN?")
	Print str1
	
	viClose(instr)
	
	viClose(defaultRM)

	return numErrors
End


Function/S QueryVISA(resourceName, cmd)
	String resourceName, cmd
	
	String outstr
	Variable defaultRM, instr
	
	Variable status
	Variable numErrors = 0

	status = viOpenDefaultRM(defaultRM)
	//	Printf "DefaultRM=%d\r", defaultRM
	if (status < 0)
		ReportVISAError("viOpenDefaultRM", defaultRM, status)
		numErrors += 1
	endif

	status = viOpen(defaultRM, resourceName, 0, 0, instr)
	//	Printf "instr=%d\r", instr
	if (status < 0)
		ReportVISAError("viOpen", instr, status)
		numErrors += 1
	endif
	
	outstr = VISABinQuery(instr, cmd)

	viClose(instr)
	
	viClose(defaultRM)
	return	outstr
//	return numErrors
End

// Prints out result

Function PrintQueryVISA(resourceName, cmd)
	String resourceName, cmd
	
	String outstr
	
	Variable defaultRM, instr
	
	Variable status
	Variable numErrors = 0

	status = viOpenDefaultRM(defaultRM)
	//	Printf "DefaultRM=%d\r", defaultRM
	if (status < 0)
		ReportVISAError("viOpenDefaultRM", defaultRM, status)
		numErrors += 1
	endif

	status = viOpen(defaultRM, resourceName, 0, 0, instr)
	//	Printf "instr=%d\r", instr
	if (status < 0)
		ReportVISAError("viOpen", instr, status)
		numErrors += 1
	endif
	
	outstr = VISABinQuery(instr, cmd)
	
	Print outstr

	viClose(instr)
	
	viClose(defaultRM)
	return numErrors
End

Function WriteVISA(resourceName,cmd)
	String resourceName, cmd
	
	Variable defaultRM, instr
	
	Variable status
	Variable numErrors = 0

	status = viOpenDefaultRM(defaultRM)
	if (status < 0)
		ReportVISAError("viOpenDefaultRM", defaultRM, status)
		numErrors += 1
	endif

	status = viOpen(defaultRM, resourceName, 0, 0, instr)
	//	Printf "instr=%d\r", instr
	if (status < 0)
		ReportVISAError("viOpen", instr, status)
		numErrors += 1
	endif
	
	VISAWrite/Q instr, cmd
	
	if(V_flag == 0)
		ReportVISAError("viWrite", instr, V_status)
		numErrors += 1
	endif

	viClose(instr)
	
	viClose(defaultRM)

	return numErrors
End

Function/S VISABinQuery(instr, query)
	Variable instr
	String query
	
	String retstr = ""
	
	String str = "viWrite: " + query
	VISAWrite instr, query

	if(V_flag == 0)
		ReportVISAError(str, instr, V_status)
		return retstr
	endif

	str = "viRead: " + query

	// WARNING: will not work if return string is larger than 1KB
	VISAReadBinary/S=1024 instr, retstr

//	VISARead instr, retstr
	if(V_flag == 0)
		if(V_status == -1)
			Print "Error trying to call viRead!\r"
		else
			ReportVISAError(str, instr, V_status)
		endif
	endif

	// Remove newline
	retstr = RemoveEnding(retstr)
		
	return retstr
End


Function/S VISAQuery(instr, query)
	Variable instr
	String query
	
	String retstr = ""
	
	String str = "viWrite: " + query
	VISAWrite instr, query

	if(V_flag == 0)
		ReportVISAError(str, instr, V_status)
		return retstr
	endif

	str = "viRead: " + query

	VISARead/T="\n" instr, retstr
//	VISARead instr, retstr
	if(V_flag == 0)
		if(V_status == -1)
			Print "Error trying to call viRead!\r"
		else
			ReportVISAError(str, instr, V_status)
		endif
	endif
	
	return retstr
End

Function/S GetYoko7651(resourceName)  // DOESN'T WORK
	String resourceName
	Variable voltage, autorange
	
	Variable defaultRM, status, instr
	
	status = viOpenDefaultRM(defaultRM)
	if(status < 0)
		ReportVISAError("viOpenDefaultRM", instr, status)
		viClose(defaultRM)
		Abort
	endif
	
	status = viOpen(defaultRM, resourceName, 0, 0, instr)
	if(status < 0)
		ReportVISAError("viOpen", instr, status)
		viClose(instr)
		viClose(defaultRM)
		Abort
	endif

	String progstr, response

	progstr = "OSE"
	response = VISABinQuery(instr, progstr)
	
	progstr = "Remote Clear"
	status = viGpibControlREN(instr,0)
	AbortOnValue status<0,99

	viClose(instr)
	viClose(defaultRM)

	return response
End

Function Yoko7651SetRange(resourceName, output,[current])
	String resourceName
	Variable output, current
	
	if(ParamIsDefault(current))
		current = 0
	endif
	
	Variable defaultRM, status, instr
	
	status = viOpenDefaultRM(defaultRM)
	if(status < 0)
		ReportVISAError("viOpenDefaultRM", instr, status)
		viClose(defaultRM)
		Abort
	endif
	
	status = viOpen(defaultRM, resourceName, 0, 0, instr)
	if(status < 0)
		ReportVISAError("viOpen", instr, status)
		viClose(instr)
		viClose(defaultRM)
		Abort
	endif

	String progstr

// Does not change output range.  This should be done manually.
//
//	if(current)
//		progstr = "F5E"
//	else
//		progstr = "F1E"
//	endif
	
//	VISAWrite instr, progstr
//	status = V_status
//	AbortOnValue V_flag==0,1
//
	if(current)
		if(output <= 1e-3)
			progstr = "R4E"
		elseif(output <= 1e-2)
			progstr = "R5E"
		elseif(output <= 1e-1)
			progstr = "R6E"
		else
			Print "ERROR: Current out of range"
			status = -1
		endif
	else
		if(output <= 1e-2)
			progstr = "R2E"
		elseif(output <= 1e-1)
			progstr = "R3E"
		elseif(output <= 1)
			progstr = "R4E"
		elseif(output <= 10)
			progstr = "R5E"
		elseif(output <= 30)
			progstr = "R6E"
		else
			Print "ERROR: Voltage out of range"
			status = -1
		endif
	endif
	
	if(status != -1)
		VISAWrite instr, progstr
		status = V_status
		AbortOnValue V_flag==0,1
	endif
	
	progstr = "Remote Clear"
	status = viGpibControlREN(instr,0)
	AbortOnValue status<0,99

	viClose(instr)
	viClose(defaultRM)

	return status
End

Function Yoko7651SetOutput(resourceName, output, [autorange])
	String resourceName
	Variable output, autorange
	
	if(ParamIsDefault(autorange))
		autorange = 0
	endif
	
	Variable defaultRM, status, instr
	
	status = viOpenDefaultRM(defaultRM)
	if(status < 0)
		ReportVISAError("viOpenDefaultRM", instr, status)
		viClose(defaultRM)
		Abort
	endif
	
	status = viOpen(defaultRM, resourceName, 0, 0, instr)
	if(status < 0)
		ReportVISAError("viOpen", instr, status)
		viClose(instr)
		viClose(defaultRM)
		Abort
	endif

	String progstr

	if(autorange)
		progstr = "SA" + num2str(output) + "EO1E"
	else
		progstr = "S" + num2str(output) + "EO1E"
	endif
	
	VISAWrite instr, progstr
	status = V_status
	AbortOnValue V_flag==0,1
	
	progstr = "Remote Clear"
	status = viGpibControlREN(instr,0)
	AbortOnValue status<0,99

	viClose(instr)
	viClose(defaultRM)

	return status
End

Function/S ykQuery(resourceName, query)
	String resourceName
	String query
	
	String response = "", tmpstr = ""

	String progstr
	Variable defaultRM, instr,status

	status = viOpenDefaultRM(defaultRM)
	if(status < 0)
		ReportVISAError("viOpenDefaultRM", instr, status)
		viClose(defaultRM)
		Abort
	endif
	
	status = viOpen(defaultRM, resourceName, 0, 0, instr)
	if(status < 0)
		ReportVISAError("viOpen", instr, status)
		viClose(instr)
		viClose(defaultRM)
		Abort
	endif
	
	progstr = "COMMUNICATE:HEADER OFF"
	VISAWrite instr, progstr
	status = V_status
	AbortOnValue V_flag==0,1
	
	response = VISABinQuery(instr,query)
	AbortOnValue strlen(response) == 0,2
	
	progstr = "COMMUNICATE:HEADER ON"
	VISAWrite instr, progstr
	status = V_status
	AbortOnValue V_flag==0,2
		
	progstr = "Remote Clear"
	status = viGpibControlREN(instr,0)
	AbortOnValue status<0,99

	viClose(instr)
	viClose(defaultRM)
	
	return response
End




Function ykWrite(resourceName, command)
	String resourceName
	String command
	
	Variable defaultRM, instr,status

	String progstr

	status = viOpenDefaultRM(defaultRM)
	if(status < 0)
		ReportVISAError("viOpenDefaultRM", instr, status)
		viClose(defaultRM)
		Abort
	endif
	
	status = viOpen(defaultRM, resourceName, 0, 0, instr)
	if(status < 0)
		ReportVISAError("viOpen", instr, status)
		viClose(instr)
		viClose(defaultRM)
		Abort
	endif
	
//	progstr = "COMMUNICATE:HEADER OFF"
//	VISAWrite instr, progstr
//	status = V_status
//	AbortOnValue V_flag==0,1
//	
	VISAWrite instr, command
	status = V_status
	AbortOnValue V_flag== 0,1.5
	
//	progstr = "COMMUNICATE:HEADER ON"
//	VISAWrite instr, progstr
//	status = V_status
//	AbortOnValue V_flag==0,2
//		
	progstr = "Remote Clear"
	status = viGpibControlREN(instr,0)
	AbortOnValue status<0,99

	viClose(instr)
	viClose(defaultRM)
	
	return status
End


Function ykGetTrace(resourceName, traceNum, wName, [record])
	String resourceName                     	// Resource name for Yokogawa DL750 (i.e. "GPIB0::2::INSTR")
	Variable traceNum                        	// Trace number to transfer
	String wName                               	// Name of wave
	Variable record
	
	if(ParamIsDefault(record))
		record = 0
	else
		record = -abs(record)
	endif
	
	Variable defaultRM, instr
	Variable status, err=0, retCnt
	String progstr, tmpstr

	//Make/O wtmp

	Variable dt, dV, offset, npnts, trigpos, bitlength, minrecord, nr

	status = viOpenDefaultRM(defaultRM)
	if(status < 0)
		ReportVISAError("viOpenDefaultRM", instr, status)
		viClose(defaultRM)
		return -1
	endif
	
	status = viOpen(defaultRM, resourceName, 0, 0, instr)
	if(status < 0)
		ReportVISAError("viOpen", instr, status)
		viClose(instr)
		viClose(defaultRM)
		return -1
	endif
	
	try
		progstr = "*CLS"
		VISAWrite instr, progstr
		status = V_status
		AbortOnValue V_flag==0,0
	
		progstr = "COMMUNICATE:HEADER OFF"
		VISAWrite instr, progstr
		status = V_status
		AbortOnValue V_flag==0,1
	
		progstr = "TIMEBASE:SRATE?"
		dt = str2num(VISABinQuery(instr,progstr))
		AbortOnValue (numtype(dt)==2),2
		dt = 1/dt
		
		if(traceNum < kNumAcq)
			progstr = "WAVEFORM:TRACE " + num2str(traceNum)
		elseif(traceNum < kNumAcq + kNumDSP)
			progstr = "WAVEFORM:TRACE DSP" + num2str(traceNum-kNumAcq)
		else
			progstr = "WAVEFORM:TRACE MATH" + num2str(traceNum-kNumAcq-kNumDSP)
		endif
		VISAWrite instr, progstr
		status = V_status
		AbortOnValue V_flag==0,21
		
		progstr = "WAVEFORM:RECord? MINimum"
		minrecord = str2num(VISABinQuery(instr,progstr))
		AbortOnValue (numtype(minrecord)==2),22
		
		if(record < minrecord)
			record = 0
		endif		
				
		progstr = "WAVEFORM:RECORD " + num2str(record)
		VISAWrite instr, progstr
		status = V_status
		AbortOnValue V_flag==0,23
		
		progstr = "WAVEFORM:RANGE?"
		dV = str2num(VISABinQuery(instr,progstr))
		AbortOnValue (numtype(dV)==2),3

		progstr = "WAVEFORM:OFFSET?"
		offset = str2num(VISABinQuery(instr,progstr))
		AbortOnValue (numtype(offset)==2),4
		
		progstr = ":STOP"
		VISAWrite instr, progstr
		status = V_status
		AbortOnValue V_flag==0,5
		
		progstr = "WAVeform:LENGth?"
//		progstr = "ACQ:RLEN?"
		npnts = str2num(VISABinQuery(instr,progstr))
		AbortOnValue (numtype(npnts)==2),6
		
		progstr = "WAVeform:TRIGGER?"
		trigpos = str2num(VISABinQuery(instr,progstr))
		AbortOnValue (numtype(trigpos)==2),7
		
		progstr = "WAVEFORM:FORMAT WORD"
		VISAWrite instr, progstr
		status = V_status
		AbortOnValue V_flag==0,8
		
		progstr = "WAVeform:BITS?"
		bitlength = str2num(VISABinQuery(instr,progstr))
		AbortOnValue (numtype(bitlength)==2),81
		bitlength = (bitlength == 16 ? 0x10 : 0x08)         // Whether 16 bit or 8 bit channel
		
		progstr = "WAVEFORM:BYTEORDER MSBFIRST"
		VISAWrite instr, progstr
		status = V_status
		AbortOnValue V_flag==0,9
		
		if(npnts > kChunkSize)
			printf "Transferring"
		endif
		
		Make/FREE/N=(kChunkSize+1) wtmp
		Make/FREE/N=(npnts) w
		Variable i, n = floor(npnts/kChunkSize), m, numtrans
	//	Make/FREE/N=(30) tmpw
		
		For(i=0;i<=n;i+=1)			// transfer 1e5 points at a time to work around the VISA timeout for larger transfers
			if(npnts > kChunkSize)
				printf "."
			endif
			
			m = min(npnts,(i+1)*kChunkSize)-1
			sprintf progstr, "WAVEFORM:START %d;:WAVEFORM:END %d", i*kChunkSize, m
			VISAWrite instr, progstr
			status = V_status
			AbortOnValue V_flag==0,61

			progstr = "WAVEFORM:SEND?"
			VISAWrite instr, progstr
			status = V_status
			AbortOnValue V_flag==0,10
				
			progstr = "viRead: data length"
			//status = viRead(instr, tmpstr, 2, retCnt)
			status = viRead(instr, tmpstr,11, retCnt)
			AbortOnValue status<0,11

//			Print tmpstr
//			nr = str2num(tmpstr[1])
//			progstr = "viRead: data length"
//			status = viRead(instr, tmpstr, nr, retCnt)
//			AbortOnValue status<0,12
//			numtrans = str2num(tmpstr)/2			// divide by two since each word is two bytes
//			Print numtrans
				
			VISAReadBinaryWave/Y={offset,dV/2400}/TYPE=(bitlength) instr, wtmp	// for scaling factors, see pg 6-178 of DL750 Communications Reference Manual
			progstr = "VISAReadBinaryWave"
			status = V_status
			AbortOnValue V_flag==0,13
//			Print "***", status			
			if(status!=0)
				ReportVISAError(progstr, instr, status)
			endif

			AbortOnValue status<0,14
			
			w[i*kChunkSize,m] = wtmp[p-i*kChunkSize]
		EndFor
		
	//	Print toc-tic, status, V_flag
		
		
		if(npnts > kChunkSize) 
			printf "\r"
		endif
		
		Duplicate/O w, $wName
		WAVE w = $wName

		SetScale/P x,-trigpos*dt,dt,"s",w
		SetScale d,0,0,"V",w

		progstr = "COMMUNICATE:HEADER ON"
		VISAWrite instr, progstr
		status = V_status
		AbortOnValue V_flag==0,20
		
		progstr = "*CLS"
		VISAWrite instr, progstr
		status = V_status
		AbortOnValue V_flag==0,90

		progstr = "Remote Clear"
		status = viGpibControlREN(instr,0)
		AbortOnValue status<0,99
	catch
		//		Print progstr
		ReportVISAError(progstr, instr, status)
		err = -1
	endtry
	
	//Make/N=(npnts) $wName	
	//WAVE w = $wName
	
	//w = dV / 9.375			// for scaling factors, see pg 6-178 of DL750 Communications Reference Manual
	//w += offset
	
	// time to acquire is <2 sec per 100K points

	//debug	Print dV, dt, offset, npnts

	//KillWaves wtmp
	viClose(instr)
	viClose(defaultRM)

	return err
End

Function LoopAndMeasureYokoDL750(loopw,setFunc,dataName)
	WAVE loopw
	FUNCREF protoSetFunc setFunc
	String dataName

	Variable status

	DFREF saveDFR = GetDataFolderDFR()
	
	NewDataFolder/O/S root:Data
	NewDataFolder/O/S root:Data:$dataName
	
	Variable i, n = numpnts(loopw)
	
	For(i=0;i<n;i+=1)
		
		if(setFunc(loopw[i]))
			Print "ERROR: setFunc failed!"
			status = -1
			break
		endif
		
		
		
	EndFor
		
	SetDataFolder saveDFR

	return status
End

Function protoSetFunc(value)
	Variable value
	
	// Set output of device to desired value
	Variable status = 0 // success
	
	return status
End

Function ykAvailableTraces(resourceName, w)
	String resourceName                     	// Resource name for Yokogawa DL750 (i.e. "GPIB0::2::INSTR")
	WAVE w                               	// Name of wave
	
	Variable defaultRM, instr
	Variable status, err=0
	String progstr, response

	Variable nchan = 0, i, channelOn=0

	status = viOpenDefaultRM(defaultRM)
	if(status < 0)
		ReportVISAError("viOpenDefaultRM", instr, status)
		viClose(defaultRM)
		return -1
	endif
	
	status = viOpen(defaultRM, resourceName, 0, 0, instr)
	if(status < 0)
		ReportVISAError("viOpen", instr, status)
		viClose(instr)
		viClose(defaultRM)
		return -1
	endif
	
	//	status = viGpibSendIFC(instr)
	//	if(status < 0)
	//		ReportVISAError("viGpibSendIFC", instr, status)
	//		viClose(instr)
	//		viClose(defaultRM)
	//		return -1
	//	endif
	
	Redimension/N=(kNumChan) w
	
	try
		progstr = "COMMUNICATE:HEADER OFF"
		VISAWrite instr, progstr
		status = V_status
		AbortOnValue V_flag==0,1
	
		
	
		for(i=1;i<=kNumAcq;i+=1)		// Check acquisition channels
			sprintf progstr, "CHANNEL%d:MODULE?", i
			response = VISABinQuery(instr,progstr)
//			AbortOnValue V_flag==0,1
			
			if(stringmatch(response,"NOMODULE"))
				continue
			endif
			
			sprintf progstr, "CHANNEL%d:DISPLAY?", i
			//Print progstr
			channelOn = str2num(VISABinQuery(instr,progstr))
			AbortOnValue (numtype(channelOn)==2),6
		
			if(channelOn!=0)
				w[nchan] = i
				nchan+=1
			endif
		endfor			
		
		for(i=1;i<=kNumDSP;i+=1)		// Check DSP channels
			sprintf progstr, "DSP%d:DISPLAY?", i
			//Print progstr
			channelOn = str2num(VISABinQuery(instr,progstr))
			AbortOnValue (numtype(channelOn)==2),6
		
			if(channelOn!=0)
				w[nchan] = i+kNumAcq
				nchan+=1
			endif
		endfor			

		for(i=1;i<=kNumMath;i+=1)		// Check MATH channels
			sprintf progstr, "MATH%d:DISPLAY?", i
			//Print progstr
			channelOn = str2num(VISABinQuery(instr,progstr))
			AbortOnValue (numtype(channelOn)==2),6
		
			if(channelOn!=0)
				w[nchan] = i+kNumAcq+kNumDSP
				nchan+=1
			endif
		endfor			

		progstr = "COMMUNICATE:HEADER ON"
		VISAWrite instr, progstr
		status = V_status
		AbortOnValue V_flag==0,98
	
		progstr = "Remote Clear"
		status = viGpibControlREN(instr,0)
		AbortOnValue status<0,99
	catch
		ReportVISAError(progstr, instr, status)
		err = -1
	endtry
	
	Redimension/N=(nchan) w

	viClose(instr)
	viClose(defaultRM)

	return err
End

Function ykBusy()
	return str2num(ykQuery(kAddress,"STAT:COND?"))
End

Function ykReady()
	return !(stringmatch(ykQuery(kAddress,"ACQuire:MODE?"),"NORM") && stringmatch(ykQuery(kAddress,"TRIG:MODE?"),"NORM"))
End

Function AddCommandsToItx(fileName,pathStr,cmdStr)
	String fileName, pathStr, cmdStr
	
	String outStr
		
	Variable refnum
	Open/Z=1/A refnum as (pathStr + ":" + fileName)

	outStr = ReplaceString("\r",cmdStr,"\rX ")
	
//	printf "%s", outStr[0,1001]
//	Print strlen(outStr)
	fprintf refnum, "X %s", outStr[0,900]
	fprintf refnum, "%s", outStr[901,2000]

	Close refnum
		
	return V_flag
End

Function LoadDataButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	String savDF= GetDataFolder(1) // Save current DF for restore.
	SetDataFolder root:Packages:GQ:ykDL750
	
	Variable err = 0, status
	
	WAVE tw = tracew
	WAVE/T tS = traceSettings
	SVAR address, basename, filePathStr
	NVAR scaleVar
	NVAR fileIndex, graphStyle //, showLegend
	Duplicate/O tw, twOld
	
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			String wName, pathName, cmdStr = ""
			String progstr, xymode, tmpstr
			Variable xaxis
								
			// Obtain all available traces on scope
			ykAvailableTraces(address,tw)
			Variable i, n = numpnts(tw)
					
			// Check to see if default graph "ykGraph" is already open for graphing trace data
			if(strlen(WinList("ykGraph",";","")) == 0)
				Display/N=ykGraph
				DoWindow/T ykGraph, "DL750 Graph"
			else
				DoWindow/F ykGraph
			endif

			// Clear Graph
			ClearGraph("ykGraph")
			
			// Create path for new data				
			sprintf pathName, "root:Data:%s%04d", basename, fileIndex
			fileIndex = CheckDataPath(pathName,basename,fileIndex)
			NewDataFolder/O/S $pathName

			// Get data for each available trace from scope
			For(i=0;i<n;i+=1)
				sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
				ykGetTrace(address, tw[i], wName)
			EndFor

			// Check to see if we are in XY mode
			progstr = "XY1:MODE?"
			xymode = ykQuery(address,progstr) 

			if(strlen(xymode) == 0)
				Print "ERROR GPIB communication!"
				err = -1
			endif
			
			// Plot data
			strswitch(xymode)
				case "TY_XY":
				case "XY":					// Plot data as Y vs X
					progstr = "XY1:XTRACE?"
					tmpstr = ykQuery(address,progstr)
					xaxis = str2num(tmpstr)
					
					if(strlen(tmpstr) == 0)
						Print "ERROR GPIB communication!"
						err = -1
						break
					endif
					
					sprintf wName, "%s%04d_%s", basename, fileIndex, tS[xaxis-1][2]
					WAVE wx = $wName

					for(i=0;i<n;i+=1)
						if(tw[i] == xaxis)
							continue
						endif
						sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
						WAVE wy = $wName
						AppendToGraph/W=ykGraph wy vs wx
					endfor
					break
				case "TY":					// Plot data as Y vs T
					for(i=0;i<n;i+=1)
						sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
						WAVE w = $wName
						AppendToGraph/W=ykGraph w
					endfor
					break
				default:
					Print "ERROR: GPIB response ", xymode, "invalid"
					err = -1
					break
			endswitch
			
			cmdStr = "NewDataFolder/O root:Data\r"
			cmdStr += "NewDataFolder/O " + pathName + "\r"
			
			Notebook ykPanel#Trace_Note selection={startOfFile,endOfFile}
			GetSelection notebook, ykPanel#Trace_Note, 3			
			
			// Scale and process data
			For(i=0;i<n;i+=1)
				sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
				WAVE w = $wName
				Execute/Q "root:Packages:GQ:ykDL750:scaleVar = " + tS[tw[i]-1][0]
				w /= scaleVar
				Note w, "Scaling: " + tS[tw[i]-1][0]
				Note w, S_Selection
				//Print tS[tw[i]-1][1]
				SetScale d,0,0,tS[tw[i]-1][1],w
//				Rename w, $(wName)
				//Duplicate/O tS, traceSettings
				cmdStr += "MoveWave " + wName + ", " + pathName + ":\r"
			EndFor
				
			// Make default graph style "ykGraphStyle"			
			if(graphStyle)
				Execute/Q "ykGraphStyle()"
			endif
			
			// Put legend if checkbox selected
			ControlInfo/W=ykPanel ShowLegend
			Legend/C/W=ykGraph/V=(V_value)/N=ykGraphLegend/A=RT/F=0/B=1/E=0

			// Save window recreation macro
			sprintf wName, "%s%04d", basename, fileIndex
			String/G winrec = WinRecreation("ykGraph",0)
			winrec = ReplaceString("ykGraph",winrec,wName)
			winrec = ReplaceString("DL750 Graph",winrec,wName)

			// Save data file
			if(numpnts(w) < kChunkSize)
				Save/O/J/B/T/W WaveList("*",";","TEXT:0") as (filePathStr + ":" + wName + ".itx")
			else
				NewPath/O/Q dataFolderPath, filePathStr
				Save/O/B/C/W/P=dataFolderPath WaveList("*",";","TEXT:0")
			endif
			
			// Save command
			// appendText
			tmpstr = RemoveListItem(0,winrec,"\r")
			tmpstr = RemoveListItem(0,tmpstr,"\r")
			tmpstr = RemoveListItem(ItemsInList(tmpstr,"\r")-1,tmpstr,"\r")
			tmpstr = ReplaceString("String",tmpstr,"String/G")
			tmpstr = ReplaceString("fldrSav0",tmpstr,"root:fldrSav0")
			tmpstr = ReplaceString("\t",tmpstr,"")
			cmdStr += tmpstr
			AddCommandsToItx(wName+".itx",filePathStr,cmdStr)
			//Print cmdStr
			
			// Increase fileIndex
			fileIndex += 1			
			break
		case -1: // control being killed
			break
	endswitch
	
	SetDataFolder savDF // Restore current DF.

	return err
End

Function ApplyGraph()
	Variable ind
	WAVE M_colors = M_colors
	For(ind=0;ind<8;ind+=1)
		ModifyGraph rgb[ind] = (M_colors[ind][0],M_colors[ind][1],M_colors[ind][2])
	EndFor
End

Function ShowLegendCheckProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			Legend/C/W=ykGraph/V=(checked)/N=ykGraphLegend/A=RT/F=0/B=1/E=0
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Window ykPanel() : Panel
   Initialize_ykDL750(0)
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(610,127,952,624) as "Yokogawa DL750 Scope"
	ShowTools/A
	SetDrawLayer UserBack
	SetDrawEnv fstyle= 1
	DrawText 64,25,"Data Information"
	SetDrawEnv fstyle= 1
	DrawText 43,138,"Trace Information"
	SetDrawEnv fstyle= 1
	DrawText 82,241,"Wave Note"
	Button LoadData,pos={14.00,442.00},size={218.00,45.00},proc=LoadDataButtonProc,title="Load Data"
	Button LoadData,fStyle=1
	SetVariable Basename,pos={8.00,414.00},size={60.00,18.00},title=" "
	SetVariable Basename,help={"Base name of data series"}
	SetVariable Basename,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:basename
	SetVariable FileIndex,pos={133.00,414.00},size={45.00,18.00},title=" "
	SetVariable FileIndex,format="%04d"
	SetVariable FileIndex,limits={0,9999,1},value= root:Packages:GQ:ykDL750:fileIndex
	SetVariable FilePath,pos={7.00,388.00},size={318.00,18.00},title=" "
	SetVariable FilePath,help={"File path to save data"}
	SetVariable FilePath,value= root:Packages:GQ:ykDL750:filePathStr
	CheckBox GraphStyle,pos={16.00,364.00},size={80.00,15.00},title="Graph Style?"
	CheckBox GraphStyle,variable= root:Packages:GQ:ykDL750:graphStyle,side= 1
	CheckBox ShowLegend,pos={8.00,344.00},size={88.00,15.00},proc=ShowLegendCheckProc,title="Show legend?"
	CheckBox ShowLegend,variable= root:Packages:GQ:ykDL750:showLegend,side= 1
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Edit/W=(9,33,328,216)/HOST=#  traceSettings.ld
	ModifyTable frameStyle= 5,format(Point)=1,width(traceSettings.l)=36,alignment(traceSettings.d)=1
	ModifyTable width(traceSettings.d)=82,width[2]=96,width[3]=32,width[4]=44
	ModifyTable showParts=0x65
	ModifyTable statsArea=85
	SetDataFolder fldrSav0
	RenameWindow #,T0
	SetActiveSubwindow ##
	NewNotebook /F=0 /N=Trace_Note /W=(7,245,328,329) /HOST=# 
	Notebook kwTopWin, defaultTab=20, autoSave= 1
	Notebook kwTopWin font="Arial", fSize=10, fStyle=0, textRGB=(0,0,0)
	Notebook kwTopWin, zdata= "GaqDU%ejN7!Z)%D?io>lbN?PWL]d_/WWX="
	Notebook kwTopWin, zdataEnd= 1
	RenameWindow #,Trace_Note
	SetActiveSubwindow ##
EndMacro

Window SweepScopePanel1() : Panel
	Initialize_ykDL750(0)
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1348,74,1598,619) as "Sweep and Scope"
	SetDrawLayer UserBack
	SetDrawEnv fstyle= 1
	DrawText 64,25,"Data Information"
	Button StartButton,pos={14,442},size={102,41},proc=ControlSweepProc,title="Start"
	Button StartButton,fStyle=1,fColor=(40960,65280,16384)
	SetVariable Basename,pos={8,257},size={31,16},title=" "
	SetVariable Basename,help={"Base name of data series"}
	SetVariable Basename,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:basename
	SetVariable FileIndex,pos={42,257},size={45,16},title=" ",format="%04d"
	SetVariable FileIndex,limits={0,9999,1},value= root:Packages:GQ:ykDL750:fileIndex
	SetVariable FilePath,pos={7,234},size={228,16},title=" "
	SetVariable FilePath,help={"File path to save data"}
	SetVariable FilePath,value= root:Packages:GQ:ykDL750:filePathStr
	SetVariable yoko7651address,pos={6,291},size={118,16},title=" Yoko 7651 GPIB"
	SetVariable yoko7651address,format="%2d"
	SetVariable yoko7651address,limits={0,21,1},value= root:Packages:GQ:ykDL750:yoko7651address
	CheckBox yoko7651currentsrc,pos={134,291},size={95,14},proc=UpdateCurrentSrcProc,title="Source Current?"
	CheckBox yoko7651currentsrc,variable= root:Packages:GQ:ykDL750:yoko7651currentsrc
	SetVariable yoko7651initialvalue,pos={11,321},size={106,16},proc=UpdateSweepParamsProc,title="Initial Value"
	SetVariable yoko7651initialvalue,format="%g"
	SetVariable yoko7651initialvalue,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651initialvalue
	SetVariable yoko7651finalvalue,pos={13,341},size={104,16},proc=UpdateSweepParamsProc,title="Final Value"
	SetVariable yoko7651finalvalue,format="%g"
	SetVariable yoko7651finalvalue,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651finalvalue
	SetVariable yoko7651step,pos={43,362},size={74,16},proc=UpdateSweepParamsProc,title="Step"
	SetVariable yoko7651step,format="%.3g"
	SetVariable yoko7651step,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651step
	SetVariable yoko7651numpnts,pos={26,383},size={92,16},proc=UpdateSweepParamsProc,title="# Points"
	SetVariable yoko7651numpnts,format="%d"
	SetVariable yoko7651numpnts,limits={1,inf,0},value= root:Packages:GQ:ykDL750:yoko7651numpnts
	Button PauseButton,pos={15,494},size={102,41},proc=ControlSweepProc,title="Pause"
	Button PauseButton,fStyle=1
	SetVariable yoko7651currpnt,pos={7,417},size={114,16},title="Current Point"
	SetVariable yoko7651currpnt,format="%g",fStyle=1
	SetVariable yoko7651currpnt,limits={1,inf,0},value= root:Packages:GQ:ykDL750:yoko7651currpnt,noedit= 1
	GroupBox EditTraceSettings,pos={4,285},size={232,127}
	PopupMenu ScaleAxis,pos={132,452},size={90,21},title="Scale Axis:"
	PopupMenu ScaleAxis,help={"Choose the channel number for scaling the data axes (0 for no scaling)"}
	PopupMenu ScaleAxis,mode=1,popvalue="0",value= #"\"0;\"+WaveToString(root:Packages:GQ:ykDL750:tracew)"
	CheckBox TransposeCheck,pos={133,480},size={74,14},title="Transpose?"
	CheckBox TransposeCheck,help={"Fast axis is acquired along x-direction. Check to transpose final matrix."}
	CheckBox TransposeCheck,value= 1,side= 1
	Button SwapButton,pos={124,329},size={50,20},proc=SwapButtonProc,title="Swap"
	CheckBox NoSweepCheck,pos={134,309},size={68,14},proc=NoSweepCheckProc,title="No Sweep"
	CheckBox NoSweepCheck,value= 0
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Edit/W=(9,33,220,216)/HOST=#  traceSettings.ld
	ModifyTable frameStyle= 5,format(Point)=1,width(traceSettings.l)=28,alignment(traceSettings.d)=1
	ModifyTable width(traceSettings.d)=82,width[2]=96,width[3]=32,width[4]=44
	ModifyTable showParts=0x45
	ModifyTable statsArea=85
	SetDataFolder fldrSav0
	RenameWindow #,T0
	SetActiveSubwindow ##
EndMacro

Window SweepScopePanel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(973,74,1272,671) as "Sweep and Scope"
	SetDrawLayer UserBack
	SetDrawEnv fstyle= 1
	DrawText 64,25,"Data Information"
	Button StartButton,pos={12.00,475.00},size={102.00,41.00},proc=ControlSweepProc,title="Start"
	Button StartButton,fStyle=1,fColor=(40960,65280,16384)
	SetVariable Basename,pos={8.00,257.00},size={31.00,18.00},title=" "
	SetVariable Basename,help={"Base name of data series"}
	SetVariable Basename,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:basename
	SetVariable FileIndex,pos={42.00,257.00},size={45.00,18.00},title=" "
	SetVariable FileIndex,format="%04d"
	SetVariable FileIndex,limits={0,9999,1},value= root:Packages:GQ:ykDL750:fileIndex
	SetVariable FilePath,pos={7.00,234.00},size={228.00,18.00},title=" "
	SetVariable FilePath,help={"File path to save data"}
	SetVariable FilePath,value= root:Packages:GQ:ykDL750:filePathStr
	SetVariable yoko7651address,pos={6.00,291.00},size={135.00,18.00},title=" Yoko 7651 GPIB"
	SetVariable yoko7651address,format="%2d"
	SetVariable yoko7651address,limits={0,21,1},value= root:Packages:GQ:ykDL750:yoko7651address
	CheckBox yoko7651currentsrc,pos={171.00,291.00},size={99.00,15.00},proc=UpdateCurrentSrcProc,title="Source Current?"
	CheckBox yoko7651currentsrc,variable= root:Packages:GQ:ykDL750:yoko7651currentsrc
	SetVariable yoko7651initialvalue,pos={11.00,321.00},size={124.00,18.00},proc=UpdateSweepParamsProc,title="Initial Value"
	SetVariable yoko7651initialvalue,format="%g"
	SetVariable yoko7651initialvalue,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651initialvalue
	SetVariable yoko7651finalvalue,pos={13.00,341.00},size={119.00,18.00},proc=UpdateSweepParamsProc,title="Final Value"
	SetVariable yoko7651finalvalue,format="%g"
	SetVariable yoko7651finalvalue,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651finalvalue
	SetVariable yoko7651step,pos={43.00,362.00},size={93.00,18.00},proc=UpdateSweepParamsProc,title="Step"
	SetVariable yoko7651step,format="%.3g"
	SetVariable yoko7651step,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651step
	SetVariable yoko7651numpnts,pos={26.00,383.00},size={109.00,18.00},proc=UpdateSweepParamsProc,title="# Points"
	SetVariable yoko7651numpnts,format="%d"
	SetVariable yoko7651numpnts,limits={1,inf,0},value= root:Packages:GQ:ykDL750:yoko7651numpnts
	Button PauseButton,pos={13.00,532.00},size={102.00,41.00},proc=ControlSweepProc,title="Pause"
	Button PauseButton,fStyle=1
	SetVariable yoko7651currpnt,pos={7.00,417.00},size={134.00,18.00},title="Current Point"
	SetVariable yoko7651currpnt,format="%g",fStyle=1
	SetVariable yoko7651currpnt,limits={1,inf,0},value= root:Packages:GQ:ykDL750:yoko7651currpnt,noedit= 1
	GroupBox EditTraceSettings,pos={4.00,282.00},size={290.00,130.00}
	PopupMenu ScaleAxis,pos={132.00,452.00},size={85.00,19.00},title="Scale Axis:"
	PopupMenu ScaleAxis,help={"Choose the channel number for scaling the data axes (0 for no scaling)"}
	PopupMenu ScaleAxis,mode=1,popvalue="0",value= #"\"0;\"+WaveToString(root:Packages:GQ:ykDL750:tracew)"
	CheckBox TransposeCheck,pos={133.00,480.00},size={73.00,15.00},title="Transpose?"
	CheckBox TransposeCheck,help={"Fast axis is acquired along x-direction. Check to transpose final matrix."}
	CheckBox TransposeCheck,value= 1,side= 1
	Button SwapButton,pos={170.00,329.00},size={50.00,20.00},proc=SwapButtonProc,title="Swap"
	CheckBox NoSweepCheck,pos={171.00,309.00},size={68.00,15.00},proc=NoSweepCheckProc,title="No Sweep"
	CheckBox NoSweepCheck,value= 0
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Edit/W=(9,33,282,216)/HOST=#  traceSettings.ld
	ModifyTable frameStyle= 5,format(Point)=1,width(traceSettings.l)=38,alignment(traceSettings.d)=1
	ModifyTable width(traceSettings.d)=82,width[2]=150,width[3]=35,width[4]=44
	ModifyTable showParts=0x45
	ModifyTable statsArea=85
	SetDataFolder fldrSav0
	RenameWindow #,T0
	SetActiveSubwindow ##
EndMacro

Window SweepScopePanel_Twice() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1668,74,1986,671) as "Sweep and Scope Twice"
	ShowTools/A
	SetDrawLayer UserBack
	SetDrawEnv fstyle= 1
	DrawText 64,25,"Data Information"
	Button StartButton,pos={13.00,545.00},size={102.00,41.00},proc=ControlSweepProc_Twice,title="Start"
	Button StartButton,fStyle=1,fColor=(40960,65280,16384)
	SetVariable Basename,pos={8.00,257.00},size={31.00,17.00},title=" "
	SetVariable Basename,help={"Base name of data series"}
	SetVariable Basename,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:basename
	SetVariable FileIndex,pos={42.00,257.00},size={45.00,17.00},title=" "
	SetVariable FileIndex,format="%04d"
	SetVariable FileIndex,limits={0,9999,1},value= root:Packages:GQ:ykDL750:fileIndex
	SetVariable FilePath,pos={7.00,234.00},size={228.00,17.00},title=" "
	SetVariable FilePath,help={"File path to save data"}
	SetVariable FilePath,value= root:Packages:GQ:ykDL750:filePathStr
	SetVariable yoko7651address,pos={6.00,291.00},size={135.00,17.00},title=" Yoko 7651 GPIB"
	SetVariable yoko7651address,format="%2d"
	SetVariable yoko7651address,limits={0,21,1},value= root:Packages:GQ:ykDL750:yoko7651address
	CheckBox yoko7651currentsrc,pos={171.00,291.00},size={104.00,14.00},proc=UpdateCurrentSrcProc,title="Source Current?"
	CheckBox yoko7651currentsrc,variable= root:Packages:GQ:ykDL750:yoko7651currentsrc
	SetVariable yoko7651initialvalue,pos={18.00,316.00},size={124.00,17.00},proc=UpdateSweepParamsProc,title="Initial Value"
	SetVariable yoko7651initialvalue,format="%g"
	SetVariable yoko7651initialvalue,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651initialvalue
	SetVariable yoko7651finalvalue,pos={24.00,338.00},size={119.00,17.00},proc=UpdateSweepParamsProc,title="Final Value"
	SetVariable yoko7651finalvalue,format="%g"
	SetVariable yoko7651finalvalue,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651finalvalue
	SetVariable yoko7651step,pos={54.00,361.00},size={93.00,17.00},proc=UpdateSweepParamsProc,title="Step"
	SetVariable yoko7651step,format="%.3g"
	SetVariable yoko7651step,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651step,noedit= 1
	SetVariable yoko7651numpnts,pos={53.00,496.00},size={109.00,17.00},proc=UpdateSweepParamsProc,title="# Points"
	SetVariable yoko7651numpnts,format="%d"
	SetVariable yoko7651numpnts,limits={1,inf,0},value= root:Packages:GQ:ykDL750:yoko7651numpnts
	Button PauseButton,pos={170.00,546.00},size={102.00,41.00},proc=ControlSweepProc,title="Pause"
	Button PauseButton,fStyle=1
	SetVariable yoko7651currpnt,pos={16.00,518.00},size={134.00,17.00},title="Current Point"
	SetVariable yoko7651currpnt,format="%g",fStyle=1
	SetVariable yoko7651currpnt,limits={1,inf,0},value= root:Packages:GQ:ykDL750:yoko7651currpnt,noedit= 1
	CheckBox TransposeCheck,pos={179.00,512.00},size={77.00,14.00},title="Transpose?"
	CheckBox TransposeCheck,help={"Fast axis is acquired along x-direction. Check to transpose final matrix."}
	CheckBox TransposeCheck,value= 1,side= 1
	Button SwapButton,pos={177.00,332.00},size={50.00,20.00},proc=SwapButtonProc,title="Swap"
	SetVariable yoko7651address2,pos={11.00,392.00},size={135.00,17.00},title=" Yoko 7651 GPIB"
	SetVariable yoko7651address2,format="%2d"
	SetVariable yoko7651address2,limits={0,21,1},value= root:Packages:GQ:ykDL750:yoko7651address2
	CheckBox yoko7651currentsrc2,pos={168.00,399.00},size={104.00,14.00},proc=UpdateCurrentSrcProc,title="Source Current?"
	CheckBox yoko7651currentsrc2,variable= root:Packages:GQ:ykDL750:yoko7651currentsrc2
	SetVariable yoko7651initialvalue2,pos={16.00,412.00},size={124.00,17.00},proc=UpdateSweepParamsProc,title="Initial Value"
	SetVariable yoko7651initialvalue2,format="%g"
	SetVariable yoko7651initialvalue2,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651initialvalue2
	SetVariable yoko7651finalvalue2,pos={22.00,433.00},size={119.00,17.00},proc=UpdateSweepParamsProc,title="Final Value"
	SetVariable yoko7651finalvalue2,format="%g"
	SetVariable yoko7651finalvalue2,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651finalvalue2
	SetVariable yoko7651step2,pos={55.00,453.00},size={93.00,17.00},proc=UpdateSweepParamsProc,title="Step"
	SetVariable yoko7651step2,format="%.3g"
	SetVariable yoko7651step2,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651step2,noedit= 1
	GroupBox EditTraceSettings2,pos={5.00,387.00},size={291.00,104.00}
	Button SwapButton3,pos={184.00,441.00},size={50.00,20.00},proc=SwapButtonProc,title="Swap"
	GroupBox EditTraceSettings4,pos={4.00,280.00},size={291.00,104.00}
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Edit/W=(9,33,282,216)/HOST=#  traceSettings.ld
	ModifyTable frameStyle= 5,format(Point)=1,width(traceSettings.l)=38,alignment(traceSettings.d)=1
	ModifyTable width(traceSettings.d)=82,width[2]=150,width[3]=35,width[4]=44
	ModifyTable showParts=0x45
	ModifyTable statsArea=85
	SetDataFolder fldrSav0
	RenameWindow #,T0
	SetActiveSubwindow ##
EndMacro

Window SweepScopeFluxPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1162,227,1613,1218) as "Sweep Scope Flux"
	SetDrawLayer UserBack
	SetDrawEnv fstyle= 1
	DrawText 64,25,"Data Information"
	Button StartButton,pos={10.00,552.00},size={102.00,41.00},proc=ControlSweepProc_Twice,title="Start"
	Button StartButton,fStyle=1,fColor=(40960,65280,16384)
	SetVariable Basename,pos={8.00,257.00},size={78.00,17.00},title=" "
	SetVariable Basename,help={"Base name of data series"}
	SetVariable Basename,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:basename
	SetVariable FileIndex,pos={133.00,257.00},size={71.00,17.00},title=" "
	SetVariable FileIndex,format="%04d"
	SetVariable FileIndex,limits={0,9999,1},value= root:Packages:GQ:ykDL750:fileIndex
	SetVariable FilePath,pos={7.00,234.00},size={314.00,17.00},title=" "
	SetVariable FilePath,help={"File path to save data"}
	SetVariable FilePath,value= root:Packages:GQ:ykDL750:filePathStr
	SetVariable yoko7651address,pos={6.00,291.00},size={135.00,17.00},title=" Yoko 7651 GPIB"
	SetVariable yoko7651address,format="%2d"
	SetVariable yoko7651address,limits={0,21,1},value= root:Packages:GQ:ykDL750:yoko7651address
	CheckBox yoko7651currentsrc,pos={171.00,291.00},size={104.00,14.00},proc=UpdateCurrentSrcProc,title="Source Current?"
	CheckBox yoko7651currentsrc,variable= root:Packages:GQ:ykDL750:yoko7651currentsrc
	SetVariable yoko7651initialvalue,pos={149.00,319.00},size={87.00,17.00},proc=UpdateFluxSweepParamsProc,title="A:"
	SetVariable yoko7651initialvalue,format="%g",valueColor=(65535,32768,32768)
	SetVariable yoko7651initialvalue,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651initialvalue,noedit= 1
	SetVariable yoko7651finalvalue,pos={148.00,338.00},size={88.00,17.00},proc=UpdateFluxSweepParamsProc,title="A:"
	SetVariable yoko7651finalvalue,format="%g",valueColor=(65535,32768,32768)
	SetVariable yoko7651finalvalue,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651finalvalue,noedit= 1
	SetVariable yoko7651step,pos={134.00,358.00},size={101.00,17.00},proc=UpdateFluxSweepParamsProc,title="Step"
	SetVariable yoko7651step,format="%.3g",valueColor=(65535,32768,32768)
	SetVariable yoko7651step,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651step,noedit= 1
	Button PauseButton,pos={327.00,553.00},size={102.00,41.00},proc=ControlSweepProc,title="Pause"
	Button PauseButton,fStyle=1
	SetVariable yoko7651currpnt,pos={141.00,573.00},size={148.00,17.00},title="Current Point"
	SetVariable yoko7651currpnt,format="%g",fStyle=1,valueColor=(65535,32768,32768)
	SetVariable yoko7651currpnt,limits={1,inf,0},value= root:Packages:GQ:ykDL750:yoko7651currpnt,noedit= 1
	CheckBox TransposeCheck,pos={27.00,517.00},size={77.00,14.00},title="Transpose?"
	CheckBox TransposeCheck,help={"Fast axis is acquired along x-direction. Check to transpose final matrix."}
	CheckBox TransposeCheck,value= 1,side= 1
	Button SwapButton,pos={293.00,287.00},size={50.00,20.00},proc=SwapButtonProc,title="Swap"
	SetVariable yoko7651address2,pos={11.00,392.00},size={135.00,17.00},title=" Yoko 7651 GPIB"
	SetVariable yoko7651address2,format="%2d"
	SetVariable yoko7651address2,limits={0,21,1},value= root:Packages:GQ:ykDL750:yoko7651address2
	CheckBox yoko7651currentsrc2,pos={169.00,397.00},size={104.00,14.00},proc=UpdateCurrentSrcProc,title="Source Current?"
	CheckBox yoko7651currentsrc2,variable= root:Packages:GQ:ykDL750:yoko7651currentsrc2
	GroupBox EditTraceSettings2,pos={5.00,388.00},size={440.00,122.00}
	Button SwapButton3,pos={297.00,394.00},size={50.00,20.00},proc=SwapButtonProc,title="Swap"
	GroupBox EditTraceSettings4,pos={4.00,280.00},size={442.00,105.00}
	SetVariable phizero1,pos={250.00,318.00},size={187.00,17.00},proc=UpdateFluxSweepParamsProc,title="Zero Phase Current"
	SetVariable phizero1,format="%g"
	SetVariable phizero1,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651phizero1
	SetVariable phi2pi1,pos={260.00,337.00},size={177.00,17.00},proc=UpdateFluxSweepParamsProc,title="2Pi Phase Current"
	SetVariable phi2pi1,format="%g"
	SetVariable phi2pi1,limits={1e-009,inf,0},value= root:Packages:GQ:ykDL750:yoko7651phi2pi1
	SetVariable phi01,pos={15.00,319.00},size={124.00,17.00},proc=UpdateFluxSweepParamsProc,title="Initial Phase"
	SetVariable phi01,format="%g"
	SetVariable phi01,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651phi01
	SetVariable phi11,pos={21.00,338.00},size={117.00,17.00},proc=UpdateFluxSweepParamsProc,title="Final Phase"
	SetVariable phi11,format="%g"
	SetVariable phi11,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651phi11
	SetVariable yoko7651initialvalue2,pos={151.00,424.00},size={87.00,17.00},proc=UpdateFluxSweepParamsProc,title="A:"
	SetVariable yoko7651initialvalue2,format="%g",valueColor=(65535,32768,32768)
	SetVariable yoko7651initialvalue2,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651initialvalue2,noedit= 1
	SetVariable yoko7651finalvalue2,pos={150.00,443.00},size={88.00,17.00},proc=UpdateFluxSweepParamsProc,title="A:"
	SetVariable yoko7651finalvalue2,format="%g",valueColor=(65535,32768,32768)
	SetVariable yoko7651finalvalue2,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651finalvalue2,noedit= 1
	SetVariable yoko7651step2,pos={137.00,463.00},size={101.00,17.00},proc=UpdateFluxSweepParamsProc,title="Step"
	SetVariable yoko7651step2,format="%.3g",valueColor=(65535,32768,32768)
	SetVariable yoko7651step2,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651step2,noedit= 1
	SetVariable phi02,pos={18.00,426.00},size={124.00,17.00},proc=UpdateFluxSweepParamsProc,title="Initial Phase"
	SetVariable phi02,format="%g"
	SetVariable phi02,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651phi02
	SetVariable phi12,pos={24.00,445.00},size={117.00,17.00},proc=UpdateFluxSweepParamsProc,title="Final Phase"
	SetVariable phi12,format="%g"
	SetVariable phi12,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651phi12
	SetVariable phizero2,pos={245.00,426.00},size={190.00,17.00},proc=UpdateFluxSweepParamsProc,title="Zero Phase Current"
	SetVariable phizero2,format="%g"
	SetVariable phizero2,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651phizero2
	SetVariable phi2pi2,pos={291.00,445.00},size={144.00,17.00},proc=UpdateFluxSweepParamsProc,title="2Pi Current"
	SetVariable phi2pi2,format="%g"
	SetVariable phi2pi2,limits={1e-009,inf,0},value= root:Packages:GQ:ykDL750:yoko7651phi2pi2
	SetVariable sinamp,pos={269.00,486.00},size={168.00,17.00},proc=UpdateFluxSweepParamsProc,title="Sine Amplitude"
	SetVariable sinamp,format="%g"
	SetVariable sinamp,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651sinamp
	SetVariable gamma12,pos={271.00,466.00},size={165.00,17.00},proc=UpdateFluxSweepParamsProc,title="Mutual gamma"
	SetVariable gamma12,format="%g"
	SetVariable gamma12,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651gamma12
	SetVariable yoko7651numpnts,pos={179.00,549.00},size={109.00,17.00},proc=UpdateFluxSweepParamsProc,title="# Points"
	SetVariable yoko7651numpnts,format="%d"
	SetVariable yoko7651numpnts,limits={1,inf,0},value= root:Packages:GQ:ykDL750:yoko7651numpnts
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Edit/W=(9,33,438,216)/HOST=#  traceSettings.ld
	ModifyTable frameStyle= 5,format(Point)=1,width(traceSettings.l)=38,alignment(traceSettings.d)=1
	ModifyTable width(traceSettings.d)=82,width[2]=250,width[3,4]=50
	ModifyTable showParts=0x45
	ModifyTable statsArea=85
	SetDataFolder fldrSav0
	RenameWindow #,T0
	SetActiveSubwindow ##
	String fldrSav1= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Display/W=(8,607,443,983)/HOST=#  yoko7651loopwave2 vs yoko7651loopwave
	SetDataFolder fldrSav1
	Label left "Current Yoko 2 (\\U)"
	Label bottom "Current Yoko 1 (\\U)"
	ModifyGraph mode=3
	ModifyGraph marker=19
	ModifyGraph rgb=(39168,0,15616)
	ModifyGraph msize=1
	RenameWindow #,G0
	SetActiveSubwindow ##
EndMacro


Window SweepScopePanel_R() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(2074,387,2373,984) as "Sweep and Scope RIGOL"
	SetDrawLayer UserBack
	SetDrawEnv fstyle= 1
	DrawText 64,25,"Data Information"
	Button StartButton,pos={13.00,475.00},size={102.00,41.00},proc=ControlRigolSweepProc,title="Start"
	Button StartButton,fStyle=1,fColor=(40960,65280,16384)
	SetVariable Basename,pos={8.00,257.00},size={31.00,17.00},title=" "
	SetVariable Basename,help={"Base name of data series"}
	SetVariable Basename,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:basename
	SetVariable FileIndex,pos={42.00,257.00},size={99.00,17.00},title=" "
	SetVariable FileIndex,format="%04d"
	SetVariable FileIndex,limits={0,9999,1},value= root:Packages:GQ:ykDL750:fileIndex
	SetVariable FilePath,pos={7.00,233.00},size={280.00,17.00},title=" "
	SetVariable FilePath,help={"File path to save data"}
	SetVariable FilePath,value= root:Packages:GQ:ykDL750:filePathStr
	SetVariable rigolch1amp,pos={8.00,327.00},size={124.00,17.00},proc=UpdateRigolSweepParamsProc,title="Amplitude CH1"
	SetVariable rigolch1amp,format="%g"
	SetVariable rigolch1amp,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:rigolCH1amp
	SetVariable rigolnumpnts,pos={27.00,378.00},size={109.00,17.00},proc=UpdateRigolSweepParamsProc,title="# Points"
	SetVariable rigolnumpnts,format="%d"
	SetVariable rigolnumpnts,limits={1,inf,0},value= root:Packages:GQ:ykDL750:rigolnumpnts
	Button PauseButton,pos={13.00,532.00},size={102.00,41.00},proc=ControlSweepProc,title="Pause"
	Button PauseButton,fStyle=1
	SetVariable rigolcurrpnt,pos={7.00,417.00},size={134.00,17.00},title="Current Point"
	SetVariable rigolcurrpnt,format="%g",fStyle=1
	SetVariable rigolcurrpnt,limits={1,inf,0},value= root:Packages:GQ:ykDL750:rigolcurrpnt,noedit= 1
	SetVariable rigolch2amp,pos={9.00,349.00},size={124.00,17.00},proc=UpdateRigolSweepParamsProc,title="Amplitude CH2"
	SetVariable rigolch2amp,format="%g"
	SetVariable rigolch2amp,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:rigolCH2amp
	SetVariable rigoladdress,pos={7.00,299.00},size={280.00,17.00},title=" Rigol Address"
	SetVariable rigoladdress,help={"File path to save data"}
	SetVariable rigoladdress,value= root:Packages:GQ:ykDL750:rigoladdress
	SetWindow kwTopWin,userdata(ACL_desktopNum)=  "11"
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Edit/W=(9,33,282,216)/HOST=#  traceSettings.ld
	ModifyTable frameStyle= 5,format(Point)=1,width(traceSettings.l)=38,alignment(traceSettings.d)=1
	ModifyTable width(traceSettings.d)=82,width[2]=150,width[3]=35,width[4]=44
	ModifyTable showParts=0x45
	ModifyTable statsArea=85
	SetDataFolder fldrSav0
	RenameWindow #,T0
	SetActiveSubwindow ##
EndMacro



Function UpdateSweepParamsProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	DFREF saveDFR = GetDataFolderDFR()
	
	Variable status = 0
	
	SetDataFolder root:Packages:GQ:ykDL750
	
	NVAR initval = yoko7651initialvalue
	NVAR finalval = yoko7651finalvalue
	NVAR step = yoko7651step
	NVAR npnts = yoko7651numpnts
	
	WAVE loopwave = yoko7651loopwave

	NVAR initval2 = yoko7651initialvalue2
	NVAR finalval2 = yoko7651finalvalue2
	NVAR step2 = yoko7651step2
	
	WAVE loopwave2 = yoko7651loopwave2

	if(initval == finalval)
		initval = 0
		finalval = 3e-4
		status = -1
	elseif(step == 0)
		step = 1e-6
		status = -1
	else

		switch( sva.eventCode )
			case 1: // mouse up
			case 2: // Enter key
			case 3: // Live update
				strswitch(sva.vName)
					case "yoko7651initialvalue":
					case "yoko7651finalvalue":
					case "yoko7651step":
						npnts = round(1+abs((finalval-initval)/step))
						step = (finalval - initval)/(npnts-1)
						break
					case "yoko7651numpnts":
						step = (finalval - initval)/(npnts-1)
						step2 = (finalval2 - initval2)/(npnts-1)
						break
					case "yoko7651initialvalue2":
					case "yoko7651finalvalue2":
						step2 = (finalval2 - initval2)/(npnts-1)
						break
					default:
						break
				endswitch				
				break
			case -1: // control being killed
				break
		endswitch
	
		Redimension/N=(npnts) loopwave, loopwave2
//		ControlInfo yoko7651currentsrc
//		if(V_Value)
//			SetScale d,0,0,"A",loopwave
//		else
//			SetScale d,0,0,"V",loopwave
//		endif
//			
		loopwave = initval + p*step
		loopwave2 = initval2 + p*step2
	endif

	SetDataFolder saveDFR

	return status
End

Function UpdateFluxSweepParamsProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	DFREF saveDFR = GetDataFolderDFR()
	
	Variable status = 0
	Variable currperiod1, currperiod2
	
	SetDataFolder root:Packages:GQ:ykDL750
	
	// Linear part
	
	NVAR initval = yoko7651initialvalue
	NVAR finalval = yoko7651finalvalue
	NVAR step = yoko7651step
	NVAR npnts = yoko7651numpnts

	NVAR initval2 = yoko7651initialvalue2
	NVAR finalval2 = yoko7651finalvalue2
	NVAR step2 = yoko7651step2


	// Conversions from current to flux

	NVAR phizero1 = yoko7651phizero1		// Bias current offset for zero flux
	NVAR phi2pi1 = yoko7651phi2pi1  		// Bias current equivalent to one flux quantum
		
	NVAR phizero2 = yoko7651phizero2
	NVAR phi2pi2 = yoko7651phi2pi2

	// Mutual coupling
	// Flux compensation only applied to yoko 2
	// Usually should be negative
	NVAR gamma12 = yoko7651gamma12 	// Mutual coupling: slope given by (compensating phase shift for yoko2)/(phase shift yoko1)
	

	// Sinusoidal part, leave sinamp = 0 for no nonlinearity
	// Amplitude for sinusoidal part, always applied to yoko 7651 2
	// Units of current!
	NVAR sinamp = yoko7651sinamp
	

	// Desired sweep values, in units of 2pi flux
	NVAR phi01 = yoko7651phi01
	NVAR phi11 = yoko7651phi11
	
	NVAR phi02 = yoko7651phi02
	NVAR phi12 = yoko7651phi12


	// Waves containing current sweep values	

	WAVE loopwave = yoko7651loopwave
	WAVE loopwave2 = yoko7651loopwave2

	
	// some error correction	
	if(phi2pi1 == 0)
		phi2pi1 = 10e-6
		status = -1
	endif
	
	if(phi2pi2 == 0)
		phi2pi1 = 10e-6
		status = -1
	endif
	
	if(npnts == 0)
		npnts = 8
		status = -1
	endif


	Redimension/N=(npnts) loopwave, loopwave2
	
	// Process changes to parameters
	
	step = (phi11-phi01)/(npnts-1)*phi2pi1
	initval = phizero1 + phi01*phi2pi1
	finalval = phizero1 + phi11*phi2pi1
			
	// Yoko 2 linear part 
	step2 = (phi12-phi02)/(npnts-1)*phi2pi2
	initval2 = phizero2 + phi02*phi2pi2
	finalval2 = phizero2 + phi12*phi2pi2
			
	// Yoko 2 flux compensation
	step2 += gamma12*(phi11-phi01)/(npnts-1)*phi2pi1
	initval2 += gamma12*phi01/(npnts-1)*phi2pi1
	finalval2 += gamma12*phi11/(npnts-1)*phi2pi1
			
	ControlInfo yoko7651currentsrc
	if(V_Value)
		SetScale d,0,0,"A",loopwave, loopwave2
		
	else
		SetScale d,0,0,"V",loopwave, loopwave2
	endif
			
	loopwave = initval + p*step

	// Old linear sweep:
	// loopwave2 = initval2 + p*step2
	
	// New sinusoidal sweep:
	//loopwave2 = initval2 + p*step2 + sinamp*sin(2*pi*(p*step2+initval2-phizero2)/phi2pi2)
	loopwave2 = initval2 + p*step2 + sinamp*sin(2*pi*(initval + p*step)/phi2pi2)
	//loopwave2 = initval2 + p*step2 + sinamp*sin(2*pi*(initval + p*step)/phi2pi1)
	
	SetDataFolder saveDFR

	return status
End

Function UpdateRigolSweepParamsProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	DFREF saveDFR = GetDataFolderDFR()
	
	Variable status = 0
	
	SetDataFolder root:Packages:GQ:ykDL750
	
	NVAR ch1amp = rigolCH1amp
	NVAR ch2amp = rigolCH2amp
	NVAR npnts = rigolnumpnts
	
	WAVE loopwave = rigolloopwave

	if(npnts <= 0)
		npnts = 0
	else
		switch( sva.eventCode )
			case 1: // mouse up
			case 2: // Enter key
			case 3: // Live update
			case -1: // control being killed
				break
		endswitch
	
		Redimension/N=(npnts,2) loopwave
		loopwave[][0] = ch1amp*cos(p/npnts*2*pi)
		loopwave[][1] = ch2amp*sin(p/npnts*2*pi)
//		ControlInfo yoko7651currentsrc
//		if(V_Value)
//			SetScale d,0,0,"A",loopwave
//		else
//			SetScale d,0,0,"V",loopwave
//		endif
//			
//		loopwave = initval + p*step
	endif

	SetDataFolder saveDFR

	return status
End


Function TransferAllRecords([scaleChannel,maxRecord])
	Variable scaleChannel
	Variable maxRecord
	
	DFREF saveDFR = GetDataFolderDFR()
	
	SetDataFolder root:Packages:GQ:ykDL750
	
	WAVE tw = tracew
	WAVE/T tS = traceSettings
	
	SVAR address, basename, filePathStr
	NVAR fileIndex, scaleVar
	
	String pathName, wName
	
	Variable numrecords = abs(str2num(ykQuery(address,"WAV:REC? MIN")))+1
	Variable numiv = str2num(ykQuery(address,"WAV:LENG?"))
	
	if(!ParamIsDefault(maxRecord))
		Print "Maximum number of records:", numrecords

		if(maxRecord <= numrecords)
			numrecords = maxRecord
		endif
	endif
	
	Print "Records to transfer:", numrecords
	
	ykAvailableTraces(address,tw)
	 
	Variable numtraces = numpnts(tw)
	
	sprintf pathName, "root:Data:%s%04d", basename, fileIndex
	fileIndex = CheckDataPath(pathName,basename,fileIndex)
	NewDataFolder/O $pathName
	NewPath/O/Q dataFolderPath, filePathStr
	
	Variable i,j
	
	Notebook ykPanel#Trace_Note selection={startOfFile,endOfFile}
	GetSelection notebook, ykPanel#Trace_Note, 3		

	printf "Transferring"
	for(i = 0;i<numrecords;i+=1)
		for(j = 0;j<numtraces;j+=1)
			ykGetTrace(address, tw[j],"tmpwave",record=i)
			WAVE tmpw = tmpwave
			Execute/Q "root:Packages:GQ:ykDL750:scaleVar = " + tS[tw[j]-1][0]  // scaling could be performed more efficiently
			FastOp tmpw = (1/scalevar)*tmpw

			sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[j]-1][2]

			if(i==0)
				Make/O/N=(numiv,numrecords) $wName
				CopyScales tmpw, $wName
				SetScale d,0,0,tS[tw[j]-1][1],$wName
				
				if(!ParamIsDefault(scaleChannel))
					Note $wName, "Scaling Axis: " + tS[scaleChannel-1][2]
				endif
				
				Note $wName, S_Selection
			endif

			WAVE w = $wName
			w[][i] = tmpw[p]
		endfor
		
		if(mod(i,10) == 0)
			printf "."
		endif
	endfor
	printf "\r"
	
	KillWaves/Z tmpwave

	if(!ParamIsDefault(scaleChannel))
		for(j = 0;j<numtraces;j+=1)
			sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[j]-1][2]
			WAVE wy = $wName
			sprintf wName, "%s%04d_%s", basename, fileIndex, tS[scaleChannel-1][2]
			WAVE wx = $wName
			
			ScaleWave(wx,wy)
		endfor
	endif
	
	sprintf pathName,"%s%04d.pxp", basename, fileIndex
	SaveData /O/Q/R/P=dataFolderPath pathName
	
	fileindex += 1
	
	SetDataFolder saveDFR
	
	return 0
End

Function SaveDataFolder(basename,fileIndex,dataFolderPathStr)
	String basename
	Variable fileIndex
	String dataFolderPathStr
	
	String pathName
	
	sprintf pathName, "root:Data:%s%04d", basename, fileIndex
	
	SaveData /O/Q/R/P=$dataFolderPathStr pathName
	
	return V_Flag
End

Function ScaleWave(wx,wy)
	WAVE wx
	WAVE wy
	
	Make/FREE/N=(DimSize(wx,0)) wcut
	
	wcut = wx
	
	StatsLinearRegression wcut
	WAVE wlinreg = W_StatsLinearRegression
	Variable offset = wlinreg[1]
	Variable	dx = wlinreg[2]
	KillWaves/Z W_StatsLinearRegression
						
	SetScale/P x,offset,dx,WaveUnits(wx,-1),wy
	
	return 0
End

//
//Function NewDataFolderAndPath(basename,fileIndex,filePathStr)
//	String basename
//	Variable fileIndex
//	String filePathStr
//
//	String pathName
//
//	
//	return pathName
//End


Function ControlSweepProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	DFREF saveDFR = GetDataFolderDFR()
	
	Variable err = 0, status
	
	SetDataFolder root:Packages:GQ:ykDL750

	SVAR address, basename, filePathStr
	NVAR fileIndex
	
	WAVE tw = tracew
	WAVE/T tS = traceSettings
	
	WAVE loopwave = yoko7651loopwave

	NVAR initval = yoko7651initialvalue
	NVAR finalval = yoko7651finalvalue
	NVAR step = yoko7651step
	NVAR npnts = yoko7651numpnts
	NVAR currentsrc = yoko7651currentsrc
	NVAR currpnt = yoko7651currpnt
	NVAR yoko7651address

	String pathName, wName
	Variable i, n
	Variable numiv, trigpos, dt
					
	switch( ba.eventCode )
		case 2: // mouse up 
			// click code here
			
			strswitch(ba.ctrlName)
				case "StartButton":
					// Check if yoko is not in Normal mode
					if(!ykReady())
						DoAlert 0, "Must set Yoko DL750 to AVERAGE or SINGLE mode!"
						err = -1
						break
					endif
					
					ControlInfo/W=SweepScopePanel NoSweepCheck
					if(V_Value==0)
						// DO NOT Put Yoko7651 in proper range
						//Yoko7651SetRange("GPIB0::"+num2str(yoko7651address)+"::INSTR", wavemax(loopwave),current=currentsrc)
			
						// Set Initial Value for Yoko7651
						Yoko7651SetOutput(kGPIBn+"::"+num2str(yoko7651address)+"::INSTR", loopwave[currpnt])
					endif


					// Prepare waves to be written
					//
					//
					ykAvailableTraces(address,tw)
//					numiv = str2num(ykQuery(address,"WAV:LENG?"))
					numiv = str2num(ykQuery(address,"ACQ:RLEN?")) 
					trigpos = str2num(ykQuery(address,"WAVeform:TRIGGER?"))
					dt = 1/str2num(ykQuery(address,"TIMEBASE:SRATE?"))
					n = numpnts(tw)				

					// Path for new data
					sprintf pathName, "root:Data:%s%04d", basename, fileIndex
					fileIndex = CheckDataPath(pathName,basename,fileIndex)
					NewDataFolder/O $pathName
					NewPath/O/Q dataFolderPath, filePathStr
				
					// Set current sweep point to 0
					currpnt = 0
			
					// Create empy data
					ControlInfo/W=SweepScopePanel ScaleAxis
					Variable scaleIndex = str2num(S_Value)
					ControlInfo/W=SweepScopePanel NoSweepCheck
					
																
					for(i = 0;i<n;i+=1)
						sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
						Make/O/N=(numiv,numpnts(loopwave)) $(pathName+":"+wName)
						WAVE w = $(pathName+":"+wName)
						w = NaN
						SetScale d,0,0,tS[tw[i]-1][1],w
						SetScale/P x,-trigpos*dt,dt,"s",w
						
						if(V_Value == 0)
							SetScale/P y,loopwave[0],loopwave[1]-loopwave[0],WaveUnits(loopwave,-1),w
						endif
						
						Note w, "Scaling: " + tS[tw[i]-1][0]
						Note w, "Scaling Axis: " + ts[scaleIndex-1][2]
					endfor
			
					// Launch background process
					CtrlNamedBackground SweepScope, period=60, proc=SweepScopeTask
					CtrlNamedBackground SweepScope, start
			
					// Change enable state of buttons
					Button $ba.ctrlName,rename=StopButton,title="Stop",fColor=(65280,0,0),win=SweepScopePanel
					//Button StartSweep,disable=2,win=SweepScopePanel
					//Button StopSweep,disable=0,win=SweepScopePanel,fColor=(40960,65280,16384)
			
					break
				case "StopButton":
				
					// Stop background process
					CtrlNamedBackground SweepScope, stop
			
					// Redimension data and transpose if requested
					ControlInfo/W=SweepScopePanel TransposeCheck
										
					n = numpnts(tw)
					i = 0
					
					sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
					sprintf pathName, "root:Data:%s%04d", basename, fileIndex
					WAVE w = $(pathName+":"+wName)
					
					numiv = DimSize(w,0)
					
					for(i = 0;i<n;i+=1)
						sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
						WAVE w = $(pathName+":"+wName)
						Redimension/N=(numiv,currpnt-1) w
						
						if(V_Value)
							MatrixTranspose w
						endif
					endfor
			
					// set currpnt to NaN
					currpnt = NaN
	
					// Change buttons to default state
					Button $ba.ctrlName,rename=StartButton,title="Start",fColor=(40960,65280,16384),win=SweepScopePanel
					Button PauseButton,valueColor=(0,0,0),win=SweepScopePanel
					
					// Save data file
					sprintf wName, "%s%04d", basename, fileIndex
					sprintf pathName, "root:Data:%s%04d", basename, fileIndex

					SetDataFolder $pathName
					
//					// Kill temporary wave if it exists
//					KillWaves/Z tmpwave
//					
					Save/O/B/C/W/P=dataFolderPath WaveList("*",";","DIMS:2,TEXT:0")

					// Increment file Index
					fileIndex += 1
						
					break
				case "PauseButton":
					if(numtype(currpnt) == 0)
						CtrlNamedBackground SweepScope, status
						if(NumberByKey("RUN",S_Info))
							CtrlNamedBackground SweepScope, stop
							Button PauseButton,valueColor=(65280,0,0),win=SweepScopePanel
						else
							CtrlNamedBackground SweepScope, start
							Button PauseButton,valueColor=(0,0,0),win=SweepScopePanel							
						endif
					endif
					break							
			endswitch
			
			break
		case -1: // control being killed
			break
	endswitch

	SetDataFolder saveDFR // Restore current DF.

	return err
End



Function RigolSetOutputs(rigoladdress, ch1amp, ch2amp)
	String rigoladdress
	Variable ch1amp, ch2amp
	
	String cmd, cmdstr, outstr
	Variable numerr = 0, inv1, inv2

	// Check polarities
	cmdstr = ":OUTP1:POL?"
	inv1 = StringMatch(QueryVISA(rigoladdress,cmdstr),"INVERTED*")

	cmdstr = ":OUTP2:POL?"
	inv2 = StringMatch(QueryVISA(rigoladdress,cmdstr),"INVERTED*")



	// Treat CH1

	if(ch1amp < 0) // smallest possible output amplitude is 2mV
		cmd = "LOW"
		
		if(!inv1) // need to invert output CH1
			cmdstr = ":SOUR1:VOLT:LOW -0.002;:SOUR1:VOLT:HIGH 0"
			numerr +=	WRITEVisa(rigoladdress,cmdstr)

			cmdstr = ":OUTP1:POL INV"
			numerr +=	WRITEVisa(rigoladdress,cmdstr)
		endif
	else
		cmd = "HIGH"

		if(inv1) // need to uninvert output CH1
			cmdstr = ":SOUR1:VOLT:HIGH 0.002;:SOUR1:VOLT:LOW 0"
			numerr +=	WRITEVisa(rigoladdress,cmdstr)

			cmdstr = ":OUTP1:POL NORM"
			numerr +=	WRITEVisa(rigoladdress,cmdstr)
		endif

	endif 
	
	cmdstr = ":SOUR1:VOLT:"+cmd+" "+num2str(ch1amp)
	numerr = WRITEVisa(rigoladdress,cmdstr)
	
	
	// Treat CH2
	
	if(ch2amp < 0) // smallest possible output amplitude is 2mV
		cmd = "LOW"
		
		if(!inv2) // need to invert output CH1
			cmdstr = ":SOUR2:VOLT:LOW -0.002;:SOUR2:VOLT:HIGH 0"
			numerr +=	WRITEVisa(rigoladdress,cmdstr)

			cmdstr = ":OUTP2:POL INV"
			numerr +=	WRITEVisa(rigoladdress,cmdstr)
		endif
	else
		cmd = "HIGH"

		if(inv2) // need to uninvert output CH1
			cmdstr = ":SOUR2:VOLT:HIGH 0.002;:SOUR2:VOLT:LOW 0"
			numerr +=	WRITEVisa(rigoladdress,cmdstr)

			cmdstr = ":OUTP2:POL NORM"
			numerr +=	WRITEVisa(rigoladdress,cmdstr)
		endif

	endif
	
	cmdstr = ":SOUR2:VOLT:"+cmd+" "+num2str(ch2amp)
	numerr += WRITEVisa(rigoladdress,cmdstr)
	
	return numerr
End
	
	
// Function to sweep Rigol high values
Function ControlRigolSweepProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	DFREF saveDFR = GetDataFolderDFR()
	
	Variable err = 0, status
	
	SetDataFolder root:Packages:GQ:ykDL750 

	SVAR address, basename, filePathStr
	NVAR fileIndex
	
	WAVE tw = tracew
	WAVE/T tS = traceSettings
	 
	WAVE loopwave = rigolloopwave

	NVAR npnts = rigolnumpnts
	NVAR currpnt = rigolcurrpnt
	SVAR rigoladdress

	String pathName, wName
	Variable i, n
	Variable numiv, trigpos, dt
					
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			
			strswitch(ba.ctrlName)
				case "StartButton":
					// Check if yoko is not in Normal mode
					if(!ykReady())
						DoAlert 0, "Must set Yoko DL750 to AVERAGE or SINGLE mode!"
						err = -1
						break
					endif
					
					// ControlInfo/W=SweepScopePanel_R NoSweepCheck
					// if(V_Value==0)
					// 	// Set Initial Value for Yoko7651
					// 	RigolSetOutputs(rigoladdress, loopwave[currpnt][0], loopwave[currpnt][1])
					// endif


					// Prepare waves to be written
					//
					//
					ykAvailableTraces(address,tw)
//					numiv = str2num(ykQuery(address,"WAV:LENG?"))
					numiv = str2num(ykQuery(address,"ACQ:RLEN?")) 
					trigpos = str2num(ykQuery(address,"WAVeform:TRIGGER?"))
					dt = 1/str2num(ykQuery(address,"TIMEBASE:SRATE?"))
					n = numpnts(tw)				

					// Path for new data
					sprintf pathName, "root:Data:%s%04d", basename, fileIndex
					fileIndex = CheckDataPath(pathName,basename,fileIndex)
					NewDataFolder/O $pathName
					NewPath/O/Q dataFolderPath, filePathStr
				
					// Set current sweep point to 0
					currpnt = 0
			
					// Create empy data
					ControlInfo/W=SweepScopePanel_R ScaleAxis
					Variable scaleIndex = str2num(S_Value)
																										
					for(i = 0;i<n;i+=1)
						sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
						Make/O/N=(numiv,numpnts(loopwave)) $(pathName+":"+wName)
						WAVE w = $(pathName+":"+wName)
						w = NaN
						SetScale d,0,0,tS[tw[i]-1][1],w
						SetScale/P x,-trigpos*dt,dt,"s",w
						
						SetScale/P y,loopwave[0],loopwave[1]-loopwave[0],WaveUnits(loopwave,-1),w
						
						Note w, "Scaling: " + tS[tw[i]-1][0]
						Note w, "Scaling Axis: " + ts[scaleIndex-1][2]
					endfor
			
					// Launch background process
					CtrlNamedBackground SweepScope_R, period=60, proc=SweepScopeTask_R
					CtrlNamedBackground SweepScope_R, start
			
					// Change enable state of buttons
					Button $ba.ctrlName,rename=StopButton,title="Stop",fColor=(65280,0,0),win=SweepScopePanel_R
			
					break
				case "StopButton":
				
					// Stop background process
					CtrlNamedBackground SweepScope_R stop
			
					n = numpnts(tw)
					i = 0
					
					sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
					sprintf pathName, "root:Data:%s%04d", basename, fileIndex
					WAVE w = $(pathName+":"+wName)
					
					numiv = DimSize(w,0)
					
					for(i = 0;i<n;i+=1)
						sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
						WAVE w = $(pathName+":"+wName)
						Redimension/N=(numiv,currpnt-1) w
					endfor
			
					// set currpnt to NaN
					currpnt = NaN
	
					// Change buttons to default state
					Button $ba.ctrlName,rename=StartButton,title="Start",fColor=(40960,65280,16384),win=SweepScopePanel_R
					Button PauseButton,valueColor=(0,0,0),win=SweepScopePanel_R
					
					// Save data file
					sprintf wName, "%s%04d", basename, fileIndex
					sprintf pathName, "root:Data:%s%04d", basename, fileIndex

					SetDataFolder $pathName
					
//					// Kill temporary wave if it exists
//					KillWaves/Z tmpwave
//					
					Save/O/B/C/W/P=dataFolderPath WaveList("*",";","DIMS:2,TEXT:0")

					// Increment file Index
					fileIndex += 1
						
					break
				case "PauseButton":
					if(numtype(currpnt) == 0)
						CtrlNamedBackground SweepScope_R, status
						if(NumberByKey("RUN",S_Info))
							CtrlNamedBackground SweepScope_R, stop
							Button PauseButton,valueColor=(65280,0,0),win=SweepScopePanel_R
						else
							CtrlNamedBackground SweepScope_R, start
							Button PauseButton,valueColor=(0,0,0),win=SweepScopePanel_R					
						endif
					endif
					break							
			endswitch
			
			break
		case -1: // control being killed
			break
	endswitch

	SetDataFolder saveDFR // Restore current DF.

	return err
End

Function SweepScopeTask_R(s)
	STRUCT WMBackgroundStruct &s
	
	if(ykBusy())
		return 0
	endif

	NVAR currpnt = root:Packages:GQ:ykDL750:rigolcurrpnt
	NVAR ch1amp = root:Packages:GQ:ykDL750:rigolCH1amp
	NVAR ch2amp = root:Packages:GQ:ykDL750:rigolCH2amp
	NVAR npnts = root:Packages:GQ:ykDL750:rigolnumpnts
	SVAR rigoladdress = root:Packages:GQ:ykDL750:rigoladdress
	WAVE lw = root:Packages:GQ:ykDL750:rigolloopwave
	
	SVAR address = root:Packages:GQ:ykDL750:address
	WAVE tw = root:Packages:GQ:ykDL750:tracew
	WAVE/T tS = root:Packages:GQ:ykDL750:traceSettings
	
	SVAR baseName = root:Packages:GQ:ykDL750:basename
	NVAR fileIndex = root:Packages:GQ:ykDL750:fileIndex
	
	String pathName
	sprintf pathName, "root:Data:%s%04d", basename, fileIndex

	String wName
	Variable i, n = numpnts(tw)

	// Change Rigol outputs for current point in sweep
	RigolSetOutputs(rigoladdress,lw[currpnt][0],lw[currpnt][1])
		
	if(currpnt==0)  // First time run, nothing to do

	else  // transfer data
		for(i = 0;i<n;i+=1)
			sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
			ykGetTrace(address, tw[i],"tmpwave")
			WAVE tmpw = tmpwave
			WAVE w = $(pathName+":"+wName)
			// Vestiges of scaling from previous versions:
			//			Execute/Q "root:Packages:GQ:ykDL750:scaleVar = " + tS[tw[i]-1][0]  // scaling could be performed more efficiently
			//			FastOp tmpw = (1/scalevar)*tmpw
			w[][currpnt-1] = tmpw[p]
		endFor
		
		KillWaves/Z tmpwave
	endif
 
	if(currpnt == npnts) // Last run, will end
		// Change enable state of buttons
		Button StopButton,rename=StartButton,title="Start",fColor=(40960,65280,16384),win=SweepScopePanel_R
			
		// set currpnt to NaN
		currpnt = NaN
		
		// Save data file
		SVAR filePathStr = root:Packages:GQ:ykDL750:filePathStr
		DFREF saveDFR = GetDataFolderDFR()
		sprintf wName, "%s%04d", basename, fileIndex
		SetDataFolder $pathName
		Save/O/B/C/W/P=dataFolderPath WaveList("*",";","DIMS:2,TEXT:0")
		
		for(i = 0;i<n;i+=1)
			sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
			WAVE w = $(pathName+":"+wName)
		endfor

		// Increment file Index
		fileIndex += 1
			
//		KillWaves/Z tmpwave
//	
		SetDataFolder saveDFR

		return 1
	endif
	
	// Restart data acquisition
	ykWrite(kAddress,":START")
	
	// Increment counter
	currpnt += 1
	
	return 0
End



// Function to sweep both Yoko 7651

Function ControlSweepProc_Twice(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	DFREF saveDFR = GetDataFolderDFR()
	
	Variable err = 0, status
	
	SetDataFolder root:Packages:GQ:ykDL750

	SVAR address, basename, filePathStr
	NVAR fileIndex
	
	WAVE tw = tracew
	WAVE/T tS = traceSettings
	
	// Parameters for YOKO 7651 1

	WAVE loopwave = yoko7651loopwave

	NVAR initval = yoko7651initialvalue
	NVAR finalval = yoko7651finalvalue
	NVAR step = yoko7651step
	NVAR npnts = yoko7651numpnts
	NVAR currentsrc = yoko7651currentsrc
	NVAR currpnt = yoko7651currpnt
	NVAR yoko7651address

	// Parameters for YOKO 7651 2

	WAVE loopwave2 = yoko7651loopwave2

	NVAR initval = yoko7651initialvalue2
	NVAR finalval = yoko7651finalvalue2
	NVAR step = yoko7651step2
	NVAR currentsrc = yoko7651currentsrc2
	NVAR yoko7651address2

	String pathName, wName
	Variable i, n
	Variable numiv, trigpos, dt
					
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			
			strswitch(ba.ctrlName)
				case "StartButton":
					// Check if yoko is not in Normal mode
					if(!ykReady())
						DoAlert 0, "Must set Yoko DL750 to AVERAGE or SINGLE mode!"
						err = -1
						break
					endif
					
					// DO NOT Put Yoko7651 in proper range
					//Yoko7651SetRange("GPIB0::"+num2str(yoko7651address)+"::INSTR", wavemax(loopwave),current=currentsrc)
			
					// Set Initial Value for Yoko7651
					Yoko7651SetOutput(kGPIBn+"::"+num2str(yoko7651address)+"::INSTR", loopwave[currpnt])
					Yoko7651SetOutput(kGPIBn+"::"+num2str(yoko7651address2)+"::INSTR", loopwave2[currpnt])
					

					// Prepare waves to be written
					//
					//
					ykAvailableTraces(address,tw)
//					numiv = str2num(ykQuery(address,"WAV:LENG?"))
					numiv = str2num(ykQuery(address,"ACQ:RLEN?")) 
					trigpos = str2num(ykQuery(address,"WAVeform:TRIGGER?"))
					dt = 1/str2num(ykQuery(address,"TIMEBASE:SRATE?"))
					n = numpnts(tw)				

					// Path for new data
					sprintf pathName, "root:Data:%s%04d", basename, fileIndex
					fileIndex = CheckDataPath(pathName,basename,fileIndex)
					NewDataFolder/O $pathName
					NewPath/O/Q dataFolderPath, filePathStr
				
					// Set current sweep point to 0
					currpnt = 0
			
					// Create empy data				
					for(i = 0;i<n;i+=1)
						sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
						Make/O/N=(numiv,numpnts(loopwave)) $(pathName+":"+wName)
						WAVE w = $(pathName+":"+wName)
						w = NaN
						SetScale d,0,0,tS[tw[i]-1][1],w
						SetScale/P x,-trigpos*dt,dt,"s",w
						
						SetScale/P y,loopwave[0],loopwave[1]-loopwave[0],WaveUnits(loopwave,-1),w
							
						Note w, "Scaling: " + tS[tw[i]-1][0]
						Note w, "Scaling Axis: " + ts[0][2]
					endfor
			
					// Launch background process
					CtrlNamedBackground SweepScopeTwice, period=60, proc=SweepScopeTwiceTask
					CtrlNamedBackground SweepScopeTwice, start
			
					// Change enable state of buttons
					Button $ba.ctrlName,rename=StopButton,title="Stop",fColor=(65280,0,0),win=SweepScopePanel_Twice
					//Button StartSweep,disable=2,win=SweepScopePanel
					//Button StopSweep,disable=0,win=SweepScopePanel,fColor=(40960,65280,16384)
			
					break
				case "StopButton":
				
					// Stop background process
					CtrlNamedBackground SweepScopeTwice, stop
			
					// Redimension data and transpose if requested
					ControlInfo/W=SweepScopePanel_Twice TransposeCheck
										
					n = numpnts(tw)
					i = 0
					
					sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
					sprintf pathName, "root:Data:%s%04d", basename, fileIndex
					WAVE w = $(pathName+":"+wName)
					
					numiv = DimSize(w,0)
					
					for(i = 0;i<n;i+=1)
						sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
						WAVE w = $(pathName+":"+wName)
						Redimension/N=(numiv,currpnt-1) w
						
						if(V_Value)
							MatrixTranspose w
						endif
					endfor
			
					// set currpnt to NaN
					currpnt = NaN
	
					// Change buttons to default state
					Button $ba.ctrlName,rename=StartButton,title="Start",fColor=(40960,65280,16384),win=SweepScopePanel_Twice
					Button PauseButton,valueColor=(0,0,0),win=SweepScopePanel_Twice
					
					// Save data file
					sprintf wName, "%s%04d", basename, fileIndex
					sprintf pathName, "root:Data:%s%04d", basename, fileIndex

					SetDataFolder $pathName
					
//					// Kill temporary wave if it exists
//					KillWaves/Z tmpwave
//					
					Save/O/B/C/W/P=dataFolderPath WaveList("*",";","DIMS:2,TEXT:0")

					// Increment file Index
					fileIndex += 1
						
					break
				case "PauseButton":
					if(numtype(currpnt) == 0)
						CtrlNamedBackground SweepScopeTwice, status
						if(NumberByKey("RUN",S_Info))
							CtrlNamedBackground SweepScopeTwice, stop
							Button PauseButton,valueColor=(65280,0,0),win=SweepScopePanel_Twice
						else
							CtrlNamedBackground SweepScopeTwice, start
							Button PauseButton,valueColor=(0,0,0),win=SweepScopePanel_Twice
						endif
					endif
					break							
			endswitch
			
			break
		case -1: // control being killed
			break
	endswitch

	SetDataFolder saveDFR // Restore current DF.

	return err
End

Function SwapButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	switch( ba.eventCode )
		case 2: // mouse up
			DFREF saveDFR = GetDataFolderDFR()

			SetDataFolder root:Packages:GQ:ykDL750
		
			NVAR initval = yoko7651initialvalue
			NVAR finalval = yoko7651finalvalue
			NVAR step = yoko7651step
			NVAR npnts = yoko7651numpnts

			NVAR initval2 = yoko7651initialvalue2
			NVAR finalval2 = yoko7651finalvalue2
			NVAR step2 = yoko7651step2

			Variable tmpval
			
			tmpval = initval
			initval = finalval
			finalval = tmpval
			step *= -1

			tmpval = initval2
			initval2 = finalval2
			finalval2 = tmpval
			step2 *= -1

			// click code here
			break
		case -1: // control being killed
			break
	endswitch
	
	SetDataFolder saveDFR // Restore current DF.

	return 0
End


Function CheckDataPath(pathName,basename,fileIndex)
	String &pathName
	String basename
	Variable fileIndex

	sprintf pathName, "root:Data:%s%04d", basename, fileIndex

	if(DataFolderExists(pathName))
		DoAlert 1, "Overwrite existing data?"
		if(V_flag == 1)
			NewDataFolder/O/S $pathName
			KillWaves/A/Z
			KillVariables/A/Z
		else
			do
				fileIndex += 1
				sprintf pathName, "root:Data:%s%04d", basename, fileIndex
			while(DataFolderExists(pathName))
		endif
	endif
			
	return fileIndex
End

Function SweepScopeTask(s)
	STRUCT WMBackgroundStruct &s
	
	if(ykBusy())
		return 0
	endif
	
	NVAR currpnt = root:Packages:GQ:ykDL750:yoko7651currpnt
	NVAR initval = root:Packages:GQ:ykDL750:yoko7651initialvalue
	NVAR step = root:Packages:GQ:ykDL750:yoko7651step
	NVAR npnts = root:Packages:GQ:ykDL750:yoko7651numpnts
	NVAR currentsrc = root:Packages:GQ:ykDL750:yoko7651currentsrc
	NVAR yoko7651address = root:Packages:GQ:ykDL750:yoko7651address
	NVAR scalevar = 	root:Packages:GQ:ykDL750:scaleVar
	WAVE lw = root:Packages:GQ:ykDL750:yoko7651loopwave
	
	SVAR address = root:Packages:GQ:ykDL750:address
	WAVE tw = root:Packages:GQ:ykDL750:tracew
	WAVE/T tS = root:Packages:GQ:ykDL750:traceSettings
	
	ControlInfo/W=SweepScopePanel NoSweepCheck
	if(V_Value==0)
		Yoko7651SetOutput(kGPIBn+"::"+num2str(yoko7651address)+"::INSTR", lw[currpnt])
	endif

	SVAR baseName = root:Packages:GQ:ykDL750:basename
	NVAR fileIndex = root:Packages:GQ:ykDL750:fileIndex
	
	String pathName
	sprintf pathName, "root:Data:%s%04d", basename, fileIndex

	String wName
	Variable i, n = numpnts(tw)

	if(currpnt==0)  // First time run, nothing to do

	else  // transfer data
		for(i = 0;i<n;i+=1)
			sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
			ykGetTrace(address, tw[i],"tmpwave")
			WAVE tmpw = tmpwave
			WAVE w = $(pathName+":"+wName)
			Execute/Q "root:Packages:GQ:ykDL750:scaleVar = " + tS[tw[i]-1][0]  // scaling could be performed more efficiently
			FastOp tmpw = (1/scalevar)*tmpw
			w[][currpnt-1] = tmpw[p]
		endFor
		
		KillWaves/Z tmpwave
	endif
 
	if(currpnt == 1) // set proper scaling of data if requested
		ControlInfo/W=SweepScopePanel ScaleAxis
		
		if(V_Value > 1) // Axis scaling not requested if "0", first element, is chosen in ScaleAxis popupmenu
			Variable dx, offset
			Variable scaleAxisIndex = str2num(S_Value)
			
			sprintf wName, "%s%04d_%s", basename, fileIndex, tS[scaleAxisIndex-1][2]
			WAVE w = $(pathName+":"+wName)
			Make/FREE/N=(DimSize(w,0)) wcut = w[p][0]
			StatsLinearRegression wcut
			WAVE wlinreg = W_StatsLinearRegression
			offset = wlinreg[1]
			dx = wlinreg[2]
			KillWaves/Z W_StatsLinearRegression
						
			for(i = 0;i<n;i+=1)
				sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
				WAVE w = $(pathName+":"+wName)
				SetScale/P x,offset,dx,tS[scaleAxisIndex-1][1],w
			endFor
		endif
	endif
	
	if(currpnt == npnts) // Last run, will end
		
		// Change enable state of buttons
		Button StopButton,rename=StartButton,title="Start",fColor=(40960,65280,16384),win=SweepScopePanel
			
		// set currpnt to NaN
		currpnt = NaN
		
		// Save data file
		SVAR filePathStr = root:Packages:GQ:ykDL750:filePathStr
		DFREF saveDFR = GetDataFolderDFR()
		sprintf wName, "%s%04d", basename, fileIndex
		SetDataFolder $pathName
		Save/O/B/C/W/P=dataFolderPath WaveList("*",";","DIMS:2,TEXT:0")
		
		// Transpose if requested
		ControlInfo/W=SweepScopePanel TransposeCheck
		
		for(i = 0;i<n;i+=1)
			sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
			WAVE w = $(pathName+":"+wName)
				
			if(V_Value)
				MatrixTranspose w
			endif
		endfor

		// Increment file Index
		fileIndex += 1
			
//		KillWaves/Z tmpwave
//	
		SetDataFolder saveDFR

		return 1
	endif
	
	// Restart data acquisition
	ykWrite(kAddress,":START")
	
	// Increment counter
	currpnt += 1
	
	return 0
End

Function SweepScopeTwiceTask(s)
	STRUCT WMBackgroundStruct &s
	
	if(ykBusy())
		return 0
	endif
	
	NVAR currpnt = root:Packages:GQ:ykDL750:yoko7651currpnt
	NVAR npnts = root:Packages:GQ:ykDL750:yoko7651numpnts

	NVAR initval = root:Packages:GQ:ykDL750:yoko7651initialvalue
	NVAR step = root:Packages:GQ:ykDL750:yoko7651step
	NVAR currentsrc = root:Packages:GQ:ykDL750:yoko7651currentsrc
	NVAR yoko7651address = root:Packages:GQ:ykDL750:yoko7651address
	WAVE lw = root:Packages:GQ:ykDL750:yoko7651loopwave
	
	NVAR initval2 = root:Packages:GQ:ykDL750:yoko7651initialvalue2
	NVAR step2 = root:Packages:GQ:ykDL750:yoko7651step2
	NVAR currentsrc2 = root:Packages:GQ:ykDL750:yoko7651currentsrc2
	NVAR yoko7651address2 = root:Packages:GQ:ykDL750:yoko7651address2
	WAVE lw2 = root:Packages:GQ:ykDL750:yoko7651loopwave2
	
	
	NVAR scalevar = 	root:Packages:GQ:ykDL750:scaleVar
	SVAR address = root:Packages:GQ:ykDL750:address
	WAVE tw = root:Packages:GQ:ykDL750:tracew
	WAVE/T tS = root:Packages:GQ:ykDL750:traceSettings
	
	Yoko7651SetOutput(kGPIBn+"::"+num2str(yoko7651address)+"::INSTR", lw[currpnt])
	Yoko7651SetOutput(kGPIBn+"::"+num2str(yoko7651address2)+"::INSTR", lw2[currpnt])

	SVAR baseName = root:Packages:GQ:ykDL750:basename
	NVAR fileIndex = root:Packages:GQ:ykDL750:fileIndex
	
	String pathName
	sprintf pathName, "root:Data:%s%04d", basename, fileIndex

	String wName
	Variable i, n = numpnts(tw)

	if(currpnt==0)  // First time run, nothing to do

	else  // transfer data
		for(i = 0;i<n;i+=1)
			sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
			ykGetTrace(address, tw[i],"tmpwave")
			WAVE tmpw = tmpwave
			WAVE w = $(pathName+":"+wName)
			Execute/Q "root:Packages:GQ:ykDL750:scaleVar = " + tS[tw[i]-1][0]  // scaling could be performed more efficiently
			FastOp tmpw = (1/scalevar)*tmpw
			w[][currpnt-1] = tmpw[p]
		endFor
		
		KillWaves/Z tmpwave
	endif

		if(currpnt == npnts) // Last run, will end
		
		// Change enable state of buttons
		Button StopButton,rename=StartButton,title="Start",fColor=(40960,65280,16384),win=SweepScopePanel_Twice
			
		// set currpnt to NaN
		currpnt = NaN
		
		// Save data file
		SVAR filePathStr = root:Packages:GQ:ykDL750:filePathStr
		DFREF saveDFR = GetDataFolderDFR()
		sprintf wName, "%s%04d", basename, fileIndex
		SetDataFolder $pathName
		Save/O/B/C/W/P=dataFolderPath WaveList("*",";","DIMS:2,TEXT:0")
		
		// Transpose if requested
		ControlInfo/W=SweepScopePanel_Twice TransposeCheck
		
		for(i = 0;i<n;i+=1)
			sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
			WAVE w = $(pathName+":"+wName)
				
			if(V_Value)
				MatrixTranspose w
			endif
		endfor

		// Increment file Index
		fileIndex += 1
			
//		KillWaves/Z tmpwave
//	
		SetDataFolder saveDFR

		return 1
	endif
	
	// Restart data acquisition
	ykWrite(kAddress,":START")
	
	// Increment counter
	currpnt += 1
	
	return 0
End



Function NoSweepCheckProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			
			checked *= 2
			
			SetVariable yoko7651initialvalue win=SweepScopePanel,disable=checked
			SetVariable yoko7651finalvalue win=SweepScopePanel,disable=checked
			SetVariable yoko7651step win=SweepScopePanel,disable=checked
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function UpdateCurrentSrcProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	DFREF saveDFR = GetDataFolderDFR()
	
	Variable status = 0
	
	SetDataFolder root:Packages:GQ:ykDL750
	
	WAVE loopwave = yoko7651loopwave

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			if(checked)
				SetScale d,0,0,"A",loopwave
			else
				SetScale d,0,0,"V",loopwave
			endif

			break
		case -1: // control being killed
			break
	endswitch

	SetDataFolder saveDFR

	return 0
End

Proc ykGraphStyle() : GraphStyle
	PauseUpdate; Silent 1
	
	ModifyGraph/Z rgb[0]=(0,26112,39168),rgb[1]=(0,0,0),rgb[2]=(0,0,63222),rgb[3]=(39168,0,15616)
	ModifyGraph/Z rgb[4]=(65280,16384,55552),rgb[5]=(0,37008,0),rgb[6]=(65280,0,0),rgb[7]=(59367,49344,0)
	ModifyGraph/Z grid=2
	ModifyGraph/Z gridRGB=(56576,56576,56576)
	ModifyGraph/Z mode=2
	//	ModifyGraph/Z margin(right)=0
	//	ModifyGraph/Z lblMargin=5
	Label/Z left "\\Z12\\f01\\U"
	Label/Z bottom "\\Z12\\f01\\U"
	//Legend/W=ykGraph/C/V=1/N=ykGraphLegend/A=RT/F=0/B=1/E=0
EndMacro

#endif