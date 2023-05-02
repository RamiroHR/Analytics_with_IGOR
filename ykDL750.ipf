#pragma rtGlobals=1		// Use modern global access method.
#pragma version=1.23		
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
// CHANGES:
// 1.23 finished modifying flux sweep
// 1.22 started modifying flux sweep to include cross-inductance
// 1.21 added SweepScopeFlux to sweep flux in a better way
// 1.20 added "Initialize_ykDL750(0)" to ykPanel in order to force config directory creation

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
//StrConstant kAddress = "GPIB0::6::INSTR"
//Constant kNumDSP = 0

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
	
//	viStatusDesc(session, status, desc)
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
		String/G filePathStr = "U:shared", basename = "TST"
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

		// Prepare YOKO 7651 1 (Coil)

		Variable/G yoko7651address = 16
		Variable/G yoko7651currentsrc = 1
		Variable/G yoko7651initialvalue = 0
		Variable/G yoko7651finalvalue = 3e-4
		Variable/G yoko7651step = 1e-6

		Make/O/N=1 yoko7651currvalue = NaN
	
		Variable/G yoko7651numpnts = 301
		Variable/G yoko7651currpnt = NaN
	
		Make/O/N=(yoko7651numpnts) yoko7651loopwave = yoko7651initialvalue + p*yoko7651step
		SetScale d,0,0,"A",yoko7651loopwave
		

		// Prepare YOKO 7651 2 (Gradio)

		Variable/G yoko7651address2 = 5
		Variable/G yoko7651currentsrc2 = 1
		Variable/G yoko7651initialvalue2 = 0
		Variable/G yoko7651finalvalue2 = 3e-4
		Variable/G yoko7651step2 = 1e-6

		Make/O/N=1 yoko7651currvalue2 = NaN

		Make/O/N=(yoko7651numpnts) yoko7651loopwave2 = yoko7651initialvalue2 + p*yoko7651step2
		SetScale d,0,0,"A",yoko7651loopwave2
	
	
		// For Flux Sweeps
		Variable/G yoko7651phi01 = 0
		Variable/G yoko7651phi11 = 1

		Variable/G yoko7651phi02 = 0
		Variable/G yoko7651phi12 = 1

		Variable/G yoko7651slowramptime = 10
		Variable/G yoko7651fastramptime = 1		
		
		// For flux compensation:
		//
		
		// beta_L
		Make/O/N=(2) fluxbetaL = 0	

		// mutual coupling
		Make/O/N=(2,2) fluxgamma = p == q
				
		// inductance matrix
		Make/O/N=(2,2) Lmatrix	
		
		// zero phase current
		Make/O/N=2 currentzero = 0
		
		
		
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

Function Yoko7651SetRange(resourceName, output,[current])  // OLD needs updating 1.7.2020
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

Function Yoko7651SetOutputOld(resourceName, output, [autorange])
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

Function Yoko7651Busy(resourceName)
	String resourceName
	
	String outstr = QueryVISA(resourceName,"OC")
	
	// see pg 6-33 of manual
	// BIT : DESCRIPTION : 0 : 1
	// 5 : output on? : off : on
	// 4 : output unstable? : normal : unstable
	// 3 : command error? : OK : error
	// 2 : program executing? : no : yes
	// 1 : programming? : no : yes
	
	if (str2num(outstr[5,strlen(outstr)]) & 15) // check if executing program, in process of programming, output unstable, or cmd error
		return 1
	else
		return 0
	endif
End

Function Yoko7651RampOutput(resourceName, output, [duration, autorange])
	String resourceName
	Variable output
	Variable duration, autorange
	
	Variable status
	
	if(Yoko7651Busy(resourceName))
		Abort "Yoko 7651 " + resourceName + " Busy!"
	endif
	
	String outputstr, progstr = ""
	
	if(ParamIsDefault(autorange))
		autorange = 0
	endif
	
	if(ParamIsDefault(duration))
		duration = 1.0
	endif
	
	if(autorange)
		outputstr = "SA"
	else
		outputstr = "S"
	endif
	
	progstr += "PRS" // start programming
	//progstr += outputstr + num2str(initval) // program initial value
	//progstr += outputstr + num2str(finalval) // program final value
	progstr += outputstr + num2str(output) // program output value
	progstr += "PRE" // end programming
	progstr += "SW" + num2str(duration) // sweep time
	progstr += "PI" + num2str(duration) // sweep time
	progstr += "M1" // run once
	progstr += "RU2" // run once
	
	status	= WriteVISA(resourceName,progstr)
	
	return status
End

Function Yoko7651SetOutput(resourceName, output, [autorange])
	String resourceName
	Variable output, autorange

	if(Yoko7651Busy(resourceName))
		Abort "Yoko 7651 " + resourceName + " Busy!"
	endif

	Variable status
	
	if(ParamIsDefault(autorange))
		autorange = 0
	endif
	
	String progstr

	if(autorange)
		progstr = "SA" + num2str(output) + "EO1E"
	else
		progstr = "S" + num2str(output) + "EO1E"
	endif
		
	status = WriteVISA(resourceName,progstr)
	
	return status
End

Function Yoko7651ReadOutput(resourceName)
	String resourceName

	String outstr = QueryVISA(resourceName,"H0OD")

	return str2num(outstr)
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


// v1.23
//
// Procedures to sweep flux when there are multiple loops and flux sources
//
//
//

Window FluxCalibPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1026,344,1326,544) as "Flux Calibration Panel"
EndMacro


Function SweepFluxButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	Variable err = 0, status = 0
	
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			DFREF saveDFR = GetDataFolderDFR()

			SetDataFolder root:Packages:GQ:ykDL750


			SVAR address, basename, filePathStr
			NVAR fileIndex
	
			WAVE tw = tracew
			WAVE/T tS = traceSettings
	

			NVAR npnts = yoko7651numpnts
			NVAR currpnt = yoko7651currpnt
			
			// Desired sweep values, in units of 2pi flux
			NVAR phi01 = yoko7651phi01
			NVAR phi02 = yoko7651phi02
			
			NVAR phi11 = yoko7651phi11
			NVAR phi12 = yoko7651phi12

			// Parameters for YOKO 7651 1
			WAVE loopwave = yoko7651loopwave
			WAVE yoko1val = yoko7651currvalue
			NVAR yoko7651address
			
			String yoko1address = kGPIBn+"::"+num2str(yoko7651address)+"::INSTR"

			// Parameters for YOKO 7651 1
			WAVE loopwave2 = yoko7651loopwave2
			WAVE yoko2val = yoko7651currvalue2
			NVAR yoko7651address2
			
			String yoko2address = kGPIBn+"::"+num2str(yoko7651address2)+"::INSTR"
			

			NVAR slowramptime = yoko7651slowramptime
			NVAR fastramptime = yoko7651fastramptime		


			// Flux calibration parameters
			WAVE czero = currentzero
			WAVE betaL = fluxbetaL
			WAVE mutual = fluxgamma
			WAVE L = Lmatrix
			
			WAVE phi = phi
			
			Variable tmp
			
			strswitch(ba.ctrlName)
				case "SwapPhaseButton":
					tmp = phi01
					phi01 = phi02
					phi02 = tmp
			
					tmp = phi11
					phi11 = phi12
					phi12 = tmp

				case "UpdateButton":
					Redimension/N=(npnts) loopwave, loopwave2
					Make/O/N=(2,npnts) phi, current
					phi[0][] = phi01+(phi02-phi01)*q/(npnts-1)
					phi[1][] = phi11+(phi12-phi11)*q/(npnts-1)

					CalcCurrent(phi,current,L,mutual,betaL)
			
					loopwave = current[0][p]+czero[0]
					loopwave2 = current[1][p]+czero[1]

					yoko1val[0] = Yoko7651ReadOutput(yoko1address)
					yoko2val[0] = Yoko7651ReadOutput(yoko2address)

					// Add arrow
					GetWindow kwTopWin wtitle
					if(strlen(S_Value) == 0)
						GetWindow kwTopWin activeSW
					else
						GetWindow kwTopWin activeSW
						S_Value += "#G0"	
					endif
					
					SetDrawLayer/W=$(S_Value) UserFront
					DrawAction/W=$(S_Value)/L=UserFront delete, begininsert
					SetDrawEnv/W=$(S_Value) gstart, xcoord= bottom,ycoord= left,linefgc= (52428,52428,52428,49151),dash= 3,arrow= 1,arrowlen= 20
					DrawLine/W=$(S_Value) loopwave[0],loopwave2[0],loopwave[numpnts(loopwave)-1],loopwave2[numpnts(loopwave)-1]
					//DrawLine/W=#G0 yoko1val,yoko2val,loopwave[0],loopwave2[0]
					SetDrawEnv/W=$(S_Value) gstop
					DrawAction/W=$(S_Value)/L=UserFront endinsert
					break
				case "RampToStartButton":
					Yoko7651RampOutput(yoko1address,loopwave[0],duration=slowramptime)
					Yoko7651RampOutput(yoko2address,loopwave2[0],duration=slowramptime)
					
					yoko1val[0] = Yoko7651ReadOutput(yoko1address)
					yoko2val[0] = Yoko7651ReadOutput(yoko2address)
					break
					
				case "StartButton":
					String pathName, wName
					Variable i, n
					Variable numiv, trigpos, dt

				
					//Check if DL750 is not in Normal mode
					if(!ykReady())
						DoAlert 0, "Must set Yoko DL750 to AVERAGE or SINGLE mode!"
						Print "Abort DL750"
						err = -1
						break
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
					for(i = 0;i<n;i+=1)
						sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
						Make/O/N=(numiv,numpnts(loopwave)) $(pathName+":"+wName)
						WAVE w = $(pathName+":"+wName)
						w = NaN
						SetScale d,0,0,tS[tw[i]-1][1],w
						SetScale/P x,-trigpos*dt,dt,"s",w

						Note w, "Scaling: " + tS[tw[i]-1][0]
					endfor
					
					// Copy phi and current sweep source waves
					Duplicate/O phi $(pathName+":phi")
					Duplicate/O loopwave $(pathName+":loopwave")
					Duplicate/O loopwave2 $(pathName+":loopwave2")

					// Check Yoko 7651 outputs
					// Need fancy check since Yoko output has finite digits
					
					tmp = Yoko7651ReadOutput(yoko1address)
					if(!AlmostSame(tmp,loopwave[currpnt]))
						DoAlert 1, "Yoko 7651 1 (Coil) output is " + num2str(tmp) + " and requested output is " + num2str(loopwave[currpnt]) + ". Continue anyway?"
						if(V_flag == 2)
							Print "Abort Yoko 7651 1"
							err = -1
							break
						endif
					endif
					
					tmp = Yoko7651ReadOutput(yoko2address)
					if(!AlmostSame(tmp,loopwave2[currpnt]))
						DoAlert 1, "Yoko 7651 1 (Gradio) output is " + num2str(tmp) + " and requested output is " + num2str(loopwave2[currpnt]) + ". Continue anyway?"
						if(V_flag == 2)
							Print "Abort Yoko 7651 2"
							err = -1
							break
						endif
					endif
					
					// Launch background process
					CtrlNamedBackground SweepScopeFlux, period=60, proc=SweepFluxTask
					CtrlNamedBackground SweepScopeFlux, start
			
					// Change enable state of buttons
					Button $ba.ctrlName,rename=StopButton,title="Stop",fColor=(65280,0,0),win=#
			
					break
				
				case "StopButton":
				
					// Stop background process
					CtrlNamedBackground SweepScopeFlux, stop
			
					// Redimension data and transpose if requested
					ControlInfo/W=# TransposeCheck
										
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
					Button $ba.ctrlName,rename=StartButton,title="Start",fColor=(40960,65280,16384),win=#
					Button PauseButton,valueColor=(0,0,0),win=#
					
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
						CtrlNamedBackground SweepScopeFlux, status
						if(NumberByKey("RUN",S_Info))
							CtrlNamedBackground SweepScopeFlux, stop
							Button PauseButton,valueColor=(65280,0,0),win=#
						else
							CtrlNamedBackground SweepScopeFlux, start
							Button PauseButton,valueColor=(0,0,0),win=#
						endif
					endif
					break	
			endswitch
					
			SetDataFolder saveDFR

			break
		case -1: // control being killed
			break
	endswitch

	
	return 0
End

Function FluxScreening(w, xW, yW)
	// function to determine value of phase corresponding to given currents
	// see help for FindRoot
	Wave w, xW, yW
	
	DFREF saveDFR = GetDataFolderDFR()

	SetDataFolder root:Packages:GQ:ykDL750
	
	// Flux calibration parameters
	WAVE czero = currentzero
	WAVE betaL = fluxbetaL
	WAVE mutual = fluxgamma
	WAVE L = Lmatrix

	SetDataFolder saveDFR
	
	Make/FREE/N=2 curroffset = w[p][0]
		
	MatrixOp/O yW = inv(L) x (xW + 1/2/pi * mutual x diagonal(betaL) x sin(2*pi*xW)) + czero - curroffset
End


Function CurrentPhase(phi, currents)
	Wave phi, currents
	
	FindRoots /X=phi FluxScreening, currents
End

Function TestFindCurrentPhase(phi)
	Wave phi
	
	Make/FREE/N=2 actualCurrents, current
	
	DFREF saveDFR = GetDataFolderDFR()

	SetDataFolder root:Packages:GQ:ykDL750
	
	// Desired sweep values, in units of 2pi flux
	NVAR phi01 = yoko7651phi01
	NVAR phi02 = yoko7651phi02
			
	NVAR phi11 = yoko7651phi11
	NVAR phi12 = yoko7651phi12

	// Parameters for YOKO 7651 1
	WAVE loopwave = yoko7651loopwave
	NVAR yoko7651address
			
	String yoko1address = kGPIBn+"::"+num2str(yoko7651address)+"::INSTR"

	// Parameters for YOKO 7651 1
	WAVE loopwave2 = yoko7651loopwave2
	NVAR yoko7651address2
			
	String yoko2address = kGPIBn+"::"+num2str(yoko7651address2)+"::INSTR"
			
	// Flux calibration parameters
	WAVE czero = currentzero
	WAVE betaL = fluxbetaL
	WAVE mutual = fluxgamma
	WAVE L = Lmatrix

	SetDataFolder saveDFR
	
	actualCurrents[0] = Yoko7651ReadOutput(yoko1address)
	actualCurrents[1] = Yoko7651ReadOutput(yoko2address)
	Print actualCurrents
		
	MatrixOp/FREE guessphi = L x (actualCurrents - czero)  // Guess solution for betaL = 0
	Print "guessphi:", guessphi
	Print "phi0:", phi01, phi02
	Print "phi1:", phi11, phi12
	Redimension/N=(2,2) actualCurrents
	FindRoots /X=guessphi FluxScreening, actualCurrents
	
End

Function AlmostSame(v1, v2, [digits,percentage])
	Variable v1, v2, digits, percentage
	
	if(ParamIsDefault(digits) && ParamIsDefault(percentage))
		digits = 2
	endif
	
	String s1 
	String s2
	
	
	sprintf s1, "%.*e", digits, v1
	sprintf s2, "%.*e", digits, v2
		
	return StringMatch(s1,s2)		
End


Function CalcCurrent(phi,current,L,mutual,betaL)
	//
	// Calculate coil currents necessary to produce desired superconducting phases
	//
	// For our purposes, m = 2 (superconducting phase differences of two loops)
	//
	Wave phi			// m x n phase matrix (m independent loops, n values of phase)
	Wave current		// m x n current matrix (m loops, n values of current)
	Wave L 			    // m x m symmetric inductance matrix such that phi = L x current
	Wave mutual         // m x m mutual coupling matrix (should have 1 along diagonal)
	Wave betaL          // length m vector with betaL for each loop

	Variable m = DimSize(phi,0)
	Variable n = DimSize(phi,1)

	//Make/FREE/N=(n) phiout
	// Our definition is Phi = L I
	MatrixOp/O current = inv(L) x (phi + 1/2/pi * mutual x diagonal(betaL) x sin(2*pi*phi))
End

Function SweepFluxTask(s)
	STRUCT WMBackgroundStruct &s
	
	if(ykBusy())
		return 0
	endif
	
	NVAR currpnt = root:Packages:GQ:ykDL750:yoko7651currpnt
	NVAR npnts = root:Packages:GQ:ykDL750:yoko7651numpnts

	NVAR yoko7651address = root:Packages:GQ:ykDL750:yoko7651address
	WAVE lw = root:Packages:GQ:ykDL750:yoko7651loopwave
	
	NVAR yoko7651address2 = root:Packages:GQ:ykDL750:yoko7651address2
	WAVE lw2 = root:Packages:GQ:ykDL750:yoko7651loopwave2
	
	NVAR scalevar = 	root:Packages:GQ:ykDL750:scaleVar
	SVAR address = root:Packages:GQ:ykDL750:address
	WAVE tw = root:Packages:GQ:ykDL750:tracew
	WAVE/T tS = root:Packages:GQ:ykDL750:traceSettings
	
	SVAR baseName = root:Packages:GQ:ykDL750:basename
	NVAR fileIndex = root:Packages:GQ:ykDL750:fileIndex


	// Parameters for YOKO 7651 1
	WAVE yoko1val = root:Packages:GQ:ykDL750:yoko7651currvalue
	String yoko1address = kGPIBn+"::"+num2str(yoko7651address)+"::INSTR"

	// Parameters for YOKO 7651 2
	WAVE yoko2val = root:Packages:GQ:ykDL750:yoko7651currvalue2
	String yoko2address = kGPIBn+"::"+num2str(yoko7651address2)+"::INSTR"

	NVAR fastramptime = root:Packages:GQ:ykDL750:yoko7651fastramptime

	if(Yoko7651Busy(yoko1address) || Yoko7651Busy(yoko2address))
		yoko1val[0] = Yoko7651ReadOutput(yoko1address)
		yoko2val[0] = Yoko7651ReadOutput(yoko2address)
		return 0
	endif
	
	String pathName
	sprintf pathName, "root:Data:%s%04d", basename, fileIndex

	String wName
	Variable i, n = numpnts(tw)

	yoko1val[0] = Yoko7651ReadOutput(yoko1address)
	yoko2val[0] = Yoko7651ReadOutput(yoko2address)

	if(currpnt == npnts) // Last run, will end
		for(i = 0;i<n;i+=1) // transfer data
			sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
			ykGetTrace(address, tw[i],"tmpwave")
			WAVE tmpw = tmpwave
			WAVE w = $(pathName+":"+wName)
			Execute/Q "root:Packages:GQ:ykDL750:scaleVar = " + tS[tw[i]-1][0]  // scaling could be performed more efficiently
			FastOp tmpw = (1/scalevar)*tmpw
			w[][currpnt-1] = tmpw[p]
		endFor
		
		KillWaves/Z tmpwave

		// Change enable state of buttons
		Button StopButton,rename=StartButton,title="Start",fColor=(40960,65280,16384),win=SweepScopeFluxPanel
			
		// set currpnt to NaN
		currpnt = NaN
		
		// Save data file
		SVAR filePathStr = root:Packages:GQ:ykDL750:filePathStr
		DFREF saveDFR = GetDataFolderDFR()
		sprintf wName, "%s%04d", basename, fileIndex
		SetDataFolder $pathName
		Save/O/B/C/W/P=dataFolderPath WaveList("*",";","DIMS:2,TEXT:0")
		
		// Transpose if requested
		ControlInfo/W=SweepScopeFluxPanel TransposeCheck
		
		for(i = 0;i<n;i+=1)
			sprintf wName, "%s%04d_%s", basename, fileIndex, tS[tw[i]-1][2]
			WAVE w = $(pathName+":"+wName)
				
			if(V_Value)
				MatrixTranspose w
			endif
		endfor

		// Increment file Index
		fileIndex += 1
			
		SetDataFolder saveDFR

		return 1
	else	
		yoko1val[0] = Yoko7651ReadOutput(yoko1address)
		yoko2val[0] = Yoko7651ReadOutput(yoko2address)

		if(currpnt>0)  // If first time run, do nothing, else transfer data
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

		// Change outputs
		Yoko7651RampOutput(yoko1address, lw[currpnt], duration=fastramptime)
		Yoko7651RampOutput(yoko2address, lw2[currpnt], duration=fastramptime)
	endif

	// Restart data acquisition
	ykWrite(kAddress,":START")
	
	// Increment counter
	currpnt += 1
	
	return 0
End



Window SweepScopeFluxPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1154,53,2150,632) as "Sweep Scope Flux"
	SetDrawLayer UserBack
	SetDrawEnv fstyle= 1
	DrawText 162,23,"Data Information"
	SetDrawEnv fsize= 20
	DrawText 13,264,"φ\\B1"
	SetDrawEnv fsize= 20
	DrawText 10,352,"φ\\B2\\M"
	DrawLine 3,472,177,472
	DrawLine 449,15,449,560
	DrawLine 5,224,434,224
	DrawText 270,361,"Mutual Coupling Matrix"
	DrawText 366,260,"β\\BL\\M vector"
	DrawLine 187,225,187,472
	DrawText 280,475,"Inductance Matrix"
	DrawLine 464,494,978,494
	DrawText 229,258,"Ref φ Current"
	DrawLine 0.00016818183567887,0.00213454727547406,3.25545174585367e-005,0.00999674492602187
	Button StartButton,pos={741.00,530.00},size={102.00,41.00},proc=SweepFluxButtonProc,title="Start"
	Button StartButton,fStyle=1,fColor=(40960,65280,16384)
	SetVariable Basename,pos={743.00,503.00},size={78.00,18.00},title=" "
	SetVariable Basename,help={"Base name of data series"}
	SetVariable Basename,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:basename
	SetVariable FileIndex,pos={848.00,504.00},size={71.00,18.00},title=" "
	SetVariable FileIndex,format="%04d"
	SetVariable FileIndex,limits={0,9999,1},value= root:Packages:GQ:ykDL750:fileIndex
	SetVariable FilePath,pos={530.00,504.00},size={200.00,18.00},title=" "
	SetVariable FilePath,help={"File path to save data"}
	SetVariable FilePath,value= root:Packages:GQ:ykDL750:filePathStr
	SetVariable yoko7651address,pos={25.00,288.00},size={135.00,18.00},title=" Yoko Coil GPIB"
	SetVariable yoko7651address,format="%2d"
	SetVariable yoko7651address,limits={0,21,1},value= root:Packages:GQ:ykDL750:yoko7651address
	Button PauseButton,pos={865.00,531.00},size={102.00,41.00},proc=SweepFluxButtonProc,title="Pause"
	Button PauseButton,fStyle=1
	SetVariable yoko7651currpnt,pos={490.00,553.00},size={143.00,18.00},title="Current Point"
	SetVariable yoko7651currpnt,format="%g",fStyle=1,valueColor=(65535,32768,32768)
	SetVariable yoko7651currpnt,limits={1,inf,0},value= root:Packages:GQ:ykDL750:yoko7651currpnt,noedit= 1
	CheckBox TransposeCheck,pos={646.00,533.00},size={74.00,15.00},title="Transpose?"
	CheckBox TransposeCheck,help={"Fast axis is acquired along x-direction. Check to transpose final matrix."}
	CheckBox TransposeCheck,value= 1,side= 1
	SetVariable yoko7651address2,pos={10.00,369.00},size={147.00,18.00},title=" Yoko Gradio GPIB"
	SetVariable yoko7651address2,format="%2d"
	SetVariable yoko7651address2,limits={0,21,1},value= root:Packages:GQ:ykDL750:yoko7651address2
	SetVariable phi01,pos={42.00,235.00},size={124.00,18.00},title="Initial Phase"
	SetVariable phi01,format="%g"
	SetVariable phi01,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651phi01
	SetVariable phi11,pos={46.00,254.00},size={120.00,18.00},title="Final Phase"
	SetVariable phi11,format="%g"
	SetVariable phi11,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651phi02
	SetVariable phi02,pos={42.00,322.00},size={124.00,18.00},title="Initial Phase"
	SetVariable phi02,format="%g"
	SetVariable phi02,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651phi11
	SetVariable phi12,pos={46.00,341.00},size={120.00,18.00},title="Final Phase"
	SetVariable phi12,format="%g"
	SetVariable phi12,limits={-inf,inf,0},value= root:Packages:GQ:ykDL750:yoko7651phi12
	SetVariable yoko7651numpnts,pos={522.00,531.00},size={109.00,18.00},title="# Points"
	SetVariable yoko7651numpnts,format="%d"
	SetVariable yoko7651numpnts,limits={1,inf,0},value= root:Packages:GQ:ykDL750:yoko7651numpnts
	Button UpdateButton,pos={8.00,494.00},size={203.00,61.00},proc=SweepFluxButtonProc,title="Update"
	Button UpdateButton,fStyle=0,fColor=(1,3,39321),valueColor=(61166,61166,61166)
	Button SwapPhaseButton,pos={20.00,407.00},size={124.00,32.00},proc=SweepFluxButtonProc,title="Swap Phase"
	Button RampToStartButton,pos={478.00,454.00},size={114.00,32.00},proc=SweepFluxButtonProc,title="Ramp to start"
	SetVariable SlowRampTime,pos={604.00,462.00},size={122.00,18.00},title="Ramp time (s) "
	SetVariable SlowRampTime,limits={1,inf,0},value= root:Packages:GQ:ykDL750:yoko7651slowramptime
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Edit/W=(21,31,416,214)/HOST=#  traceSettings.ld
	ModifyTable frameStyle= 5,format(Point)=1,width(traceSettings.l)=38,alignment(traceSettings.d)=1
	ModifyTable width(traceSettings.d)=82,width[2]=250,width[3,4]=50
	ModifyTable showParts=0x45
	ModifyTable statsArea=85
	SetDataFolder fldrSav0
	RenameWindow #,T0
	SetActiveSubwindow ##
	String fldrSav1= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Display/W=(470,11,980,448)/HOST=#  yoko7651loopwave2 vs yoko7651loopwave
	AppendToGraph yoko7651currvalue2 vs yoko7651currvalue
	SetDataFolder fldrSav1
	ModifyGraph mode=3
	ModifyGraph marker=19
	ModifyGraph rgb(yoko7651loopwave2)=(39168,0,15616)
	ModifyGraph msize(yoko7651loopwave2)=1
	ModifyGraph zero=2
	Label left "Current Yoko 2 (\\U)"
	Label bottom "Current Yoko 1 (\\U)"
	SetAxis/A/N=1 left
	SetAxis/A/N=1 bottom
	SetDrawLayer UserFront
	SetDrawEnv xcoord= bottom,ycoord= left,linefgc= (52428,52428,52428,49151),dash= 3,arrow= 1,arrowlen= 20
	SetDrawEnv gstart
	DrawLine 0.000259752909187227,0.0022149751894176,9.50247194850817e-005,0.0121985031291842
	SetDrawEnv gstop
	RenameWindow #,G0
	SetActiveSubwindow ##
	String fldrSav2= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Edit/W=(250,364,413,453)/HOST=#  fluxgamma
	ModifyTable format(Point)=1,width(Point)=17,width(fluxgamma)=70
	ModifyTable showParts=0x4d
	ModifyTable statsArea=19
	SetDataFolder fldrSav2
	RenameWindow #,T1
	SetActiveSubwindow ##
	String fldrSav3= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Edit/W=(349,265,442,334)/HOST=#  fluxbetaL
	ModifyTable format(Point)=1,width(Point)=17,width(fluxbetaL)=70
	ModifyTable showParts=0x4d
	ModifyTable statsArea=17
	SetDataFolder fldrSav3
	RenameWindow #,T2
	SetActiveSubwindow ##
	String fldrSav4= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Edit/W=(249,481,412,571)/HOST=#  Lmatrix
	ModifyTable format(Point)=1,width(Point)=17,width(Lmatrix)=70
	ModifyTable showParts=0x4d
	ModifyTable statsArea=19
	SetDataFolder fldrSav4
	RenameWindow #,T3
	SetActiveSubwindow ##
	String fldrSav5= GetDataFolder(1)
	SetDataFolder root:Packages:GQ:ykDL750:
	Edit/W=(207,263,331,332)/HOST=#  currentzero
	ModifyTable format(Point)=1,width(Point)=17,width(currentzero)=100
	ModifyTable showParts=0x4d
	ModifyTable statsArea=18
	SetDataFolder fldrSav5
	RenameWindow #,T4
	SetActiveSubwindow ##
EndMacro



#endif