#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.




Menu "JJ Spectro"

	"Initialize", InitializeConstants()
	"Get good gain mismatch", GetGoodGainMismatch()
	"Correct tilt", CorrectTilt()
	"-"
	"Fit RN at room temperature", FitRN()
	"Correct for a series resistance", SeriesR()
	"-"
	"Extract switching current/7", SuperIsw()
	"Auto XYZ3Matrix/8", SuperXYZ()
	"-"
	"Get R, L and C for an in-loop mode", GetRLC()
	"-"
	"Get a nice IV plot", NiceIV()
	"Make a nice IV from ykdl750 screen",MakeNicePlotFromYkdl()
	"Add legend for 0 and π flux",Legend0Pi()
	"Get a nice map", NiceMap()
	"Add a frequency axis on the axis opposite to the voltage", AddFreqAxis()
	"Add a top axis with flux (for maps)", TopAxisFlux()
	"Add a frequency axis on the axis opposite to the current",AddFreqAxisCurrent()
	"Add a Sensibility axis on the axis opposite to the current",AddSensAxisCurrent()
End


Function InitializeConstants()
//To do when loading this procedure
	Variable/G root:Data:GainMismatch = 0
	Variable/G root:Data:io = 0
	Variable/G root:Data:al = 1
	
	String savDF= GetDataFolder(1) // Save current DF for restore.
	
	//Load more colors
	NewDataFolder/O/S root:Colors
	LoadWave/Q/P=Igor "S:analysis:igor:ColorsForJJ.ibw"
	Wave W_colors
	
	Make/O/N=(161*3) GreenBrown
	Make/O/N=(101*3) Inferno
	Make/O/N=(101*3) Magma
	Make/O/N=(101*3) Plasm
	Make/O/N=(101*3) Viridis
	Make/O/N=(101*3) Seismic
	Make/O/N=(101*3) RdBu
	Make/O/N=(101*3)	 PuOr
	
	GreenBrown = W_colors[p]
	Inferno = W_colors[p+3*161]
	Magma = W_colors[p+3*161+3*101]
	Plasm = W_colors[p+3*161+2*3*101]
	Viridis = W_colors[p+3*161+3*3*101]	
	Seismic = W_colors[p+3*161+4*3*101]	
	RdBu = W_colors[p+3*161+5*3*101]
	PuOr = W_colors[p+3*161+6*3*101]
	
	Redimension/N=(161,3) GreenBrown
	Redimension/N=(101,3) Inferno,Magma,Plasm,Viridis,Seismic,RdBu,PuOr
	
	KillWaves W_colors
	
	SetDataFolder savDF // Restore current DF.
End

Function FitRN()
	String Traces = TraceNameList("",";",1+4)
	String yname = StringFromList(0,Traces)
	Wave/Z wy = TraceNameToWaveRef("",yname)	
	Wave/Z wx = XWaveRefFromTrace("",yname)
	
	Make/FREE/D/O/N=2 fitc	
	CurveFit/Q line kwCWave=fitc wy /X=wx
	Variable/G Rn = 1/fitc[1]

	Legend/C/N=text0/J/F=0/A=LT "\\s(" + yname + ")\tBefore cooldown\r\tR\\BN\\M = " + num2str(Rn) + " Ω"

	Modifygraph grid=1
	ModifyGraph mirror=2	
	Label bottom "Voltage (\\U)"
	Label left "Current (\\U)"
	
End


Function GetGoodGainMismatch()
//Calculates the gain mismatch between V and X amps
	Variable/G root:Data:GainMismatch
	
	String Traces = TraceNameList("",";",1+4)
	String yname = StringFromList(0,Traces)
	Wave/Z wy = TraceNameToWaveRef("",yname)	
	Wave/Z wx = XWaveRefFromTrace("",yname)
	
	String corname = yname + "_cor"
	Duplicate/O wy Ijcor
	
	ShowInfo
	
	NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursors"
	DoWindow/C tmp_PauseforCursors
	DrawText 21,20,"Adjust the cursors and then"
	DrawText 21,40,"Click Continue."
	Button button0,pos={80,58},size={92,20},title="Continue"
	Button button0,proc=UserCursorAdjust_ContButtonProc
	
	String graphname = WinName(0,1)
	PauseForUser tmp_PauseforCursors,$graphname
	
	
	NVAR gainm=root:Data:GainMismatch
	Variable gm=gainm
	
	Make/FREE/D/O/N=2 fitc	
	CurveFit/Q line kwCWave=fitc wy[pcsr(A),pcsr(B)] /X=wx
	
	//gm = (vcsr(B)-vcsr(A))/(hcsr(B)-hcsr(A))
	gm = fitc[1]
	gainm = gm
	Ijcor = wy - wx*gm
	
	Duplicate/O Ijcor $corname
	KillWaves Ijcor
	
	if(strsearch(Traces,corname,0)<0)	
		AppendToGraph $corname	 vs wx		
		ModifyGraph rgb($corname)=(0,0,65535)
	endif

End

Function UserCursorAdjust_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K tmp_PauseforCursors // Kill panel
End

Function CorrectTilt()
//Correct gain mismatch between the X and V amps
//If you used GetGoodGainMismatch before
	Variable/G root:Data:GainMismatch
	NVAR gainm=root:Data:GainMismatch
	Variable gm=gainm
	
	String Traces = TraceNameList("",";",1+4)
	String yname = StringFromList(0,Traces)
	Wave/Z wy = TraceNameToWaveRef("",yname)	
	Wave/Z wx = XWaveRefFromTrace("",yname)
	
	String corname = yname + "_cor"
	Duplicate/O wy Ijcor

	
	if(gm==0)
		Prompt gm, "Gain mismatch: "
		DoPrompt "Initialization of gain mismatch", gm
	endif	
	gainm=gm
	
	Ijcor = wy - wx*gm
	
	
	Duplicate/O Ijcor $corname
	KillWaves Ijcor
	
	if(strsearch(Traces,corname,0)<0)	
		AppendToGraph $corname	 vs wx		
		ModifyGraph rgb($corname)=(0,0,65535)	
	endif
	
End




Function SeriesR()
//Corrects a series resistance
	WAVE/Z wI = CsrWaveRef(A)
	WAVE/Z wIb = CsrWaveRef(B)
	if(!WaveExists(wI) || !WaveExists(wIb))
		print("Put both cursors on a trace")
		Return 0
	endif
	
	WAVE/Z wV = CsrXWaveRef(A)
	
	Variable Rs=(hcsr(B)-hcsr(A))/(vcsr(B)-vcsr(A))
	Print "Rs = " + num2str(Rs) + " Ω"
	
	Duplicate/FREE/O wV Vcor
	
	Variable pA=pcsr(A)
	Variable pB=pcsr(B)
	
	if(pA>pB)
		pA=pcsr(B)
		pB=pcsr(A)
	endif
	
	Variable i=0
	For(i=pA;i<=pB;i++)
		Vcor[i]=Vcor[pA]
	Endfor
	
	Duplicate/O Vcor wV
End

Function PlotIVs(V2D,I2D,[nIVs,sym])
//Plot 1D IVs of 2D maps
//nIVs is the number of IVs you want to plot (default is number of flux values). They are linearly spaced
//sym = 1 symmetrizes IVs (default is 0)
//This also creates a wave FluxW containing the currents at which the IVs are taken
	WAVE V2D,I2D
	Variable nIVs
	Variable sym
	Variable ncol = DimSize(V2D,0)
	if(ParamIsDefault(sym))
		sym = 0
	endif
	if(ParamIsDefault(nIVs))
		nIVs = DimSize(V2D,0)
	endif
	Make/FREE/T/O/N=(nIVs) FluxWave
	Make/O/N=(nIVs) FluxW
	Setscale d 0,0,"A",FluxW
	
	Variable npts = DimSize(V2D,1)
	Variable i = 0
	String V1D = ""
	String I1D = ""
	For(i=0;i<nIVs;i++)
		V1D = "Vj_" + num2str(i)
		I1D = "Ij_" + num2str(i)
		
		Make/FREE/O/N=(npts) Vj
		SetScale d 0,0,"V",Vj
		Make/FREE/O/N=(npts) Ij
		SetScale d 0,0,"A",Ij
		
		Vj = V2D[floor(i*(ncol-1)/(nIVs-1))][p]
		Ij = I2D[floor(i*(ncol-1)/(nIVs-1))][p]
		
		FluxW[i] = pnt2x(V2D,floor(i*(ncol-1)/(nIVs-1)))
		FluxWave[i] = num2str(floor(1e6*FluxW[i]))
		
		if(sym==1)
			WaveStats/Q Vj
			Variable xmin = V_min
			Variable xmax = V_max
			WaveStats/Q Ij
			Variable ymin = V_min
			Variable ymax = V_max

			Symmetrize(Vj,Ij,1024,xmin,xmax,ymin,ymax)
		endif
		
		Duplicate/O Vj $V1D
		Duplicate/O Ij $I1D
		
		if(i==0)
			Display $I1D vs $V1D
			Modifygraph grid=1
			TextBox/C/N=text0/J/F=0/A=MC "\\s(" + I1D + ") I\Bflux\M = " + FluxWave[i] + " µA"
		else
			Appendtograph $I1D vs $V1D
			AppendText/N=text0 "\\s('" + I1D + "') I\Bflux\M = " + FluxWave[i] + " µA"
		endif
	Endfor
	
	ColorTab2Wave Rainbow	// creates M_colors
	Wave M_colors
	
	Variable numRows= DimSize(M_colors,0)
	Variable red, green, blue
	Variable index
	for(i=0; i<nIVs; i+=1)
		index = round(i/(nIVs-1) * (numRows-1))	// spread entire color range over all traces.
		ModifyGraph rgb[i]=(M_colors[index][0], M_colors[index][1], M_colors[index][2])
	endfor
	
	KillWaves M_colors	
End


Function GetRLC()
//Get R,L and C for an in-loop mode
	WAVE/Z Iw = CsrWaveRef(A)
	WAVE/Z Iwb = CsrWaveRef(B)
	if(!WaveExists(Iw) || !WaveExists(Iwb))
		print("Put both cursors on a trace")
		Return 0
	endif
	SetDataFolder GetWavesDataFolder(Iw,1)
	
	WAVE/Z Vw = CsrXWaveRef(A)
	
	Variable ptmin=min(pcsr(A),pcsr(B))
	Variable ptmax=max(pcsr(A),pcsr(B))
	
	Make/FREE/D/N=3/O W_coef
	W_coef = {2*vcsr(B),0.01,hcsr(B)}
	FuncFit/Q FitPeakLowCoupling W_coef Iw[pcsr(A),pcsr(B)] /X=Vw
	String fitname = NameOfWave(Iw)+"_fit"
	Make/FREE/N=500 fit
	Setscale/I x,min(hcsr(A),hcsr(B)),max(hcsr(A),hcsr(B)),fit
	fit = FitPeakLowCoupling(W_coef,x)
	
	Duplicate/O fit $fitname
	
	
	String Traces = TraceNameList("",";",1+4)
	
	if(strsearch(Traces,fitname,0)<0)	
		AppendToGraph $fitname	
		ModifyGraph lstyle($fitname)=3,lsize($fitname)=2,rgb($fitname)=(1,12815,52428)	
	endif
	
	
	Variable gam=W_coef[0]
	Variable bet=W_coef[1]
	Variable V0=W_coef[2]
	
	NVAR/Z Icm = Icmax
	Variable Ic = Icm
	if(!NVAR_Exists(Icm))
		Prompt Ic, "Maximal critical current: "
		DoPrompt "Critical current", Ic
	endif
	
	Variable R,L,C
	L = bet*3.291059e-16/(Ic/2)
	C = 2*(3.291059e-16)^2/(V0^2*L)
	R = gam*V0/(Ic/2)
	
	print "Γ = " + num2str(gam)
	print "β = " + num2str(bet)
	print "V0 = " + num2str(V0*1e6) + " µV"
	
	print "L = " + num2str(L*1e12) + " pH"
	print "C = " + num2str(C*1e15) + " fF"
	print "R = " + num2str(R) + " Ω"
End

Function FitPeakLowCoupling(w,V) : FitFunc
	Wave w
	Variable V

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(V) = V/V0/(2*g*(4/b^2*(1-(V/V0)^2)^2+(V/V0)^2/g^2)) + i0
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ V
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = g
	//CurveFitDialog/ w[1] = b
	//CurveFitDialog/ w[2] = V0
	//CurveFitDialog/ w[3] = i0

	return V/w[2]/(2*w[0]*(4/w[1]^2*(1-(V/w[2])^2)^2+(V/w[2])^2/w[0]^2)) + w[3]
End
?

Function FitPeakOffLoop(w,V) : FitFunc
	Wave w
	Variable V
	
		
	// w[0] = Ic
	// w[1] = r
	// w[2] = Z0 = sqrt(L/C)
	// w[3] = Vr
	// w[4] = ioff
	
	//Return 1/2*w[0]^2/(V*w[1])*1/(1+(2*w[2]/w[1])^2*(1-V/w[3])^2)+w[4]
	Return 1/2*w[0]/(1+w[1]^2*(1-V/w[2])^2) + w[3]
End



Function ToBeIntegratedqp(pw,nrj)
	Wave pw
	Variable nrj
	Variable omega=pw[0]
	Variable Delta=pw[1]
	return 1/(khbar*omega)*abs(E*(E+khbar*omega)+Delta^2)/(sqrt(E^2-Delta^2)*sqrt((E+khbar*omega)^2-Delta^2))
End

Function Iqp(V,alpha,Delta)
//Returns the real part of the impedance of a superconductor of gap Delta at the frequency 2eV/h
//Not actually the real part, but something proportional to it
	Variable V,alpha,Delta
	Delta*=kElectronCharge
	Variable omega=V*2*kElectronCharge/khbar
	Make/O/FREE/N=2 pw={omega,Delta}
	return (khbar*omega>2*Delta)*alpha/V*Integrate1D(ToBeIntegratedqp,Delta-khbar*omega,-Delta,0,1,pw)
End

Function qpexcitation(w,V) : FitFunc
//Fit function for the real part of impedance of a superconductor
	Wave w
	Variable V

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(V) = Iqp(V,alpha,Delta) + ioff
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ V
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = alpha
	//CurveFitDialog/ w[1] = Delta
	//CurveFitDialog/ w[2] = ioff

	return Iqp(V,w[0],w[1]) + w[2]
End


Function FitPeakLowCouplingQP(w,V) : FitFunc
//Fit function for the real part of impedance of a superconductor superimposed with a peak
	Wave w
	Variable V

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(V) = V/V0/(2*g*(4/b^2*(1-(V/V0)^2)^2+(V/V0)^2/g^2)) + i0
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ V
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = g
	//CurveFitDialog/ w[1] = b
	//CurveFitDialog/ w[2] = V0
	//CurveFitDialog/ w[3] = i0
	//CurveFitDialog/ w[4] = alpha
	//CurveFitDialog/ w[5] = Delta

	return V/w[2]/(2*w[0]*(4/w[1]^2*(1-(V/w[2])^2)^2+(V/w[2])^2/w[0]^2)) + Iqp(V,w[4],w[5]) + w[3]
End



Function RemoveHalfPoints(w1,w2)
//Remove half points of waves w1 and w2
	Wave w1,w2
	Variable n=numpnts(w1)
	Make/FREE/N=(floor(n/2)) wd1,wd2
	
	Variable i
	
	For(i=0;i<floor(n/2);i++)
		wd1[i]=(w1[2*i]+w1[2*i+1])/2
		wd2[i]=(w2[2*i]+w2[2*i+1])/2
	Endfor
	
	Duplicate/O wd1 w1
	Duplicate/O wd2 w2	
End


Function TrueMod(x)
//Modulo function between -pi and pi
	Variable x
	Variable mo = mod(x,2*pi)
	if(mo<=pi && mo>=-pi)
		Return mo
	else
		if(mo>0)
			Return mo-2*pi
		else
			Return mo+2*pi
		endif
	endif
	
End


Function IcLagr(coefw,px) : FitFunc
//FitFunction for CPR of an asymmetrical SQUID with 
//different inductances in both arms
//based on a Lagrangian multiplicator
	Wave coefw
	Variable px
	// w[0] = i1
	// w[1] = i2
	// w[2] = betal1
	// w[3] = betal2
	// w[4] = phoff
	// w[5] = Dph
	Make/FREE/D/O/N=5 cw={coefw[0],coefw[1],coefw[2],coefw[3],TrueMod((px-coefw[4])/coefw[5])}
	Make/FREE/D/O/N=3 ttl = {pi/2,pi/2,0}
	
	FindRoots/Q/X=ttl L1,cw,L2,cw,L3,cw
	Return coefw[0]*sin(ttl[0])+coefw[1]*sin(ttl[1])
End

Function IcLagrsameL(coefw,px) : FitFunc
//Same as IcLagr, with the same inductance in both arms of the SQUID
	Wave coefw
	Variable px
	// w[0] = i1
	// w[1] = i2
	// w[2] = betal1
	// w[3] = phoff
	// w[4] = Dph
	Make/FREE/D/O/N=5 cw={coefw[0],coefw[1],coefw[2],coefw[2]*coefw[1]/coefw[0],TrueMod((px-coefw[3])/coefw[4])}
	Make/FREE/D/O/N=3 ttl = {pi/2,pi/2,0}
	
	FindRoots/Q/X=ttl L1,cw,L2,cw,L3,cw
	Return coefw[0]*sin(ttl[0])+coefw[1]*sin(ttl[1])
End

Function IcLagrsameLm(w,px) : FitFunc
//Same as IcLagr, with the same inductance in both arms of the SQUID
//Fits the negative switching current
	Wave w
	Variable px

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ // w[0] = i1
	//CurveFitDialog/ 
	//CurveFitDialog/ // w[1] = i2
	//CurveFitDialog/ 
	//CurveFitDialog/ // w[2] = betal1
	//CurveFitDialog/ 
	//CurveFitDialog/ // w[3] = phoff
	//CurveFitDialog/ 
	//CurveFitDialog/ // w[4] = Dph
	//CurveFitDialog/ 
	//CurveFitDialog/ Make/FREE/D/O/N=5 cw={coefw_0,coefw_1,-coefw_2,-coefw_2*coefw_1/coefw_0,TrueMod((px-coefw_3)/coefw_4)}
	//CurveFitDialog/ Make/FREE/D/O/N=3 ttl = {pi/2,pi/2,0}
	//CurveFitDialog/ 
	//CurveFitDialog/ FindRoots/Q/X=ttl L1,cw,L2,cw,L3,cw
	//CurveFitDialog/ f(px) = -(coefw_0*sin(ttl[0])+coefw_1*sin(ttl[1]))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ px
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = coefw_0
	//CurveFitDialog/ w[1] = coefw_1
	//CurveFitDialog/ w[2] = coefw_2
	//CurveFitDialog/ w[3] = coefw_3
	//CurveFitDialog/ w[4] = coefw_4

	// w[0] = i1
	
	// w[1] = i2
	
	// w[2] = betal1
	
	// w[3] = phoff
	
	// w[4] = Dph
	
	Make/FREE/D/O/N=5 cw={w[0],w[1],-w[2],-w[2]*w[1]/w[0],TrueMod((px-w[3])/w[4])}
	Make/FREE/D/O/N=3 ttl = {pi/2,pi/2,0}
	
	FindRoots/Q/X=ttl L1,cw,L2,cw,L3,cw
	return -(w[0]*sin(ttl[0])+w[1]*sin(ttl[1]))
End


Function L1(cw,t1,t2,lambda)
	Wave cw
	Variable t1,t2,lambda
	Variable b1 = cw[2]
	Variable b2 = cw[3]

	Return -cos(t1)*(t1+t2)/2+b1*sin(t1)*cos(t1)+lambda*(1+b1*cos(t1))
End

Function L2(cw,t1,t2,lambda)
	Wave cw
	Variable t1,t2,lambda
	Variable gam = cw[1]/cw[0]
	Variable b1 = cw[2]
	Variable b2 = cw[3]

	Return -gam*cos(t2)*(t1+t2)/2+gam*b2*sin(t2)*cos(t2)-lambda*(1+b2*cos(t2))
End

Function L3(cw,t1,t2,lambda)
	Wave cw
	Variable t1,t2,lambda
	Variable b1 = cw[2]
	Variable b2 = cw[3]
	Variable phix = cw[4]

	Return phix+t1+b1*sin(t1)-t2-b2*sin(t2)
End



Function ToBeSolvedPlas(wp,phi)
	Wave wp
	//wp[0] = betal
	//wp[1] = phix
	Variable phi
	
	Return wp[1]-phi-wp[0]*sin(phi)
End

Function Plasma(phix,betal,wp0)
	Variable phix,betal,wp0
	Make/O/Free/N=2 wp={betal,phix}
	FindRoots/Q ToBeSolvedPlas,wp
	
	//Return wp0*sqrt(abs(sin(V_root)))
	Return wp0*sqrt(1+betal*cos(V_root))
End



Function GetVmax(wV,wI,wout)
//Returns a 1D wave with voltages at which the current is maximal (wV and wI are 2D waves)
//wout doesn't need to be scaled
	Wave wV,wI,wout
	Variable nout = DimSize(wV,0)
	Variable npts = DimSize(wV,1)
	Variable i
	Redimension/N=(nout) wout
	Setscale/P x,DimOffset(wV,0),DimDelta(wV,0),wout
	
	For(i=0;i<nout;i++)
		Make/O/FREE/N=(npts) Ij,Vj
		Vj = wV[i][p]
		Ij = wI[i][p]
		//Smooth/B 10, Vj
		//Smooth/B 10, Ij
		WaveStats/Q Ij
		wout[i]=Vj[V_maxloc]		
	EndFor
End

Function GetImax(wV,wI,wout)
//Returns a 1D wave with the maximal currents of slices of wI (wV and wI are 2D waves)
//wout doesn't need to be scaled
	Wave wV,wI,wout
	Variable nout = DimSize(wV,0)
	Variable npts = DimSize(wV,1)
	Variable i
	Redimension/N=(nout) wout
	Setscale/P x,DimOffset(wV,0),DimDelta(wV,0),wout
	
	For(i=0;i<nout;i++)
		Make/O/FREE/N=(npts) Ij,Vj
		Vj = wV[i][p]
		Ij = wI[i][p]
		//Smooth/B 10, Vj
		//Smooth/B 10, Ij
		WaveStats/Q Ij
		wout[i]=V_max	
	EndFor
End


Function GetIsw(wV,wI,wout)
//Returns a 1D wave with the switching currents (max of current for voltage below 5 µV)
//(wV and wI are 2D waves)
//wout doesn't need to be scaled
	Wave wV,wI,wout
	Variable nout = DimSize(wV,0)
	Variable npts = DimSize(wV,1)
	Variable i
	Redimension/N=(nout) wout
	Setscale/P x,DimOffset(wV,0),DimDelta(wV,0),"A",wout
	SetScale d 0,0,"A", wout
	
	For(i=0;i<nout;i++)
		Make/O/FREE/N=(npts) Ij,Vj
		Vj = wV[i][p]
		Ij = wI[i][p]*(abs(Vj[p])<5e-6)
		//Smooth/B 10, Vj
		//Smooth/B 10, Ij
		WaveStats/Q Ij
		wout[i]=V_max	
	EndFor
End

Function GetIswm(wV,wI,wout)
//Returns a 1D wave with the switching currents (max of current for voltage below 5 µV)
//(wV and wI are 2D waves)
//wout doesn't need to be scaled
	Wave wV,wI,wout
	Variable nout = DimSize(wV,0)
	Variable npts = DimSize(wV,1)
	Variable i
	Redimension/N=(nout) wout
	Setscale/P x,DimOffset(wV,0),DimDelta(wV,0),wout
	
	For(i=0;i<nout;i++)
		Make/O/FREE/N=(npts) Ij,Vj
		Vj = wV[i][p]
		Ij = wI[i][p]*(abs(Vj[p])<5e-6)
		//Smooth/B 10, Vj
		//Smooth/B 10, Ij
		WaveStats/Q Ij
		wout[i]=V_min	
	EndFor
End

Function SuperIsw()
//Extract switching current from a color plot of current (not map)
//It's better to have it flatenned beforehand
//You will have to choose the corresponding voltage wave
//In the end, you will get a plot of Iswitching vs Iflux
	String list= ImageNameList("",";")
	String imagePlot = StringFromList(0, list)
	if (strlen(imagePlot))	// one image
		String info,Iwave,df
		info=ImageInfo("",imagePlot,0)
		Iwave=StringByKey("ZWAVE",info)
		df=StringByKey("ZWAVEDF",info)
	else	// no images
		Abort "No Image in graph"
	endif
	
	SetDataFolder df
	DFREF dfr = GetDataFolderDFR()	
	
	Variable index = 0
	String ListOfWaves = ""
	do
		Wave/Z w = WaveRefIndexedDFR(dfr, index)
		if (!WaveExists(w))
			break
		endif
		index += 1
		ListOfWaves += NameOfWave(w) + ";"
	while(1)
	String Vwave
	
	Prompt Vwave,"Name of the V wave:",popup,ListOfWaves
	DoPrompt "Extracting switching current...",Vwave
	
	String Iswwave = UniqueName(CleanupName(Iwave + "_sw",1),1,0)
	
	Make $Iswwave


	GetIsw($Vwave,$Iwave,$Iswwave)

	Display $Iswwave
	Modifygraph grid=1
	SetAxis left 0,*
	Label left "I\\Bswitching\\M (\\U)"
	Label bottom "I\\Bflux\\M (\\U)"
	ModifyGraph mirror=2
End

Function PlasmaRF(w,ifl) : FitFunc
//Fit function for the resonance of a RF-SQUID
	Wave w
	Variable ifl

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ ifl
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = ioffset
	//CurveFitDialog/ w[1] = itoflux
	//CurveFitDialog/ w[2] = bl
	//CurveFitDialog/ w[3] = wp0


	return Plasma((ifl-w[0])/w[1],w[2],w[3]/sqrt(w[2]))
End


Function GetIp(wV,wI,wout)
//Similar to GetImax but fits each slice with a Lorentzian
//Returns a 1D wave with the currents at the top of the Lorentzian fits of slices of wI
//wV and wI are 2D waves
//wout doesn't need to be scaled
	Wave wV,wI,wout
	Variable nout = DimSize(wV,0)
	Variable npts = DimSize(wV,1)
	Variable i
	Redimension/N=(nout) wout
	Setscale/P x,DimOffset(wV,0),DimDelta(wV,0),wout
	
	
	For(i=0;i<nout;i++)
		Make/O/FREE/N=(npts) Ij,Vj
		Vj = wV[i][p]
		Ij = wI[i][p]
		Make/D/O/FREE/N=4 W_coef
		CurveFit/Q lor kwCWave=W_coef Ij /X=Vj
		//wout[i] = W_coef[0]+W_coef[1]/W_coef[3]
		wout[i] = W_coef[1]/W_coef[3]
	EndFor
End

Function GetVp(wV,wI,wout)
//Similar to GetVmax but fits each slice with a Lorentzian
//Returns a 1D wave with the voltage at the top of the Lorentzian fits of slices of wI
//wV and wI are 2D waves
//wout doesn't need to be scaled
	Wave wV,wI,wout
	Variable nout = DimSize(wV,0)
	Variable npts = DimSize(wV,1)
	Variable i
	Redimension/N=(nout) wout
	Setscale/P x,DimOffset(wV,0),DimDelta(wV,0),wout	
	
	For(i=0;i<nout;i++)
		Make/O/FREE/N=(npts) Ij,Vj
		Vj = wV[i][p]
		Ij = wI[i][p]
		Make/D/O/FREE/N=4 W_coef
		CurveFit/Q lor kwCWave=W_coef Ij /X=Vj
		//wout[i] = W_coef[0]+W_coef[1]/W_coef[3]
		wout[i] = W_coef[2]
	EndFor
End

Function FindPeaks(wV,wI,wVout,wIout)
//Smooths (Box of 30 points) Ij vs Vj slices and fits them with Lorentzians
//Return both V and I of peaks
//Return NaN when the fit doesn't converge well
	Wave wV,wI,wVout,wIout
	Variable nout = DimSize(wV,0)
	Variable npts = DimSize(wV,1)
	Variable i
	Redimension/N=(nout) wVout
	Redimension/N=(nout) wIout
	Setscale/P x,DimOffset(wV,0),DimDelta(wV,0),"A",wVout
	SetScale d 0,0,"V",wVout
	Setscale/P x,DimOffset(wV,0),DimDelta(wV,0),"A",wIout
	SetScale d 0,0,"A",wIout
	
	For(i=0;i<nout;i++)
		Make/O/FREE/N=(npts) Ij,Vj
		Vj = wV[i][p]
		Ij = wI[i][p]
		Make/D/O/N=4 W_coef		
		Make/D/O/N=4 W_sigma
		Smooth/B 20, Vj
		Smooth/B 20, Ij
		CurveFit/Q lor kwCWave=W_coef Ij /X=Vj
		//wout[i] = W_coef[0]+W_coef[1]/W_coef[3]
		Make/D/O/N=4 W_err
		W_err = abs(W_sigma/W_coef)
		if ((sum(W_err)==Nan) || (sum(W_err)>1) || W_coef[1]<0 || W_coef[2]<0 || W_coef[3]<0)
			wVout[i] = Nan
			wIout[i] = Nan
		else
			wVout[i] = W_coef[2]
			wIout[i] = W_coef[1]/W_coef[3]
		endif
	EndFor
End

Function FindPeaksMax(wV,wI,wVout,wIout)
//Find max of Ij vs Vj slices
//Return both V and I of peaks
	Wave wV,wI,wVout,wIout
	Variable nout = DimSize(wV,0)
	Variable npts = DimSize(wV,1)
	Variable i
	Redimension/N=(nout) wVout
	Redimension/N=(nout) wIout
	Setscale/P x,DimOffset(wV,0),DimDelta(wV,0),"A",wVout
	SetScale d 0,0,"V",wVout
	Setscale/P x,DimOffset(wV,0),DimDelta(wV,0),"A",wIout
	SetScale d 0,0,"A",wIout
	
	
	For(i=0;i<nout;i++)
		Make/O/FREE/N=(npts) Ij,Vj
		Vj = wV[i][p]
		Ij = wI[i][p]
		WaveStats/Q Ij
		wIout[i]=V_max
		wVout[i]=Vj[V_maxloc]
	EndFor
End

Proc FindPeaksExtra() : GraphMarquee
Silent 1;PauseUpdate
	String list= ImageNameList("",";")
	if( strlen(StringFromList(1, list)) )	// two or more images, comment this test out to make it work with the first of several images
		Abort "CopyImageSubset works only with a single image in a graph"
	endif
	String imagePlot = StringFromList(0, list)
	if (strlen(imagePlot))	// one image
		String info,vaxis,haxis,image,df,xwave,ywave
		Variable i
		info=ImageInfo("",imagePlot,0)
		vaxis=StringByKey("YAXIS",info)
		haxis=StringByKey("XAXIS",info)
		image=StringByKey("ZWAVE",info)
		df=StringByKey("ZWAVEDF",info)
		xwave=StringByKey("XWAVE",info)
		ywave=StringByKey("YWAVE",info)
		if( strlen(xwave)+strlen(ywave) )
			Abort "CopyImageSubset does not work on images with X or Y waves"
		endif
		String modInfo= WMGetRECREATIONFromInfo(info)
		info= AxisInfo("",haxis)
		String axisFlags= StringByKey("AXFLAG",info)
		info= AxisInfo("",vaxis)
		axisFlags+=StringByKey("AXFLAG",info)
		String winStyle= WinRecreation("",1)
		Variable swapxy= strsearch(winStyle,"swapXY=1",0) >= 0
		GetMarquee/K $haxis,$vaxis
		if( V_Flag != 1)
			Abort "CopyImageSubset() requires a marquee in the target graph"
		endif
		// copy the marqueed image subset into a clean, liberal, unique name.
		String copy= UniqueName(CleanupName(image+"Copy",1),1,0)
		CopyMatrix(df+PossiblyQuoteName(image),copy,V_left,V_right,V_top,V_bottom,swapxy)
		// make a graph just like the source graph
		copy=PossiblyQuoteName(copy)
		Preferences 1
		Display
		String cmd="AppendImage"+axisFlags+" "+copy
		Execute cmd
		if( strlen(modInfo) )
			// modInfo has ALL the recreation pieces separated by semicolons, not commas
			modInfo= RemoveEnding(ReplaceString(";", modInfo, ","),",")
			cmd="ModifyImage "+copy+" "+modInfo
			Execute cmd 
		endif
		ApplyStyleMacro(winStyle)
		
		DoAutoSizeImage(0,-1)			// -1 for flip vert flag new. Means don't change current setting.
	else	// no images
		Abort "No Image in graph"
	endif
	
	SetDataFolder df
	GetVname(V_left,V_right,V_top,V_bottom,swapxy,copy)
End

Function GetVName(V_left,V_right,V_top,V_bottom,swapxy,copyI)
	Variable V_left,V_right,V_top,V_bottom,swapxy
	String copyI
	DFREF dfr = GetDataFolderDFR()	
	Variable index = 0
	String ListOfWaves = ""
	do
		Wave/Z w = WaveRefIndexedDFR(dfr, index)
		if (!WaveExists(w))
			break
		endif
		index += 1
		ListOfWaves += NameOfWave(w) + ";"
	while(1)
	String imageV
	
	Prompt imageV,"Name of the V wave:",popup,ListOfWaves
	DoPrompt "Extracting peak profile",imageV
	
	String copyV= UniqueName(CleanupName(imageV+"Copy",1),1,0)
	CopyMatrixf(PossiblyQuoteName(imageV),copyV,V_left,V_right,V_top,V_bottom,swapxy)
	copyV=PossiblyQuoteName(copyV)

	String Vpeak = copyV + "_peak"
	String Ipeak	= copyI + "_peak"
	Make $Vpeak
	Make $Ipeak

	String lormax = ""
	String lorormax = "Lorentzian fit;Find maximum of current"
	Prompt lormax,"Fit with Lorentzians or find max of current?",popup,lorormax
	DoPrompt "Extracting peak profile...",lormax

	if(CmpStr(lormax,"Lorentzian fit")==0)
		FindPeaks($copyV,$copyI,$Vpeak,$Ipeak)
	else
		FindPeaksMax($copyV,$copyI,$Vpeak,$Ipeak)
	endif	
	
	Display $Vpeak
	Modifygraph grid=1
	Label left "Position of the peak (\\U)"
	Label bottom "I\\Bflux\\M (\\U)"
	ModifyGraph mirror(bottom)=2
	
	NewFreeAxis/R rightax
	Label rightax "Frequency (\\U)"
	ModifyFreeAxis rightax, master=left, hook=VoltToHz
	ModifyGraph tick(rightax)=0,lblPos(rightax)=50,freePos(rightax)={0,kwFraction}	
	
	Display $Ipeak
	Modifygraph grid=1
	Label left "Height of the peak (\\U)"
	Label bottom "I\\Bflux\\M (\\U)"
	ModifyGraph mirror=2
End



Function CopyMatrixf(inmatrix,outmatrix,left,right,top,bottom,swapxy)
	String inmatrix,outmatrix
	Variable left,right,top,bottom,swapxy
	
	if( swapxy )
		Duplicate/O/R=(bottom,top)(left,right) $inmatrix,$outmatrix
	else
		Duplicate/O/R=(left,right)(bottom,top) $inmatrix,$outmatrix
	endif
End


Function NiceIV()
	ModifyGraph width=340.157,height=226.772
	Modifygraph grid=1
	Label left "Current (\\U)"
	Label bottom "Voltage (\\U)"
	ModifyGraph mirror=2
	Modifygraph mode=0,lsize=1
End

Function Legend0Pi()	
	String traceList = TraceNameList(WinName(0,1),";",5)
	
	Variable trace0 = 0
	Prompt trace0,"Name of wave at zero flux:",popup,traceList
	DoPrompt "Coloring traces",trace0
	
	Variable tracepi = 0
	Prompt tracepi,"Name of wave at pi flux:",popup,traceList
	DoPrompt "Coloring traces",tracepi
	
	
	if(trace0 != 0 && tracepi != 0)
		ModifyGraph rgb[trace0-1]=(65535,0,0)
		ModifyGraph rgb[tracepi-1]=(0,16019,65535)
	
		Legend/C/N=text0/J/F=0/A=LT "\\s(" + StringFromList(trace0-1,traceList) + ") φ = 0\r\\s(" + StringFromList(tracepi-1,traceList) + ") φ = π"
	endif
End

Function MakeNicePlotFromYkdl()
	String Traces = TraceNameList("",";",5)
	String yname = StringFromList(0,Traces)
	Wave/Z wy = TraceNameToWaveRef("",yname)	
	Wave/Z wx = XWaveRefFromTrace("",yname)

	Display wy vs wx
	NiceIV()
End


Function VoltToHz(info)
	STRUCT WMAxisHookStruct &info
	GetAxis/Q/W=$info.win $info.mastName	// get master (bottom) axis' range in V_min, V_Max
	Variable minF =  V_min*2*1.6e-19/6.63e-34
	Variable maxF =  V_max*2*1.6e-19/6.63e-34
	info.min = minF	// new min for free axis
	info.max= maxF	// new max for free axis
	info.units = "Hz"
	return 0
End

Function AmpereToHz(info)
	STRUCT WMAxisHookStruct &info
	GetAxis/Q/W=$info.win $info.mastName	// get master (bottom) axis' range in V_min, V_Max
	Variable minF =  V_min/(2*1.6e-19)
	Variable maxF =  V_max/(2*1.6e-19)
	info.min = minF	// new min for free axis
	info.max= maxF	// new max for free axis
	info.units = "Hz"
	return 0
End

//Modification by RHR 29-11-2019:
Function AmpereToSens(info)
	STRUCT WMAxisHookStruct &info
	GetAxis/Q/W=$info.win $info.mastName	// get master (bottom) axis' range in V_min, V_Max
	Variable minF =  sqrt(V_min/(2*1.6e-19))
	Variable maxF =  sqrt(V_max/(2*1.6e-19))
	info.min = minF	// new min for free axis
	info.max= maxF	// new max for free axis
	info.units = "Hz"
	return 0
End

Function NiceMap()
	ModifyGraph width=226.772*2,height=226.772*2
	ModifyGraph margin(right)=150
	Label left "Voltage (\\U)"
	Label bottom "Coil current (\\U)"
	ModifyGraph mirror=2
	ColorScale/C/A=RC/E/N=text0 "Current (\\U)"
	//ColorScale/C/N=text0 "Current (\\U)"	
	
	String plotwave = StringFromList(0, ImageNameList("",";"))
	
	Wavestats/Q $plotwave
	If(V_min<-V_max/3)
	//if this condition is true, this map is prone to be a positive/negative map
	//otherwise, it is more likely to be positive only
		//ModifyImage $plotwave ctab= {*,*,:Colors:RedWhiteBlue,0}
		ModifyImage $plotwave ctab= {-V_max,V_max,root:Colors:GreenBrown,0}
	else
		ModifyImage $plotwave ctab= {0,V_max,root:Colors:Plasm,1}
	endif	
End

Function PlotAll2DWaves()
	DFREF dfr = GetDataFolderDFR()	
	
	Variable index = 0
	String ListOfWaves = ""
	do
		Wave/Z w = WaveRefIndexedDFR(dfr, index)
		if (!WaveExists(w))
			break
		endif
		index += 1
		Display
		Appendimage w
		ModifyImage $NameOfWave(w) ctab= {*,*,root:Colors:GreenBrown,0}
	while(1)
End

Function AddFreqAxis()
	String ImageList = ImageNameList("",";")
	
	If(CmpStr(ImageList,""))
	//In that case, there is an image
	//Axis is to be added on the right
		Modifygraph mirror(left)=0
		ModifyGraph mirror(bottom)=2,axisOnTop=0
		NewFreeAxis/R rightax
		Label rightax "Frequency (\\U)"
		ModifyFreeAxis rightax, master=left, hook=VoltToHz
		ModifyGraph tick(rightax)=0,lblPos(rightax)=50,freePos(rightax)={0,kwFraction}
	
	else	
		//In that case, there is no image. Only 1D traces
		//Axis is to be added at the top
		ModifyGraph mirror(bottom)=0
		ModifyGraph mirror(left)=2,axisOnTop=0	
		NewFreeAxis/T topax
		Label topax "Frequency (\\U)"
		ModifyFreeAxis topax, master=bottom, hook=VoltToHz
		ModifyGraph tick(topax)=0,lblPos(topax)=50,freePos(topax)={0,kwFraction}	
	Endif
End

Function AddFreqAxisCurrent()
	ModifyGraph mirror(left)=0
	ModifyGraph mirror(bottom)=2,axisOnTop=0	
	NewFreeAxis/R rightax
	Label rightax "Frequency (\\U)"
	ModifyFreeAxis rightax, master=left, hook=AmpereToHz
	ModifyGraph tick(rightax)=0,lblPos(rightax)=50,freePos(rightax)={0,kwFraction}	
End

//by Ramiro: 
//AmpereToSens(info) function
Function AddSensAxisCurrent()
	ModifyGraph mirror(left)=0
	ModifyGraph mirror(bottom)=2,axisOnTop=0	
	NewFreeAxis/R rightax
	Label rightax "Sensibility (\\U)"
	ModifyFreeAxis rightax, master=left, hook=AmpereToSens
	ModifyGraph tick(rightax)=0,lblPos(rightax)=50,freePos(rightax)={0,kwFraction}	
End

Function CurrentToFlux(info)
	STRUCT WMAxisHookStruct &info
	Variable Ioff = 0
	Variable alpha = 1
	NVAR Iof = root:Data:io
	NVAR alph = root:Data:al
	Ioff = Iof
	alpha = alph
	GetAxis/Q/W=$info.win $info.mastName	// get master (bottom) axis' range in V_min, V_Max
	Variable minF =  (V_min-Ioff)*alpha
	Variable maxF =  (V_max-Ioff)*alpha
	info.min = minF	// new min for free axis
	info.max= maxF	// new max for free axis
	info.units = ""
	return 0
End

Function TopAxisFlux()
//Adds a top axis showing flux for a map where bottom axis is current
	ShowInfo
	
	NewPanel/K=2 /W=(187,368,437,531) as "Pause for Cursors"
	DoWindow/C tmp_PauseforCursors2
	DrawText 71,40,"Adjust the cursors to"
	DrawText 101,60,"-π and π"
	DrawText 61,80,"and then click Continue."
	Button button0,pos={80,110},size={92,20},title="Continue"
	Button button0,proc=UserCursorAdjust
	
	String graphname = WinName(0,1)
	PauseForUser tmp_PauseforCursors2,$graphname
	

	if(CmpStr(CsrInfo(A),"") && CmpStr(CsrInfo(B),""))
	//CmpStr compares strings and returns 0 if they are equal
		Variable Ipi = max(xcsr(A),xcsr(B))
		Variable Impi = min(xcsr(A),xcsr(B))
	
		Variable/G root:Data:io = (Impi+Ipi)/2
		Variable/G root:Data:al = 1/(Ipi-Impi)
		
		Cursor/M/K A //Hides cursor A
		Cursor/M/K B //Hides cursor B
		HideInfo
	
		ModifyGraph mirror(bottom)=0
		NewFreeAxis/O/T topax
		Label topax "Flux in the SQUID (Φ\B0\M)"
		ModifyFreeAxis topax, master=bottom, hook=CurrentToFlux
		ModifyGraph tick(topax)=0,lblPos(topax)=50,freePos(topax)={0,kwFraction}
		ModifyGraph manTick(topax)={0,0.5,0,1},manMinor(topax)={0,50}
	endif
End

Function UserCursorAdjust(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K tmp_PauseforCursors2 // Kill panel
End



Function SuperXYZ()
//Performs XYZ3Matrix
	String list= ImageNameList("",";")
	String imagePlot = StringFromList(0, list)
	if (strlen(imagePlot))	// one image
		String info,Iwave,df
		info=ImageInfo("",imagePlot,0)
		Iwave=StringByKey("ZWAVE",info)
		df=StringByKey("ZWAVEDF",info)
	else	// no images
		Abort "No Image in graph"
	endif
	
	SetDataFolder df
	DFREF dfr = GetDataFolderDFR()	
	
	Variable index = 0
	String ListOfWaves = ""
	do
		Wave/Z w = WaveRefIndexedDFR(dfr, index)
		if (!WaveExists(w))
			break
		endif
		index += 1
		ListOfWaves += NameOfWave(w) + ";"
	while(1)
	String Vwave
	
	Prompt Vwave,"Name of the V wave",popup,ListOfWaves
	DoPrompt "Auto XYZ3Matrix 1/2",Vwave

	
	Duplicate/FREE/O $Iwave Icoil
	Icoil=x

	String Map = UniqueName(CleanupName(df + "_map",1),1,0)
	
	String listny="128;256;512;1024"
	
	Variable nypts
	Prompt nypts,"Number of points in the y direction",popup,listny
	DoPrompt "Auto XYZ3Matrix 2/2",nypts
	
	nypts=2^(6+nypts)
		
	XYZ3Matrix(Icoil,$Vwave,$Iwave,Map,nx=Dimsize(Icoil,0),ny=nypts)
	WaveStats/Q $Map

	Display
	Appendimage $Map
	
	
	If(V_min<-V_max/3) 
	//if this condition is true, this map is prone to be a positive/negative map
	//otherwise, it is more likely to be positive only
		//ModifyImage $Map ctab= {-V_max,V_max,RedWhiteBlue,0}
		ModifyImage $Map ctab= {-V_max,V_max,root:Colors:GreenBrown,0}
	else
		ModifyImage $Map ctab= {0,V_max,root:Colors:Plasm,1}
	endif

	Label left "Voltage (\\U)"
	Label bottom "Flux current (\\U)"
	ModifyGraph mirror=2	
End


Function GetRvst(Vw,Iw,Rname)
//Creates a wave with resistance
//Name of this wave is Rname
	Wave Vw,Iw
	String Rname
	Variable i,npts,nt
	npts = Dimsize(Vw,1)
	nt = DimSize(Vw,0)
	
	Make/FREE/O/N=(nt) Rtemp
	SetScale d 0,0,"Ω", Rtemp
	
	For(i=0;i<nt;i++)
		Make/O/FREE/N=(npts) Vtmp,Itmp
		Vtmp = Vw[i][p]
		Itmp = Iw[i][p]
		
		Make/O/FREE/N=2 coefs
		
		CurveFit/Q line,kwCWave=coefs,Vtmp /X=Itmp
		
		Rtemp[i] = coefs[1]
	EndFor
	
	Wavestats/Q Rtemp
	
	Duplicate/O Rtemp $Rname
	Display $Rname
	Modifygraph grid=1
	Label left "Resistance (\\U)"
	SetAxis left 0,1.1*V_max
	ModifyGraph mirror=2
	Label bottom "Time (a.u.)"
End

Function GetRvsTemp(Vw,Iw,Tw,Rname,Tname)
//Creates a wave with resistance and one with temperature
//Vw and Iw are 2D waves acquired while Temperature is changing and recorded in Tw
//Works only if Tw was recorded on the oscilloscope as a 2D wave
//Names of these wave are Rname and Tname
	Wave Vw,Iw,Tw
	String Rname,Tname
	Variable i,npts,nt
	npts = Dimsize(Vw,1)
	nt = DimSize(Vw,0)
	
	Make/FREE/O/N=(nt) Rtemp
	SetScale d 0,0,"Ω", Rtemp
	Make/FREE/O/N=(nt) Ttemp
	SetScale d 0,0,"K", Ttemp
	
	For(i=0;i<nt;i++)
		Make/O/FREE/N=(npts) Vtmp,Itmp,Ttmp
		Vtmp = Vw[i][p]
		Itmp = Iw[i][p]
		Ttmp = Tw[i][p]
		
		Make/O/FREE/N=2 coefs		
		CurveFit/Q line,kwCWave=coefs,Vtmp /X=Itmp		
		Rtemp[i] = coefs[1]
		
		Wavestats/Q Ttmp
		Ttemp[i] = V_avg
	EndFor
	
	Wavestats/Q Rtemp
	
	Duplicate/O Rtemp $Rname
	Duplicate/O Ttemp $Tname
	Display $Rname vs $Tname
	Modifygraph grid=1
	Label left "Resistance (\\U)"
	SetAxis left 0,1.1*V_max
	ModifyGraph mirror=2
	Label bottom "Temperature (\\U)"
End

Function GetRvsTemp1(Vw,Iw,Tw,Rname,Tname)
//Creates a wave with resistance and one with temperature
//Vw and Iw are 2D waves acquired while Temperature is changing and recorded in Tw
//Here Tw is a 1D wave starting and stopping at the time as Vw and Iw, but possibly has a different number of points
//Names of the returned waves are Rname and Tname
	Wave Vw,Iw,Tw
	String Rname,Tname
	Variable i,npts,nt,ntT,startT,endT
	npts = Dimsize(Vw,1)
	nt = DimSize(Vw,0)
	startT = DimOffset(Vw,0)
	endT = DimOffset(Vw,0) + (nt-1)*DimDelta(Vw,0)
	
	Make/FREE/O/N=(nt) Rtemp
	SetScale d 0,0,"Ω", Rtemp
	Make/FREE/O/N=(nt) Ttemp
	SetScale d 0,0,"K", Ttemp
	Duplicate/FREE/O Tw Tw2
	Setscale/I x,startT,endT,Tw2
	
	For(i=0;i<nt;i++)
		Make/O/FREE/N=(npts) Vtmp,Itmp
		Vtmp = Vw[i][p]
		Itmp = Iw[i][p]
		Variable tim = IndexToScale(Vw,i,0)
		Ttemp[i] = Tw2(tim)
		
		Make/O/FREE/N=2 coefs		
		CurveFit/Q line,kwCWave=coefs,Vtmp /X=Itmp		
		Rtemp[i] = coefs[1]
	EndFor
	
	Wavestats/Q Rtemp
	
	Duplicate/O Rtemp $Rname
	Duplicate/O Ttemp $Tname
	Display $Rname vs $Tname
	Modifygraph grid=1
	Label left "Resistance (\\U)"
	SetAxis left 0,1.1*V_max
	ModifyGraph mirror=2
	Label bottom "Temperature (\\U)"
End

