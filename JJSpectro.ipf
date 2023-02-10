#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Menu "JJ Spectro"

	"Correct for a series resistance/7", SeriesR()
	"-"
	"Plot 1D IVs", PlotIVs()
	"-"
	"Get R, L and C for a in-loop mode/8", GetRLC()
End

Function SeriesR()
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

Function PlotIVs(V2D,I2D,[FluxWave,nIVs])
	WAVE V2D,I2D
	WAVE/T FluxWave
	Variable nIVs
	Variable ncol = DimSize(V2D,0)
	if(ParamIsDefault(nIVs))
		nIVs = DimSize(V2D,0)
	endif
	if(ParamIsDefault(FluxWave))
		Make/FREE/O/N=(nIVs) FluxW
		FluxW = p
		Make/FREE/T/O/N=(nIVs) FluxWave
		FluxWave = num2str(FluxW)
	endif
	Variable npts = DimSize(V2D,1)
	Variable i = 0
	String V1D = ""
	String I1D = ""
	For(i=0;i<nIVs;i++)
		V1D = "Vj_" + FluxWave[i]
		I1D = "Ij_" + FluxWave[i]
		
		Make/FREE/O/N=(npts) Vj
		SetScale d 0,0,"V",Vj
		Make/FREE/O/N=(npts) Ij
		SetScale d 0,0,"A",Ij
		
		Vj = V2D[floor(i*ncol/nIVs)][p]
		Ij = I2D[floor(i*ncol/nIVs)][p]
		
		Duplicate/O Vj $V1D
		Duplicate/O Ij $I1D
		
		if(i==0)
			Display $I1D vs $V1D
			Modifygraph grid=1
			TextBox/C/N=text0/J/F=0/A=MC "\\s(" + I1D + ") φ = " + FluxWave[i] 
		else
			Appendtograph $I1D vs $V1D
			AppendText/N=text0 "\\s('" + I1D + "') φ = " + FluxWave[i] 
		endif
	Endfor
	//ShowKBColorizePanel()
End


Function GetRLC()
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
	//CurveFitDialog/ f(V) = V/V0/(2*g*(4/b^2*(1-(V/V0)^2)^2+(V/V0)^2/g^2))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ V
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = g
	//CurveFitDialog/ w[1] = b
	//CurveFitDialog/ w[2] = V0

	return V/w[2]/(2*w[0]*(4/w[1]^2*(1-(V/w[2])^2)^2+(V/w[2])^2/w[0]^2))
End

Function FitPeakHighCoupling(w,V,i) : FitFunc
//Not to be used with Analysis:Curve Fitting. Use it with these command lines.
//Make/D/O coefs={1,0.1,100e-6}
//Duplicate Vr Vfit
//Duplicate ir ifit
//FuncFit/ODR=3 FitPeakHighCoupling, coefs /X={Vr,ir}/XD={Vfit,ifit}

	Wave w
	Variable V
	Variable i

	//CurveFitDialog/
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = g
	//CurveFitDialog/ w[1] = b
	//CurveFitDialog/ w[2] = V0
	
	return -1+1*(2*(1-(V/w[2])^2)/w[1]*sqrt(2*(i)*w[0]*w[2]/V)/(BesselJ(0,sqrt(2*(i)*w[0]*w[2]/V))-BesselJ(2,sqrt(2*(i)*w[0]*w[2]/V))))^2+((V/w[2])/w[0]*sqrt(2*(i)*w[0]*w[2]/V)/(BesselJ(0,sqrt(2*(i)*w[0]*w[2]/V))+BesselJ(2,sqrt(2*(i)*w[0]*w[2]/V))))^2
End


Function ToBeSolvedResPeak(w,del)
	Wave w
	Variable del
	Variable gam=w[0]
	Variable bet=w[1]
	Variable v=w[2]
	
	Return -1+1*(2*(1-v^2)/bet*del/(BesselJ(0,del)-BesselJ(2,del)))^2+(v/gam*del/(BesselJ(0,del)+BesselJ(2,del)))^2
End

Function PeakHighCoupling(gam,bet,x)
	Variable gam,bet,x
	Make/FREE/N=3 param={gam,bet,x}
	
	FindRoots/Q/L=0/H=1.84 ToBeSolvedResPeak,param
	Variable del=V_root
	
	Return del^2*x/(2*gam)
End

Function PeakHighCoupling2(gam,bet,x)
	Variable gam,bet,x
	Make/FREE/N=3 param={gam,bet,x}
	
	FindRoots/Q/L=1.84/H=3.83 ToBeSolvedPeak,param
	
	Variable del=V_root
	
	Return del^2*x/(2*gam)
End

Function PeakHighCoupling3(gam,bet,x)
	Variable gam,bet,x
	Make/FREE/N=3 param={gam,bet,x}
	
	FindRoots/Q/L=1.84/H=3.83 ToBeSolvedPeak,param
	
	Variable del=V_root
	
	if(V_numRoots==2)
		del=V_Root2
	endif
	
	Return del^2*x/(2*gam)
End

Function ToBeIntegratedqp(pw,nrj)
	Wave pw
	Variable nrj
	Variable omega=pw[0]
	Variable Delta=pw[1]
	return 1/(khbar*omega)*abs(E*(E+khbar*omega)+Delta^2)/(sqrt(E^2-Delta^2)*sqrt((E+khbar*omega)^2-Delta^2))
End

Function Iqp(V,alpha,Delta)
	Variable V,alpha,Delta
	Delta*=kElectronCharge
	Variable omega=V*2*kElectronCharge/khbar
	Make/O/FREE/N=2 pw={omega,Delta}
	return (khbar*omega>2*Delta)*alpha/V*Integrate1D(ToBeIntegratedqp,Delta-khbar*omega,-Delta,0,1,pw)
End
Function qpexcitation(w,V) : FitFunc
	Wave w
	Variable V

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(V) = alpha*V+Delta
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ V
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = alpha
	//CurveFitDialog/ w[1] = Delta

	return Iqp(V,w[0],w[1])
End



Function RemoveHalfPoints(w1,w2)
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