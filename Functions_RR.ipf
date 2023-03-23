#pragma TextEncoding = "Windows-1252"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Function:
// EllipseCirc(LongAxis,ShortAxis)
// seriesC(capa1,capa2)
// freqLC(indL,capa)
// ExtractLayer(wave3D)
// xtract1Layer(wave3D,NL)
// xtractZoom(wave2D)
// AreaJJ_v3(Map,zmin)
// AreaJJ_v4(Map,zmin)
// MeanValues(dName1,wavesNum)
// PlotLinesProfiles(wave_ij,wave_vj,n1,n2,dn,icoil)
// AddRetrapingToMap(w_up,w_down,extraN,wNewName)
// AddRetrapingToMap_Max(w_up,w_down,v0,dV,Nh)
// AlignPeak(wV,wI,Vmin,Vmax,Voffset)
// GetV_when_Iths(wV,wI,wout,Itsh)
// HideRetraping(w_up,wNewName,dVhide,V0,dv,dN)
// HideSwitching(wI,wV,wNewName,V0,dv,dN)
// AlignAtMaxV(wV,Voffset)
// CutRangeVJ_v0(waveI,waveV,Vi,Vf)
// CutRangeVJ_v1(waveI,waveV,Vi,Vf,RetSlope)




Function EllipseCirc(LongAxis,ShortAxis)
	Variable LongAxis,ShortAxis
	Variable a,b,P1,h,P2
	
	a=LongAxis*0.5
	b=ShortAxis*0.5
	
	//Perimeter aprox #1:
	P1 = pi*( 3*(a+b) - sqrt((3*a+b)*(a+3*b)) )
	
	//Perimeter aprox #2:
	
	h=(a-b)^2/(a+b)^2
	P2=pi*(a+b)*(1+ 3*h/(10+sqrt(4-3*h)) )
	
	print "Per_1 =", P1
	print "Per_2 =", P2
	
	return P1
End

Function seriesC(capa1,capa2)
	Variable capa1,capa2
	Variable capa

	capa = capa1*capa2/(capa1+capa2)
	
	print "SerieCapa =",capa,"Farad"
	return capa
End

Function freqLC(indL,capa)
	Variable indL,capa
	Variable freq

	freq = 1/(2*pi*sqrt(indL*capa))/1e9
	
	print "f =",freq,"GHz"
	print "Vj =",freq*2.08,"µV"
	return freq*1e9
End


Function ExtractLayer(wave3D)
	wave wave3D
	
	variable NL
	NL=Dimsize(wave3D,2)
	
	variable i
	for(i=0; i<NL; i+=1)
		string Lname = "Lyer"+num2str(i)
		print Lname
		duplicate/R=[][][i,0] wave3D, $Lname
		redimension/N=(-1,-1,0,0) $Lname	
	
	endfor
End

Function xtract1Layer(wave3D,NL)
	wave wave3D
	
	
	variable NL
	//NL=Dimsize(wave3D,2)
	
	//variable i
	//for(i=0; i<NL; i+=1)
	string Lname = NameOfWave(wave3D)+"_L"+num2str(NL)
	print Lname
	duplicate/R=[][][NL][0] wave3D, $Lname
	redimension/N=(-1,-1,0,0) $Lname	
	
	display;appendimage $Lname
	Legend/C/N=text0/J Lname
	
	//endfor
End


Function xtractZoom(wave2D)
	wave wave2D

	string Lname = NameOfWave(wave2D)+"_c"
	//print Lname
	//print min(pcsr(A),pcsr(B))
	duplicate/R=[min(pcsr(A),pcsr(B)),max(pcsr(A),pcsr(B))][min(qcsr(A),qcsr(B)),max(qcsr(A),qcsr(B))] wave2D, $Lname
	
	display;appendimage $Lname
	Legend/C/N=text0/J Lname+"\r Zoom"
	
	//endfor
End


function AreaJJ_v3(Map,zmin)
	wave Map
	variable zmin
	variable dx,dy
	dx= DimDelta(Map,0)
	dy= DimDelta(Map,1)
	
	variable Npxls,N0,N1
	Npxls=0
	N0=Dimsize(Map,0)
	N1=Dimsize(Map,1)

	//string newmap = NameOfWave(Map)+"Area"
	duplicate Map, $(NameOfWave(Map)+"Area")
	wave newmap = $(NameOfWave(Map)+"Area")
	
	variable i,j
	for(i=0;i<N0;i+=1)
		for(j=0;j<N1;j+=1)
			if (Map[i][j] > zmin) 
				Npxls += 1
				newmap[i][j] = map[i][j]
			else	
				newmap[i][j] = nan
			endif	
		endfor
	endfor	
	
	variable A=0.0
	A=Npxls*(dx*dy)
	
	//return A
	print A/1e-12,"µm^2"
	
	display;appendimage Map
	appendimage newmap
	ModifyImage $NameOfWave(newmap) ctab= {*,*,YellowHot,1}
	Legend/C/N=text0/J NameOfWave(Map)+"\r zmin= "+num2str(zmin) +"\r dx= "+num2str(dx)+"\r dy= "+num2str(dy)+"\r Area= "+num2str(A/1e-12)+" µm\S2\M"
	
End

function AreaJJ_v4(Map,zmin)
//measure the area of the JJ in Map (2D-wave picture taken with AFM)
//Analyse the data only within the square, with oposite corners defined by the position of the cursos
//Creates a new map copy only with data within the defined square and plots in color the measured AreaJJ.
//The area of the JJ is defined if the zvalue of the 2D-wave is higher than "zmin"

	wave Map
	variable zmin
	variable dx,dy
	dx= DimDelta(Map,0)
	dy= DimDelta(Map,1)

	variable Npxls,N0,N1
	Npxls=0
	N0=Dimsize(Map,0)
	N1=Dimsize(Map,1)
	
	if (numtype(pcsr(A))==2)
		print "THERE IS NO CURSOR ON TOP GRAPH"
		return 0
	endif
	
	variable Nxmin,Nxmax,Nymin,Nymax
	Nxmin=min(pcsr(A),pcsr(B))
	Nxmax=max(pcsr(A),pcsr(B))
	Nymin=min(qcsr(A),qcsr(B))
	Nymax=max(qcsr(A),qcsr(B))
	
	//string newmap = NameOfWave(Map)+"Area"
	string StrName
	//StrName = NameOfWave(Map)+"Area"
	variable k
	k=0
	StrName = NameOfWave(Map) +"_A"+num2str(k)	
	do
		if (exists(StrName) == 0)
			duplicate Map, $StrName
			wave newmap = $StrName
			break 		
		else
			k=k+1
			StrName = NameOfWave(Map) +"_A"+num2str(k)
		endif
	while (1)
	///R=[min(pcsr(A),pcsr(B)),max(pcsr(A),pcsr(B))][min(qcsr(A),qcsr(B)),max(qcsr(A),qcsr(B))]
	
	variable i,j
	for(i=0;i<N0;i+=1)
		for(j=0;j<N1;j+=1)
				
			if ((i>Nxmin) && (i<Nxmax))
				if ((j>Nymin) && (j<Nymax))
					if (Map[i][j] > zmin) 
						Npxls += 1
						newmap[i][j] = map[i][j]
					else	
						newmap[i][j] = nan
					endif
				else
					newmap[i][j] = nan	
				endif
			else
				newmap[i][j] = nan	
			endif
		
		endfor
	endfor	
	
	variable A=0.0
	A=Npxls*(dx*dy)
	
	//return A
	print A/1e-12,"µm^2"
	
	display;appendimage Map
	appendimage newmap
	ModifyImage $NameOfWave(newmap) ctab= {*,*,YellowHot,1}
	Legend/C/N=text0/J NameOfWave(newmap)+"\r zmin= "+num2str(zmin) +"\r dx= "+num2str(dx)+"\r dy= "+num2str(dy)+"\r Area= "+num2str(A/1e-12)+" µm\S2\M"	
	
End


Function MeanValues(dName1,wavesNum)//,dName2
// It computes the mean value of many waves located in different folders.
// It store each individual result in a new wave of a single element in the corresponding folder.
// The waves {W0001_ij,W0002_ij,W0003_ij} are given for example as: 
// 	dName1='W00', wavesNum={'01','02','03'}
// Finally it plots the set of point Avg_ij vs Avg_Vj
	
	String dName1
	wave/T wavesNum
	//String dName2
	
	Variable numWaves = numpnts(wavesNum)
	Variable i
	for(i=0; i<numWaves; i+=1)	
		PauseUpdate; Silent 1		// building window...
		String fldrSav0= GetDataFolder(1)
		String dfName
		String dName
		sprintf dName, dName1 +"%s" ,wavesNum[i]
		sprintf dfName,"root:Data:%s:",dName
		SetDataFolder dfName
		
		string wVj, wIj, AvgVj, AvgIj
		sprintf wIj, dName1 +"%s"+"_ij",wavesNum[i]
		sprintf wVj, dName1 +"%s"+"_vj",wavesNum[i]
		
		sprintf AvgVj, wVj +"_Avg"	
		sprintf AvgIj, wIj +"_Avg"	
		
		wave dataVj = $wVj
		wave dataIj = $wIj
		make/N=1 $AvgVj=mean(dataVj)
		make/N=1 $AvgIj=mean(dataIj)
				
		if (i == 0)		
			Display $AvgIj vs $AvgVj
		else
			appendtograph $AvgIj vs $AvgVj
		endif
		
		setscale d,0,0,"V",$AvgVj
		setscale d,0,0,"A",$AvgIj
		
		ModifyGraph mode($AvgIj)=3,marker($AvgIj)=28,msize($AvgIj)=5,mrkThick($AvgIj)=1.5		
	endfor
	
	modifygraph grid=1
	ModifyGraph gridRGB=(52428,52428,52428)
	legend
	
	return 0
End


function PlotLinesProfiles(wave_ij,wave_vj,n1,n2,dn,icoil)
//plot a function

	wave wave_ij,wave_vj
	variable n1,n2,dn
	wave icoil
	
	variable k
	for(k=n1;k<n2;k+=dn)
		if(k==n1)
			display wave_ij[k][] vs wave_vj[k][]
		else
			appendtograph wave_ij[k][] vs wave_vj[k][]
		endif
	endfor
	
	k=n1
	string Leg
	Leg = "\\s("+NameOfWave(wave_ij)+") k="+num2str(k)+"; I= "+ num2str(icoil[k][0]/1e-6) +" µA"
	for(k=n1+dn;k<n2;k+=dn)
		Leg=Leg+ "\r"+"\\s("+NameOfWave(wave_ij)+") k="+num2str(k)+"; I= "+ num2str(icoil[k][0]/1e-6) +" µA"
	endfor
	ModifyGraph margin(right)=100
	Legend/C/N=text0/J Leg
	
	return 0
end


Function AddRetrapingToMap(w_up,w_down,extraN,wNewName)
// find the position (i,j) where there is a NAN in the wave w_up
// replace the NaNs, and the following "extraN" points, in w_up by the values
// of w_down at the same positions: (i,j),(i,j+1),(i,j+2),...,(i,j+extraN). 
 
	wave w_up,w_down
	variable extraN
	string wNewName
	
	//string w_updw=Nameofwave(w_up)+"_UpDown"
	duplicate/O w_up $wNewName 
	Wave wUD = $wNewName
	
	variable i,j,n
	for(i=0;i<dimsize(w_up,0);i+=1)
		for(j=0;j<dimsize(w_up,1);j+=1)
			if(numtype(w_up[i][j]) == 2)
				
				wUD[i][j] = w_down[i][j]
				
				for(n=1;n<extraN;n+=1)  // to replace just a section
				//for(n=0;n<(dimsize(w_up,1)-j);n+=1) //to replace up to the end
					//wUD[i][j+n] = 0.5*(w_down[i][j+n]+w_up[i][j+n])
					wUD[i][j+n] = w_down[i][j+n]
				endfor				
				//print "col=",i,"row=",j
			
			endif
		endfor
	endfor
	
	return 0
End


Function AddRetrapingToMap_Max(w_up,w_down,v0,dV,Nh)
// find the position (i,jm) where there is a max current in the wave w_up
// within the range (v0-dv;v0+dv) for each IV.
// replace the w_up[i][j] values by the w_down[i][j] for j>jm+jh
// where jh is the point of the voltage dVh beyond the mawimum of the peak (5 µV range for example)
// all copie to the newwave. preserves the original waves.
 
	wave w_up,w_down
	variable v0,dv,Nh//,dVh
	
	string wNewname=Nameofwave(w_up)+"_UD"
	duplicate/O w_up $wNewName 
	Wave wUD = $wNewName
	
	variable Nx=dimsize(w_up,0)
	variable Ny=dimsize(w_up,1)
	Make/Free/O/N=(Ny) auxI	
	Make/Free/O/N=(Ny) auxV		
	auxV = dimoffset(w_up,1)+p*dimdelta(w_up,1)
	
	variable i,j,P1,P2,jm
	for(i=0;i<Nx;i+=1)
		
		auxI = w_up[i][p]*(abs(auxV[p]-V0)<dv)
		wavestats/Q auxI
		
		//findlevel/Q/P auxV, auxV[V_maxRowLoc]+dVh
		
		for(j=V_maxRowLoc+1;j<Ny;j+=1)
//			wUD[i][j] = w_down[i][j]
			wUD[i][j+Nh] = w_down[i][j+Nh] //if we want start to replace 5 point before the numerical max
												  //This in case there are many points near the max of the peak.
		endfor
	endfor
End


Function AlignPeak(wV,wI,Vmin,Vmax,Voffset)
// the function will find the point position where wI-vs-wV is maximum within the range [Vmin;Vmax] in wV.
// create a duplicated wave from wV which will be aligned (IV by IV, as a flattening)
// All Iv will be aligned such that the maximum find in [Vmin;Vmax] occurs at Voffset
// Basically it requires that a known peak has to have a fix position, at Voffset, for every single IV.
	wave wV,wI
	variable Vmin, Vmax, Voffset
	
	variable Nx= dimsize(wV,0)
	variable Ny= dimsize(wV,1)
	
	string A_wV = NameofWave(wV)+"_A"
	Duplicate/O wV $A_wV
	Wave AwV = $A_wV
	
	Make/O/N=(Ny) auxV
	Make/O/N=(Ny) auxI
	
	variable i,j,Pmin,Pmax
	for(i=0;i<Nx;i+=1)
		//print i
		auxV = wV[i][p]
		auxI = wI[i][p]
		//wAux2 = wV[200][p]
		Findlevel/EDGE=1/P/Q auxV, Vmin ;Pmin=round(V_levelX)
		Findlevel/EDGE=1/P/Q auxV, Vmax ;Pmax=round(V_levelX)
		//print Pmin, Pmax
		wavestats/Q/R=[Pmin,Pmax] auxI	
		
		awV[i][] += -wV[i][V_maxloc] + Voffset	
	endfor
End

Function GetV_when_Iths(wV,wI,wout,Itsh)
//Similar to GetVamx
//Returns a 1D wave with voltages at which the current is equal to Itsh (wV and wI are 2D waves)
//It gets a one point per slide. 
//wout doesn't need to be scaled
	Wave wV,wI,wout
	variable Itsh
	Variable nout = DimSize(wV,0)
	Variable npts = DimSize(wV,1)
	Variable i
	Redimension/N=(nout) wout
	Setscale/P x,DimOffset(wV,0),DimDelta(wV,0),wout
	
	For(i=0;i<nout;i++)
		Make/O/FREE/N=(npts) Ij,Vj
		Vj = wV[i][p]
		Ij = wI[i][p]

		Findlevel/Q/P Ij, Itsh
		//WaveStats/Q Ij
		wout[i]=Vj[round(V_LevelX)]		
	EndFor
End


Function HideRetraping(w_up,wNewName,dVhide,V0,dv,dN)
// dVhide,V0,dv all three negative for retraping.
// find position of max current in w_up in range V0+/-dv for its volatge dimmension
// replace all values by nan in the points on the retraping region.

	wave w_up
	string wNewName
	variable dVhide,v0,dv,dN
	
	duplicate/O w_up $wNewName 
	Wave wHide = $wNewName

	variable Nx=dimsize(w_up,0)
	variable Ny=dimsize(w_up,1)
	Make/Free/O/N=(Ny) auxI	
	Make/Free/O/N=(Ny) auxV		
	auxV = dimoffset(w_up,1)+p*dimdelta(w_up,1)
	
	variable i,j,P1,val
	for(i=0;i<Nx;i+=1)
		
		auxI = w_up[i][p]*(abs(auxV[p]-V0)<abs(dv))
		wavestats/Q auxI
		val=auxV[V_maxRowLoc]+dVhide
		Findlevel/Q/P auxV,val
		P1=V_LevelX
		
		for(j=P1;j<V_maxRowLoc-dN;j+=1) //plots
			wHide[i][j]=NaN
		endfor
	
	endfor
	
	return 0
End


Function HideSwitching(wI,wV,wNewName,V0,dv,dN)
// dVhide,V0,dv all three negative for retraping.
// find position of max current in w_up in range V0+/-dv for its volatge dimmension
// replace all values by nan in the points on the switching region.

	wave wI,wV
	string wNewName
	variable v0,dv,dN

	duplicate/O wI $wNewName 
	Wave wHide = $wNewName

	variable Nx=dimsize(wI,0)
	variable Ny=dimsize(wI,1)
	Make/Free/O/N=(Ny) auxI	
	Make/Free/O/N=(Ny) auxV		
	//auxV = dimoffset(w_up,1)+p*dimdelta(w_up,1)

	variable i,j,P1,val
	for(i=0;i<Nx;i+=1)
		auxV = wV[i][p]
		auxI = wI[i][p]*(abs(auxV[p]-V0)<abs(dv))
		wavestats/Q/P auxI
		//val=auxV[V_maxRowLoc]+dVhide
		//Findlevel/Q/P auxV,val
		//P1=V_LevelX

		for(j=V_maxRowLoc;j<V_maxRowLoc+dN;j+=1)
			wHide[i][j]=NaN
		endfor

	endfor

	return 0
End


Function AlignAtMaxV(wV,Voffset)
// find the max Voltage. it should be at the bottom of the gap, before the oscillations
// Align all IV at this point to be set at the value Voffset
// for the Gap, Voffset=396 µV for example, but it depends on sample.
  
	wave wV
	variable Voffset
	
	variable Nx= dimsize(wV,0)
	variable Ny= dimsize(wV,1)
	
	string A_wV = NameofWave(wV)+"_A"
	Duplicate/O wV $A_wV
	Wave AwV = $A_wV
	
	Make/O/N=(Ny) auxV
	
	variable i,j
	for(i=0;i<Nx;i+=1)
	
		auxV = wV[i][p]	
		wavestats/Q auxV
		
		awV[i][] += -wV[i][V_maxRowLoc] + Voffset	
	endfor
End


Function CutRangeVJ_v0(waveI,waveV,Vi,Vf)
//2020.03.15.RR: This function create a new set of waves
//copied from waveI and waveV only in the voltage range [Vi,Vf] 

	wave waveI,waveV
	variable Vi,Vf
	
	variable Ni,Nf
	Findlevel/Q/P/R=[0,10010] waveV,Vi
	Ni=V_levelX
	Findlevel/Q/P/R=[0,10010] waveV,Vf
	Nf=V_levelX
		
	string Xwi = NameofWave(waveI)+"_X"
	Duplicate/O/R=[Ni,Nf] waveI $Xwi
	string Xwv = NameofWave(waveV)+"_X"
	Duplicate/O/R=[Ni,Nf] waveV $Xwv

End


Function CutRangeVJ_v1(waveI,waveV,Vi,Vf,RetSlope) //,outI,outV
//2020.03.18.RR: This function creates a new set of wave for an I-V curve
//copied from waveI and waveV only in the voltage range [Vi,Vf]
//It combines switching (from positive Vj) with the retrapping (at negative Vj)
//Uses the value of the slope of the load line (RetSlope) in rettraping to skip the hidden region. 
//Select the case with minus sign if retraping is at negative voltages.
	wave waveI,waveV
	variable Vi,Vf,RetSlope

//find the point number of voltages values Vi, Vf.	
	variable Np,Ni_up,Nf_up
	Np=Dimsize(waveI,0)
	Findlevel/Q/P/R=[0,Np] waveV,Vi
	Ni_up=V_levelX
	Findlevel/Q/P/R=[0,Np] waveV,Vf
	Nf_up=V_levelX
	Wavestats/Q/R=[Ni_up,Nf_up] waveI
	Nf_up=V_maxRowLoc

//copie the waves, invert the I-V to work in rettraping as it was a switching event.
	Duplicate/FREE waveI, wI_ret
	Duplicate/FREE waveV, wV_ret

///-------- TO SELECT THE CASE -----------------
	//to invert the sense of the sweept if retraping if Up  & Down are both at positive voltages
		//wI_ret = waveI[20019-p]	
		//wV_ret = waveV[20019-p]	
	//to invert the sense of the sweept is retraping if Up  & Down are at oposite sides of the IV.
	wI_ret = -waveI[Np-1-p]	
	wV_ret = -waveV[Np-1-p]	
	//search for the max I_value
	variable Ni_ret,Nf_ret
	Findlevel/Q/P/R=[0,Np] wV_ret,Vi
	Ni_ret=V_levelX
	Findlevel/Q/P/R=[0,Np] wV_ret,Vf
	Nf_ret=V_levelX
	Wavestats/Q/R=[Ni_ret,Nf_ret] wI_ret
	
	//over-estimated the Landing voltage after small-switching
	Variable Imax,Vmax,Vlanding
	Imax = V_max
	Vmax = wV_ret[V_maxRowLoc]
	Vlanding = Vmax + (-Imax)/RetSlope
	
	//findlevel the point of Landing_Vj (this is Ni_ret)
	Findlevel/Q/P/R=[0,Np] wV_ret,Vlanding
	Ni_ret=V_levelX
	
// duplicate waveI,waveV in range Ni_up, Nf_up
	Duplicate/FREE/R=[Ni_up,Nf_up] waveI wI_upX
	Duplicate/FREE/R=[Ni_up,Nf_up] waveV wV_upX
// duplicate wI_ret,wV_ret in range Ni_ret, Nf_ret
	Duplicate/FREE/R=[Ni_ret,Nf_ret] wI_ret wI_retX
	Duplicate/FREE/R=[Ni_ret,Nf_ret] wV_ret wV_retX	
// concatenate both parts.
	string UDwI = NameofWave(waveI)+"_UD"
	Concatenate/O {wI_upX,wI_retX}, $UDwI 
	string UDwV = NameofWave(waveV)+"_UD"
	Concatenate/O {wV_upX,wV_retX}, $UDwV
	
End
