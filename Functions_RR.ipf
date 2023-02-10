#pragma TextEncoding = "Windows-1252"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

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


//
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
	duplicate/R=[][][NL,0] wave3D, $Lname
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
	Legend/C/N=text0/J Lname+"_zoom"
	
	//endfor
End


function AreaJJ_v3(Map,zmin,dx)
	wave Map
	variable zmin,dx


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
	
	variable A
	A=Npxls*(dx)^2
	
	//return A
	print A/1e-12,"µm^2"
	
	display;appendimage Map
	appendimage newmap
	ModifyImage $NameOfWave(newmap) ctab= {*,*,YellowHot,1}
	Legend/C/N=text0/J NameOfWave(Map)+"\r zmin= "+num2str(zmin)+"\r Area= "+num2str(A/1e-12)+" µm\S2\M" 
	
End
