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