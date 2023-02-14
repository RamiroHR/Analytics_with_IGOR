#pragma TextEncoding = "Windows-1252"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// PSD calculate Power Spectral Density
Function PSD(w, segment_size)
	WAVE w // real wave to calculate PSD
	Variable segment_size // size for each segment when calculating the FFT
						// preferably a power of 2
	
	DSPPeriodogram /NOR=((segment_size)^2/2) /SEGN={(segment_size),(segment_size/2)}  w
	
	WAVE w_psd = W_Periodogram

	// Fix DC component
	w_psd[0] /= 2
	// Fix Nyquist component
	w_psd[numpnts(w_psd)-1] /= 2
	
	// Change normalization:
	//  Default is that sum(w_psd') == rms(w)^2
	// For spectral density, want
	//  sum(w_psd)*deltax(w_psd) == rms(w)^2
	
	w_psd /= deltax(w_psd)
	
	// Fix units
	SetScale d,0,0,WaveUnits(w,-1)+"\S2\M/Hz", w_psd
	
	Duplicate/O w_psd, $(NameOfWave(w)+"_psd")
	
	killwaves w_psd //remove useless w_periodogram
End


//function that returns rms value of an input wave
function rms(w)
	wave w
	duplicate w wtmp
	wtmp=wtmp^2
	variable rms=sum(wtmp)
	rms = rms/numpnts(wtmp)
	killwaves wtmp
	return sqrt(rms)
end



//use this to compute the psd, the rms, and the uncorrelated data
// name_pref = "BF14_"
// name_list = wave of string like ("0001, 0003, 0005, ...")
// name_suf = "C1_"
// common folder of measured data is root:data:$(name_pref+name_list[i])
// then each measurement is in $(name_pref+'name_list[i]+name_suf+"j"
//    full path example : root:Data:BF14_0007:BF14_0007_C1_2 
//       where i=0007 and j=2
Function noise12(name_pref, name_list, name_suf, name_out_suf)
	string name_pref, name_suf, name_out_suf
	wave/t name_list

	// make a new folder that will contain the psd and rms and uncorrelated arrays
	string output_name=name_pref+name_suf+name_out_suf
	newdatafolder/o root:$output_name
	print "create folder : root:"+output_name
	
	//go to this folder
	setdatafolder $"root:"+output_name
	
	//initiate the arrays 	
	make/n=(65537,12) $("Psd_map")
	wave psd_map= $("Psd_map")
	
	make/n=12 $("Rms_map")
	wave rms_map= $("Rms_map")
	
	make/n=(65537,6) $("Psd_uncorr_map")
	wave psd_uncorr_map=$("Psd_uncorr_map")
	
	variable i,j //i is used to select the appropriate data folder. j is used to select both dataset in the folder
	for(i=0; i<6; i+=1)
		string input_dir = "root:data:"+name_pref+name_list[i] //adress of the folder to go to
		print("looking at "+input_dir)
		for(j=1; j<=2; j+=1)
			//define which wave will be analysed here
			string input_wave = input_dir+":"+name_pref+name_list[i]+"_"+name_suf+num2str(2*i+j)
			print "treating data : "+input_wave
			// create reference to this wave
			wave wtmp = $input_wave

			//then do psd and add it to the psd map
			psd(wtmp,2^17)
			wave wtmp_psd = $(input_wave+"_psd_")
			psd_map[][2*i+j-1] = wtmp_psd
			
			//now compute rms, and add it to rms map
			rms_map[][2*i+j-1] = rms(wtmp)
		endfor
		//now compute the uncorellated noise of the two data j=1 and j=2
		// compute difference of the two waves
		string input_wave1 = input_dir+":"+name_pref+name_list[i]+"_"+name_suf+num2str(2*i+1)
		wave wtmp_diff = $input_wave1
		wtmp_diff -= wtmp
		//take psd of it, then store it in the map
		psd(wtmp_diff, 2^15)
		wave wtmp_diff_psd = $(input_wave1+"_psd_")
		psd_uncorr_map[][i]=wtmp_diff_psd		
	endfor
	print "fini"
End

// what still needs to be done :
//  - make a version of noise12 for quad (different file name, and just one file per folder)
//  - setup a graph style or something to automatically plot the psd maps
//  - plot the rms automatically too



function testnoise()
	string name_pref = "BF14_"
	make/n=6/o/t name_list
	name_list= {"0001","0003","0005","0007","0009","0011"}
	string name_suf="C1_"
	string name_out_suf="pt_on"
	
	noise12(name_pref, name_list, name_suf, name_out_suf)
end



//////////////////////////////////////////////////////////////////
// New version of noise12, which uses name_data and name_folder
//    name_data :
//    name_folder :
//////////////////////////////////////////////////////////////////

//use this to compute the psd, the rms, and the uncorrelated data
// name_pref = "N1"
//
// name_folder = wave of string for the different folder name
// name_data = wave of string for the datas contained in each folder (2 per folder)
//  so name_folder is like {"0001","0002" ...
//   and name_data is like {"1","2","5","6", ...
//
// name_suf = "K"
// common folder of measured data is root:data:$(name_pref+name_list[i])
// then each measurement is in $(name_pref+'name_folder[i]+name_suf+name_data[2*i+j]
//    full path example : root:Data:N10001:N10001_K7
//       where i=0001 and j=2
//
//exemplo : noise10("N3",name_folder,name_data,"","PT ON")
//
Function noise10(name_pref, name_folder,name_data, name_suf, name_out_suf)
	string name_pref, name_suf, name_out_suf
	wave/t name_folder, name_data

	// make a new folder that will contain the psd and rms and uncorrelated arrays
	string output_name=name_pref+name_suf+name_out_suf
	newdatafolder/o root:$output_name
	print "create folder : root:"+output_name
	
	//go to this folder
	setdatafolder $"root:"+output_name
	
	//initiate the arrays 	
	make/n=(16385,10) $("Psd_map")
	wave psd_map= $("Psd_map")
	
	make/n=10 $("Rms_map")
	wave rms_map= $("Rms_map")
	
//	make/n=(65537,6) $("Psd_uncorr_map")
//	wave psd_uncorr_map=$("Psd_uncorr_map")
	
	variable i,j //i is used to select the appropriate data folder. j is used to select both dataset in the folder
	for(i=0; i<5; i+=1)
		string input_dir = "root:data:"+name_pref+name_folder[i] //adress of the folder to go to
		print("looking at "+input_dir)
		for(j=0; j<=1; j+=1)
			//define which wave will be analysed here
			string input_wave = input_dir+":"+name_pref+name_folder[i]+"_"+name_suf+name_data[2*i+j]
			print "treating data : "+input_wave
			// create reference to this wave
			wave wtmp = $input_wave

			//then do psd and add it to the psd map
			psd(wtmp,2^15)
			wave wtmp_psd = $(name_pref+name_folder[i]+"_"+name_suf+name_data[2*i+j]+"_psd")
			print(input_wave+"_psd")
			psd_map[][2*i+j] = wtmp_psd[p]
			
			//now compute rms, and add it to rms map
			print(rms(wtmp))
			rms_map[2*i+j] = rms(wtmp)
		endfor
		//now compute the uncorellated noise of the two data j=1 and j=2
		// compute difference of the two waves
//		string input_wave1 = input_dir+":"+name_pref+name_list[i]+"_"+name_suf+num2str(2*i+1)
//		wave wtmp_diff = $input_wave1
//		wtmp_diff -= wtmp
		//take psd of it, then store it in the map
//		psd(wtmp_diff, 2^17)
//		wave wtmp_diff_psd = $(input_wave1+"_psd_")
//		psd_uncorr_map[][i]=wtmp_diff_psd		
	endfor
	print "fini"
End


function testnoise_pton()
	string name_pref = "N2"
	make/n=5/o/t name_folder
	name_folder= {"0001","0003","0005","0007","0009"}
	string name_suf=""
	make/n=10/o/t name_data
	name_data= {"3","4","5","6","7","8","9","10","11","12"}
	string name_out_suf="PT_On"
	
	noise10(name_pref, name_folder,name_data, name_suf, name_out_suf)
end

function testnoise_ptoff()
	string name_pref = "N2"
	make/n=5/o/t name_folder
	name_folder= {"0002","0004","0006","0008","0010"}
	string name_suf=""
	make/n=10/o/t name_data
	name_data= {"3","4","5","6","7","8","9","10","11","12"}
	string name_out_suf="PT_Off"
	
	noise10(name_pref, name_folder,name_data, name_suf, name_out_suf)
end


//command to make psd map pretty :
function plot_psdmap()
	ModifyImage psd_map_50K_off ctab= {*,*,Grays,1},log=1,ctabAutoscale=1,lookup= $""
	ModifyGraph log(bottom)=1
	SetAxis bottom 10,*
end