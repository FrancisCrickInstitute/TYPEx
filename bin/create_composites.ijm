#@String panel
#@String rootDir
#@String compositeDir
#@String tumour
#@String stroma
#@String immune
#@string auxStroma
#@String dna
#@Integer nrCategs


minImgSize = 100;
// Choose whether to analyse all images ("all") or specific e.g. "P1_TMA003_R_20190619-roi_12", or "P1_TMA004_L_20190619-roi_13"
image="all";
// Choose whether to apply directional filter ("direct") or "none"
directional='none';

// Choose the NextFlow run iteration
panck_threshold = 15;

// DNA
dnaList = split(dna, "=|,");
dnaMarkers = "";
for (m = 1; m < dnaList.length; m++) {
	dnaMarkers = dnaMarkers + "image" + m + "=[Median of " + dnaList[m] + "] ";
}
dnaMarkers = dnaMarkers + "image" + dnaList.length + "=[-- None --]";
print(dnaMarkers)

// Tumour
tumourList = split(tumour, "=|,");
tumourMarkers = "";
for (m = 1; m < tumourList.length; m++) {
    tumourMarkers = tumourMarkers + "image" + m + "=[Median of " + tumourList[m] + "] ";
}
tumourMarkers = tumourMarkers + "image" + tumourList.length + "=[-- None --]";
print(tumourMarkers)

// Immune
immuneList = split(immune, "=|,");
immuneMarkers = "";
for (m = 1; m < immuneList.length; m++) {
    immuneMarkers = immuneMarkers + "image" + m + "=[Median of " + immuneList[m] + "] ";
}
immuneMarkers = immuneMarkers + "image" + immuneList.length + "=[-- None --]";
print(immuneMarkers)

// Stroma
stromaList = split(stroma, "=|,");
stromaMarkers = "";
for (m = 1; m < stromaList.length; m++) {
	  print(stromaList[m]);
    stromaMarkers = stromaMarkers + "image" + m + "=[Median of " + stromaList[m] + "] ";
}

stromaMrgMarkers = stromaMarkers + "image" + stromaList.length + "=auxStromaMrg image" +
    ( stromaList.length + 1 ) + "=[-- None --]";
stromaMarkers = stromaMarkers + "image" + (stromaList.length) + "=[-- None --]";
print(stromaMarkers)


// AUX Stroma
auxStromaList = split(auxStroma, "=|,");
auxStromaMarkers = "";
if(auxStromaList.length > 1) {
	for (m = 1; m < auxStromaList.length; m++) {
    	auxStromaMarkers = auxStromaMarkers + "image" + m + "=[Median of " + auxStromaList[m] + "] ";
	}
	auxStromaMarkers = auxStromaMarkers + "image" + auxStromaList.length + "=[-- None --]";
	print(auxStromaMarkers);
}
nonstroma='image1=Tumour image2=SUM_DNA image3=[-- None --]';

setBatchMode(true);
runs=getFileList(rootDir);

print("Number of images found in " + rootDir + " " + runs.length);
for (k = 0; k< runs.length;  k++)	{
	// print('RUNL ', runs[k]);

	runDir=rootDir + runs[k];
	roiList=getFileList(runDir);
	print('Run dir', runDir, roiList.length);
		
	for (j = 0; j < roiList.length; j++) {

		tma=replace(runs[k],  "/", "");
		roi=replace(roiList[j], "/", "");
		
		print(tma);
		imgName=tma + "-" + roi;
		print('IMG name:  ', imgName);
		fileOut = tma + "-" + roi + ".tiff";
		print(fileOut);
		if(File.exists(fileOut)) {
			print('File exists', fileOut);
			continue;	
		}
		if(image != "all" && imgName != image) {
			print(image, 'skipping');
			continue;
		}	
		stacks=getFileList(runDir+roiList[j]);
		print("Stack dir: " + runDir + roiList[j]);
		if(stacks.length == 0) {
			print('Stack size is 0. Skipping');
			continue;
		}
		for (m = 0; m < stacks.length; m++)	{
	
			if(! endsWith(stacks[m], "full_stack/")) 
				continue;
			
			print('Stack: ', stacks[m]);
			
			imgDir=runDir + roiList[j] + stacks[m];
			if(stacks.length == 0) {
				print("ERROR: empty stack " + stacks[m] + "\n");
				exit;
			}
			if(File.isDirectory(imgDir)) {
			
				imgList=getFileList(imgDir);

				for ( i=0; i<imgList.length; i++ ) {
					
					if(indexOf(tumour, imgList[i]) == -1 &&
						indexOf(immune, imgList[i]) == -1 &&
						indexOf(stroma, imgList[i]) == -1 &&
						indexOf(auxStroma, imgList[i]) == -1 &&
						indexOf(dna, imgList[i]) == -1) {
							continue;
						}
						
					print("Opening:" + imgList[i]);		
					open(imgDir + imgList[i]);
					// very small images willl not work with Ilastik
					width = getWidth;
					height = getHeight;
					print(imgList[i], width, height, minImgSize);
					if(width < minImgSize || height < minImgSize){
							run("Close All");
							continue;
					}
					isTumourMarker = "false";
					for(m = 0; m < tumourList.length; m++) {
                        if(imgList[i] == tumourList[m]) {
							isTumourMarker = 'true';
							break;
						}
					}
					if(isTumourMarker == "false") {
						run("Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=3");
						run("Median (3D)");
					}
				  }
				  autoAdjust();
				}
			}

			if(isOpen(tumourList[1]) == 0)
				continue;
			print(imgName);

			if(width < minImgSize || height < minImgSize) {
							run("Close All");
							continue;
			}
			// DNA Sum
			run("Concatenate...", " title=DNA open " + dnaMarkers);
			run("Z Project...", "projection=[Sum Slices]");
			run("Remove Outliers", "block_radius_x=2 block_radius_y=2 standard_deviations=3");
			autoAdjust();
			run("Enhance Contrast", "saturated=0.35");
		//	saveAs("PNG",tma + "-" + roi +  "dna.png");

			print(tumourList.length);
			// Tumour
			print("Tumour markers");
			for (tind = 1; tind < tumourList.length; tind++) {

			    print(tumourList[tind]);
				selectWindow(tumourList[tind]);
				getMinAndMax(min, max);
				print(tumourList[tind]);

				if(max < panck_threshold) {

					run("Morphological Filters", "operation=[Dilation] element=Square radius=1");
					run("Directional Filtering", "type=Max operation=Median line=2 direction=32");
					run("Morphological Filters", "operation=[Erosion] element=Square radius=1");
				} else {

					run("Duplicate...", "title=TumourDupl");
					run("Morphological Filters", "operation=[Dilation] element=Square radius=20");
					rename("ProcessNoise");
					run("Calculator Plus", "i1=TumourDupl i2=ProcessNoise operation=[Divide: i2 = (i1/i2) x k1 + k2] k1=1 k2=0 create");
					close("ProcessNoise");
					close("TumourDupl");
					//run("Directional Filtering", "type=Max operation=Median line=1 direction=32");
				}
				close(tumourList[tind]);
				run("Median (3D)");
				rename("Median of " + tumourList[tind]);
			}

			if(tumourList.length > 2) {

				print("Concatenate tumour");
				print(tumourMarkers);
				run("Concatenate...", " title=Test open " + tumourMarkers);
				run("Z Project...", "projection=[Sum Slices]");
				rename('Tumour');
				autoAdjust();
				run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
			} else {
				selectWindow("Median of " + tumourList[1]);
				autoAdjust();
				print("Enhance");
				run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
				rename('Tumour');
			}
			// Stroma
			if(auxStromaMarkers != "") {
			
				// Concatenation works for more than one marker
				if(auxStromaList.length > 2) {
					run("Concatenate...", "  title=Nonsp open " + auxStromaMarkers);
					run("Z Project...", "projection=[Sum Slices]");
				} else {
					print(auxStromaList[1]);
					
					if(isOpen(auxStromaList[1]) == 0) {
							print(auxStromaList[1], "windows is not open");
							run("Close All");
							continue;
					}
					print(auxStromaList[1], "window is open");
					selectWindow(auxStromaList[1]);
				}
				
				run("Remove Outliers", "block_radius_x=2 block_radius_y=2 standard_deviations=3");
				autoAdjust();
				run("Median (3D)");
				autoAdjust();
				run("Enhance Contrast", "saturated=0.35");
				rename('auxStroma');
			}
		
			// Immune sum
			// Concatenation works for more than one marker
			print("Concatenate immune");
			if(immuneList.length > 2) {
				run("Concatenate...", "  title=Immune open " + immuneMarkers);
				run("Z Project...", "projection=[Sum Slices]");
			} else {
				selectWindow(immuneList[1]);
			}
			run("Remove Outliers", "block_radius_x=2 block_radius_y=2 standard_deviations=3");
			print('Immune');
			autoAdjust();
			run("Enhance Contrast", "saturated=0.35");
           
			run("Concatenate...", "  title=NonStroma keep open " + nonstroma); 
			run("Z Project...", "projection=[Sum Slices]");
			print(stromaMarkers);
			if(auxStromaMarkers != "") {
				// Max to generate a tumour mask - considering nuclei are not in it
				run("Median (3D)");
				setAutoThreshold("Triangle dark");
				run("Convert to Mask");
				run("Morphological Filters", "operation=Closing element=Square radius=1");
				run("Invert");
				run("32-bit");
				run("Calculator Plus", "i1=Median-Closing i2=Median-Closing operation=[Divide: i2 = (i1/i2) x k1 + k2] k1=1 k2=0");
				run("Calculator Plus", "i1=auxStroma i2=Median-Closing operation=[Multiply: i2 = (i1*i2) x k1 + k2] k1=1 k2=0 create");
				selectWindow("Result");
				run("Median (3D)");
				rename('auxStromaMrg');
				autoAdjust();
				print(stromaMrgMarkers);
				run("Concatenate...", "  title=Stroma open " + stromaMrgMarkers);
				run("Z Project...", "projection=[Sum Slices]");
			} else if(stromaList.length > 1)	{
				if(stromaList.length > 2) {
				 run("Concatenate...", "  title=Stroma open " + stromaMarkers);
	 				run("Z Project...", "projection=[Sum Slices]");
				} else {
					print(stromaList[1]);
					selectWindow(stromaList[1]);
				}
			}

			run("Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=3");
			run("Median (3D)");
			autoAdjust();
			rename('Stroma_Only');
			// saveAs("PNG", tma + "-" + roi + "stroma.png");

			run("Concatenate...", "  title=StromaStack open image1=SUM_Immune image2=[Stroma_Only] image3=[-- None --]");
			run("Z Project...", "projection=[Sum Slices]");
			rename("Stroma_merge");
			autoAdjust();
//			saveAs("PNG", tma + "-" + roi + "stroma_merg.png");

			run("Enhance Contrast", "saturated=0.35");
			run("Merge Channels...", "c3=SUM_DNA c7=Tumour c6=Stroma_merge create");
			run("RGB Color");
			saveAs("Tiff", fileOut);
			print(fileOut);
			run("Close All");

	}
}

function autoAdjust() {
	// creating a new function - to include auto brightness and contrast in macro
	
        /*
         * rewriting “Auto contrast adjustment” button of “Brightness/Contrast”
         * Damien Guimond
         * 20120516
         * Acknowledgements: Kota Miura
        */
        /*
         * Edits by 
         * Mihaela Angelova
        */
          
		  aggregateMax=5000;
		  aggregateMin=15;
        noiseLevel=3;
        AUTO_THRESHOLD = 5000;
        pctPixels=10;

		 getRawStatistics(pixcount);
		 limit = pixcount/pctPixels;
		 threshold = pixcount/AUTO_THRESHOLD;
		 nBins = 256; 
		
		 getHistogram(values, histA, nBins);
		
		 i = -1;
		 found = false;
		 do {
		    counts = histA[++i];
		    if (counts > limit) counts = 0;
		    found = counts > threshold;
		 
		 } while ((!found) && (i < histA.length-1))
		 hmin = values[i];
		 i = histA.length;
		 do {
			counts = histA[--i];
		    if (counts > limit) counts = 0;
		    found = counts > threshold;
		 } while ((!found) && (i > 0))
		 hmax = values[i];
		 // the maximum value
		
		if(hmax > aggregateMax) { // || hmax == 0) {
			hmin=noiseLevel;
			hmax=15;
		}
		if(hmin > aggregateMin) hmin=noiseLevel;
		if(hmin < 0) hmin=0;
		if (hmax > noiseLevel) setMinAndMax(hmin, hmax); 
		 print(hmin, hmax);
		// run(“Apply LUT”);
}
