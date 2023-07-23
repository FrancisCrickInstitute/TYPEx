#@String panel
#@String rootDir
#@String compositeDir
#@String tumour
#@String stroma
#@String immune
#@string auxStroma
#@String dna
#@Integer nrCategs


// Choose whether to analyse all images ("all") or specific e.g. "P1_TMA003_R_20190619-roi_12", or "P1_TMA004_L_20190619-roi_13"
image="all";
// Choose whether to apply directional filter ("direct") or "none"
directional='none';

// Choose the NextFlow run iteration
f = File.open('panck_preprocessed_images.txt');
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
    stromaMarkers = stromaMarkers + "image" + m + "=[Median of " + stromaList[m] + "] ";
}
stromaMarkers = stromaMarkers + "image" + stromaList.length + "=auxStromaMrg image" + 
	( stromaList.length + 1 ) + "=[-- None --]";
print(stromaMarkers)


// AUX Stroma
auxStromaList = split(auxStroma, "=|,");
auxStromaMarkers = "";
if(auxStromaList.length > 1) {
	for (m = 1; m < auxStromaList.length; m++) {
    	auxStromaMarkers = auxStromaMarkers + "image" + m + "=[Median of " + auxStromaList[m] + "] ";
	}
	auxStromaMarkers = auxStromaMarkers + "image" + auxStromaList.length + "=[-- None --]";
	print(auxStromaMarkers)
}
nonstroma='image1=[Tumour,image2=SUM_DNA image3=[-- None --]';

setBatchMode(true);
runs=getFileList(rootDir);

print("Number of images found in " + rootDir + " " + runs.length);
for (k = 0; k< runs.length;  k++)	{
	print('RUNL ', runs[k]);

	runDir=rootDir + runs[k];
	roiList=getFileList(runDir);
	print('Run dir', runDir, roiList.length);
		
	for (j = 0; j < roiList.length; j++) {

		tma=replace(runs[k],  "/", "");
		roi=replace(roiList[j], "/", "");

		imgName=tma + "-" + roi;
		print('IMG name:  ', imgName);
		if(image != "all" && imgName != image) {
			print(image, 'skipping');
			continue;
		}	
		stacks=getFileList(runDir+roiList[j]);
		print("Stack dir: " + runDir + roiList[j]);
		if(stacks.length==0) {
			print('Stack size is 0. Skipping');
			continue;
		}
		newPanel=0;
		for (m = 0; m < stacks.length; m++)	{
	
			if(! endsWith(stacks[m], "full_stack/")) continue;
			
			print('Stack: ', stacks[m]);
			
			imgDir=runDir + roiList[j] + stacks[m];
			if(stacks.length==0) {
				print("ERROR: empty stack " + stacks[m] + "\n");
				exit;
			}
			if(stacks[m] != "full_stack/")
				continue;
			fileOut = tma + "-" + roi + ".tiff";
		    if(File.exists(fileOut)) {
				print('File exists', fileOut);
				 continue;	
			}
			if(File.isDirectory(imgDir)) {
				imgList=getFileList(imgDir);

				for ( i=0; i<imgList.length; i++ ) {
					if(imgList[i] == '131Xe.tiff' || imgList[i] == '80ArAr.tiff' || imgList[i] == '100Ru_ruthenium.tiff' || 
						imgList[i] == '100Ru.tiff' || imgList[i] == '134Xe.tiff' || 
						imgList[i] == '141Pr_CD38.tiff' || imgList[i] == '143Nd_MCT4.tiff' ||
						imgList[i] == '148Nd_KIR2DL3.tiff' || imgList[i] == '150Nd_PDL1.tiff' ||
						imgList[i] == '163Dy_CD103.tiff' || imgList[i] == '173Yb_MHCII.tiff' || 
						imgList[i] ==  '153Eu_LAG3.tiff' || imgList[i] == '172Yb_CD56.tiff' ||
						imgList[i] == '154Sm_TIM3.tiff' || imgList[i] == '155Gd_IDO.tiff' ||
						imgList[i] == '176Yb_CAIX.tiff'  || imgList[i] == '165Ho_PD1.tiff' ||
						imgList[i] == '166Er_CLEC9a.tiff' || imgList[i] == '167Er_GZMB.tiff' || 
						imgList[i] == '160Gd_VISTA.tiff' || imgList[i] == '168Er_CD73.tiff' ||
						imgList[i] == '168Er_Ki67.tiff' || imgList[i] == '176Yb_TCF1.tiff' || 
						imgList[i] == '173Yb_B2M.tiff' || imgList[i] == '174Yb_pSTAT1.tiff' || 
						imgList[i] == '175Lu_CD25.tiff' || imgList[i] == '166Er_CD45RA.tiff' ||
						imgList[i] == '158Gd_CXCL12.tiff' || imgList[i] == '159Tb_CXCR4.tiff' ||
						imgList[i] == '161Dy_CD39.tiff' || imgList[i] == '160Gd_GITR.tiff' ||
						imgList[i] == '155Gd_FOXP3.tiff' || imgList[i] == '149Sm_GATA3.tiff' ||
						imgList[i] == '144Sm_CD57.tiff' || imgList[i] == '171Yb_CD27.tiff' ||
						imgList[i] == '142Nd_CCR7.tiff' || imgList[i] == '172Yb_CASP3.tiff' ||
						imgList[i] == '145Nd_CTLA4.tiff' || imgList[i] == '146Nd_FAP1.tiff' ||
						imgList[i] == '147Sm_CXCR6.tiff' || imgList[i] == '148Nd_ICOS.tiff')
						
						continue;
						
					open( imgDir + imgList[i] );
				   	if(imgList[i] == '142Nd_CAM52.tiff' || imgList[i] == '142Nd_CAM52Nd142Di.tiff')
						newPanel=1;
					print(imgList[i]);
					autoAdjust();

					if(imgList[i] != '164Dy_panCK.tiff') { 
						run("Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=3");
						run("Median (3D)");
					}
				}
			}
			print(newPanel);
			// DNA Sum
			run("Concatenate...", " title=DNA open " + dnaMarkers);
			run("Z Project...", "projection=[Sum Slices]");
			run("Remove Outliers", "block_radius_x=2 block_radius_y=2 standard_deviations=3");
			autoAdjust();
			run("Enhance Contrast", "saturated=0.35");

			if(newPanel) {
				run("Concatenate...", " title=Test open " + tumourMarkers);
				run("Z Project...", "projection=[Sum Slices]");
				rename('Tumour');
				autoAdjust();
			} else {
				selectWindow("164Dy_panCK.tiff");
				autoAdjust();
				getMinAndMax(min, max);
				print(imgName, 'Tumour', min, max);
				if(max < panck_threshold) {
					open(compositeDir + 'panckf_' + imgName + '.tif');
					print('panckf_' + imgName  + '.tif');
				} else {
					run("Median (3D)");
				}
				run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
				singleMarkerImg="Tumour";
				rename(singleMarkerImg);
			}

			// Stroma
			if(panel == 'p1') {
				run("Concatenate...", "  title=Nonsp open " + auxStromaMarkers);
				run("Z Project...", "projection=[Sum Slices]");
				run("Remove Outliers", "block_radius_x=2 block_radius_y=2 standard_deviations=3");
				autoAdjust();
				run("Median (3D)");
				autoAdjust();
				run("Enhance Contrast", "saturated=0.35");
				rename('auxStroma');
				print('Aux Stroma');
			} else {
				selectWindow(auxStroma);
				rename('auxStroma');
			}
		
			// Immune sum
			run("Concatenate...", "  title=Immune open " + immuneMarkers);
			run("Z Project...", "projection=[Sum Slices]");
			run("Remove Outliers", "block_radius_x=2 block_radius_y=2 standard_deviations=3");
			print('Immune');
			autoAdjust();
			run("Enhance Contrast", "saturated=0.35");
           
			run("Concatenate...", "  title=NonStroma keep open " + nonstroma); 
			run("Z Project...", "projection=[Sum Slices]");

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
			print(stromaMarkers);

            run("Concatenate...", "  title=Stroma open " + stromaMarkers);
            run("Z Project...", "projection=[Sum Slices]");

			run("Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=3");
			run("Median (3D)");
			autoAdjust();
			rename('Stroma_Only');

			run("Concatenate...", "  title=StromaStack open image1=SUM_Immune image2=[Stroma_Only] image3=[-- None --]");
			run("Z Project...", "projection=[Sum Slices]");
			rename("Stroma_merge");
			autoAdjust();

			run("Enhance Contrast", "saturated=0.35");
			run("Merge Channels...", "c3=SUM_DNA c7=Tumour c6=Stroma_merge create");
			run("RGB Color");
			saveAs("Tiff", fileOut);
			print(fileOut);
			run("Close All");

		}
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


