#@String panel
#@String runPattern
#@Boolean mcd
#@String rootDir
#@String outDir

// Choose whether to analyse all images ("all") or specific e.g. "P1_TMA003_R_20190619-roi_12", or "P1_TMA004_L_20190619-roi_13"
image="all";
// Choose whether to apply directional filter ("direct") or "none"
directional='none';
// Choose the NextFlow run iteration

panck_thres=150

print(rootDir);
if(File.isDirectory(outDir) == 0)
	File.makeDirectory(outDir);

setBatchMode(true);
runs=getFileList(rootDir);
print(runs.length);
for (k=0; k < runs.length; k++)	{
	print(runs[k]);
	if(startsWith(runs[k], "run")) {
		runDir=rootDir + runs[k] + "results/imctools/";
		slideDirList=getFileList(runDir);
		print(runDir, slideDirList.length);
		for (h = 0; h < slideDirList.length; h++) {
			roiList=getFileList(runDir + slideDirList[h]);
			// if(startsWith(slideDirList[h], 'Pano')) continue;
			print(roiList.length);
			for (j = 0; j < roiList.length; j++) {
				tma=replace(slideDirList[h],  "/", "");
				roi=replace(roiList[j], "/", "");
				imgName=tma + "-" + roi;
				print(imgName);
				
				if(image != "all" && imgName != image)
					continue;
				stacks=getFileList(runDir+slideDirList[h]+ roiList[j]);
				if(stacks.length==0)
					continue;
				fileOut=outDir + 'panckf_' + imgName;
				if(File.exists(fileOut)) continue;
				newPanel=0;
				for (m = 0; m < stacks.length; m++)	{
					print(stacks[m]);
					imgDir=runDir + slideDirList[h] + roiList[j] + stacks[m];
					if(stacks.length==0) {
						print("ERROR: empty stack " + stacks[m] + "\n");
						exit;
					}
					if(stacks[m] != "full_stack/")
						continue;
					fileOut = outDir + tma + "-" + roi + ".tiff";
					print(fileOut);
					
					print('IMG dir: ', imgDir);
					if(File.isDirectory(imgDir)) {
						imgList=getFileList(imgDir);
						for ( i=0; i<imgList.length; i++ ) {
						   	if(imgList[i] == '142Nd_CAM52.tiff' || imgList[i] == '142Nd_CAM52Nd142Di.tiff')
							newPanel=1;
							
						}
					}
					if(newPanel == 1) continue;
					open(imgDir + "164Dy_panCK.tiff");
					// run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
					run("Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=3");

					autoAdjust();
					selectWindow("164Dy_panCK.tiff");
					getMinAndMax(min, max);
					print('Tumour', min, max);
					if(max > panck_thres && max < 5000 && imgName != 'P1_TMA006_R_20190619-roi_19' &&
						imgName != 'P1_TMA004_L_20190619-roi_2' && imgName != 'P1_TMA006_R_20190619-roi_14' &&
						imgName != 'P1_TMA003_R_20190619-roi_9' && imgName != 'P1_TMA003_R_20190619-roi_9' &&
						imgName != 'P1_TMA003_L_20190619_correct-roi_7' && imgName != 'P1_TMA003_L_20190619_correct-roi_6' &&
						imgName != 'P1_TMA003_L_20190619_correct-roi_19' && imgName != 'P1_TMA002_L_20190619-roi_9' &&
						imgName != 'P1_TMA002_L_20190619-roi_18' && imgName != 'P1_TMA002_L_20190619-roi_13' && imgName != 'P1_TMA001_R_20190508-roi_10' &&
						imgName != 'P1_TMA001_L_20190508-roi_5' && imgName != 'P1_TMA001_L_20190508-roi_19' && imgName != 'P1_TMA_REC_20190508-roi_7' ) {
						run("Close All");
						continue;
					}
					print('Additional  filtering on the tumour intensities');
					run("Morphological Filters", "operation=[Dilation] element=Square radius=1");
					run("Directional Filtering", "type=Max operation=Median line=5 direction=32");
					run("Morphological Filters", "operation=[Erosion] element=Square radius=1");
					run('Median (3D)');
//					run("Median...", "radius=1");
			
					run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
					singleMarkerImg="Tumour";
					rename(singleMarkerImg);
					saveAs('Tiff', outDir + 'panckf_' + imgName);
					run("Close All");
				}
			}
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
          * Lubna Ahmad
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
		
		if(hmax > aggregateMax && hmax == 0) {
			hmin=noiseLevel;
			hmax=10;
		}
		if(hmin > aggregateMin) hmin=noiseLevel;
		if(hmin < 0) hmin=0;
		if (hmax > noiseLevel) setMinAndMax(hmin, hmax); 
		 print(hmin, hmax);
		// run(“Apply LUT”);
}



