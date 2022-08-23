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

panck_thres=60
print(rootDir);
if(File.isDirectory(outDir) == 0)
	File.makeDirectory(outDir);

setBatchMode(true);
runs=getFileList(rootDir);
print(runs.length);
for (k=0; k < runs.length; k++)	{
//	print(runs[k]);
	if(startsWith(runs[k], "run")) {
		runDir=rootDir + runs[k] + "results/imctools/";
		slideDirList=getFileList(runDir);
//		print(runDir, slideDirList.length);
		for (h = 0; h < slideDirList.length; h++) {
			roiList=getFileList(runDir + slideDirList[h]);
			// if(startsWith(slideDirList[h], 'Pano')) continue;
//			print(roiList.length);
			for (j = 0; j < roiList.length; j++) {
				tma=replace(slideDirList[h],  "/", "");
				roi=replace(roiList[j], "/", "");
				imgName=tma + "-" + roi;
			    //if(imgName != 'P2_TMA007_L_20190619-roi_6') continue;
//				print(imgName);
				if(image != "all" && imgName != image)
					continue;
				stacks=getFileList(runDir+slideDirList[h]+ roiList[j]);
//				print(stacks.length);
				if(stacks.length==0) continue;
				fileOut=outDir + 'panckf_' + imgName + '.tif'; 
//				print(fileOut);
				if(File.exists(fileOut)) {
					print('File exists', fileOut);
					continue;
				}

				newPanel=0;
				for (m = 0; m < stacks.length; m++)	{
//					print(stacks[m]);
					imgDir=runDir + slideDirList[h] + roiList[j] + stacks[m];
					if(stacks.length==0) {
						print("ERROR: empty stack " + stacks[m] + "\n");
						exit;
					}
					if(stacks[m] != "full_stack/")
						continue;
					//fileOut = outDir + tma + "-" + roi + ".tiff";
					//if(File.exists(fileOut)) continue;	
					if(File.isDirectory(imgDir)) {
						imgList=getFileList(imgDir);
						for ( i=0; i<imgList.length; i++ ) {
						   	if(imgList[i] == '142Nd_CAM52.tiff' || imgList[i] == '142Nd_CAM52Nd142Di.tiff')
							newPanel=1;
							
						}
					}
//					print('New panel', newPanel);
					if(newPanel == 1) continue;
					open(imgDir + "164Dy_panCK.tiff");
					// run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
					run("Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=2");

					autoAdjust();
					selectWindow("164Dy_panCK.tiff");
					getMinAndMax(min, max);
					print(imgName, 'Tumour', min, max);
					
					if(max > panck_thres && max < 5000 && imgName != 'P2_TMA_REC_20190508-roi_7' &&
                            imgName != 'P2_TMA002_R_20190619-roi_8' && imgName != 'P2_TMA004_R_20190619-roi_20' &&
                            imgName != 'P2_TMA005_R_20190619-roi_5' && imgName != 'P2_TMA006_20190619_L-roi_15' &&
                            imgName != 'P2_TMA006_R_20190619-roi_13' && imgName != 'P2_TMA006_R_20190619-roi_14' && 
							imgName != 'P2_TMA004_L_20190619-roi_20' &&
                            imgName != 'P2_TMA006_20190619_L-roi_15' && imgName != 'P2_TMA006_R_20190619-roi_1' &&  
                            imgName != 'P2_TMA002_L_20190619-roi_19 ' && imgName != 'P2_TMA005_R_20190619-roi_14' &&
                            imgName != 'P2_TMA002_R_20190619-roi_7' && imgName != 'P2_TMA004_L_20190619-roi_26' &&
                            imgName != 'P2_TMA004_L_20190619-roi_6' && imgName != 'P2_TMA004_R_20190619-roi_13' &&
                            imgName != 'P2_TMA006_R_20190619-roi_11' && imgName != 'P2_TMA007_L_20190619-roi_12' &&
							imgName != 'P2_TMA_REC_20190508-roi_8' && imgName != 'P2_TMA002_L_20190619-roi_3' &&
                            imgName != 'P2_TMA004_L_20190619-roi_7' && imgName != 'P2_TMA004_R_20190619-roi_19' &&
                            imgName != 'P2_TMA005_R_20190619-roi_4' && imgName != 'P2_TMA005_R_20190619-roi_5' &&
                            imgName != 'P2_TMA006_R_20190619-roi_15' && imgName != 'P2_TMA007_L_20190619-roi_12'  && 
							imgName != "P2_TMA001_L_20190508-roi_4" && imgName != "P2_TMA001_L_20190508-roi_5" && 
							imgName != "P2_TMA001_L_20190508-roi_7" && imgName != "P2_TMA001_R_20190508-roi_14" &&
							imgName != "P2_TMA005_L_20190619-roi_16" && imgName != "P2_TMA005_R_20190619-roi_13" &&
							imgName != "P2_TMA006_R_20190619-roi_4" && imgName != "P2_TMA001_R_20190508-roi_10" &&
							imgName != "P2_TMA002_L_20190619-roi_7" && imgName != "P2_TMA003_L_20190806-roi_4" && 
							imgName != "P2_TMA004_L_20190619-roi_8" && imgName != "P2_TMA005_L_20190619-roi_14" &&
							imgName != "P2_TMA005_R_20190619-roi_18" && imgName != "P2_TMA006_20190619_L-roi_10" &&
							imgName != "P2_TMA006_R_20190619-roi_18" && imgName != "P2_TMA006_R_20190619-roi_8" && 
							imgName != "P2_TMA007_L_20190619-roi_7" && imgName != "P2_TMA001_L_20190508-roi_14" &&
							imgName != "P2_TMA002_R_20190619-roi_15" && imgName != "P2_TMA002_R_20190619-roi_2" && 
							imgName != "P2_TMA002_R_20190619-roi_3" && imgName != "P2_TMA002_R_20190619-roi_7" && 
							imgName != "P2_TMA002_R_20190619-roi_8" && imgName != "P2_TMA003_L_20190806-roi_8" && 
							imgName != "P2_TMA004_L_20190619-roi_13" && imgName != "P2_TMA004_L_20190619-roi_20" && 
							imgName != "P2_TMA004_L_20190619-roi_21" && imgName != "P2_TMA004_L_20190619-roi_6" && 
							imgName != "P2_TMA004_R_20190619-roi_6" && imgName != "P2_TMA005_L_20190619-roi_15" &&
							imgName != "P2_TMA005_L_20190619-roi_4" && imgName != "P2_TMA005_R_20190619-roi_14" &&
							imgName != "P2_TMA005_R_20190619-roi_16" && imgName != "P2_TMA005_R_20190619-roi_17" &&
							imgName != "P2_TMA005_R_20190619-roi_2" && imgName != "P2_TMA006_20190619_L-roi_15" &&
							imgName != "P2_TMA006_R_20190619-roi_11" && imgName != "P2_TMA006_R_20190619-roi_17" &&
							imgName != "P2_TMA006_R_20190619-roi_7" && imgName != "P2_TMA007_R_20190619-roi_5" && 
							imgName != "P2_TMA007_R_20190619-roi_6"  && imgName != "P2_TMA002_R_20190619-roi_7" && 
							imgName != "P2_TMA002_R_20190619-roi_8" && imgName != "P2_TMA003_R_20190619-roi_2" && 
							imgName != "P2_TMA004_L_20190619-roi_20" && imgName != "P2_TMA004_R_20190619-roi_10" && 
							imgName != "P2_TMA005_L_20190619-roi_15" && imgName != "P2_TMA005_L_20190619-roi_4" && 
							imgName != "P2_TMA005_R_20190619-roi_17" && imgName != "P2_TMA005_R_20190619-roi_2" && 
							imgName != "P2_TMA006_R_20190619-roi_4" && imgName != "P2_TMA006_R_20190619-roi_7" && 
							imgName != "P2_TMA005_R_20190619-roi_18" && imgName != 'P2_TMA001_L_20190508-roi_13' &&
							imgName != 'P2_TMA003_L_20190806-roi_10' && imgName != 'P2_TMA003_L_20190806-roi_14' && 
							imgName != 'P2_TMA004_L_20190619-roi_5' && imgName != 'P2_TMA004_L_20190619-roi_7' && 
							imgName != 'P2_TMA004_L_20190619-roi_8' && imgName != 'P2_TMA004_L_20190619-roi_9' && 
							imgName != 'P2_TMA004_L_20190619-roi_14' && imgName != 'P2_TMA004_L_20190619-roi_26' &&
							imgName != 'P2_TMA004_R_20190619-roi_4' && imgName != 'P2_TMA004_R_20190619-roi_6' &&
							imgName != 'P2_TMA004_R_20190619-roi_2' && imgName != 'P2_TMA005_L_20190619-roi_4' &&
							imgName != 'P2_TMA005_R_20190619-roi_5' && imgName != 'P2_TMA005_R_20190619-roi_6' && 
							imgName != 'P2_TMA005_R_20190619-roi_11' && imgName != 'P2_TMA005_R_20190619-roi_9' && 
							imgName != 'P2_TMA005_R_20190619-roi_4' &&  imgName != 'P2_TMA005_R_20190619-roi_20' && 
							 imgName != 'P2_TMA005_R_20190619-roi_24' && imgName != 'P2_TMA006_R_20190619-roi_9' &&
							imgName != 'P2_TMA006_R_20190619-roi_13' && imgName != 'P2_TMA007_L_20190619-roi_6' && 
							imgName != 'P2_TMA002_R_20190619-roi_15' && imgName != 'P2_TMA002_R_20190619-roi_3' &&
                        imgName != 'P2_TMA004_R_20190619-roi_6' && imgName != 'P2_TMA005_R_20190619-roi_4' &&
                        imgName != 'P2_TMA005_R_20190619-roi_14' && imgName != 'P2_TMA005_R_20190619-roi_5' &&
                        imgName !=  "P2_TMA001_L_20190508-roi_14" && imgName != "P2_TMA001_R_20190508-roi_1" &&
                        imgName !=  "P2_TMA002_R_20190619-roi_2" && imgName != "P2_TMA002_R_20190619-roi_8" &&
                        imgName !=  "P2_TMA003_L_20190806-roi_8" && imgName != "P2_TMA004_L_20190619-roi_20" &&
                        imgName !=  "P2_TMA004_L_20190619-roi_21" && imgName != "P2_TMA004_R_20190619-roi_10" &&
                        imgName !=  "P2_TMA005_L_20190619-roi_15" && imgName != "P2_TMA005_L_20190619-roi_4" &&
                        imgName !=  "P2_TMA005_R_20190619-roi_16" && imgName != "P2_TMA005_R_20190619-roi_17" &&
                        imgName !=   "P2_TMA005_R_20190619-roi_2" && imgName != 'P2_TMA004_L_20190619-roi_21' &&
                        imgName != 'P2_TMA006_R_20190619-roi_15' && imgName != 'P2_TMA006_R_20190619-roi_17' &&
                        imgName != 'P2_TMA006_R_20190619-roi_7' && imgName != 'P2_TMA007_R_20190619-roi_5' &&
                        imgName != 'P2_TMA007_R_20190619-roi_6' && imgName != 'P2_TMA005_R_20190619-roi_13' &&
                        imgName != 'P2_TMA001_L_20190508-roi_5' && imgName != 'P2_TMA006_R_20190619-roi_4') { 

						run("Close All");
						continue;
					}
				//	print('Additional  filtering on the tumour intensities');
					if(imgName == 'P2_TMA002_R_20190619-roi_15' || imgName == 'P2_TMA002_R_20190619-roi_3' || 
						imgName == 'P2_TMA004_R_20190619-roi_6' || imgName == 'P2_TMA005_R_20190619-roi_4' || 
						imgName == 'P2_TMA005_R_20190619-roi_14' || imgName == 'P2_TMA005_R_20190619-roi_5' ||
						imgName ==	"P2_TMA001_L_20190508-roi_14" || imgName == "P2_TMA001_R_20190508-roi_1" ||
						imgName ==	"P2_TMA002_R_20190619-roi_2" || imgName == "P2_TMA002_R_20190619-roi_8" ||
						imgName ==	"P2_TMA003_L_20190806-roi_8" || imgName == "P2_TMA004_L_20190619-roi_20" ||
						imgName ==	"P2_TMA004_L_20190619-roi_21" || imgName == "P2_TMA004_R_20190619-roi_10" ||
						imgName ==	"P2_TMA005_L_20190619-roi_15" || imgName == "P2_TMA005_L_20190619-roi_4" ||
						imgName ==	"P2_TMA005_R_20190619-roi_16" || imgName == "P2_TMA005_R_20190619-roi_17" ||
						imgName ==	 "P2_TMA005_R_20190619-roi_2" || imgName == 'P2_TMA004_L_20190619-roi_21' ||
						imgName == 'P2_TMA006_R_20190619-roi_15' || imgName == 'P2_TMA006_R_20190619-roi_17' ||
						imgName == 'P2_TMA006_R_20190619-roi_7' || imgName == 'P2_TMA007_R_20190619-roi_5' ||
						imgName == 'P2_TMA007_R_20190619-roi_6' || imgName == 'P2_TMA005_R_20190619-roi_13' ||
						imgName == 'P2_TMA001_L_20190508-roi_5' || imgName == 'P2_TMA006_R_20190619-roi_4') {
						run("Morphological Filters", "operation=[Dilation] element=Square radius=1");
						run("Directional Filtering", "type=Max operation=Median line=2 direction=32");
						run("Morphological Filters", "operation=[Erosion] element=Square radius=1");
					} else {
						run("Directional Filtering", "type=Max operation=Median line=5 direction=32");
					}
					run('Median (3D)');
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



