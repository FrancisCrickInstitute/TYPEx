#@String rootDir

// Choose whether to analyse all images ("all") or specific e.g. "P1_TMA003_R_20190619-roi_12", or "P1_TMA004_L_20190619-roi_13"
image="all";
// Choose whether to apply directional filter ("direct") or "none"
directional='none';
// Choose the NextFlow run iteration

panck_thres=15

setBatchMode(true);
runs=getFileList(rootDir);
print(runs.length + ' images in ' + rootDir);

f = File.open('panck_range.txt');

for (k=0; k < runs.length; k++)	{
	
	runDir=rootDir + runs[k];
	slideDirList=getFileList(runDir);

	for (h = 0; h < slideDirList.length; h++) {
		roiList=getFileList(runDir + slideDirList[h]);

		for (j = 0; j < roiList.length; j++) {
		
			tma=replace(slideDirList[h],  "/", "");
			roi=replace(roiList[j], "/", "");
			imgName=tma + "-" + roi;

			if(image != "all" && imgName != image)
				continue;
			stacks=getFileList(runDir+slideDirList[h]+ roiList[j]);

			if(stacks.length==0) continue;
			fileOut='panckf_' + imgName + '.tif'; 

			if(File.exists(fileOut)) {
				print('File exists', fileOut);
				continue;
			}

			newPanel=0;
			for (m = 0; m < stacks.length; m++)	{
				if(! endsWith(stacks[m], "full_stack/")) continue;
				
				imgDir=runDir + slideDirList[h] + roiList[j] + stacks[m];
				print(imgDir);
				if(stacks.length==0) {
					print("ERROR: empty stack " + stacks[m] + "\n");
					exit;
				}
				if(stacks[m] != "full_stack/")
					continue;
				
				if(File.isDirectory(imgDir)) {
					imgList=getFileList(imgDir);
					for ( i=0; i<imgList.length; i++ ) {
					   	if(imgList[i] == '142Nd_CAM52.tiff' || imgList[i] == '142Nd_CAM52Nd142Di.tiff')
						newPanel=1;
						
					}
				}

				if(newPanel == 1) continue;
				open(imgDir + "164Dy_panCK.tiff");
				selectWindow("164Dy_panCK.tiff");

				autoAdjust();
				getMinAndMax(min, max);
				print(imgName, 'Tumour', min, max);
				print(f, imgName + "\t" + 'Tumour' + '\t' + min + '\t' +  max);			

			//	print('Additional  filtering on the tumour intensities');
				if(max < panck_thres) {
					run("Morphological Filters", "operation=[Dilation] element=Square radius=1");
					run("Directional Filtering", "type=Max operation=Median line=2 direction=32");
					run("Morphological Filters", "operation=[Erosion] element=Square radius=1");
				} else {
					run("Directional Filtering", "type=Max operation=Median line=1 direction=32");
				}
				run('Median (3D)');
				run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
				singleMarkerImg="Tumour";
				rename(singleMarkerImg);
				saveAs('Tiff', 'panckf_' + imgName);
				run("Close All");
			}
		}
	}
}

File.close(f);

function autoAdjust() {
	// creating a new function - to include auto brightness and contrast in macro
	
        /*
         * Edited from “Auto contrast adjustment” button of “Brightness/Contrast”
         * Damien Guimond
         * 20120516
         * Acknowledgements: Kota Miura
         */
	aggregateMax=5000;
	aggregateMin=15;
        noiseLevel=2;
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

