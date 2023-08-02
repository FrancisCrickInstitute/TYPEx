#@String rootDir
#@String tumour

image="all";
directional='none';
panck_thres=15
f = File.open('panck_preprocessed_images.txt');

setBatchMode(true);
runs = getFileList(rootDir);
print(runs.length + ' images in ' + rootDir);

print(tumour);
tumourMarkers = split(tumour, "=|,");
print(tumourMarkers[0]);
print(tumourMarkers.length);
panckID = "";
for (m = 1; m < tumourMarkers.length; m++) {
	index = indexOf(tumourMarkers[m], "panCK");
	print(tumourMarkers[m]);
	if(index != -1)
		panckID = tumourMarkers[m];
}
if(panckID == "") {
	print("No panck tumour marker"); 
	exit;
}
print(panckID);

for (k=0; k < runs.length; k++)	{
	
	runDir = rootDir + runs[k];
	print(runDir);
	roiList = getFileList(runDir);

	for (j = 0; j < roiList.length; j++) {
		
			tma = replace(runs[k],  "/", "");
			roi = replace(roiList[j], "/", "");
			imgName = tma + "-" + roi;
			print(imgName);
			print(imgName);

			if(image != "all" && imgName != image)
				continue;
			stacks = getFileList(runDir + roiList[j]);

			if(stacks.length == 0) {
				print("Stack is empty");
				continue;
			}
			fileOut = 'panckf_' + imgName + '.tif'; 

			if(File.exists(fileOut)) {
				print('File exists', fileOut);
				continue;
			}

			for (m = 0; m < stacks.length; m++)	{
				print(stacks[m]);
				if(! endsWith(stacks[m], "full_stack/")) 
                	continue;
				
				imgDir = runDir + roiList[j] + "/" + stacks[m];
				print('imgDir');
				print(imgDir);
				if(stacks.length == 0) {
					print("ERROR: empty stack " + stacks[m] + "\n");
					exit;
				}
				
				open(imgDir + panckID);
				selectWindow(panckID);

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

File.close(f);

function autoAdjust() {
	// creating a new function - to include auto brightness and contrast in macro
	
        /*
         * Edited from “Auto contrast adjustment” button of “Brightness/Contrast”
         * Damien Guimond
         * 20120516
         * Acknowledgements: Kota Miura
		 * Edits: Mihaela A.
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

