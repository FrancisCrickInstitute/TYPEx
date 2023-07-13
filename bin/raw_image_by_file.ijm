#@String inDir
#@String posFile

setBatchMode("hide");

selectedImgs = File.openAsString(posFile);
lines = split(selectedImgs, '\n');

for (i=0; i < lines.length; i++) {
	
	line = split(lines[i], "\t");
	imageID = replace(line[0], "-", "/") + "/";
	
	print(line[0]);
	marker=line[1];
	marker = format_marker(marker);
	fileID = marker + ".tiff";

	imagePaths=getFiles(inDir, imageID + "full_stack/", fileID);
	if(imagePaths.length == 0) continue;

	open(imagePaths[0]);
	imageInfo=marker + ".." + line[0] + ".." + line[2];
	rename(imageInfo);
	getMinAndMax(min, max);
	autoAdjust();
	saveAs("PNG", imageInfo + ".png");
	run("Close All");
}

function getFiles(dir, imgID, extension) {

	print(dir, imgID, extension);
	result = newArray();
	imgDir=dir + imgID;
	
	list = getFileList(imgDir);
	print(list.length, imgDir);
	for (i=0; i<list.length; i++) {
		 	showProgress(i, list.length);
		    if (endsWith(list[i], "/")) {
		       files=getFiles("" + imgDir + list[i], extension);
		       result=Array.concat(result, files);
		    } else if(endsWith(list[i], extension)) {
		    	newFile=imgDir + list[i];
		    	print(newFile);
		    	result=Array.concat(result, newFile);
		    }
		}
	return result;
}

function automatedBrightnessAdjustment(image, limit, nBins, minThreshold, minThreshold)	{

	selectWindow(image);
	getRawStatistics(nPixels, mean, min, max);
	getHistogram(values, counts, nBins);
	cumulCount=0;
	max=0;
	for (i = 0; i < nBins; i++) {
		cumulCount = cumulCount + counts[i];
		print(nPixels - cumulCount, (nPixels - cumulCount)/nPixels, values[i]);
		if(values[i] < minThreshold) continue;
		if(cumulCount <= limit || min == 0) {
			min=values[i];
		}
		if(nPixels - cumulCount > limit) {
			max=values[i];
		}
	}
	if(max > maxThreshold) setMinAndMax(min, max);
	if(max < maxThreshold) setMinAndMax(min, maxThreshold);
	print(min, max, image);
}


function autoAdjust() { 
	/* 
	 * rewriting "Auto contrast adjustment" button of "Brightness/Contrast" 
	 * Damien Guimond 
	 * 20120516 
	 * Acknowledgements: Kota Miura 
	 */ 

	aggregateMin=50;
	aggregateMax=5000;
    AUTO_THRESHOLD = 5000;
    getRawStatistics(pixcount);
    limit = pixcount/10; 
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
	if(hmin > aggregateMin || hmax > aggregateMax) {
		hmin = 0; hmax = 10;
	}
	if(hmin > 3) hmin=3;
	setMinAndMax(hmin, hmax); 
	print(hmin, hmax); 
	// run("Apply LUT");  
} 

function format_marker(name) {
	if(name == "Vimentin") return "vimentin";
	if(name == "alphaSMA") return "aSMA";
	if(name == "casp3") return "CASP3";
	return name;
}

