#@String inDir
#@String posFile
#@String outDir

setBatchMode("hide");

if(File.isDirectory(outDir) == 0)
	File.makeDirectory(outDir); 
print("Opening " + posFile);
selectedImgs = File.openAsString(posFile);
lines = split(selectedImgs, '\n');

for (i=0; i < lines.length; i++) {
	
	print(i);
	line = split(lines[i], "\t");
	marker=line[1];
	marker = format_marker(marker);
	imagePaths = getFiles(inDir, line[0], marker, ".tiff");
	print(imagePaths.length);
	if(imagePaths.length == 0) 
		continue;
	if(imagePaths.length > 1) {
		print("Multiple paths found");
		print(line[0]);
		for(j = 0; j < imagePaths.length; j++) {
			print(imagePaths[j]);
			imgID = replace(line[0], "-", "//");
			lineIndex = indexOf(imagePaths[j], imgID);
			print(lineIndex + " " + imgID + " " + imagePaths[j]);
			if(lineIndex > -1) {
				pathIndex = j;
				break;
			} 
			
		}
	} else {
		pathIndex = 0;
	}
	open(imagePaths[pathIndex]);
	imageInfo=marker + ".." + line[0] + ".." + line[2];
	rename(imageInfo);
	getMinAndMax(min, max);
	autoAdjust();
	saveAs("PNG", outDir + "/" + imageInfo + ".png");
	run("Close All");
}

function getFiles(dir, imgID, marker, extension) {

	print(dir, imgID, marker, extension);
	result = newArray();

	list = getFileList(dir);
	for (i = 0; i < list.length; i++) {

		showProgress(i, list.length);
		dirName = replace(list[i], "\\/", "");
		if (endsWith(list[i], "/")) {

			dirIndex = indexOf(imgID, dirName);			
			if(endsWith(dirName, "full_stack"))
				dirIndex = 0;
				
			if(dirIndex > -1) {
				files = getFiles("" + dir + "/" + list[i], imgID, marker, extension);
				result = Array.concat(result, files);
			}
		} else if(endsWith(list[i], extension)) {
			index = indexOf(list[i], marker + ".");
			if(index == -1) continue;
			result = Array.concat(result, dir + list[i]);
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

