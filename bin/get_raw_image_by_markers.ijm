panel='p1';
segRun="publication";
inDir="/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/nextflow/p1/publication/";
outDir="/camp/lab/swantonc/working/angelom/analyses/imc/test_pub/output/single_channel/";
posFile="data/naive.txt";

//inDir='/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/nextflow/p2/publication/';
// outDir='/camp/lab/swantonc/working/angelom/analyses/imc/test_pub_p2/output/single_channel/';
// posFile='data/over_under_qc_images.txt';
print(panel);
print(posFile);


if(File.isDirectory(outDir) == 0) File.makeDirectory(outDir);
selectedImgs=File.openAsString(posFile);
lines=split(selectedImgs, '\n');
selection=0;
print(lines.length);
for (i=0; i < lines.length; i++) {
	
	line=split(lines[i], "\t");
	print(lines[i]);
	imageID=replace(line[selection], "-", "/") + "/";
	
	print(line[selection]);
	marker=line[1];
	marker=format_marker(marker);
	fileID=marker + ".tiff";

	imagePaths=getFiles(inDir, imageID + "full_stack/", fileID);
	print(imagePaths.length);
	if(imagePaths.length == 0) continue;
	print(imagePaths[0]);
	//run("Bio-Formats Importer", "open=" + imagePaths[0] + " autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

	open(imagePaths[0]);
	
	imageInfo=marker + ".." + line[selection] + ".." + line[2];
	rename(imageInfo);
	getMinAndMax(min, max);
	// automatedBrightnessAdjustment(imageInfo, pixelLimit, nBins, minThreshold, maxThreshold);
	autoAdjust();		
	getMinAndMax(min, max);
	print(min, max);
	imageInfo=imageInfo + ".." + min + ".." + max;
	rename(imageInfo);
	saveAs("PNG", outDir + imageInfo + ".png");
	run("Close All");
}

function getFiles(dir, imgID, extension) {
	
	runs=getFileList(dir);
	result=newArray();
	print(runs.length, dir);
	for(j = 0; j < runs.length; j++) {
		print(runs[j]);
		imgDir=dir + runs[j]+ '/results/imctools/' + imgID;
		
		list = getFileList(imgDir);
		for (i=0; i<list.length; i++) {
			 	showProgress(i, list.length);
			    if (endsWith(list[i], "/")) {
			       files=getFiles(""+imgDir+list[i], extension);
			       result=Array.concat(result, files);
			    } else if(endsWith(list[i], extension)) {
			    	newFile=imgDir + list[i];
			    	print(newFile);
			    	result=Array.concat(result, newFile);
			    }
 		}
	}
	return result;
}


function automatedBrightnessAdjustment(image, limit, nBins, minThreshold, maxThreshold)	{

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

