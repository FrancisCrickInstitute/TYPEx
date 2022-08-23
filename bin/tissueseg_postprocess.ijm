#@String wDir

run("Set Measurements...", "area centroid perimeter fit shape area_fraction limit display redirect=None decimal=3");

minArea=100;
iter=2;

// INPUT SETTINGS
probDir=wDir + '/composites/';
outDirArea = "area_info/";
outDirMasks = "masks/";
outDirLabels = "labels/";
outDirOverlay = "overlay/";

if(File.exists(outDirLabels) == 0) {
	File.makeDirectory(outDirArea);
	File.makeDirectory(outDirMasks);
	File.makeDirectory(outDirOverlay);
	File.makeDirectory(outDirLabels);
}
print(probDir);

slideDirList = getFileList(probDir);
print(slideDirList.length);
setBatchMode(true);
run("ROI Manager...");

for (i = 0; i < slideDirList.length; i++) {

	if(! endsWith(slideDirList[i], "_Probabilities.h5")) continue;
		imgName = replace(slideDirList[i], "_Probabilities.h5", "");
		//if(imgName != 'P2_TMA002_R_20190619-roi_3') continue;
		print(imgName);
		fileOut=outDirOverlay + "overlay_" + imgName +".png";
		if(File.exists(fileOut)) continue;

		run("Import HDF5", "select=" + probDir+slideDirList[i] + 
			" datasetname=[/exported_data: (1747, 1756, 2) uint8] axisorder=yxc");
		run("Make Composite", "display=Color");
		run("Duplicate...", "duplicate");
		run("Auto Threshold", "method=Default white show use_stack_histogram");
		run("Make Composite", "display=Composite");
		Stack.setActiveChannels("1110");
		saveAs("PNG", outDirMasks + "raw_" + imgName);
		
		run("Stack to Images");
		rawImgNameBase=replace(imgName, "Composite_", "");
		rawImgName = wDir + '/composites/' + rawImgNameBase + ".tiff";
	
		imgTitle = getList("image.titles");
		for (j = 0; j < imgTitle.length; j++) {
			print(imgTitle[j]);
			selectWindow(imgTitle[j]);
			if (imgTitle[j] == '1') {
				outFile = "Tumour" + "_" + imgName + ".tiff";
			} else if(imgTitle[j] == '2') {
				outFile = "Stroma" + "_" + imgName + ".tiff";
			} else if(imgTitle[j] == '3') {
				outFile = "Background" + "_" + imgName + ".tiff";
			} else {
				print("Closed:", imgTitle[j]);
				close();
				continue;
			}
			setBackgroundColor('black');
			run("8-bit");
			run("Median (3D)");
			run("Options...", "iterations=2 count=4 black pad do=Dilate");
			run("Options...", "iterations=2 count=4 black pad do=Erode");
			rename(outFile);
			print(outFile);
			close(imgTitle[j]);
		}
		
		imgTitle = getList("image.titles");
		for (j = 0; j < imgTitle.length; j++) {
			print('Processing', imgTitle[j]);
			selectWindow(imgTitle[j]);
			//run("Open", "stack");
			
			// for intratumoral infiltration
			if(startsWith(imgTitle[j], "Tumour")) {
				run("Calculator Plus", "i1=Tumour_" + imgName + ".tiff i2=Stroma_" + imgName + 
					".tiff operation=[Add: i2 = (i1+i2) x k1 + k2] k1=1 k2=0 create");
				close(imgTitle[j]);
				run("Calculator Plus", "i1=Result i2=Stroma_" + imgName + 
					".tiff operation=[Subtract: i2 = (i1-i2) x k1 + k2] k1=1 k2=0 create");
				rename(imgTitle[j]);
				close('Result');
			} 
			run("Analyze Particles...", "size=" + minArea + "-Infinity pixel clear summarize");
			saveAs("Results", outDirArea + "Results_" + imgTitle[j] + ".csv");
			close('Results');
			close('Plot');
			run("Connected Components Labeling", "connectivity=4 type=[16 bits]");
			saveAs("tiff", outDirLabels + "labeled_regions_" + imgTitle[j]);
			close();
			if (startsWith(imgTitle[j], "Background")) {
				close();
				continue;
			}
		}
		selectWindow("Tumour_" + imgName + '.tiff');
		imgTitle = getList("image.titles");
        for (j = 0; j < imgTitle.length; j++) {
			print(imgTitle[j]);
		}
		run("Images to Stack");
		run("Make Composite", "display=Composite");
		Stack.setActiveChannels("1110");
		saveAs("tiff", outDirMasks + "segmentation_" + imgName);
		saveAs("PNG", outDirMasks + "segmentation_" + imgName);

		print('run RGB');
		run("RGB Color");
		rename("Filled");
		print(rawImgName, 'opening');
		open(rawImgName);
		run("Add Image...", "image=Filled" + " x=0 y=0 opacity=20");
		saveAs("PNG", outDirOverlay + "overlay_" + imgName +".png");
		run("Close All");
}

//save the result with appropriate names
selectWindow("Summary");
saveAs("Results", outDirArea+"Summary_area." + iter + ".csv");

