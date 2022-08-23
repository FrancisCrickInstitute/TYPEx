#@String wDir
#@String outDir

run("Set Measurements...", "area centroid perimeter fit shape area_fraction limit display redirect=None decimal=3");



minArea=100;
iter=2;

// INPUT SETTINGS
probDir=wDir + '/composites/';
outDirArea = outDir + "area_info/";
outDirMasks = outDir + "masks/";
outDirLabels = outDir + "labels/";
outDirOverlay = outDir + "overlay/";
outDirFract=outDir + 'fractal/';
if(File.exists(outDirFract) == 0) {
	File.makeDirectory(wDir + "segmentation/");
	File.makeDirectory(outDir);
	File.makeDirectory(outDirArea);
	File.makeDirectory(outDirMasks);
	File.makeDirectory(outDirOverlay);
	File.makeDirectory(outDirLabels);
	File.makeDirectory(outDirFract);
}
print(probDir);
slideDirList = getFileList(probDir);
print(slideDirList.length);
setBatchMode(true);
run("ROI Manager...");
for (i = 0; i < slideDirList.length; i++) {

	if(! endsWith(slideDirList[i], "_Probabilities.h5")) continue;
		imgName = replace(slideDirList[i], "_Probabilities.h5", "");
		// if(imgName != 'P1_TMA003_L_20190619_correct-roi_4') continue;
		print(imgName);
		fileOut=outDirOverlay + "overlay_" + imgName +".png";
		fdOut=outDirFract +  imgName + "_fd.csv";

		run("Import HDF5", "select=" + probDir+slideDirList[i] + " datasetname=[/exported_data: (1747, 1756, 2) uint8] axisorder=yxc");
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
			//setAutoThreshold("Default white");
			//run("Convert to Mask");
			rename(outFile);
			print(outFile);
		//	run("Median...", "radius=3");
			run("Options...", "iterations=2 count=4 black pad do=Dilate");
			run("Options...", "iterations=2 count=4 black pad do=Erode");
	    	// saveAs("tiff", outDirMasks + "test3_" +  imgTitle[j]);
		}
		
		imgTitle = getList("image.titles");
		for (j = 0; j < imgTitle.length; j++) {
			print('Processing', imgTitle[j]);
			selectWindow(imgTitle[j]);
			//run("Open", "stack");
			
		//	if(startsWith(imgTitle[j], "Tumour")) {
		//		imageCalculator("Subtract ", imgTitle[j], "Background_"+ imgName + ".tiff" );
	//			imageCalculator("Subtract ", imgTitle[j], "Stroma_"+ imgName + ".tiff" );
	//		} else if(isOpen("Tumour_" + imgName + ".tiff")) {
	//			imageCalculator("Subtract ", imgTitle[j], "Tumour_" + imgName + ".tiff" );
	//		}
			// run("Options...", "iterations=2 count=4 black pad do=Dilate");
			// run("Fill Holes");
			// run("Options...", "iterations=2 count=4 black pad do=Erode");
			//run("Close-", "stack");
			 list = getList("image.titles");
			if (list.length==0)
			     print("No image windows are open");
		 	else {
			     print("Image windows:");
			     for (im=0; im<list.length; im++)
			        print("   "+list[im]);
			}
			print("");
			run("Analyze Particles...", "size=" + minArea + "-Infinity pixel display clear summarize");
			// if(isOpen("Results")) {
			saveAs("Results", outDirArea + "Results_" + imgTitle[j] + ".csv");
			close('Results');
			close('Plot');
			// continue;
			//}
		
			//run("Invert");
			run("Connected Components Labeling", "connectivity=4 type=[16 bits]");
			saveAs("tiff", outDirLabels + "labeled_regions_" + imgTitle[j]);
			close();
			if (startsWith(imgTitle[j], "Background")) {
				close();
				continue;
			}
		}
		selectWindow("Tumour_" + imgName + '.tiff');
//		run("Duplicate...","title=T");
//		run("Options...", "iterations=5 count=1 black pad do=Dilate");
//		selectWindow("Stroma_" + imgName + '.tiff');
//		run("Duplicate...", "title=S");
//		run("Options...", "iterations=5 count=1 black pad do=Dilate");
//		imageCalculator("AND create", "T","S");
//		//saveAs('Tiff', fdOut + 'tiff');
//		getRawStatistics(nPixels, mean, min, max, std);
//		if(max > min) {
//			run("Fractal Box Count...", "box=10,20,30,50,100,150,200,250,300 black");
//			saveAs("Results", fdOut);
//		}
//		close('T');
//		close('S');
//		close('Result of T');
//		close('Plot');
//		// run('Close All');
//		// continue;
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
//		run("Bio-Formats Importer", "open="+rawImgName + " autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		run("Add Image...", "image=Filled" + " x=0 y=0 opacity=20");
		saveAs("PNG", outDirOverlay + "overlay_" + imgName +".png");
		run("Close All");
}

//save the result with appropriate names
selectWindow("Summary");
saveAs("Results", outDirArea+"Summary_area." + iter + ".csv");

