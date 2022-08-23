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
panck_threshold=60;

print(rootDir);
if(File.isDirectory(outDir) == 0)
	File.makeDirectory(outDir);

if(panel == "p1" && mcd) {
	immune_markers="image1=[Median of 152Sm_CD45.tiff] image2=[Median of 170Er_CD3.tiff] image3=[Median of 162Dy_CD8a.tiff] image4=[Median of 156Gd_CD4.tiff] image5=[-- None --]";
	auxStromaMarkers='image1=[Median of 143Nd_vimentin.tiff] image2=[Median of 169Tm_collagen1.tiff] image3=[-- None --]';
	tumour='image1=[Median of 164Dy_panCK.tiff] image2=[Median of 142Nd_CAM52.tiff]';
	//stroma_images='image1=141Pr_aSMA.tiff image2=151Eu_CD31.tiff image3=[-- None --]';
	stroma_images='image1=[Median of 141Pr_aSMA.tiff] image2=[Median of 151Eu_CD31.tiff] image3=auxStromaMrg image4=[-- None --]';
	dna='image1=[Median of 191Ir_DNA1.tiff] image2=[Median of 193Ir_DNA2.tiff]';
} else if (panel == 'p2' && mcd) {
	immune_markers="image1=[Median of 152Sm_CD45.tiff] image2=[Median of 170Er_CD3.tiff] image3=[Median of 162Dy_CD8a.tiff] image5=[Median of 169Tm_CD206.tiff] image6=[Median of 161Dy_CD20.tiff] image7=[Median of 149Sm_CD11b.tiff] image8=[Median of 158Gd_CD79a.tiff] image9=[Median of 147Sm_CD163.tiff] image10=[Median of 146Nd_CD16.tiff] image11=[Median of 159Tb_CD68.tiff] image12=[Median of 144Sm_CD14.tiff] image13=[Median of 156Gd_CD4.tiff] image14=[-- None --]";
	auxStroma='Median of 175Lu_panactin.tiff';
	stroma_images='image1=[Median of 151Eu_CD31.tiff] image2=auxStromaMrg image3=[-- None --]';
	tumour='image1=[Median of 164Dy_panCK.tiff] image2=[Median of 142Nd_CAM52.tiff]';
	dna='image1=[Median of 191Ir_DNA1.tiff] image2=[Median of 193Ir_DNA2.tiff]';
} else if(panel == "p1") {
	immune_markers="[Median of 152Sm_CD45Sm152Di.tiff] image2=[Median of 170Er_CD3Er170Di.tiff] image3=[Median of 162Dy_CD8aDy162Di.tiff] image4=[Median of 156Gd_CD4Gd156Di.tiff] image5=[-- None --]";
	auxStromaMarkers='image1=[Median of 143Nd_vimentin.tiff] image2=[Median of 169Tm_collagen1.tiff] image3=[-- None --]';
	stroma_images='image1=[Median of 141Pr_aSMA.tiff] image2=[Median of 151Eu_CD31.tiff] image3=auxStromaMrg image4=[-- None --]';
	dna='image1=[Median of 191Ir_DNA1Ir191Di.tiff] image2=[Median of 193Ir_DNA2Ir193Di.tiff]';
	tumour='image1=[Median of 164Dy_panCKDy164Di.tiff] image2=[Median of 142Nd_CAM52Nd142Di.tiff]';
} else if(panel == "p2") {
	dna='image1=[Median of 191Ir_DNA1Ir191Di.tiff] image2=[Median of 193Ir_DNA2Ir193Di.tiff]';
	immune_markers="image1=[Median of 152Sm_CD45Sm152Di.tiff] image2=[Median of 170Er_CD3Er170Di.tiff] image3=[Median of 162Dy_CD8aDy162Di.tiff] image5=[Median of 156Gd_CD4Gd156Di.tiff] image6=[Median of 161Dy_CD20Dy161Di.tiff] image7=[Median of 149Sm_CD11bSm149Di.tiff] image8=[Median of 158Gd_CD79aGd158Di.tiff] image9=[Median of 147Sm_CD163Sm147Di.tiff] image10=[Median of 146Nd_CD16Nd146Di.tiff] image11=[Median of 159Tb_CD68Tb159Di.tiff] image12=[Median of 144Sm_CD14.tiff] image13=[Median of 156Gd_CD4Gd156Di.tiff] image14=[-- None --]";
	auxStroma='Median of 175Lu_panactinLu175Di.tiff';
	tumour='image1=[Median of 164Dy_panCKDy164Di.tiff] image2=[Median of 142Nd_CAM52Nd142Di.tiff]';
	stroma_images='image1=[Median of 151Eu_CD31Eu151Di.tiff] image2=auxStromaMrg image3=[-- None --]';
}

setBatchMode(true);
runs=getFileList(rootDir);
print(runs.length);
for (k=0; k < runs.length; k++)	{
	print('RUNL ', runs[k]);
	if(startsWith(runs[k], "run")) {
		runDir=rootDir + runs[k] + "results/imctools/";
		slideDirList=getFileList(runDir);
		print('Run dir', runDir, slideDirList.length);
		for (h = 0; h < slideDirList.length; h++) {
			roiList=getFileList(runDir + slideDirList[h]);
			// if(startsWith(slideDirList[h], 'Pano')) continue;
			print('ROI list', roiList.length);
			for (j = 0; j < roiList.length; j++) {
				tma=replace(slideDirList[h],  "/", "");
				roi=replace(roiList[j], "/", "");
				imgName=tma + "-" + roi;
				print('IMG name:  ', imgName);
				//if(imgName != 'P2_TMA_REC_20190508-roi_8') continue;	
				if(image != "all" && imgName != image) {
					print(image, 'skipping');
						continue;
				}	
				stacks=getFileList(runDir+slideDirList[h]+ roiList[j]);
				if(stacks.length==0) {
					print('Stack size is 0. Skipping');
					continue;
				}
				newPanel=0;
				for (m = 0; m < stacks.length; m++)	{
					print('Stack: ', stacks[m]);
					imgDir=runDir + slideDirList[h] + roiList[j] + stacks[m];
					if(stacks.length==0) {
						print("ERROR: empty stack " + stacks[m] + "\n");
						exit;
					}
					if(stacks[m] != "full_stack/")
						continue;
					fileOut = outDir + tma + "-" + roi + ".tiff";
				    if(File.exists(fileOut)) {
						print('File exists', fileOut);
						 continue;	
					}
					print('Image dir', imgDir);
					if(File.isDirectory(imgDir)) {
						imgList=getFileList(imgDir);
						for ( i=0; i<imgList.length; i++ ) {
							open( imgDir + imgList[i] );
							run("Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=2");
						   	if(imgList[i] == '142Nd_CAM52.tiff' || imgList[i] == '142Nd_CAM52Nd142Di.tiff')
							newPanel=1;
							print(imgList[i]);
							autoAdjust();
							run("Median (3D)");
						}
					}
					print(newPanel);
					// DNA Sum
					run("Concatenate...", " title=DNA open " + dna);
					run("Z Project...", "projection=[Sum Slices]");
	//				run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
					run("Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=2");
					autoAdjust();
					run("Enhance Contrast", "saturated=0.35");

					if(newPanel) {
						run("Concatenate...", " title=Test open " + tumour);
						run("Z Project...", "projection=[Sum Slices]");
						rename('Tumour');
						autoAdjust();
					} else {
						selectWindow("164Dy_panCK.tiff");
						run("Median (3D)");
						run("Median (3D)");
						getMinAndMax(min, max);
						print(imgName, 'Tumour', min, max);
						panckf=0;
						if(max < panck_threshold || imgName == 'P2_TMA_REC_20190508-roi_7' || 
							imgName == 'P2_TMA002_R_20190619-roi_8' || imgName == 'P2_TMA004_R_20190619-roi_20' ||
							imgName == 'P2_TMA005_R_20190619-roi_5' || imgName == 'P2_TMA006_20190619_L-roi_15' ||
							imgName == 'P2_TMA006_R_20190619-roi_13' || imgName == 'P2_TMA006_R_20190619-roi_14' ||
							imgName == 'P2_TMA002_R_20190619-roi_15' ||
							imgName == 'P2_TMA006_20190619_L-roi_15' || imgName == 'P2_TMA006_R_20190619-roi_1' || 
							imgName == 'P2_TMA002_L_20190619-roi_19' || imgName == 'P2_TMA005_R_20190619-roi_14' ||
							imgName == 'P2_TMA002_R_20190619-roi_7' || imgName == 'P2_TMA004_L_20190619-roi_26' ||
							imgName == 'P2_TMA004_L_20190619-roi_6' || imgName == 'P2_TMA004_R_20190619-roi_13' ||
							imgName == 'P2_TMA006_R_20190619-roi_11' || imgName == 'P2_TMA007_L_20190619-roi_12' ||	
							imgName == 'P2_TMA_REC_20190508-roi_8' || imgName == 'P2_TMA002_L_20190619-roi_3' ||
							imgName == 'P2_TMA004_L_20190619-roi_7' || imgName == 'P2_TMA004_R_20190619-roi_19' ||
							imgName == 'P2_TMA005_R_20190619-roi_4' || imgName == 'P2_TMA005_R_20190619-roi_5' ||
							imgName == 'P2_TMA006_R_20190619-roi_15' || imgName == 'P2_TMA007_L_20190619-roi_12' ||
							imgName == "P2_TMA001_L_20190508-roi_4" || imgName == "P2_TMA001_L_20190508-roi_5" ||
                            imgName == "P2_TMA001_L_20190508-roi_7" || imgName == "P2_TMA001_R_20190508-roi_14" ||
                            imgName == "P2_TMA005_L_20190619-roi_16" || imgName == "P2_TMA005_R_20190619-roi_13" ||
                            imgName == "P2_TMA006_R_20190619-roi_4" || imgName == "P2_TMA001_R_20190508-roi_10" ||
                            imgName == "P2_TMA002_L_20190619-roi_7" || imgName == "P2_TMA003_L_20190806-roi_4" ||
                            imgName == "P2_TMA004_L_20190619-roi_8" || imgName == "P2_TMA005_L_20190619-roi_14" ||
                            imgName == "P2_TMA005_R_20190619-roi_18" || imgName == "P2_TMA006_20190619_L-roi_10" ||
                            imgName == "P2_TMA006_R_20190619-roi_18" || imgName == "P2_TMA006_R_20190619-roi_8" ||
                            imgName == "P2_TMA007_L_20190619-roi_7" || imgName == "P2_TMA001_L_20190508-roi_14" ||
                        	imgName == "P2_TMA002_R_20190619-roi_2" || imgName == "P2_TMA002_R_20190619-roi_3" ||
                            imgName == "P2_TMA002_R_20190619-roi_8" || imgName == "P2_TMA003_L_20190806-roi_8" ||
                            imgName == "P2_TMA004_L_20190619-roi_13" ||
                            imgName == "P2_TMA004_L_20190619-roi_21" || imgName == "P2_TMA004_L_20190619-roi_6" ||
                            imgName == "P2_TMA004_R_20190619-roi_6" || imgName == "P2_TMA005_L_20190619-roi_15" ||
                            imgName == "P2_TMA005_L_20190619-roi_4" || imgName == "P2_TMA005_R_20190619-roi_14" ||
                            imgName == "P2_TMA005_R_20190619-roi_16" || imgName == "P2_TMA005_R_20190619-roi_17" ||
                            imgName == "P2_TMA005_R_20190619-roi_2" || imgName == "P2_TMA006_20190619_L-roi_15" ||
                            imgName == "P2_TMA006_R_20190619-roi_11" || imgName == "P2_TMA006_R_20190619-roi_17" ||
                            imgName == "P2_TMA006_R_20190619-roi_7" || imgName == "P2_TMA007_R_20190619-roi_5" ||
                            imgName == "P2_TMA007_R_20190619-roi_6"  || imgName == "P2_TMA002_R_20190619-roi_7" || 
							imgName == "P2_TMA002_R_20190619-roi_8" || imgName == "P2_TMA003_R_20190619-roi_2" || 
							imgName == "P2_TMA004_L_20190619-roi_20" || imgName == "P2_TMA004_R_20190619-roi_10" || 
							imgName == "P2_TMA005_L_20190619-roi_15" || imgName == "P2_TMA005_L_20190619-roi_4" || 
							imgName == "P2_TMA005_R_20190619-roi_17" || imgName == "P2_TMA005_R_20190619-roi_2" || 
							imgName == "P2_TMA006_R_20190619-roi_4" || imgName == "P2_TMA006_R_20190619-roi_7" || 
							imgName == "P2_TMA005_R_20190619-roi_18" || imgName == 'P2_TMA001_L_20190508-roi_13' ||
							imgName == 'P2_TMA003_L_20190806-roi_10' || imgName == "P2_TMA003_L_20190806-roi_14" ||
							imgName == 'P2_TMA004_L_20190619-roi_5' || imgName == 'P2_TMA004_L_20190619-roi_9' ||
                            imgName == 'P2_TMA004_L_20190619-roi_14' ||
							imgName == 'P2_TMA004_R_20190619-roi_4' || imgName == 'P2_TMA004_R_20190619-roi_6' ||
                            imgName == 'P2_TMA004_R_20190619-roi_2' || imgName == 'P2_TMA005_L_20190619-roi_4' ||
							imgName == 'P2_TMA005_R_20190619-roi_5' || imgName == 'P2_TMA005_R_20190619-roi_6' ||
                            imgName == 'P2_TMA005_R_20190619-roi_11' || imgName == 'P2_TMA005_R_20190619-roi_9' ||
                            imgName == 'P2_TMA005_R_20190619-roi_4' ||  imgName == 'P2_TMA005_R_20190619-roi_20' || 
                            imgName == 'P2_TMA005_R_20190619-roi_24' || imgName == 'P2_TMA006_R_20190619-roi_9' ||
							imgName == 'P2_TMA006_R_20190619-roi_13' || imgName == 'P2_TMA007_L_20190619-roi_6' || 
						imgName == 'P2_TMA002_R_20190619-roi_15' || imgName == 'P2_TMA002_R_20190619-roi_3' ||
                        imgName == 'P2_TMA004_R_20190619-roi_6' || imgName == 'P2_TMA005_R_20190619-roi_4' ||
                        imgName == 'P2_TMA005_R_20190619-roi_14' || imgName == 'P2_TMA005_R_20190619-roi_5' ||
                        imgName ==  "P2_TMA001_L_20190508-roi_14" || imgName == "P2_TMA001_R_20190508-roi_1" ||
                        imgName ==  "P2_TMA002_R_20190619-roi_2" || imgName == "P2_TMA002_R_20190619-roi_8" ||
                        imgName ==  "P2_TMA003_L_20190806-roi_8" || imgName == "P2_TMA004_L_20190619-roi_20" ||
                        imgName ==  "P2_TMA004_L_20190619-roi_21" || imgName == "P2_TMA004_R_20190619-roi_10" ||
                        imgName ==  "P2_TMA005_L_20190619-roi_15" || imgName == "P2_TMA005_L_20190619-roi_4" ||
                        imgName ==  "P2_TMA005_R_20190619-roi_16" || imgName == "P2_TMA005_R_20190619-roi_17" ||
                        imgName ==   "P2_TMA005_R_20190619-roi_2" || imgName == 'P2_TMA004_L_20190619-roi_21' ||
                        imgName == 'P2_TMA006_R_20190619-roi_15' || imgName == 'P2_TMA006_R_20190619-roi_17' ||
                        imgName == 'P2_TMA006_R_20190619-roi_7' || imgName == 'P2_TMA007_R_20190619-roi_5' ||
                        imgName == 'P2_TMA007_R_20190619-roi_6' || imgName == 'P2_TMA005_R_20190619-roi_13' ||
                        imgName == 'P2_TMA001_L_20190508-roi_5' || imgName == 'P2_TMA006_R_20190619-roi_4' ) {
                        


//							if(imgName != 'P2_TMA004_L_20190619-roi_13' && imgName != 'P2_TMA006_R_20190619-roi_17' &&
//								imgName != 'P2_TMA004_L_20190619-roi_13' ) {
								open(outDir + 'panckf_' + imgName + '.tif');
								print(outDir + 'panckf_' + imgName  + '.tif');
								panckf=1;
//							}
						} 
						run("Median (3D)");
						autoAdjust();
						if(panckf == 0)
							run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
						singleMarkerImg="Tumour";
						rename(singleMarkerImg);
					}

					// Stroma
					if(panel == 'p1') {
						run("Concatenate...", "  title=Nonsp open " + auxStromaMarkers);
						run("Z Project...", "projection=[Sum Slices]");
	//					run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
						run("Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=2");
						autoAdjust();
//						run("Morphological Filters", "operation=Opening element=Square radius=1");
						// run("Median...", "radius=1");
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
					run("Concatenate...", "  title=Immune open " + immune_markers);
					run("Z Project...", "projection=[Sum Slices]");
	//				run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
					run("Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=2");
					print('Immune');
					autoAdjust();
					// run("Median...", "radius=1");
					run("Enhance Contrast", "saturated=0.35");
                    
                    run("Duplicate...", "title=Immune_mask");
                    setAutoThreshold("Triangle dark");
                    run("Convert to Mask");
                    run("Invert");
                    //run("Divide...", "value=255");
					//run("Macro...", "code=v=v/255");
					//Create a binary mask as Divide doesn't work
					run("32-bit");
					run("Calculator Plus", "i1=Immune_mask i2=Immune_mask operation=[Divide: i2 = (i1/i2) x k1 + k2] k1=1 k2=0");
                    run("Calculator Plus", "i1=auxStroma i2=Immune_mask operation=[Multiply: i2 = (i1*i2) x k1 + k2] k1=1 k2=0 create");
                    selectWindow("Result");
                    rename('Corrected_for_immune');

					selectWindow("Tumour");
                    run("Duplicate...", "title=Tumour_mask");
					// Max to generate a tumour mask - considering nuclei are not in it
					run("Max (3D)");
                    setAutoThreshold("Triangle dark");
                    run("Convert to Mask");
					run("Morphological Filters", "operation=Closing element=Square radius=2");
                    run("Invert");
					// run("Median (3D)");
                    // run("Divide...", "value=255.000");
					run("32-bit");
					run("Calculator Plus", "i1=Maximum-Closing i2=Maximum-Closing operation=[Divide: i2 = (i1/i2) x k1 + k2] k1=1 k2=0");
                    run("Calculator Plus", "i1=Corrected_for_immune i2=Maximum-Closing operation=[Multiply: i2 = (i1*i2) x k1 + k2] k1=1 k2=0 create");
                    selectWindow("Result");
					run("Median (3D)");
                    rename('auxStromaMrg');
					print(stroma_images);
					//list = getList("image.titles");
					//if (list.length==0)
					//	print("No non-image windows are open");
					//else {
					//	print("Non-image windows:");
					//for (i=0; i<list.length; i++)
					//	print("   "+list[i]);
  					//}
					//print("");
                    run("Concatenate...", "  title=Stroma open " + stroma_images);
                    run("Z Project...", "projection=[Sum Slices]");
      //            run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
//					run("Morphological Filters", "operation=Opening element=Square radius=1");
					run("Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=2");
					run("Median (3D)");
					autoAdjust();
					run("Enhance Contrast", "saturated=0.35");

					//run("Median...", "radius=1");

					if(imgName == 'P1_TMA_REC_20190508-roi_10') {

                            selectWindow('SUM_Stroma');
                            run("Duplicate...", " ");
                            run("Magenta");
                            saveAs('Tiff',  outDir + "stroma" + tma + "-" + roi + ".tiff");
							close();

							selectWindow('SUM_DNA');
							run("Duplicate...", " ");
                            run("Blue");
                            saveAs('Tiff',  outDir + "dna" + tma + "-" + roi + ".tiff");
							close();

							selectWindow('SUM_Immune');
							run("Duplicate...", " ");
                            run("Magenta");
                            saveAs('Tiff',  outDir + "immune" + tma + "-" + roi + ".tiff");
							close();
							
							selectWindow('Tumour');
							run("Duplicate...", " ");
							run("Yellow");
                            saveAs('Tiff',  outDir + "tumour" + tma + "-" + roi + ".tiff");
							close();
					}
	
					run("Concatenate...", "  title=StromaStack open image1=SUM_Immune image2=[Median of SUM_Stroma] image3=[-- None --]");
					run("Z Project...", "projection=[Sum Slices]");
					rename("Stroma_merge");
					autoAdjust();
					run("Merge Channels...", "c3=SUM_DNA c7=Tumour c6=Stroma_merge create");
					run("RGB Color");
					saveAs("Tiff", fileOut);
					print(fileOut);
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


