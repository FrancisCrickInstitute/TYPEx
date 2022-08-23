rootDir='output/tissue_seg/masks/';
runs=getFileList(rootDir);
for (k=0; k < runs.length; k++) {
		if(!endsWith(runs[k], 'tif')) continue;
		if(!startsWith(runs[k], 'segmentation')) continue;
		print(runs[k]);
        open(rootDir + runs[k]);
        run("Stack to Images");
        imgTitle = getList("image.titles");
        print(runs[k], imgTitle.length);
		run("Close All");
}
