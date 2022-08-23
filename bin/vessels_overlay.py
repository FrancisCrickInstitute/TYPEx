#!/usr/bin/env python
import os
import glob
import numpy as np
import pandas as pd
import re
from matplotlib import pyplot
import argparse

parser = argparse.ArgumentParser(description='Get average marker expression per tissue category')
parser.add_argument('--run', metavar='r', default="publication", help='NextFlow Run')
parser.add_argument('--panel', default="p2", help='Antibody panel')
parser.add_argument('--analysisDir', metavar='wDir', default="output", help='NextFlow Run')

args = parser.parse_args()

maskRegEx="labeled_regions_(.*)_" + args.panel.upper() + "_(.*).tiff"

# Input file settings
feature='LocationCenter'
centerX="LocationCenter_X"
centerY="LocationCenter_Y"
imageID="imagename"

# Files
segDir=os.path.join(args.analysisDir, "tissue_seg/")
maskDir=os.path.join(segDir, "labels")
cellObjDir=os.path.join(args.analysisDir, "features", feature)
outDir=os.path.join(segDir, "cell_info")
if(not os.path.isdir(outDir)):
    os.mkdir(outDir)

cellObjFile=args.panel + "_LocationCenter_" + args.run + ".csv"
cellFrame=pd.read_csv(os.path.join(cellObjDir, cellObjFile), sep =",")

maskList=glob.glob(os.path.join(maskDir, maskRegEx.replace("(.*)", "*")))
fileOut=os.path.join(outDir, args.panel + "_regional_cell_info_" + args.run + ".csv")
print(len(maskList), " masks listed")
regionInfo = []
for mask in maskList:
	region=re.sub(maskRegEx, "\\1", os.path.basename(mask))
	imagename=re.sub(maskRegEx, "\\2", os.path.basename(mask))
	print(region, imagename)
	tiff = pyplot.imread(mask)

	dataTmp=cellFrame[cellFrame[imageID] == args.panel.upper() + "_" + imagename]
	print(dataTmp[centerX])
	regions=[tiff[int(np.round(y)), int(np.round(x))] for x, y in zip(dataTmp[centerX], dataTmp[centerY])]
	d = {'centerX':dataTmp[centerX], 'centerY':dataTmp[centerY],
         'ObjectNumber':dataTmp['ObjectNumber'],
         'imagename':imagename,
         "region":region,
         "regionID":regions}
	d = pd.DataFrame(d)
	regionInfo.append(d)

summary = pd.concat(regionInfo, sort = True)
summary.to_csv(fileOut, index=False)
print('Output saved in ', fileOut)
