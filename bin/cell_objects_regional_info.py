#!/usr/bin/env python
import os
import glob
import numpy as np
import pandas as pd
import re
from matplotlib import pyplot
import argparse

parser = argparse.ArgumentParser(description='Allocate cell objects to tissue areas based on binary masks')
parser.add_argument('--panel', default="p2",
    help='Antibody panel')
parser.add_argument('--maskDir', metavar='wDir',
    default="tissue_seg/labels", help='NextFlow Run')
parser.add_argument('--tissueAreaDir', metavar='wDir',
    default="tissue_seg/labels", help='NextFlow Run')
parser.add_argument('--cellObjFile', metavar='coords',
    help='NextFlow Run')
parser.add_argument('--maskRegEx', metavar='coords',
    help='NextFlow Run')
parser.add_argument('--regionType', default = 'regional',
    metavar='coords', help='Tissue category')
parser.add_argument('--imageRegEx', default = '\\2',
    metavar='coords', help='RegEx index for imageID')
parser.add_argument('--regionRegEx', default = None,
    metavar='coords', help='RegEx index for region type')
parser.add_argument('--outDir', default = None,
    metavar='out', help='RegEx index for region type')
    
args = parser.parse_args()

outDir = os.path.join(args.outDir, "cell_info")
if(not os.path.isdir(outDir)):
    print("Creating ", outDir)
    os.mkdir(outDir)

cellFrame=pd.read_csv(args.cellObjFile, sep =",")
print(cellFrame.columns.values.tolist())

# if not os.path.exists(args.maskDir) and args.regionType == 'regional':
if args.regionType == 'regional':
    args.maskDir = os.path.join(args.tissueAreaDir, 'labels')

maskList=glob.glob(os.path.join(args.maskDir,
    args.maskRegEx.replace(".*", "*").
                   replace("(", "").
                   replace(")", "")))
print(args.maskRegEx.replace(".*", "*").
                   replace("(", "").
                   replace(")", ""))
fileOut=os.path.join(outDir, args.panel + "_" + args.regionType + "_cell_info.csv")
print(len(maskList), " masks listed in dir ", 
    args.maskDir, ' with regex ', args.maskRegEx)

if(len(maskList) > 0):
    regionInfo = []
    for mask in maskList:
        
        if(args.regionRegEx is not None):
            region=re.sub(args.maskRegEx, args.regionRegEx, os.path.basename(mask))
        else:
            region=args.regionType
        
        imagename=re.sub(args.maskRegEx, args.imageRegEx, os.path.basename(mask))
        tiff = pyplot.imread(mask)

        imgIndex=cellFrame['imagename'] == imagename

        dataTmp=cellFrame[imgIndex]
        
        print(region, imagename)
        print(dataTmp.shape)
        if(dataTmp.shape[0] == 0):
            continue
        
        print(max(dataTmp["LocationCenter_X"]))
        print(max(dataTmp["LocationCenter_Y"]))
        print(tiff.shape)
        print(region, imagename)
        # Check the orientation of the mask
        if(max(dataTmp["LocationCenter_X"]) <= tiff.shape[0] and max(dataTmp["LocationCenter_Y"]) <= tiff.shape[1]):
            regions=[
                tiff[int(np.round(x)), int(np.round(y))]
                    for x, y in zip(dataTmp['LocationCenter_X'],
                                    dataTmp['LocationCenter_Y'])
            ]
        elif(max(dataTmp["LocationCenter_X"]) <= tiff.shape[1] and max(dataTmp["LocationCenter_Y"]) <= tiff.shape[0]):
	        regions=[
    	        tiff[int(np.round(y)), int(np.round(x))]
        	        for x, y in zip(dataTmp['LocationCenter_X'], 
            	                    dataTmp['LocationCenter_Y'])
        	]
        d = {'centerX':dataTmp['LocationCenter_X'],
             'centerY':dataTmp['LocationCenter_Y'],
             'ObjectNumber':dataTmp['ObjectNumber'],
             'imagename':imagename,
             "region":region,
             "regionID":regions
        }
        d = pd.DataFrame(d)
        regionInfo.append(d)

    if(len(regionInfo) == 0):
        print('ERROR: No data to append. Check the regular expression ' +
             ' for extracting the imagename from the input files. Exiting')
    else:
        print(len(regionInfo))
        summary = pd.concat(regionInfo, sort = True)
        summary.to_csv(fileOut, index=False)
        print('Output saved in ', fileOut)
