#!/usr/bin/env python

import os
import sys
import glob
import numpy as np
import re
import skimage as sk
import pandas as pd
from scipy import ndimage
from skimage import measure
from skimage import io

# INPUT

# SETTINGS
if(len(sys.argv) < 3):
	sys.exit("Usage: python %s <panel> <run>" % sys.argv[0])
infile=sys.argv[1]
tiffFile=sys.argv[2]
imageid=sys.argv[3]
outDir=sys.argv[4]

#cohort=sys.argv[3]

#inDir = os.path.join("/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx", cohort, "imc/outputs/nextflow")
#outDir = os.path.join("/camp/home/angelom/labwd/analyses/imc/typing/major/ilastik_stack/csm/", run, panel, cohort)

if not(os.path.exists(outDir) and os.path.isdir(outDir)):
	os.makedirs(outDir, exist_ok=True)


print(infile)
frame = pd.read_csv(infile)
# imageid = file.replace(inDir, "").replace(filePattern, "").replace("/", "_")
coordinates = []
marker = re.sub("raw_([^_]+)_*mask.tiff", "\\1", os.path.basename(tiffFile))
fileOut = os.path.join(outDir, ".".join([imageid, marker, "ratios.csv"]))
tiff = io.imread(tiffFile)
# labeled_img = sk.color.label2rgb(tiff, bg_label=0)
# labeled_mask, num_labels = ndimage.label(labeled_img

objects = sk.measure.regionprops(tiff)
for obj in objects:
	diff = np.sqrt(abs(frame['Location_Center_Y'] - obj.centroid[0]) ** 2 +
			   abs(frame['Location_Center_X'] - obj.centroid[1]) ** 2)
	idx = np.where(diff == min(diff))[0][0]
	d = np.array([[obj.centroid[1], obj.centroid[0],
		frame['Location_Center_X'][idx],
		frame['Location_Center_Y'][idx],
		frame['ObjectNumber'][idx],
	marker, np.round(min(diff))]])
	d = pd.DataFrame(d, columns = [ 'X', "Y", 'Location_Center_X', 'Location_Center_Y',
		'ObjectNumber', 'marker', "NN"])
	coordinates.append(d)
if(len(coordinates) > 0):
	coordinates = pd.concat(coordinates, sort = False)
	coordinates.to_csv(fileOut)

