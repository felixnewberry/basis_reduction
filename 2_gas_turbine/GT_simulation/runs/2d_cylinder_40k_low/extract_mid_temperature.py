#!/usr/bin/python

import sys
import os
import numpy as np
from paraview.simple import *
from paraview.servermanager import *
from paraview.vtk import *


if len(sys.argv) != 4:
  sys.exit('Usage: pvpython {} run_identifier pht_file.pht'.format(os.path.basename(__file__)))
else:
  runID = sys.argv[1]
  dataFileName = sys.argv[2]
  lastTimeStep = sys.argv[3]

pluginPath = '/usr/local/paraview/5.3.0/lib/paraview-5.3/plugins/libPhastaSyncIOReader.so'
#pluginPath = '/usr/local/paraview/5.4.1/lib/paraview-5.4/plugins/libPhastaSyncIOReader.so'
LoadPlugin(pluginPath, remote=False)

print 'Loading ParaView file... {}'.format(dataFileName)
#phtFile = PhastaReader(FileName=dataFileName)
phtFile = PhastaSyncIOReader(FileName=dataFileName)
#phtFile = PhastaReader(FileName=flowPosix_extract.pht)
print 'Loaded.'

### GENERATED COMMANDS FROM PARAVIEW ###

###############################################################################
# Extract temperature at end of domain
###############################################################################
# might have to figure this stuff out anew...

# breaks on this:  UpdatePipeline
#phtFile.UpdatePipeline()

#### disable automatic camera reset on 'Show'
#paraview.simple._DisableFirstRenderCameraReset()

print 'flowPosix ready to go'
# create a new 'Merge Blocks'
mergeBlocks1 = MergeBlocks(Input=phtFile)

###
# Merge blocs for ease of saving later. All filters start from merge blocks from now
###

print 'mergeBlocks ready'

###
## create a new 'Slice' for mid line at 0.2. 
###
slice1 = Slice(Input=mergeBlocks1)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

## init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.2, 0.0, 0.0]
slice1.SliceType.Normal = [0.2, 0.0, 0.0]


###
## create a new 'Slice' to make 2d. 
###
slice2 = Slice(Input=slice1)
slice2.SliceType = 'Plane'
slice2.SliceOffsetValues = [0.0]

## init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [0.0, 0.0, 0.001]
slice2.SliceType.Normal = [0.0, 0.0, 1.0]


print 'saving'

SaveData('mid_temp_%d.csv' % int(lastTimeStep), proxy=slice2) #lastTimeStep... 
