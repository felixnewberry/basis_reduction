#!/usr/bin/python

import sys
import os
import numpy as np
from paraview.simple import *
from paraview.servermanager import *
from paraview.vtk import *


if len(sys.argv) != 5:
  sys.exit('Usage: pvpython {} run_identifier pht_file.pht'.format(os.path.basename(__file__)))
else:
  runID = sys.argv[1]
  dataFileName = sys.argv[2]
  lastTimeStep = sys.argv[3]
  Colorbar_on = sys.argv[4]

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
# Save image of cylinder temp field
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

#SaveData('cylinder_temp_%d.csv' % int(lastTimeStep), proxy=slice1) #lastTimeStep... 
rv = GetActiveViewOrCreate('RenderView')
rv.Background = [1, 1, 1]

"""Set camera properties. """
rv.ViewSize = [1920, 1080]
rv.CameraPosition = [0.15, -0.05, 1]
rv.CameraFocalPoint = [0.15, -0.05, 0.0] # also reasonable could be 0.03
rv.CameraViewAngle = 30.0
rv.CameraViewUp = [0, 1, 0]



# Put a colorbar in if desired


aaaaaaaaaaaa
# set scalar coloring
ColorBy(flowPosix_extract_1.phts, ('POINTS', 'temp'))
#ColorBy(temp, ('SURFACE', field))
set_rgb_colormap(field, minVal=285, maxVal=400, pr1max=pr1max)
show_colorbar(d, rv, field, numberFormat=strFmt, numTicks=ticks)
position_colorbar(d, rv, fet scalar coloring
ColorBy(flowPosix_extract_1phtsDisplay, ('POINTS', 'temp'))
 
# rescale color and/or opacity maps used to include current data range
flowPosix_extract_1phtsDisplay.RescaleTransferFunctionToDataRange(True, False)
 
# show color bar/color legend
flowPosix_extract_1phtsDisplay.SetScalarBarVisibility(renderView1, True)
 
# get color transfer function/color map for 'temp'
tempLUT = GetColorTransferFunction('temp')

# get opacity transfer function/opacity map for 'temp'
tempPWF = GetOpacityTransferFunction('temp') 
 
# Set scale to be same between runs
# Rescale transfer function
tempPWF.RescaleTransferFunction(285.0, 400.0)

# Rescale transfer function
tempLUT.RescaleTransferFunction(285.0, 400.0)


eld, [0.9, 0.05], 0.4)
HideUnusedScalarBars()
fname = data_path(dataDir, 'spanslice_{0}_slice{1}_{2}.png'.format(field,i,xstr))
SaveScreenshot(fname, view=rv)

