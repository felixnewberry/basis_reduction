from __future__ import division # force floating-point division

import sys
import os
import math
import argparse
from paraview.simple import *

def data_path(dataDir, filename):
  return os.path.join(dataDir, filename)

def format_coord_for_filename(coord):
  return '{0:.3f}'.format(coord).replace('.','p')

def check_fields(phtFile):
  """Ensure we have access to the right point data files in the pht file"""
  requiredPointData = {
    'dwal': 'distance to wall',
    'u': 'vector of velocity components',
    'T': 'temperature'
  }
  if not all(key in phtFile.PointData.keys() for key in requiredPointData.keys()):
    print 'error; required fields are...'
    l = len(max(requiredPointData.keys(), key=len))
    for key, val in requiredPointData.iteritems():
      print '  ' + key.rjust(l), '=', val
    err('ParaView file does not contain all of the required requisite point data fields (see above)')

def print_fields(proxy):
  """Print the fields that proxy has access to; useful for debugging"""
  print 'Fields in', proxy, 'are', proxy.PointData.keys()

def start_over_visualization(renderView=None):
  """Hide all sources in pipeline."""
  rv = renderView
  if not rv:
    rv = GetActiveViewOrCreate('RenderView')
  for source in GetSources().values():
    Hide(source)
  rv.Background = [1, 1, 1]

# Alter these from aip, duct and duct iso to appropriate for cylinder applcation. Ie not duct at all... 

def set_view(rv, idstr):
  """Set camera properties.
  
  Keyword arguments:
  rv -- the RenderView to modify
  idstr -- the string specifying the desired view
  """
  if idstr == 'cylinder':
    rv.ViewSize = [1920, 1080]
    rv.CameraPosition = [0.13, 0, 1]
    rv.CameraFocalPoint = [0.13, 0, 0.0005] # also reasonable could be 0.03
    rv.CameraViewAngle = 30.0
    rv.CameraViewUp = [0, 1, 0]
  elif idstr == 'aip':
    rv.ViewSize = [1920, 1080]
    rv.CameraPosition = [500, 20, 40]
    rv.CameraFocalPoint = [aipx, aipy, 0.025] # also reasonable could be 0.03
    rv.CameraViewAngle = 0.02
    rv.CameraViewUp = [0, 1, 0]
  elif idstr == 'duct':
    rv.ViewSize = [1920, 1080]
    rv.CameraPosition = [rampx, 0, 1200]
    rv.CameraFocalPoint = [rampx, 0, 0]
    rv.CameraViewAngle = 0.02
    rv.CameraViewUp = [0, 1, 0]
  elif idstr == 'ductiso':
    rv.ViewSize = [1920, 1080]
    rv.CameraPosition = [1, 0.5, 1.2]
    rv.CameraFocalPoint = [0.14, 0, 0]
    rv.CameraViewAngle = 15.0
    rv.CameraViewUp = [0, 1, 0]
  elif idstr == 'separation':
    # Separation should be the 'ductiso' view, but also rotate about AIP center for movie
    pass
  else:
    raise ValueError('View {0} not supported'.format(idstr))

def set_solid_color(proxy, view):
  """Set a proxy to render as a solid color"""
  display = GetDisplayProperties(proxy, view)
  display.DiffuseColor = [1, 1, 1] # Set to white (default color, shows up as gray)
  display.ColorArrayName = [None, '']

def set_flat_lighting(proxy, view):
  """Set a proxy to render without shadows"""
  display = GetDisplayProperties(proxy, view)
  display.Specular = 0.0
  display.Ambient = 1.0
  display.Diffuse = 0.0

def set_rgb_colormap(fieldName, minVal, maxVal, pr1max):
  """Set desired field to RGB colormap with min/max values
  
  Keyword arguments:
  fieldName -- name of the field to set colormap for
  minVal -- minimum field value on colormap
  maxVal -- maximum field value on colormap
  """
  lut = GetColorTransferFunction(fieldName)
  if fieldName == 'PR':
    if pr1max:
      lut.RGBPoints = [0.0/6, 0, 0, 1, \
                       1.0/6, 0, 1, 1, \
                       2.0/6, 0, 1, 1, \
                       3.0/6, 0, 1, 0, \
                       4.0/6, 1, 1, 0, \
                       6.0/6, 1, 0, 0  ]
    else:
      lut.RGBPoints = [0.0/7, 0, 0, 1, \
                       1.0/7, 0, 1, 1, \
                       3.0/7, 0, 1, 0, \
                       4.0/7, 1, 1, 0, \
                       6.0/7, 1, 0, 0, \
                       6.5/7, 1, 0, 1, \
                       7.0/7, 0.5, 0.5, 0.5 ]
  else:
    lut.RGBPoints = [0.00, 0, 0, 1, \
                     0.25, 0, 1, 1, \
                     0.50, 0, 1, 0, \
                     0.75, 1, 1, 0, \
                     1.00, 1, 0, 0]
  lut.ColorSpace = 'RGB'
  lut.NanColor = [0.5, 0.5, 0.5]
  lut.BelowRangeColor = [0.0, 0.0, 0.0]
  lut.AboveRangeColor = [1.0, 1.0, 1.0]
  lut.UseBelowRangeColor = 0
  lut.UseAboveRangeColor = 0
  # For some reason you cannot toggle whether to use a NaN color; it must use it no matter what
  lut.ScalarRangeInitialized = 1
  lut.RescaleTransferFunction(minVal, maxVal)

  pwf = GetOpacityTransferFunction(fieldName)
  pwf.Points = [0.0, 0.0, 0.5, 0.0, \
                1.0, 1.0, 0.5, 0.0]
  pwf.ScalarRangeInitialized = 1
  pwf.RescaleTransferFunction(minVal, maxVal)

def show_colorbar(display, renderView, fieldName, numberFormat, numTicks=5):
  lut = GetColorTransferFunction(fieldName)
  bar = GetScalarBar(lut, renderView)
  bar.Title = fieldName
  bar.ComponentTitle = ''
  bar.TitleColor = [0, 0, 0]
  bar.LabelColor = [0, 0, 0]
  bar.TitleFontFamily = 'Courier'
  bar.LabelFontFamily = 'Courier'
  bar.AutomaticLabelFormat = 0
  bar.LabelFormat = numberFormat 
  bar.RangeLabelFormat = numberFormat
  #this disappeared in PV541 #bar.NumberOfLabels = numTicks
  bar.Visibility = True

def position_colorbar(display, renderView, fieldName, position, size, orient='Vertical'):
  lut = GetColorTransferFunction(fieldName)
  bar = GetScalarBar(lut, renderView)
  bar.WindowLocation = 'AnyLocation'
  bar.Position = position
  bar.ScalarBarLength = size
  bar.ScalarBarThickness = 40
  bar.Orientation = orient

# Not really applicable to cylinder. Don't call this function. It would be nice to save inputs... maybe just run number in name is most feasible. 
def construct_detail_text(idstr, mdotMain, mdotUB, mdotLB, relUB, relLB, percentPR):
  s  = '      Run ID: ' + idstr + '\n'
  s += '   PR at AIP: {0:.2f}%\n'.format(percentPR)
  s += '   Main Duct: {0:.5f} kg/s\n'.format(mdotMain)
  s += 'Upper Blower: {0:.5f} kg/s ({1:.2f}%)\n'.format(mdotUB, relUB)
  s += 'Lower Blower: {0:.5f} kg/s ({1:.2f}%)\n'.format(mdotLB, relLB)
  return s

def set_text_detail(textSource, text, renderView):
  textSource.Text = text
  display = Show(textSource, renderView)
  display.Color = [0, 0, 0]
  display.FontFamily = 'Courier'
  display.FontSize = 5
  display.WindowLocation = 'UpperLeftCorner'

def set_text_header(textSource, text, renderView, location='UpperCenter', fontsize=16):
  textSource.Text = text
  display = Show(textSource, renderView)
  display.Color = [0, 0, 0]
  display.FontFamily = 'Courier'
  display.FontSize = fontsize
  display.WindowLocation = location

def saveImgCylinder(dataDir, rv, vol, aip, slicex, slicexol, textDetail, text, textHeader, pr1max):
  """Save image of data at AIP.
  
  Keyword arguments:
  dataDir -- directory to save data
  rv -- the RenderView to use
  vol -- the clip volume that gives context to the AIP
  aip -- the slice at the AIP
  Not sure that I need aip or vol... 
  """
  start_over_visualization()
  set_view(rv, 'cylinder')

  #set_text_detail(textDetail, text, rv)

  d = Show(vol, rv)
  set_solid_color(vol, rv)

  #d = Show(aip, rv)
  set_flat_lighting(aip, rv)
  if pr1max:
    details = (('PR', (0.7, 1), 7, '%.2f'), ('Mach', (0.0, 1.0), 6, '%.1f'))
  else:
    details = (('PR', (0.7, 1.05), 8, '%.2f'), ('Mach', (0.0, 1.0), 6, '%.1f'))
  i = 0
  for x, xol in zip(slicex, slicexol):
    i += 1
    print '    x/L = {0:+.3f}'.format(xol)
    for field, bounds, ticks, strFmt in details:
      set_text_header(textHeader, 'x/L = {0:+.3f}'.format(xol), rv, 'LowerLeftCorner', fontsize=14)
      xstr = format_coord_for_filename(xol)
      vol.ClipType.Origin = [x - tinyOffset, 0, 0]
      aip.SliceType.Origin = [x, 0, 0]
      ColorBy(d, ('POINTS', field))
      set_rgb_colormap(field, minVal=bounds[0], maxVal=bounds[1], pr1max=pr1max)
      show_colorbar(d, rv, field, numberFormat=strFmt, numTicks=ticks)
      position_colorbar(d, rv, field, [0.9, 0.05], 0.4)
      HideUnusedScalarBars()
      fname = data_path(dataDir, 'spanslice_{0}_slice{1}_{2}.png'.format(field,i,xstr))
      SaveScreenshot(fname, view=rv)

def saveImgDuct(dataDir, rv, vol, side, slicez, slicezod, tinyOffset, textDetail, text, textHeader, pr1max):
  """Save image of data along a z-normal slice along the duct.
  
  Keyword arguments:
  dataDir -- directory to save data
  rv -- the RenderView to use
  vol -- the clip volume that gives context to the side slice 
  side -- the slice of the side
  slicez -- iterable of slice z-coordinates
  slicezod -- iterable of slice z/D-coordinates
  tinyOffset -- size of gap between reference volume and slice with data
  """
  start_over_visualization()

  set_text_detail(textDetail, text, rv)

  d = Show(vol, rv)
  set_solid_color(vol, rv)

  d = Show(side, rv)
  set_flat_lighting(side, rv)

  if pr1max:
    details = (('PR', (0.7, 1), 7, '%.2f'), ('Mach', (0.0, 1.0), 6, '%.1f'))
  else:
    details = (('PR', (0.7, 1.05), 8, '%.2f'), ('Mach', (0.0, 1.0), 6, '%.1f'))
  for view in ('duct', 'ductiso'):
    print '  {0}'.format(view)
    set_view(rv, view)
    i = 0
    for z, zod in zip(slicez, slicezod):
      i += 1
      print '    z/D = {0:+.3f}'.format(zod)
      set_text_header(textHeader, 'z/D = {0:+.3f}'.format(zod), rv)
      zstr = format_coord_for_filename(zod)
      vol.ClipType.Origin = [0, 0, z - tinyOffset]
      side.SliceType.Origin = [0, 0, z]
      for field, bounds, ticks, strFmt in details:
        ColorBy(d, ('POINTS', field))
        set_rgb_colormap(field, minVal=bounds[0], maxVal=bounds[1], pr1max=pr1max)
        show_colorbar(d, rv, field, numberFormat=strFmt, numTicks=ticks)
        position_colorbar(d, rv, field, [0.65, 0.03], 0.3, 'Horizontal')
        HideUnusedScalarBars()
        fname = data_path(dataDir, '{0}_{1}_slice{2}_{3}.png'.format(view,field,i,zstr))
        SaveScreenshot(fname, view=rv)

def main(runID, phtPath, overwrite, overrideDir, pr1max):

  ###############################################################################
  # Process input
  ###############################################################################
  
  prog_header('Starting post-processing script')
  
  # disable automatic camera reset on 'Show'
  paraview.simple._DisableFirstRenderCameraReset()
  
  if overrideDir:
    print 'User requested output files written to', overrideDir
    if not os.path.exists(overrideDir):
      err('Overridden directory to write files does not exist at ' + overrideDir)
  
  if not os.path.exists(phtPath):
    err('ParaView file does not exist at ' + phtPath)
  
  ###############################################################################
  # Load plugin
  ###############################################################################
  
  prog_header('Loading plugins and files')
  
  #pluginPath = '/projects/TRAssembly_2/paraview/ParaViewSyncIOReaderPlugin/build-v5.4.1/libPhastaSyncIOReader.so'
  pluginPath = '/home/skinnerr/ParaView/ParaViewSyncIOReaderPlugin/build/libPhastaSyncIOReader.so'
  print 'PHASTA SyncIO reader plugin path:', pluginPath
  print '  Loading local version...',
  LoadPlugin(pluginPath, remote=False, ns=globals())
  print 'done'
  print '  Loading remote version...',
  LoadPlugin(pluginPath, remote=True, ns=globals())
  print 'done'
  
  ###############################################################################
  # Load file and check it contains the expected fields
  ###############################################################################
  
  print 'SyncIO file path:', phtPath
  print '  Loading...',
  phtFile = PhastaSyncIOReader(FileName=phtPath)
  UpdatePipeline() # The phtFile has no data fields/keys until we update the pipeline
  print 'done'
  
  print 'Checking point data...',
  check_fields(phtFile)
  print 'done'
  
  ###############################################################################
  # Create directory to save files, or inform user if already exists
  ###############################################################################
  
  # Directory to save processed data/pictures in same directory as pht file
  if overrideDir:
    dataDir = os.path.join(overrideDir, 'pvdata_' + runID)
  else:
    dataDir = os.path.join(os.path.dirname(phtPath), 'pvdata_' + runID)
  
  if not os.path.exists(dataDir):
    os.makedirs(dataDir)
  else:
    print 'The data directory already exists...', dataDir
    inp = ''
    if not overwrite:
      while not inp in ('y', 'n'):
        inp = raw_input('  Do you wish to overwrite? (y/n) ')
      if inp == 'n':
        err('User chose to abort instead of overwriting files')
    print '  Continuing while overwriting files'
  
  print 'Saving data in', dataDir
  
  ###############################################################################
  # Set up view for saving images
  ###############################################################################
  
  prog_header('Building pipeline')
  
  prog('Creating render view')
  rv = GetActiveViewOrCreate('RenderView')
  
  ###############################################################################
  # Build the entire data processing pipeline
  ###############################################################################
  
  ###
  # Merge blocks for ease of saving data later; all filters start from mergedBlocks now
  ###
  
  prog('Merging blocks')
  mergedBlocks = MergeBlocks(Input=phtFile)
  
  change interaction mode for render view
  rv.InteractionMode = '2D' 
 
  # View Temperature 
  # set scalar coloring
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
  
  
  
  
  
  
  ###############################################################################
  # Save images
  ###############################################################################
  
  prog_header('Saving images')
  
  prog('Setting up text sources...')
  textDetail = Text()
  textHeader = Text()
  
  prog('cylinder...', hasMore=False)
  saveImgAIP(dataDir, rv, refClipAIP, sliceAIP, slicex, slicexol, textDetail, detailText, textHeader, pr1max)
  print 'done'
  
  print '  ...done'
  
  prog_header('All done!')

###############################################################################
# Parse arguments and call main()
###############################################################################
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Load ParaView .phts file for duct, and save extracted data and images.',
                                   usage='pvpython %(prog)s')
  parser.add_argument('runID', type=str,
                      help='string identifier for the run being visualized')
  parser.add_argument('phtPath', type=str,
                      help='path to .phts file')
  parser.add_argument('--overwrite', action='store_true',
                      help='force data to be overwriten if it already exists')
  parser.add_argument('--directory', type=str,
                      help='specify directory to save data (otherwise defaults to location of phts file)')
  parser.add_argument('--pr1max', action='store_true',
                      help='set pressure recovery colorbar maximum to 1')
  args = parser.parse_args()

  # Define physical and geometric parameters (all units in meters)
  Rconst = 288.294 # PHASTA's gas constant
  gconst = 1.4     # Gamma constant
  channel_width = 0.2 # Width of channel in z-direction
  slicezod = tuple([           eval(t) for t in slicezod_text])
  slicez   = tuple([duct_width*eval(t) for t in slicezod_text])
  slicexol_text = ('0', '1/10', '2/10', '3/10', '4/10', '5/10', '6/10', '7/10', '8/10', '9/10', '1', ) # X-normal slice coords
  slicexol = tuple([           eval(t) for t in slicexol_text])
  slicex   = tuple([      aipx*eval(t) for t in slicexol_text])
  tinyOffset = 0.001
  
  main(args.runID, args.phtPath, args.overwrite, args.directory, args.pr1max)

# would like to capture from -0.1 in x to positive 0.2. y -0.1 to positive 0.1
