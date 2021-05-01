# run with /home/skinnerr/ParaView_Python/pvbatch_pv[0-9]+.sh

from __future__ import division # force floating-point division

import sys
import os
import math
import argparse
from paraview.simple import *

def prog(s, hasMore=False):
  """Print update of progress."""
  progString = 'Progress:'
  if hasMore:
    print progString, s,
  else:
    print progString, s
  
def prog_header(s):
  """Print update of progress with embellishment."""
  minlen = 60
  s = '*** ' + s.center(minlen) + ' ***'
  l = len(s)
  print
  print l*'*'
  print s
  print l*'*'

def err(s):
  """Print error message and quit."""
  print 'Error:', s
  print 'ABORTED python script; press Ctrl+C to exit MPI'
  sys.exit()

def data_path(dataDir, filename):
  return os.path.join(dataDir, filename)

def format_coord_for_filename(coord):
  return '{0:.3f}'.format(coord).replace('.','p')

def check_fields(phtFile):
  """Ensure we have access to the right point data files in the pht file"""
  requiredPointData = {
    'dwal': 'distance to wall',
    'u': 'vector of velocity components',
    'p': 'pressure',
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

def set_view(rv, idstr):
  """Set camera properties.
  
  Keyword arguments:
  rv -- the RenderView to modify
  idstr -- the string specifying the desired view
  """
  if idstr == 'aip':
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

def saveImgAIP(dataDir, rv, vol, aip, slicex, slicexol, textDetail, text, textHeader, pr1max):
  """Save image of data at AIP.
  
  Keyword arguments:
  dataDir -- directory to save data
  rv -- the RenderView to use
  vol -- the clip volume that gives context to the AIP
  aip -- the slice at the AIP
  """
  start_over_visualization()
  set_view(rv, 'aip')

  set_text_detail(textDetail, text, rv)

  d = Show(vol, rv)
  set_solid_color(vol, rv)

  d = Show(aip, rv)
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
  
  ###
  # Create Kiel probe
  #   Intent: save as reference probe and use to calculate PR
  ###
  
  prog('Creating Kiel probe', hasMore=True)
  
  kiel = ProbeLocation(Input=mergedBlocks)
  kiel.ProbeType.Center = [kielx, 0, 0]
  
  tmp = kiel.GetPointDataInformation().GetArray('u')
  refumag = math.sqrt(tmp.GetRange(0)[0]**2 + tmp.GetRange(1)[0]**2 + tmp.GetRange(2)[0]**2)
  refp = kiel.GetPointDataInformation().GetArray('p').GetRange()[0]
  refT = kiel.GetPointDataInformation().GetArray('T').GetRange()[0]
  
  # Reference total pressure at Kiel probe
  refp0 = refp * (1 + refumag**2 / (2 * Rconst * refT))
  print '(Reference total pressure = {0} Pa)'.format(refp0)
  
  # Reference Mach number at Kiel probe
  refM = refumag / math.sqrt(gconst * Rconst * refT)
  print '(Reference Mach = {0})'.format(refM)
  
  ###
  # Remove contraction section and downstream duct for all subsequent processing 
  ###
  
  prog('Clipping away contraction section and downstream duct')
  
  rootClip = Clip(Input=mergedBlocks) # A handful of filters will fork off this 'root' clip
  rootClip.ClipType = 'Plane'
  rootClip.ClipType.Origin = [2*aipx, 0, 0]
  rootClip.ClipType.Normal = [-1, 0, 0]
  
  ###
  # Create clips of the full geometry to give images context
  ###
  
  prog('Creating reference volume clips')
  
  refClipDuct = Clip(Input=rootClip)
  refClipDuct.ClipType = 'Plane'
  refClipDuct.ClipType.Origin = [0, 0, 0]
  refClipDuct.ClipType.Normal = [0, 0, -1]
  
  clipContractionSection = Clip(Input=rootClip)
  clipContractionSection.ClipType = 'Plane'
  clipContractionSection.ClipType.Origin = [-aipx, 0, 0]
  clipContractionSection.ClipType.Normal = [1, 0, 0]
  
  refClipAIP = Clip(Input=clipContractionSection)
  refClipAIP.ClipType = 'Plane'
  refClipAIP.ClipType.Origin = [aipx - tinyOffset, 0, 0]
  refClipAIP.ClipType.Normal = [-1, 0, 0]
  
  ###
  # Create all the calculators
  ###
  
  prog('Creating calculators')
  
  calc1 = Calculator(Input=rootClip)
  calc1.ResultArrayName = 'rho'
  calc1.Function = 'p/({0}*T)'.format(Rconst)
  
  calc2 = Calculator(Input=calc1)
  calc2.ResultArrayName = 'PR'
  calc2.Function = 'p*(1+mag(u)^2/(2*{0}*T)) / {1}'.format(Rconst, refp0)
  
  rootCalc = Calculator(Input=calc2)
  rootCalc.ResultArrayName = 'Mach'
  rootCalc.Function = 'mag(u)/sqrt({0}*{1}*T)'.format(gconst, Rconst)
  
  ###
  # Duct span-normal slice
  #   Intent: iterate over clipDuct origin and save images
  ###
  
  prog('Creating span-normal duct slice')
  
  sliceDuct = Slice(Input=rootCalc)
  sliceDuct.SliceType.Origin = [0, 0, 0]
  sliceDuct.SliceType.Normal = [0, 0, 1]
  
  ###
  # Duct upper/lower surface data
  #   Intent: iterate over sliceDuctSurf origin and save surfTop/surfBot
  ###
  
  prog('Creating upper/lower surface data slices')
  
  wallCont = Contour(Input=rootCalc)
  wallCont.ContourBy = ['POINTS', 'dwal']
  wallCont.Isosurfaces = [1e-10]
  
  sliceDuctSurf = Slice(Input=wallCont)
  sliceDuctSurf.SliceType = 'Plane'
  sliceDuctSurf.SliceType.Origin = [0, 0, 0]
  sliceDuctSurf.SliceType.Normal = [0, 0, 1]
  
  surfTop = Clip(Input=sliceDuctSurf)
  surfTop.ClipType = 'Plane'
  surfTop.ClipType.Origin = [0, 0, 0]
  surfTop.ClipType.Normal = [0, 1, 0]
  
  surfBot = Clip(Input=sliceDuctSurf)
  surfBot.ClipType = 'Plane'
  surfBot.ClipType.Origin = [0, 0, 0]
  surfBot.ClipType.Normal = [0, -1, 0]
  
  ###
  # AIP slice
  #   Intent: save as image
  ###
  
  prog('Creating AIP slice')
  
  sliceAIP = Slice(Input=rootCalc)
  sliceAIP.SliceType = 'Plane'
  sliceAIP.SliceType.Origin = [aipx, 0, 0]
  sliceAIP.SliceType.Normal = [1, 0, 0]
  
  ###
  # AIP line data
  #   Intent: iterate over lineAIP origin and save line
  ###
  
  prog('Creating AIP line data slices')
  
  lineAIP = Slice(Input=sliceAIP)
  lineAIP.SliceType = 'Plane'
  lineAIP.SliceType.Origin = [0, 0, 0]
  lineAIP.SliceType.Normal = [0, 0, 1]
  
  ###############################################################################
  # Calculate mass flow
  ###############################################################################
  
  prog_header('Calculating mass flow and integrated PR')
  
  prog('Main duct')
  
  sliceMain = Slice(Input=rootCalc)
  sliceMain.SliceType = 'Plane'
  sliceMain.SliceType.Origin = [-0.2, 0, 0]
  sliceMain.SliceType.Normal = [1, 0, 0]
  
  calc4 = Calculator(Input=sliceMain)
  calc4.ResultArrayName = 'MassFluxX'
  calc4.Function = 'rho*u_X'
  
  IV1 = IntegrateVariables(Input=calc4)
  massDuct = IV1.GetPointDataInformation().GetArray('MassFluxX').GetRange()[0]
  
  prog('Upper blower')
  
  # Clip away contraction section
  clipContractionMassFlow = Clip(Input=rootCalc)
  clipContractionMassFlow.ClipType = 'Plane'
  clipContractionMassFlow.ClipType.Origin = [-0.3, 0, 0]
  clipContractionMassFlow.ClipType.Normal = [1, 0, 0]
  
  # Clip away ramp section
  clip2 = Clip(Input=clipContractionMassFlow)
  clip2.ClipType = 'Plane'
  clip2.ClipType.Origin = [0.09, 0, 0]
  clip2.ClipType.Normal = [-1, 0, 0]
  
  slice1 = Slice(Input=clip2)
  slice1.SliceType = 'Plane'
  slice1.SliceType.Origin = [0, 0.1, 0]
  slice1.SliceType.Normal = [0, 1, 0]
  
  calc5 = Calculator(Input=slice1)
  calc5.ResultArrayName = 'MassFluxY'
  calc5.Function = '-rho*u_Y'
  
  IV2 = IntegrateVariables(Input=calc5)
  massUB = IV2.GetPointDataInformation().GetArray('MassFluxY').GetRange()[0]
  massUBpercent = 100 * massUB / massDuct
  
  prog('Lower blower')
  
  slice2 = Slice(Input=clipContractionMassFlow)
  slice2.SliceType = 'Plane'
  slice2.SliceType.Origin = [0, -0.08, 0]
  slice2.SliceType.Normal = [0, 1, 0]
  
  calc6 = Calculator(Input=slice2)
  calc6.ResultArrayName = 'MassFluxY'
  calc6.Function = 'rho*u_Y'
  
  IV3 = IntegrateVariables(Input=calc6)
  massLB = IV3.GetPointDataInformation().GetArray('MassFluxY').GetRange()[0]
  massLBpercent = 100 * massLB / massDuct
  
  prog('AIP')
  
  calc7 = Calculator(Input=sliceAIP)
  calc7.ResultArrayName = 'Unity'
  calc7.Function = '1'
  
  IV4 = IntegrateVariables(Input=calc7)
  areaAIP = IV4.GetPointDataInformation().GetArray('Unity').GetRange()[0]
  integratedPRatAIP= IV4.GetPointDataInformation().GetArray('PR').GetRange()[0]
  bestPR = 1 * areaAIP
  percentPR = 100 * integratedPRatAIP / bestPR
  
  detailText = construct_detail_text(runID, massDuct, massUB, massLB, massUBpercent, massLBpercent, percentPR)
  print detailText
  
  ###############################################################################
  # Save data
  ###############################################################################
  
  prog_header('Saving data')
  
  prog('Saving reference probe data...', hasMore=True)
  SaveData(data_path(dataDir, 'ref.csv'), proxy=kiel)
  print 'done'
  
  for z, zod in zip(slicez, slicezod):
    zstr = format_coord_for_filename(zod)
  
    prog('Saving z/D = {0:.3f} surface data...'.format(zod), hasMore=True)
    sliceDuctSurf.SliceType.Origin = [0, 0, z]
    print 'top',
    SaveData(data_path(dataDir, 'surfTop_zoverD_' + zstr + '.csv'), proxy=surfTop)
    print 'bottom',
    SaveData(data_path(dataDir, 'surfBot_zoverD_' + zstr + '.csv'), proxy=surfBot)
    print 'done'
  
    prog('Saving z/D = {0:.3f} AIP line data...'.format(zod), hasMore=True)
    lineAIP.SliceType.Origin = [0, 0, z]
    SaveData(data_path(dataDir, 'lineAIP_zoverD_' + zstr + '.csv'), proxy=lineAIP)
    print 'done'
  
  ###############################################################################
  # Save images
  ###############################################################################
  
  prog_header('Saving images')
  
  prog('Setting up text sources...')
  textDetail = Text()
  textHeader = Text()
  
  prog('AIP...', hasMore=False)
  saveImgAIP(dataDir, rv, refClipAIP, sliceAIP, slicex, slicexol, textDetail, detailText, textHeader, pr1max)
  print 'done'
  
  prog('Duct...', hasMore=False)
  saveImgDuct(dataDir, rv, refClipDuct, sliceDuct, slicez, slicezod, tinyOffset, textDetail, detailText, textHeader, pr1max)
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
  #kielx = -0.2159     # Kiel probe x-coordinate
  kielx = -2     # Kiel probe x-coordinate
  aipx = 0.217341856  # AIP x-coordinate
  aipy = 0.0394081    # Y-coordinate at center of AIP
  rampx = aipx / 2    # X-coordinate at center of ramp
  duct_width = 0.1524 # Width of duct in z-direction
  slicezod_text = ('-1/3', '-1/4', '0', '1/4', '1/3') # Z-normal slice coords
  slicezod = tuple([           eval(t) for t in slicezod_text])
  slicez   = tuple([duct_width*eval(t) for t in slicezod_text])
  slicexol_text = ('0', '1/10', '2/10', '3/10', '4/10', '5/10', '6/10', '7/10', '8/10', '9/10', '1', ) # X-normal slice coords
  slicexol = tuple([           eval(t) for t in slicexol_text])
  slicex   = tuple([      aipx*eval(t) for t in slicexol_text])
  tinyOffset = 0.001
  
  main(args.runID, args.phtPath, args.overwrite, args.directory, args.pr1max)
