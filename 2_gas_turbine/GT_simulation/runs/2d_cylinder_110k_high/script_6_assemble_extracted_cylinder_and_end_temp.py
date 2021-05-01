#!/usr/bin/python

import sys
import os

#pip install wheel
#pip install pandas

import pandas as pd
import numpy as np
#from paraview.simple import *

# Step through all run directories, load recircData, extract numbers and save to file. 

###############################################################################
# Load data
###############################################################################

meshNum=110

directorIncrement = 1
numStart = 1
numRunDirectories = 200
numRunTimes = 1
timeStart=1000
timeStep=50
# enter points manually first: then find way to do automatically, ie check one file? maybe. 

numPointsMid=128
numPointsEnd=35
numPointsCylinder=208

#rootdir = '~/PHASTA_Sensitivity/runs/2D_geom14p51_174k/16-1-Chef/runA-IC-AFC-Scan_t'

# Initialise array to save results to
#assembledData = np.zeros((numRunDirectories, 2))
#assembledRunMid = np.zeros((numRunTimes,numPointsMid))
#assembledRunEnd = np.zeros((numRunTimes,numPointsEnd))
#assembledRunCylinder = np.zeros((numRunTimes,numPointsCylinder))

assembledRunMid = np.zeros((numRunDirectories,numPointsMid))
assembledRunEnd = np.zeros((numRunDirectories,numPointsEnd))
assembledRunCylinder = np.zeros((numRunDirectories,numPointsCylinder))
# Loop over run directories

for i_runs in range(numRunDirectories):
    runNum_1 = str(numStart+i_runs)
    runNum = runNum_1.zfill(3)
    
    print 'Extracting results from run:' ,runNum
    directory = 'Run-'+runNum
    os.chdir(directory)
    

    for i_times in range(numRunTimes):
        timeStepExtract=timeStart-i_times*timeStep
#        print 'Time:' ,i_times

        file_name_cylinder='cylinder_temp_%d.csv' % timeStepExtract
        file_name_end='end_temp_%d.csv' % timeStepExtract
        file_name_mid='mid_temp_%d.csv' % timeStepExtract
        
        df_cylinder = pd.read_csv(file_name_cylinder, sep=',') 
    
        data_cylinder= df_cylinder['temp']
        cylinder_coords_x = df_cylinder['Points:0']
        cylinder_coords_y = df_cylinder['Points:1']

        df_end = pd.read_csv(file_name_end, sep=',') 
    
        data_end= df_end['temp']
        end_coords_x = df_end['Points:0']
        end_coords_y = df_end['Points:1']

       
        df_mid = pd.read_csv(file_name_mid, sep=',') 
    
        data_mid= df_mid['temp']
        mid_coords_x = df_mid['Points:0']
        mid_coords_y = df_mid['Points:1']

        #assembledRunCylinder[i_times,:]=data_cylinder
        #assembledRunMid[i_times,:]=data_mid
        #assembledRunEnd[i_times,:]=data_end
        
    assembledRunCylinder[i_runs,:]=data_cylinder
    assembledRunEnd[i_runs,:]=data_end
    assembledRunMid[i_runs,:]=data_mid

    #np.savetxt('assembledRunCylinder_%d' %meshNum,assembledRunCylinder)
    #np.savetxt('assembledRunEnd_%d' %meshNum,assembledRunEnd)
    #np.savetxt('assembledRunMid_%d' %meshNum,assembledRunMid)
         
    os.chdir('..')


np.savetxt('assembledRunCylinder_%d' %meshNum,assembledRunCylinder)
np.savetxt('assembledRunEnd_%d' %meshNum,assembledRunEnd)
np.savetxt('assembledRunMid_%d' %meshNum,assembledRunMid)

assembledCoordsCylinder= np.zeros((2,numPointsCylinder))
assembledCoordsCylinder[0,:]=cylinder_coords_x
assembledCoordsCylinder[1,:]=cylinder_coords_y

assembledCoordsEnd= np.zeros((2,numPointsEnd))
assembledCoordsEnd[0,:]=end_coords_x
assembledCoordsEnd[1,:]=end_coords_y

assembledCoordsMid= np.zeros((2,numPointsMid))
assembledCoordsMid[0,:]=mid_coords_x
assembledCoordsMid[1,:]=mid_coords_y

np.savetxt('assembledCoordsCylinder_%d' %meshNum, assembledCoordsCylinder)
np.savetxt('assembledCoordsEnd_%d' %meshNum, assembledCoordsEnd)
np.savetxt('assembledCoordsMid_%d' %meshNum, assembledCoordsMid)
#np.savetxt('assembledData',assembledData)

