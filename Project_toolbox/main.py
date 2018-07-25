# -*- coding: utf-8 -*-
'''----------------------------------------------------------------------------------
File Name      : surfaceAdjusted.py
Author         : Mehran Ghandehari, Meg Brantley, Andrew Eaman
Organization   : University of Colorado at Boulder
Created        : Mar. 20th, 2016
Python Version : 2.7

-- Description --
In this project we want to calculate surface distance using different methods and testing in different resolutions
----------------------------------------------------------------------------------'''
#import modules
import arcpy
import numpy as np
import math
from arcpy import env
import arcpy.sa as sa

#import my modules
import surfaceAdjusted
import samples


#set the input parameters.
temp = arcpy.GetParameterAsText(0)
listDEM = temp.split(';')
arcpy.AddMessage("\nselected data " + listDEM[0] + "\n")
for i in range(len(listDEM)):
    listDEM [i] = sa.Raster(listDEM[i])

#listDEM = [sa.Raster('dem10m'), sa.Raster('dem30m'), sa.Raster('dem100m'), sa.Raster('dem1000m'), sa.Raster('dem5000m')] # list of DEMs of different resolutions
transects = arcpy.GetParameterAsText(1) # this is a polyline including a set of transects
sampleSize = arcpy.GetParameterAsText(2) # transect would be converted to sample points based upon the sample size

#import arc license and set workspace
path = arcpy.GetParameterAsText(3)
env.workspace = path
env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")

#Deleting extra fields in the transects feature class
# There should be only three fields in the attribute table of the transects feature class (FID, SHAPE, and Id)
# Id field is a unique number for each transects starting from one
try:
    discard = []
    for field in [f.name for f in arcpy.ListFields(transects)if f.type <> 'Geometry']:
        if field == 'FID' or field == 'Id':
            pass
        else:
            discard.append(field)
    arcpy.DeleteField_management(transects, discard)
except:
    arcpy.GetMessages(2)

# List of diffrent resolutions, DEM left lower points and methods
rsolutions = [int(i.meanCellWidth) for i in listDEM]
llpnts = [i.extent.lowerLeft for i in listDEM]

# list of attributes that would be added to the attribute table of the transects feature class
attributes = ['p2p', 'arcBil', 'arcTIN', 'arcNN', 'arcCls', 'clos', 'weiAvr', 'biLin', 'biQua', 'biQub']

# Get the number of transects in our feature class
result = arcpy.GetCount_management(transects)
numTransects = int(result.getOutput(0))

#Add new feilds to the the transects feature class based on all combination of resolutions and methods
fields = ['SHAPE@', 'Id']
for res in rsolutions:
    for att in attributes:
        arcpy.AddField_management(transects, att + str(res), 'FLOAT')
        fields = fields + [att + str(res)]

with arcpy.da.UpdateCursor(transects, (fields)) as uCurs:
    for row in uCurs:
        i = 0
        Id = row[1]

        # create sample points
        samplePnts = samples.generateSamples(transects, Id, sampleSize)

        for res in rsolutions:
            #'p2p' method
            dist = surfaceAdjusted.p2pFn(listDEM[rsolutions.index(res)], res, transects,Id)
            row[i + 2] = dist

            #'arc' methods
            dist = surfaceAdjusted.arcFns(listDEM[rsolutions.index(res)], res, transects,Id,sampleSize)
            row[i + 3:i + 7] = dist

            #'sampled' methods
            dist = surfaceAdjusted.sampledFuns(listDEM[rsolutions.index(res)], res, llpnts[rsolutions.index(res)], transects,Id ,sampleSize, samplePnts)
            row[i + 7:i + 12] = dist

            i+=10
                
        uCurs.updateRow(row)


            







