# -*- coding: utf-8 -*-
'''----------------------------------------------------------------------------------
File Name      : surfaceAdjusted.py
Author         : Mehran Ghandehari, Meg Brantley, Andrew Eaman
Organization   : University of Colorado at Boulder
Created        : Mar. 20th, 2016
Python Version : 2.7

-- Description --

----------------------------------------------------------------------------------'''

#import arc license and set workspace
import arcpy
import numpy as np
import math
from arcpy import env
import arcpy.sa as sa
env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")
arcpy.CheckOutExtension('3D')


# import my modules
import polyInterpolation
import inverseDistanecWeighting

#********************************************************************
def p2pFn(dem, res, transects, Id):
    '''
        Using a given DEM, resolution, Id and transects this function
        calculates the 3D pixel to pixel distance.
    '''
    expression = '"Id" = {0:d}'.format(Id)
    aCurs = arcpy.da.SearchCursor(transects,['SHAPE@', 'Id'], expression)
    row = aCurs.next()
    fid = str(row[1])
    lyr = arcpy.MakeFeatureLayer_management(transects, "shplyr")
    selectfeat = arcpy.SelectLayerByAttribute_management(lyr, "NEW_SELECTION", '"Id" = {}'.format(fid))
    mask = arcpy.sa.ExtractByMask(dem,selectfeat)
    output = "masks"+str(fid)+".shp"
    arcpy.RasterToPoint_conversion(mask, output)
    del aCurs

    pCurs = arcpy.da.SearchCursor(output, ["POINTID","SHAPE@","GRID_CODE"])
    coords = []
    Xs = []
    Ys = []
    Zs = []
    #loop through the shape field to get X, Y and Z values
    for i in pCurs:
        coords.append(i[1].getPart(0))
        Zs.append(i[2])
    pixCount = len(Zs)
    #loop through the coordinate values to get start and end points
    for i in coords:
        Xs.append(i.X)
        Ys.append(i.Y)
    xstart = Xs[0]
    xend = Xs[-1]
    ystart = Ys[0]
    yend = Ys[-1]
    #find the slope of the line (m)
    if xend!=xstart:
        m = (yend-ystart) / (xend-xstart)
    else:
        m = 0

    if  m > 0:
        #if the slope is posiive, re-order the z values from the
            #default ID order
        XYZ = zip(Xs,Ys,Zs)

        t = []
        temp = [Zs[0]]
        for i in range(len(Xs)-1):
            if Ys[i]==Ys[i+1]:
                temp.append(Zs[i+1])
            else:
                temp.reverse()
                t = t + temp
                temp = [Zs[i+1]]


        vals = zip(t,t[1:])
        c = []
        #calculate the difference in the z values
        for i in vals:
            c.append((i[0]-i[1])**2)
        d = []
        for i in c:
            d.append((res**2 + i)**0.5)
        #calculate the euclidean distance by adding the sum of the elevation distance
            #and the pixel to pixel distance
        dist = (sum(d))

    else:
        vals = zip(Zs, Zs[1:])
        c = []
        for i in vals:
            c.append((i[0]-i[1])**2)
        d = []
        for i in c:
            d.append((res**2 + i)**0.5)
        dist = (sum(d))
    del pCurs

    return dist
#********************************************************************

def arcFns(dem, res, transects, Id ,sampleSize):
#Use DEM to add surface information to transects using bilinear method
#Check ID to find measurement for one transect only

  if Id == 1:
      #create a TIN file to permit the use of linear, natural neighbor, and nearest methods
      arcpy.RasterTin_3d(dem,'out_tin' +str(res))
      #bilinear
      arcpy.CopyFeatures_management(transects, 'transects'+ str(res)+'bil')
      arcpy.AddSurfaceInformation_3d('transects'+str(res)+'bil.shp',dem, 'SURFACE_LENGTH','BILINEAR', sampleSize)
      #linear
      arcpy.CopyFeatures_management(transects, 'transects'+ str(res)+'TIN')
      arcpy.AddSurfaceInformation_3d('transects'+str(res)+'TIN.shp','out_tin' +str(res), 'SURFACE_LENGTH','LINEAR', sampleSize)
      #natural neighbors
      arcpy.CopyFeatures_management(transects, 'transects'+ str(res)+'NN')
      arcpy.AddSurfaceInformation_3d('transects'+ str(res)+'NN.shp','out_tin' +str(res), 'SURFACE_LENGTH','NATURAL_NEIGHBORS', sampleSize)
      #conflate nearest
      arcpy.CopyFeatures_management(transects, 'transects'+ str(res)+'close')
      arcpy.AddSurfaceInformation_3d('transects'+ str(res)+'close.shp','out_tin' +str(res), 'SURFACE_LENGTH','CONFLATE_NEAREST', sampleSize)
      #use cursors to explore the ouputted surface distances in the attribute table
      sCursBil = arcpy.da.SearchCursor('transects'+ str(res)+'bil.shp', 'SLENGTH')
      sCursTIN = arcpy.da.SearchCursor('transects'+str(res)+'TIN.shp', 'SLENGTH')
      sCursNN = arcpy.da.SearchCursor('transects'+ str(res)+'NN.shp', 'SLENGTH')
      sCursClose = arcpy.da.SearchCursor('transects'+ str(res)+'close.shp', 'SLENGTH')

      d1 = sCursBil.next()[0]
      d2 = sCursTIN.next()[0]
      d3 = sCursNN.next()[0]
      d4 = sCursClose.next()[0]
      #delete cursors
      del sCursBil, sCursTIN, sCursNN, sCursClose
      #return surface distance values
      return [d1, d2, d3, d4]

  else:
      #for consequential transects, we only need to explore the table, as the surface information as already been added in the previous step
      expression = '"Id" = {0:d}'.format(Id)
      with arcpy.da.SearchCursor('transects'+ str(res)+'bil.shp',['SLENGTH', 'Id'], expression) as sCursBil:
          row= sCursBil.next()
          d1 = row[0]
      with arcpy.da.SearchCursor('transects'+ str(res)+'TIN.shp',['SLENGTH', 'Id'], expression) as sCursTIN:
          row = sCursTIN.next()
          d2= row[0]
      with arcpy.da.SearchCursor('transects'+ str(res)+'NN.shp',['SLENGTH', 'Id'], expression) as sCursNN:
          row = sCursNN.next()
          d3 = row[0]
      with arcpy.da.SearchCursor('transects'+ str(res)+'close.shp',['SLENGTH', 'Id'], expression) as sCursClose:
          row = sCursClose.next()
          d4 = row[0]
      return [d1, d2, d3, d4]

#********************************************************************

def sampledFuns(dem, res, llpnts, transects, Id ,sampleSize, samplePnts):

    #this arrays would be filled based on the interpolated elavtion of each sample points
    elevClosest = np.zeros(len(samplePnts))
    elevWeiAvr = np.zeros(len(samplePnts))
    elevBiLinear = np.zeros(len(samplePnts))
    elevBiQuadratic = np.zeros(len(samplePnts))
    elevBiQubic = np.zeros(len(samplePnts))

    rasterBlock_elev = np.zeros((5,5)) # this 5 by 5 matrix would contain the elevation of 16 neighbor pixels of a sample point
    rasterBlock_x = np.zeros((5,5)) # this 5 by 5 matrix would contain the x coordinate of 16 neighbor pixels of a sample point
    rasterBlock_y = np.zeros((5,5)) # this 5 by 5 matrix would contain the y coordinate of 16 neighbor pixels of a sample point

    llpntsX  = llpnts.X
    llpntsY  = llpnts.Y

    for i in range(len(samplePnts)): # in this for loop we calculate the elevation of each sample points using different methods
        # x and y coordinate of a sample point
        x = samplePnts[i].firstPoint.X
        y = samplePnts[i].firstPoint.Y

        # calculated the left lower coordinated of the cell encompass the current sample point
        cellX = (int((x - llpntsX)/ res) * res) + llpntsX
        cellY = (int((y - llpntsY)/ res) * res) + llpntsY

        # calculated the left lower coordinate of 16 cell to be used for interpolation
        llpnts16X = cellX - (2 * res)
        llpnts16Y = cellY - (2 * res)

        # Extract the elevation of a data block (we use the RasterToNumPyArray to extract only the 16 neighbor pixels of a sample point)
        lowerLeft = arcpy.Point(llpnts16X,llpnts16Y)
        rasterBlock_elev = arcpy.RasterToNumPyArray(dem,lowerLeft,5,5).astype('float')

        # calculate the coordinate of a data block
        for ii in [4,3,2,1,0]:
            dy = (abs(ii-4)+0.5) * res
            rasterBlock_y [ii,:] = llpnts16Y + dy
        for jj in [0,1,2,3,4]:
            dx = (jj+0.5) * res
            rasterBlock_x [:,jj] = llpnts16X + dx

        # nearest neghbor
        elevClosest [i] = rasterBlock_elev[2,2]

        # weighted average of 9 surrounding pixels
        xCoor = np.array([rasterBlock_x [2,1], rasterBlock_x[1,2], rasterBlock_x[2,3], rasterBlock_x[3,2], rasterBlock_x[1,2],
                 rasterBlock_x[1,3], rasterBlock_x[2,2], rasterBlock_x[3,1], rasterBlock_x[3,3]])
        yCoor = np.array([rasterBlock_y [2,1], rasterBlock_y[1,2], rasterBlock_y[2,3], rasterBlock_y[3,2], rasterBlock_y[1,2],
                 rasterBlock_y[1,3], rasterBlock_y[2,2], rasterBlock_y[3,1], rasterBlock_y[3,3]])
        elev = np.array([rasterBlock_elev [2,1], rasterBlock_elev[1,2], rasterBlock_elev[2,3], rasterBlock_elev[3,2], rasterBlock_elev[1,2],
                 rasterBlock_elev[1,3], rasterBlock_elev[2,2], rasterBlock_elev[3,1], rasterBlock_elev[3,3]])
        elevWeiAvr [i] = inverseDistanecWeighting.IDW(x, y, xCoor, yCoor, elev, 2)

        #Bilinear interpolation
        xCoor = np.array([rasterBlock_x [2,1], rasterBlock_x[1,2], rasterBlock_x[2,3], rasterBlock_x[3,2]])
        yCoor = np.array([rasterBlock_y [2,1], rasterBlock_y[1,2], rasterBlock_y[2,3], rasterBlock_y[3,2]])
        elev = np.array([rasterBlock_elev [2,1], rasterBlock_elev[1,2], rasterBlock_elev[2,3], rasterBlock_elev[3,2]])
        elevBiLinear [i] = polyInterpolation.polyval2d(x,y,polyInterpolation.polyfit2d(xCoor, yCoor, elev, 1))

        #Biquadratic interpolation
        xCoor = np.array([rasterBlock_x [2,1], rasterBlock_x[1,2], rasterBlock_x[2,3], rasterBlock_x[3,2], rasterBlock_x[1,2],
                 rasterBlock_x[1,3], rasterBlock_x[2,2], rasterBlock_x[3,1], rasterBlock_x[3,3]])
        yCoor = np.array([rasterBlock_y [2,1], rasterBlock_y[1,2], rasterBlock_y[2,3], rasterBlock_y[3,2], rasterBlock_y[1,2],
                 rasterBlock_y[1,3], rasterBlock_y[2,2], rasterBlock_y[3,1], rasterBlock_y[3,3]])
        elev = np.array([rasterBlock_elev [2,1], rasterBlock_elev[1,2], rasterBlock_elev[2,3], rasterBlock_elev[3,2], rasterBlock_elev[1,2],
                 rasterBlock_elev[1,3], rasterBlock_elev[2,2], rasterBlock_elev[3,1], rasterBlock_elev[3,3]])
        elevBiQuadratic [i] = polyInterpolation.polyval2d(x,y,polyInterpolation.polyfit2d(xCoor, yCoor, elev, 2))

        #BiQubic interpolation
        xCoor = np.array([rasterBlock_x [2,1], rasterBlock_x[1,2], rasterBlock_x[2,3], rasterBlock_x[3,2], rasterBlock_x[1,2],
                 rasterBlock_x[1,3], rasterBlock_x[3,1], rasterBlock_x[3,3], rasterBlock_x[0,0], rasterBlock_x[0,2], rasterBlock_x[0,4],
                          rasterBlock_x[2,0], rasterBlock_x[2,2], rasterBlock_x[2,4], rasterBlock_x[4,0], rasterBlock_x[4,2], rasterBlock_x[4,4]])
        yCoor = np.array([rasterBlock_y [2,1], rasterBlock_y[1,2], rasterBlock_y[2,3], rasterBlock_y[3,2], rasterBlock_y[1,2],
                 rasterBlock_y[1,3], rasterBlock_y[3,1], rasterBlock_y[3,3], rasterBlock_y[0,0], rasterBlock_y[0,2], rasterBlock_y[0,4],
                          rasterBlock_y[2,0], rasterBlock_y[2,2], rasterBlock_y[2,4], rasterBlock_y[4,0], rasterBlock_y[4,2], rasterBlock_y[4,4]])
        elev = np.array([rasterBlock_elev [2,1], rasterBlock_elev[1,2], rasterBlock_elev[2,3], rasterBlock_elev[3,2], rasterBlock_elev[1,2],
                 rasterBlock_elev[1,3], rasterBlock_elev[3,1], rasterBlock_elev[3,3], rasterBlock_elev[0,0], rasterBlock_elev[0,2], rasterBlock_elev[0,4],
                          rasterBlock_elev[2,0], rasterBlock_elev[2,2], rasterBlock_elev[2,4], rasterBlock_elev[4,0], rasterBlock_elev[4,2], rasterBlock_elev[4,4]])
        elevBiQubic [i] = polyInterpolation.polyval2d(x,y,polyInterpolation.polyfit2d(xCoor, yCoor, elev, 3))

    # After calculating the elevation of each sample points, we can calculate the surface distance
    #Closest
    d=0
    for ii in range(len(samplePnts)-1):
        d = d + math.sqrt((samplePnts[ii].firstPoint.X-samplePnts[ii+1].firstPoint.X)**2 +(samplePnts[ii].firstPoint.Y-samplePnts[ii+1].firstPoint.Y)**2 + (elevClosest[ii]-elevClosest[ii+1])**2)
    distClosest = d

    #Weighted average
    d=0
    for ii in range(len(samplePnts)-1):
        d = d + math.sqrt((samplePnts[ii].firstPoint.X-samplePnts[ii+1].firstPoint.X)**2 +(samplePnts[ii].firstPoint.Y-samplePnts[ii+1].firstPoint.Y)**2 + (elevWeiAvr[ii]-elevWeiAvr[ii+1])**2)
    distWeiAve=d

    #BiLinear
    d=0
    for ii in range(len(samplePnts)-1):
        d = d + math.sqrt((samplePnts[ii].firstPoint.X-samplePnts[ii+1].firstPoint.X)**2 +(samplePnts[ii].firstPoint.Y-samplePnts[ii+1].firstPoint.Y)**2 + (elevBiLinear[ii]-elevBiLinear[ii+1])**2)
    distBiLinear = d

    #BiQuadratic
    d=0
    for ii in range(len(samplePnts)-1):
        d = d + math.sqrt((samplePnts[ii].firstPoint.X-samplePnts[ii+1].firstPoint.X)**2 +(samplePnts[ii].firstPoint.Y-samplePnts[ii+1].firstPoint.Y)**2 + (elevBiQuadratic[ii]-elevBiQuadratic[ii+1])**2)
    distBiQuadratic = d

    #BiQubic
    d=0
    for ii in range(len(samplePnts)-1):
        d = d + math.sqrt((samplePnts[ii].firstPoint.X-samplePnts[ii+1].firstPoint.X)**2 +(samplePnts[ii].firstPoint.Y-samplePnts[ii+1].firstPoint.Y)**2 + (elevBiQubic[ii]-elevBiQubic[ii+1])**2)
    distBiQubic = d

    dist = [distClosest, distWeiAve, distBiLinear, distBiQuadratic, distBiQubic]
    return dist
#********************************************************************











