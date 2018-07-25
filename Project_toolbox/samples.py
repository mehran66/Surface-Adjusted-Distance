#import arc license and set workspace
import arcpy

# this finction creates ragular sample points along the input polyline based on the sample size

def generateSamples(transects, Id, sampleSize):
    #ref: https://geonet.esri.com/thread/78492
    pts = []
    expression = '"Id" = {0:d}'.format(Id)
    with arcpy.da.SearchCursor(transects,['SHAPE@', 'Id'], expression) as rows:
        row = rows.next()
        pts.append(row[0].positionAlongLine(0))
        i = sampleSize
        while i < row[0].length:
            pts.append(row[0].positionAlongLine(i))
            i += sampleSize
        leng = row[0].length
        pts.append(row[0].positionAlongLine(leng))
    return pts
