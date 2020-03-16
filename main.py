#Rethinking this
#
import sys
import easygui
import os
import matplotlib.pyplot as plt
from PltParser import PltParser
from DomLine import DomLine
from shapely.geometry import Point, LineString, MultiLineString, Polygon
from PlotStage import PlotStage
from descartes import PolygonPatch

inputPath = easygui.diropenbox("Choose the directory where the files are stored", default="/home/sabastien/Documents/Globus/Completed")

outputPath = easygui.diropenbox("Choose the directory to save the output", default = "/home/sabastien/Documents/Carleton/Domino Diagrams/G_FRAC Results")
# inputPath = "TestData/"
# outputPath = inputPath
pathNum = 0

intersectList = []
while True:
	
	thisStage = PlotStage(pathNum, inputPath)
	if  len(thisStage.pltList) >0:
		thisStage.plotIsos(outputPath)
		#thisStage.plotPhase(outputPath)
		#intersectList.append(thisStage.getIntersection())

		xAx = thisStage.pltList[0].xAx
		yAx =thisStage.pltList[0].yAx
		xMin = thisStage.pltList[0].Tmin
		xMax = thisStage.pltList[0].Tmax
		yMin = thisStage.pltList[0].Pmin
		yMax = thisStage.pltList[0].Pmax
	else:
		print("Does not contain stage " + str(pathNum))
		break
	
	pathNum += 1
# for line in myPlt.domLines:
# 	x, y  = line.PTline.xy
	
# ixFig = plt.figure(figsize =(16,10))

# figname = "Isopleth Intersections"

# ixFig.suptitle(figname, fontsize = 16)

# ixAx= ixFig.add_subplot()
# ixAx.set_xlabel(xAx, fontsize = 14)
# ixAx.set_ylabel(yAx, fontsize =14)

# ixAx.set_xlim(xMin, xMax)
# ixAx.set_ylim(yMin, yMax)

# for poly in intersectList:
# 	print(list(poly))
# 	polyPatch = PolygonPatch(poly)
# 	ixAx.add_patch(polyPatch)

# 	plt.plot(x, y, color = "black", marker = None, linestyle = "-", markersize = 7, linewidth = 2)
plt.show()