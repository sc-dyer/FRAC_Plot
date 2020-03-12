#Rethinking this
#
import sys
import easygui
import os
import matplotlib.pyplot as plt
from PltParser import PltParser
from DomLine import DomLine
from shapely.geometry import Point, LineString, MultiLineString
from PlotStage import PlotStage
trialPath = "TestData/"

#myPlt = PltParser(trialPath)

firstStage = PlotStage(0, trialPath)

# for line in myPlt.domLines:
# 	x, y  = line.PTline.xy
	
	
# 	plt.plot(x, y, color = "black", marker = None, linestyle = "-", markersize = 7, linewidth = 2)
# plt.show()