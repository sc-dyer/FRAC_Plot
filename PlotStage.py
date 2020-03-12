#This program is a class that will be used to plot one stage
#of the G_FRAC run. Also will allow you to get a polygon intersection

from PltParser import PltParser
from DomLine import DomLine
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, MultiLineString
from shapely.ops import polygonize
import os
import re

GRT_CODES = ["alm","gr","py","spss"]
COLOURS = ["red","orange","blue","green"]

class PlotStage():

	def __init__(self, stageNum, fileDir):

		pltList = []
		for i in range(len(GRT_CODES)):
			suffix = "Stage{:02d}".format(stageNum) + "_" + GRT_CODES[i] + ".plt"

			for filename in os.listdir(fileDir):
				fileMatch = re.compile("[^_]+?" + suffix)
				

				if fileMatch.match(filename):
					print(filename)
					pltList.append(PltParser(os.path.join(fileDir,filename)))

			