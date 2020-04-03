#This program is a class that will be used to plot one stage
#of the G_FRAC run. Also will allow you to get a polygon intersection

from PltParser import PltParser
from DomLine import DomLine
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, MultiLineString, Polygon
from shapely.ops import polygonize, polygonize_full

from labellines import labelLine, labelLines
import os
import re


GRT_CODES = ["alm","gr","py","spss"]
COLOURS = ["red","orange","blue","green"]
GRT_LABELS = ["Almandine", "Grossular","Pyrope","Spessartine"]
MID_LINE = 2 #The position of the middle line with base 0

class PlotStage():

	def __init__(self, stageNum, fileDir):

		self.fileDir = fileDir
		self.stage = stageNum
		self.pltList = []
		self.sampleName = None
		for i in range(len(GRT_CODES)):
			suffix = "_Stage{:02d}".format(stageNum) + "_" + GRT_CODES[i] + ".plt"
			pSuff = "_Stage{:02d}".format(stageNum) + "_Phase.plt"

			for filename in os.listdir(fileDir):
				fileMatch = re.compile("[^_]+?" + suffix)
				pFileMatch = re.compile("[^_]+?" + pSuff)	

				if fileMatch.match(filename):
					self.sampleName = re.split("_", filename)[0]
					#print(filename)
					self.pltList.append(PltParser(os.path.join(fileDir,filename)))
				elif pFileMatch.match(filename):
					self.phasePlt = PltParser(os.path.join(fileDir,filename), isPhase = True)

		#Set up the isopleth intersections based on the first Plt (alm)
		if self.sampleName == None:
			return

		self.isoFig = plt.figure(figsize =(12,10))	
		figname = self.sampleName + " Stage " + str(stageNum)

		self.isoFig.suptitle(figname, fontsize = 16)

		self.isoAx = self.isoFig.add_subplot()
		self.isoAx.set_xlabel(self.pltList[0].xAx, fontsize = 14)
		self.isoAx.set_ylabel(self.pltList[0].yAx, fontsize =14)
		
		self.isoAx.set_xlim(self.pltList[0].Tmin, self.pltList[0].Tmax)
		self.isoAx.set_ylim(self.pltList[0].Pmin, self.pltList[0].Pmax)

		self.phaseFig = plt.figure(figsize =(12,10))	
		figname = self.sampleName + " Stage " + str(stageNum)

		self.phaseFig.suptitle(figname, fontsize = 16)

		self.phaseAx = self.phaseFig.add_subplot()
		self.phaseAx.set_xlabel(self.pltList[0].xAx, fontsize = 14)
		self.phaseAx.set_ylabel(self.pltList[0].yAx, fontsize =14)
		
		self.phaseAx.set_xlim(self.pltList[0].Tmin, self.pltList[0].Tmax)
		self.phaseAx.set_ylim(self.pltList[0].Pmin, self.pltList[0].Pmax)

			
	def plotIsos(self, saveDir, numError = 0):
		#Plot the isopleths and garnet in
		#numError is the amount of lines beyond the "middle line" (max of 2)
		
		#First we can plot the Garnet in curve from the first Plt
		grtIn = self.pltList[0].getLines("N")

		for line in grtIn:
			x, y = line.PTline.xy
			#print("Plotting " + line.leftSide)
			self.isoAx.plot(x, y, color = "black", marker = None, linestyle = "-", markersize = 7, linewidth = 2, label="Garnet In")
		#labelLines(self.isoAx.get_lines(), zorder = 2.5)
		
		self.midLines = []
		for i in range(len(self.pltList)):
			member = self.pltList[i]
			valNum = 0
			for j in range(len(member.domLines)):

				line = member.domLines[j]

				if j > 0:
					#Checks if the value of the isopleth has changed since the last one
					#If the value is "N" it wont work and it has reached the garnet in curve and will exit the loop
					#And iterate to the next member of pltList
					try:
						lastIso = float(member.domLines[j-1].leftSide)
						thisIso = float(line.leftSide)

						if thisIso > lastIso:
							valNum += 1
					except:
						break

				#Will plot the center line as a solid line and if the options are chosen
				#Plots +- 5% as dashed lines and +-10% as dotted lines
				if valNum >= MID_LINE - numError and valNum <= MID_LINE + numError:
					x, y = line.PTline.xy

					thisLineStyle = "-"
					if valNum == MID_LINE - 2 or valNum== MID_LINE + 2:
						thisLineStyle = ":"
					elif valNum == MID_LINE -1 or valNum == MID_LINE +1:
						thisLineStyle = "--"
					elif valNum == MID_LINE:
						self.midLines.append(line)
					#print("Plotting " + line.leftSide)
					thisLabel = GRT_LABELS[i] + " = " + line.leftSide
					self.isoAx.plot(x, y, color = COLOURS[i], marker = None, linestyle = thisLineStyle, markersize = 7, linewidth = 2, label=thisLabel)	
				
				
		self.isoAx.legend(fontsize = 14, loc = 'upper right')
		saveName = self.sampleName + "_Stage" + str(self.stage) + ".svg"
		self.isoFig.show()
		self.isoFig.savefig(os.path.join(saveDir, saveName))


	def plotPhase(self, saveDir):
		#Plots the phase diagram for this stage

		#Temporary test plotting lines
		self.phasePlt.getPolys()
		# for line in self.phasePlt.domLines:
		# 	x, y = line.PTline.xy
		# 	self.phaseAx.plot(x, y, color = "black", marker = None, linestyle = "-", markersize = 7, linewidth = 2)
		count = 0
		csvName = self.sampleName + "_Stage" + str(self.stage)+"_PhaseList.csv"
		
		try:
			csvFile = open(os.path.join(self.fileDir,csvName), 'w')
		except:
			print("Issue making file " + csvName)
			exit()
		csvFile.write("Field,Phases,Valid Polygon,Coordinates\n")
		for poly in self.phasePlt.polyList:
			# if count <13:
			poly.plotPoly(self.phaseAx, count, csvFile)
			count += 1

		for line in self.phasePlt.failedPolys:
			x, y = line.PTline.xy
			self.phaseAx.plot(x, y, color = "black", marker = None, linestyle = "-", markersize = 7, linewidth = 2)
		saveName = self.sampleName + "_Stage" + str(self.stage) + "_Phase.svg"
		csvFile.close()
		self.phaseFig.show()
		self.phaseFig.savefig(os.path.join(saveDir, saveName))


	def getIntersection(self):
		#Returns a polygon of the intersection area between the four isopleths
		#Hopefully
		pltLineStrings = []

		for line in self.midLines:
			pltLineStrings.append(line.PTline)


		