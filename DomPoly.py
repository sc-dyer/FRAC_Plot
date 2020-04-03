#Class for handling a Domino Polygon
from shapely.geometry import Point, LineString, MultiLineString, Polygon
from shapely.ops import linemerge, snap
import matplotlib.pyplot as plt
import copy


class DomPoly:

	def __init__(self, polyCoords,name):
		#Coordinates will be defined outside of this class constructor
		self.field = Polygon(polyCoords)
		self.phases = name
		self.renamePhases()
		if not self.field.is_valid:
			print(self.phases + " not a valid Poly")


	def plotPoly(self, pltIn, index, csvWriter):

		index = index+1
		print("Plotting: " + self.phases)
		print(self.field)
		print(self.field.is_valid)
		x, y = self.field.exterior.xy
		
		textX, textY = self.field.centroid.xy
		
		textX = float(textX[0])
		textY = float(textY[0])
		isValid = "No"
		pltIn.plot(x, y, color = "black", marker = None, linestyle = "-", markersize = 7, linewidth = 2)
		if self.field.is_valid:
			pltIn.text(textX, textY, str(index), fontsize = 10, color = "blue")
			isValid = "Yes"
		delimit = ","
		fieldString = str(self.field).replace("POLYGON ((","")
		fieldString = fieldString.replace("))","")
		csvWriter.write(str(index) + delimit  + self.phases + delimit +  isValid + delimit + fieldString + "\n")
	def renamePhases(self):
		#This is a long method to rename names of phases in the string

		tPhases = self.phases.replace("FLUID3","Fluid")
		tPhases = tPhases.replace("FSP", "Fsp")
		tPhases = tPhases.replace("GARNET", "Grt")
		tPhases = tPhases.replace("ILM","Ilm")
		tPhases = tPhases.replace("BIO","Bt")
		tPhases = tPhases.replace("CHLR","Chl")
		tPhases = tPhases.replace("PHNG", "WM")
		tPhases = tPhases.replace("q","Qz")
		tPhases = tPhases.replace("gph","Gr")
		tPhases = tPhases.replace("and","And")
		tPhases = tPhases.replace("ky","Ky")
		tPhases = tPhases.replace("sill","Sil")
		tPhases = tPhases.replace("STAU","St")
		tPhases = tPhases.replace("CORD","Crd")
		tPhases = tPhases.replace("cz","Czo")
		tPhases = tPhases.replace("LIQtc","Melt")
		tPhases = tPhases.replace("ru","Rt")
		self.phases = tPhases