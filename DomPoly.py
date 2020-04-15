#Class for handling a Domino Polygon
from shapely.geometry import Point, LineString, MultiLineString, Polygon, MultiPolygon
from shapely.ops import linemerge, snap
import matplotlib.pyplot as plt
import copy

PHASES = ["Fluid","Fsp","Grt","Ilm","Bt","Chl","WM","Qz","Gr","And","Ky","Sil","St","Crd","Czo","Mrg","Cld","Rt","Melt"]
FILL_OPTIONS = ["barrovian"]
class DomPoly:

	def __init__(self, polyCoords,name):
		#Coordinates will be defined outside of this class constructor
		self.field = Polygon(polyCoords)
		self.phases = name
		self.renamePhases()
		if not self.field.is_valid:
			newField = self.field.buffer(0)
			# print(newField)
			# if newField.is_valid:
			# 	self.field = newField
			# 	print(self.phases + " is a valid Poly")
			# else:
			print(self.phases + " not a valid Poly")
		else:
			print(self.phases + " is a valid Poly")
		print("\n")


	def plotPoly(self, pltIn, index, csvWriter, areaMin = 0, fillOp = None):
		
		index = index+1
		delimit = ","
		print("Plotting: " + self.phases)
		print(self.field)
		print(self.field.is_valid)
		print(self.field.is_empty)
		if not self.field.is_empty:
			if isinstance(self.field, MultiPolygon):
				fields = list(self.field)
			elif isinstance(self.field, Polygon):
				fields = [self.field]
			for poly in fields:
				x, y = poly.exterior.xy
				
				textX, textY = poly.representative_point().xy
				
				textX = float(textX[0])
				textY = float(textY[0])
				
				textColour = "blue"

				if fillOp == FILL_OPTIONS[0]: #Fill the fields according to the barrovian zones
					fillColour = "white"
					if "Melt" in self.phases:
						fillColour = "gold"
						textColour = "black"
					elif "Crd" in self.phases:
						fillColour = "magenta"
						textColour = "black"
					elif "Sil" in self.phases:
						fillColour = "blueviolet"
						textColour = "white"
					elif "Ky" in self.phases:
						fillColour = "cyan"
						textColour = "black"
					elif "St" in self.phases:
						fillColour = "yellow"
						textColour = "black"
					elif "Grt" in self.phases:
						fillColour = "red"
						textColour = "black"
					elif "Bt" in self.phases:
						fillColour = "sandybrown"
						textColour = "black"
					elif "Chl" in self.phases:
						fillColour = "green"
						textColour = "white"
					pltIn.fill(x,y,color=fillColour)

				isValid = "No"
				if poly.is_valid:	
					isValid = "Yes"

				pltIn.plot(x, y, color = "black", marker = None, linestyle = "-", markersize = 7, linewidth = 2)
				
				if poly.area >= areaMin:
					pltIn.text(textX, textY, str(index), fontsize = 10, color = textColour,horizontalalignment='center', verticalalignment='center')

					
					fieldString = str(self.field).replace("POLYGON ((","")
					fieldString = fieldString.replace("))","")
					# csvWriter.write(str(index) + delimit  + self.phases + delimit +  isValid + delimit + fieldString + "\n")
					csvWriter.write(str(index) + delimit + self.phases + delimit)
					for phase in PHASES:
						if "(2)" + phase in self.phases:
							csvWriter.write("(2)" + phase + delimit)
						elif phase in self.phases:
							csvWriter.write(phase + delimit)
						else:
							csvWriter.write(delimit)
					csvWriter.write(isValid + delimit + fieldString + "\n")


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
		tPhases = tPhases.replace("ma","Mrg")
		tPhases = tPhases.replace("CHTD","Cld")
		self.phases = tPhases