#Class for definining points in a domino line
#Has temperature, pressure, and a name
#name is not defined by default
from shapely.geometry import Point, LineString, MultiLineString
from shapely.ops import linemerge, snap
import copy
T_THRESH = 1.5
P_THRESH = 35
class DomLine: 

 	def __init__(self, name = None):

 		self.PTline = LineString() #Empty PT line where x = T and y = P
 		self.leftSide = name #if only one value to a line (e.g. isopleth) only the leftside will have a name
 		self.rightSide = None

 	def addPT(self, tIn, pIn):
 		pointList = list(self.PTline.coords)
 		pointList.append((tIn,pIn))
 		if(len(self.PTline.coords)<1):
 			self.PTline = Point(pointList)
 		else:

 			self.PTline = LineString(pointList)

 	def addLeftSide(self, name):
 		self.leftSide = name

 	def addRightSide(self, name):
 		self.rightSide = name

 	def joinLine(self, toJoin):
 	#Attempts to add the t and ps of toJoin to this DomLine
 	#BAsed off of T_THRESH and P_THRESH
 	#Returns true if join is succesful
 		
 		numPoints = len(self.PTline.coords)
 		thisPTline = copy.deepcopy(self.PTline)
 		#Check if the threshold is met for P and T then adds the values of toJoin to this domLine
 		for i in range(3):
 			joinX = toJoin.PTline.coords[0][0]
 			joinY = toJoin.PTline.coords[0][1]
 			lastX = thisPTline.coords[numPoints-1][0]
 			lastY = thisPTline.coords[numPoints-1][1]
 			
	 		if abs(joinX - lastX) < T_THRESH and abs(joinY - lastY) < P_THRESH:

	 			newCoords = list(thisPTline.coords)
	 			newCoords.extend(list(toJoin.PTline.coords))
	 			self.PTline.coords = newCoords

	 			
	 			return True
	 		if i == 0:
	 			toJoin.PTline.coords = list(toJoin.PTline.coords)[::-1]
	 			#First pass no match, invert the toJoin line
	 		if i == 1:
	 			#Second pass, still no match, invert the self PT line. This is to cover the case where the first cell in each line matches
	 			thisPTline.coords = list(thisPTline.coords)[::-1]
	 	

	 	
 		return False
