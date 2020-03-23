#Class for definining points in a domino line
#Has temperature, pressure, and a name
#name is not defined by default
from shapely.geometry import Point, LineString, MultiLineString
from shapely.ops import linemerge, snap
import copy

T_THRESH = 1.5
P_THRESH = 50
EXTRAP_RATIO = 50

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

	def extrapLine(self, atEnd = True, bothEnds = False):
		#Returns a line that is extrapolated in one direction with the same direction as the first or last segment
		#If atEnd = True it will extrapolate the last segment if not, it will extrapolate the first segment
		#This is useful for checking the intersection of reaction lines in a phase diagram
		if bothEnds:
			loopNum = 2
		else:
			loopNum = 1
			newCoords = list(copy.deepcopy(self.PTline.coords))
			for i in range(loopNum):
				if atEnd:
					p2 = self.PTline.coords[len(self.PTline.coords)-1]
					p1 = self.PTline.coords[len(self.PTline.coords)-2]
				else:
					p2 = self.PTline.coords[0]
					p1 = self.PTline.coords[1]
				#Need something to check to make sure this is not extrapolating past the limits of the axes
				if len(self.axIntersec) > 1:
					break
				elif len(self.axIntersec) >0:
					if p2 == self.axIntersec[0].coords:
						continue

				nextPtX = (p2[0]-p1[0])*EXTRAP_RATIO + p2[0]
				nextPtY = (p2[1]-p1[1])*EXTRAP_RATIO + p2[1]

				nextPt = Point(nextPtX, nextPtY)


				if atEnd:
					newCoords.extend(nextPt.coords)
				else:
					newCoords.insert(0,list(nextPt.coords)[0]) #I dont know why I have to do this but it breaks if just feeding it the coords

				atEnd = not atEnd


			return LineString(newCoords)


	def axesIntersect(self, axes):
	#this checks intersections with axes and stores these 
	#Maximum of two axes intersections
		self.axIntersec = []
		self.intersectedAx = []
		self.axInterLoc = []
		for line in axes:
			intersect = self.PTline.intersection(line)

			if isinstance(intersect, Point):
				
				self.axIntersec.append(intersect)
				self.intersectedAx.append(line)
				
				x1 = intersect.coords[0][0]
				y1 = intersect.coords[0][1]
				
				x2= self.PTline.coords[0][0]
				y2 = self.PTline.coords[0][1]
				#x3, y3 = self.PTline.xy[len(self.PTline.xy)-1]
				
				eq_thresh = 0.0001
				if abs(x1-x2) <= eq_thresh and abs(y1-y2) <= eq_thresh:
					self.axInterLoc.append(0)
				else:
					self.axInterLoc.append(len(self.PTline.coords)-1)





	def extrapIntersec(self, otherLine):
	#This method will check if the extrapolation of both this line and otherLine intersect
	#Will actually return itself PLUS the intersection coordinate added to the correct side of the line

		tryEnds = [True, True]
		lineInter = list(copy.deepcopy(self.PTline.coords))
		for i in range(len(tryEnds)):
			thisExtrap = self.extrapLine(tryEnds[0])
			

			for j in range(len(tryEnds)):
				otherExtrap = otherLine.extrapLine(tryEnds[1])
				# print("ThisLine extrapolation:")
				# print(thisExtrap)
				# print("OtherLine extrapolation:")
				# print(otherExtrap)
				intersect = thisExtrap.intersection(LineString(otherExtrap))

				if isinstance(intersect, Point):
					if tryEnds[0]:
						lineInter.extend(intersect.coords)
					else:
						lineInter = lineInter[::-1]
						lineInter.extend(intersect.coords)
					return LineString(lineInter)
				else:
					tryEnds[1] = not tryEnds[1]

			tryEnds[0] = not tryEnds[0]
		return None