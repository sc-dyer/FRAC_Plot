#Class for definining points in a domino line
#Has temperature, pressure, and a name
#name is not defined by default
from shapely.geometry import Point, LineString, MultiLineString, MultiPoint
from shapely.ops import linemerge, snap, nearest_points
import copy

T_THRESH = 1.5
P_THRESH = 10
EXTRAP_RATIO = 50
MAX_EXTRAP = 50
EQ_THRESH= 0.0001
import sys
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

	def extrapLine(self, atEnd = True, bothEnds = False, extrapRat = EXTRAP_RATIO):
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

				nextPtX = (p2[0]-p1[0])*extrapRat + p2[0]
				nextPtY = (p2[1]-p1[1])*extrapRat + p2[1]

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
				
				
				#Checking which end is the intersection
				if abs(x1-x2) <= EQ_THRESH and abs(y1-y2) <= EQ_THRESH:
					self.axInterLoc.append(0)
				else:
					self.axInterLoc.append(len(self.PTline.coords)-1)
			
			else:
				#If there is an error in the data and the line does not quite intersect...
				closestPoints = nearest_points(line,self.PTline)
				x1 = closestPoints[0].coords[0][0]
				y1 = closestPoints[0].coords[0][1]

				x2 = closestPoints[1].coords[0][0]
				y2 = closestPoints[1].coords[0][1]

				x3= self.PTline.coords[0][0]
				y3 = self.PTline.coords[0][1]
				
				x4 = self.PTline.coords[len(self.PTline.coords)-1][0]
				y4 = self.PTline.coords[len(self.PTline.coords)-1][1]
				if abs(x1-x2) <= T_THRESH and abs(y1-y2) <= P_THRESH:
					
					# print("Closest Points: ")
					# print(closestPoints[0].coords[0])
					# print(closestPoints[1].coords[0])
					if abs(x2-x3) <= EQ_THRESH and abs(y2-y3) <= EQ_THRESH:
						self.axIntersec.append(closestPoints[0])
						self.intersectedAx.append(line)
						self.axInterLoc.append(0)
					elif abs(x2-x4) <= EQ_THRESH and abs(y2-y4) <= EQ_THRESH:
						self.axIntersec.append(closestPoints[0])
						self.intersectedAx.append(line)
						self.axInterLoc.append(len(self.PTline.coords)-1)


	def extrapIntersec(self, otherLine, tMin = 0, tMax = 10000, pMin = 0, pMax = 3000000, extrapolRat = EXTRAP_RATIO):
	#This method will check if the extrapolation of both this line and otherLine intersect
	#Will actually return itself PLUS the intersection coordinate added to the correct side of the line
		smallestDistance = sys.float_info.max #Used as a flag in case of divergence
		tryEnds = [True, True]
		lineInter = list(copy.deepcopy(self.PTline.coords))
		# chkDiverge= False
		parallelCount = 1
		# if extrapolRat == MAX_EXTRAP:
		# 	#This is a case to check for lines that dont intersect because the final two points diverge
		# 	# extrapolRat = EXTRAP_RATIO
		# 	chkDiverge = True
		# 	parallelCount = 2

		for i in range(len(tryEnds)):
			thisExtrap = self.extrapLine(tryEnds[0],extrapRat = extrapolRat)
			for k in range(parallelCount):
				# if chkDiverge:
				# 	#This is a case to check for lines that dont intersect because the final two points diverge
				# 	if k == 0:
				# 		thisExtrap = thisExtrap.parallel_offset(10, 'left', join_style = 2)
				# 	else:
				# 		thisExtrap = thisExtrap.parallel_offset(10, 'right', join_style = 2)
				intersect = thisExtrap.intersection(LineString(otherLine.PTline))

				if isinstance(intersect, Point):
					#Triggers if the thisLine extrapolation intersects otherLine without extrapolation
					if not tryEnds[0]:
						lineInter = lineInter[::-1]
					return LineString(lineInter), 0

				for j in range(len(tryEnds)):
					otherExtrap = otherLine.extrapLine(tryEnds[1],extrapRat = extrapolRat)
					# if chkDiverge:
					# 	#This is a case to check for lines that dont intersect because the final two points diverge
					# 	if k == 1:
					# 		otherExtrap = otherExtrap.parallel_offset(10, 'left', join_style = 2)
					# 	else:
					# 		otherExtrap = otherExtrap.parallel_offset(10, 'right', join_style = 2)
					# print("\n")
					# print("ThisLine extrapolation:")
					# print(thisExtrap)
					# print("OtherLine extrapolation:")
					# print(otherExtrap)
					intersect = self.PTline.intersection(LineString(otherExtrap))

					if isinstance(intersect, Point):
						#Triggers if the otherline extrapolation intersects thisLine without extrapolation

						firstPoint = Point(self.PTline.coords[0])
						lastPoint = Point(self.PTline.coords[len(self.PTline.coords)-1])

						distanceFirst = intersect.distance(firstPoint)
						distanceLast = intersect.distance(lastPoint)
						print("Intersection of unextrapolated line:")
						print(distanceFirst)
						print(distanceLast)
						print("ThisLine extrapolation:")
						print(thisExtrap)
						print("OtherLine extrapolation:")
						print(otherExtrap)
						print("yay an intersection!")
						print(intersect)
						print(tryEnds[0])
						print(tryEnds[1])
						if distanceFirst > distanceLast:
							lineInter = lineInter[::-1]
						return LineString(lineInter), 0


					intersect = thisExtrap.intersection(LineString(otherExtrap))

					if isinstance(intersect, Point):
						t = intersect.coords[0][0]
						p = intersect.coords[0][1]
						print("ThisLine extrapolation:")
						print(thisExtrap)
						print("OtherLine extrapolation:")
						print(otherExtrap)
						print("yay an intersection!")
						print(intersect)
						print(tryEnds[0])
						print(tryEnds[1])

						if t< tMin or t >tMax or p <pMin or p>pMax:
							return None, smallestDistance
						# if chkDiverge:
						# 	if not tryEnds[0]:
						# 		lineInter = lineInter[::-1]
						# 	return LineString(lineInter)
						if tryEnds[0]:
							lineInter.extend(intersect.coords)
						else:
							lineInter = lineInter[::-1]
							lineInter.extend(intersect.coords)

						return LineString(lineInter), 0
					elif isinstance(intersect, MultiPoint) or isinstance(intersect, MultiLineString):
					#If the extrapolations are roughly parallel, that means you dont need an extra intersection point
					#Just add the line as normal
						print("Multiple intersections ")
						print("yay an intersection!")
						print(intersect)
						if not tryEnds[0]:
							lineInter = lineInter[::-1]
						return LineString(lineInter), 0
							

					else:
						# print("No intersect found")
						# print(intersect)
						lineDistance = thisExtrap.distance(otherExtrap)
						smallestDistance = min(smallestDistance,lineDistance)
						tryEnds[1] = not tryEnds[1]

			tryEnds[0] = not tryEnds[0]

		return None, smallestDistance