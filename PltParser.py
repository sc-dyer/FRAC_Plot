#Class for parsing out a Plt file from theriak-domino
#Will find each "reaction" in the file and initialize them as DomLines
#Can also connect lines that should be attached

from DomLine import DomLine, EQ_THRESH, MAX_EXTRAP
from DomPoly import DomPoly
import re
from shapely.geometry import Point, LineString, Polygon
from shapely.ops import polygonize
import copy
import math
import sys
GRPH_CODE = "gph" #Used to remove nonexistant lines that dont conain graphite in a C saturated rock

DIST_THRESH = 50 #This is the calue used as a threshold for linking lines. If the distance between endpoints exceeds this, then it will return multiple polys

class PltParser:

	def __init__(self, fileName, isPhase = False):
		try:
			pltFile = open(fileName, 'r')
		except:
			print("File " + fileName + " not found.")
			exit()



		rawText = pltFile.read()
		textLines = rawText.splitlines()

		self.xAx = re.split(":",textLines[0])[1] #First line of the file is x-axis, expecting "L:Temperature [C]""
		self.yAx = textLines[1] #Second line of the file is the y-axis, expecting "Pressure [bar]""
		headLine = textLines[2].split() #Header line includes Temperature and pressure range of the plot
		self.Tmin = float(headLine[0]) 
		self.Tmax = float(headLine[1])
		self.Pmin = float(headLine[2])
		self.Pmax = float(headLine[3])
		metadataLines = int(textLines[3].split()[0])
		self.metadata = []
		for i in range (4,metadataLines + 4):
			#Iterate through the metadataLines to record THERIN compo and database info
			nextLineSplit = textLines[i].split()
			thisMeta = ""
			#These conditionals are so that the text does not get split up
			for j in range(0, len(nextLineSplit)):
				if j == 4:
					thisMeta += nextLineSplit[j][1:]
				elif j > 4:
					thisMeta +=" " + nextLineSplit[j]
			self.metadata.append(thisMeta)
		
		nextLineNum = 4 + metadataLines
		#The rest of the file defines the Domino lines to plot
		#Expected sequence is as follows:
		#2  NumPoints  0  0  0  0
		#X1  Y1  3  X2  Y2  2  X3  Y3  2 ... X7  Y7  2
		#X8  Y8  2  ... X14  Y14  2  and so on
		#(X1-dx)  Y1  0  0  0LeftSide
		#(X1+dx)  Y1  0  0  0RightSide
		#The last two lines are more relevant for phase diagrams for X1-dx is X1 minus a very small number and X1+dx is X1 plus a very small number
		#For isopleths the only relevant part is in the position of "LeftSide"

		self.domLines = []
		while(nextLineNum < len(textLines)):

			numPoints = int(textLines[nextLineNum].split()[1]) #Extracting NumPoints
			if numPoints > 1:
				nextLineNum +=1
				#Now extract the points
				thisDomLine = DomLine()
				numPTLines = int(numPoints/7) #Each PT line can have a maximum of 7 PT points
				if numPoints%7 > 0:
					numPTLines += 1
				for i in range(numPTLines):
					ptLine = textLines[nextLineNum].split()
					for i in range(0,len(ptLine),3):
						#Add each temperature and pressure point to the domino line
						ti = float(ptLine[i])
						pi = float(ptLine[i+1])
						thisDomLine.addPT(ti,pi)
					nextLineNum += 1
				nextLineSplit = textLines[nextLineNum].split()
				leftsideName = ""				
				#These conditionals are so text doesnt get split up if LeftSide is a group of phases
				for i in range(0, len(nextLineSplit)):
					if i == 4:
						leftsideName += nextLineSplit[i][1:]
					elif i > 4:
						leftsideName +=" " + nextLineSplit[i]
				
				thisDomLine.addLeftSide(leftsideName)
				
				nextLineNum += 1
				nextLineSplit = textLines[nextLineNum].split()

				rightsideName = ""

				for i in range(0, len(nextLineSplit)):
					if i == 4:
						if len(nextLineSplit[i]) > 1: #Wont add a rightside if the label is empty
							rightsideName += nextLineSplit[i][1:]
					elif i > 4:
						rightsideName +=" " + nextLineSplit[i]
				if len(rightsideName) > 0:
					thisDomLine.addRightSide(rightsideName)
				nextLineNum += 1	
				twoPtFlag = True 
				#Check if the line is 2 points if the two points are close enough to ignore
				#I do this because of errors that crop up with little baby lines such as going in the wrong direction
				if numPoints ==2:
					x1 = thisDomLine.PTline.coords[0][0]
					y1 = thisDomLine.PTline.coords[0][1]
					
					x2= thisDomLine.PTline.coords[1][0]
					y2 = thisDomLine.PTline.coords[1][1]
					if abs(x1-x2) <= 0.1 and abs(y1-y2) <= 1:
						twoPtFlag = False
				if twoPtFlag and ( not isPhase or (GRPH_CODE in thisDomLine.leftSide and GRPH_CODE in thisDomLine.rightSide)):
					
					self.domLines.append(thisDomLine)
			else:
				nextLineNum += 4 
			
		self.joinLines()

	def sortLeft(self,low, high):
		#Function to sort the somLines by leftSide name in alphabetical order
		#Lets do a quick sort because thats fun
		
		if low < high:

			i = low - 1
			pivot = self.domLines[high].leftSide #Compare every element in list to the pivot
			
			for j in range(low,high):

				if(self.domLines[j].leftSide == None or self.domLines[j].leftSide < pivot):
					#Increment i by 1 if j is before pivot alphabetically
					#Then swap element at i and j

					i = i + 1
					self.domLines[i], self.domLines[j] = self.domLines[j], self.domLines[i]

			#Swap element at high and element at i + 1
			#Putting all elements bigger in volume to left of the element at high
			self.domLines[i+1], self.domLines[high] = self.domLines[high], self.domLines[i+1]

			pi = i + 1
			
			self.sortLeft(low, pi-1)
			self.sortLeft(pi+1, high)

	def sortRight(self, low, high):
		#Function to sort the somLines by leftSide name in alphabetical order
		#Lets do a quick sort because thats fun
		
		if low < high:

			i = low - 1
			pivot = self.domLines[high].rightSide #Compare every element in list to the pivot
			
			for j in range(low,high):

				if(self.domLines[j].rightSide == None or self.domLines[j].rightSide < pivot):
					#Increment i by 1 if j is before pivot alphabetically
					#Then swap element at i and j

					i = i + 1
					self.domLines[i], self.domLines[j] = self.domLines[j], self.domLines[i]

			#Swap element at high and element at i + 1
			#Putting all elements bigger in volume to left of the element at high
			self.domLines[i+1], self.domLines[high] = self.domLines[high], self.domLines[i+1]

			pi = i + 1
			
			self.sortRight(low, pi-1)
			self.sortRight(pi+1, high)

	def groupLines(self):
		#Function to group domLines of the same type into one subarray
		lineGroups = []
		self.sortRight(0,len(self.domLines)-1)
		self.sortLeft(0,len(self.domLines)-1)
		thisGroup = [self.domLines[0]]

		for i in range(1, len(self.domLines)):
			#Checks if they have the same flags and groups them then places them in the lineGroups array
			if self.domLines[i].leftSide == thisGroup[0].leftSide and self.domLines[i].rightSide == thisGroup[0].rightSide:
				thisGroup.append(self.domLines[i])
			else:
				lineGroups.append(thisGroup)
				thisGroup = [self.domLines[i]]
		lineGroups.append(thisGroup)
		return lineGroups

	def joinLines(self):
		#Joins all lines in a shared group
		#Only if they touch
		#Tests head to tail, tail to tail, and head to head for each line in the group
		lineGroups = self.groupLines()
		
		newDomLines = []
		for i in range(len(lineGroups)):

			thisGroup = copy.deepcopy(lineGroups[i])
			j = 0
			# print("For group "+ thisGroup[j].leftSide + " Length = " + str(len(thisGroup)))
			#While loops to account for changing array sizes
			while j < len(thisGroup):
				
				
				didRemove = True
				while didRemove: #Does not iterate j until it iterates the array without any matches
					k = 0
					didRemove = False

					while k<len(thisGroup):
						if k!= j:
							
							didJoin = thisGroup[j].joinLine(thisGroup[k])
							if didJoin:
								# print(str(j) + " and " + str(k) + " match")		
								thisGroup.pop(k)

								if k < j:
									j -= 1
								# print("Now " + str(j) + " and " + str(k) + " Length = " + str(len(thisGroup)))
								didRemove = True
							else:

								# print(str(j) + " and " + str(k) + " DONT match")	
								k += 1
								# print("Now " + str(j) + " and " + str(k) + " Length = " + str(len(thisGroup)))
								
						else:
							k += 1
				
				j += 1
			newDomLines.extend(thisGroup)
		self.domLines = newDomLines


	def getLines(self, left = None, right = None):
	#Returns a list of DomLines that have a left or right side that match left and right
	#keep as none if one side does not matter

		matchList = []
		for line in self.domLines:
			if right == None:
				if line.leftSide == left:
					matchList.append(line)
			elif left ==None:
				if line.rightSide == right:
					matchList.append(line)
			else:
				if line.leftSide == left and line.rightSide == right:
					matchList.append(line)

		return matchList


	def getPolys(self):
		#This will return any polygons that exist in the plt file
		#Useful more for phase diagrams

		#Join matching domLines
		self.joinLines()
		self.sortLeft(0,len(self.domLines)-1)
		for line in self.domLines:
			print(line.leftSide)
			print(line.rightSide)
			print(line.PTline)
			print("\n")
		leftAx = LineString([(self.Tmin,self.Pmin),(self.Tmin,self.Pmax)])
		rightAx = LineString([(self.Tmax,self.Pmin),(self.Tmax,self.Pmax)])
		botAx = LineString([(self.Tmin,self.Pmin),(self.Tmax,self.Pmin)])
		topAx = LineString([(self.Tmin,self.Pmax),(self.Tmax,self.Pmax)])

		#Axes used for lines that dont intersect on one side
		axes = [leftAx,rightAx,botAx,topAx]

		for line in self.domLines:
			#Find axis intersections
			line.axesIntersect(axes)

		self.polyList = []
		self.failedPolys = []
		thisField = self.domLines[0].leftSide
		checkLeft = True
		commonLines = []
		addLines = True
		print("")
		for h in range(2):
			for i in range(len(self.domLines)+1): #+1 because it wont run through the else statement and check right side if the lat entry is a different value

				#Group together all lines that have the same leftside
				if i < len(self.domLines):
					print("Checking line " + self.domLines[i].leftSide)
				if checkLeft and i<len(self.domLines) and self.domLines[i].leftSide == thisField:
					print("Matches")
					x1 = self.domLines[i].PTline.coords[0][0]
					y1 = self.domLines[i].PTline.coords[0][1]
					
					x2= self.domLines[i].PTline.coords[len(self.domLines[i].PTline.coords)-1][0]
					y2 = self.domLines[i].PTline.coords[len(self.domLines[i].PTline.coords)-1][1]
					
					if not (abs(x1-x2) <= EQ_THRESH and abs(y1-y2) <= EQ_THRESH):
						#This is a conditional for the situation where a line has the same first and last point in its sequence
						commonLines.append(self.domLines[i])
						# commonLines.append(self.domLines[i].extrapLine(bothEnds = True))
						print(self.domLines[i].leftSide)
						print(self.domLines[i].PTline)

				elif (not checkLeft) and i <len(self.domLines) and self.domLines[i].rightSide == thisField:
					#This is for the second iteration through the second forloop

					hasLeft = False
					for j in range(len(self.domLines)):
						if self.domLines[j].leftSide == thisField:
							hasLeft = True

					if hasLeft:
						print(thisField + " is already present")
						addLines = False
					else:
						print(thisField + " is not present")
						addLines = True
						x1 = self.domLines[i].PTline.coords[0][0]
						y1 = self.domLines[i].PTline.coords[0][1]
						
						x2= self.domLines[i].PTline.coords[len(self.domLines[i].PTline.coords)-1][0]
						y2 = self.domLines[i].PTline.coords[len(self.domLines[i].PTline.coords)-1][1]
						
						if not (abs(x1-x2) <= EQ_THRESH and abs(y1-y2) <= EQ_THRESH):
							#This is a conditional for the situation where a line has the same first and last point in its sequence
							commonLines.append(self.domLines[i])

							# commonLines.append(self.domLines[i].extrapLine(bothEnds = True))
							print(self.domLines[i].rightSide)
							print(self.domLines[i].PTline)

				elif addLines and (len(commonLines) > 0 or checkLeft):
					#Once a different leftSide is reached, search through the whole array for matching rightSides and add them to the commonLines
					if checkLeft:
						# print("Not a match")
						for j in range(len(self.domLines)):
							# print("Checking line " + self.domLines[j].rightSide)
							if self.domLines[j].rightSide == thisField:
								# print("Matches")
								commonLines.append(self.domLines[j])
								# commonLines.append(self.domLines[j].extrapLine(bothEnds = True))
								print(self.domLines[j].rightSide)
								print(self.domLines[j].PTline)
							# else:
							# 	print("Not a match")

					print("Number of lines = " + str(len(commonLines)))
					sortedLines = self.sortLines(commonLines)

					print("Lines in order:")
					polyPts = []
					for group in sortedLines:

						for line in group:
							print(line.PTline)
							polyPts.extend(line.PTline.coords)
						if len(polyPts) > 2:
							self.polyList.append(DomPoly(polyPts,thisField))
						polyPts = []

					print("\n")
					


					if i < len(self.domLines):
						if checkLeft:
							x1 = self.domLines[i].PTline.coords[0][0]
							y1 = self.domLines[i].PTline.coords[0][1]
							
							x2= self.domLines[i].PTline.coords[len(self.domLines[i].PTline.coords)-1][0]
							y2 = self.domLines[i].PTline.coords[len(self.domLines[i].PTline.coords)-1][1]
							thisField = self.domLines[i].leftSide
							if not (abs(x1-x2) <= EQ_THRESH and abs(y1-y2) <= EQ_THRESH):
								
								commonLines = [self.domLines[i]]
							else:
								commonLines = []
						else:
							thisField = self.domLines[i].rightSide
							commonLines = []
							hasLeft = False
							for j in range(len(self.domLines)):
								if self.domLines[j].leftSide == thisField:
									hasLeft = True

							if hasLeft:
								print(thisField + " is already present")
								addLines = False
							else:
								print(thisField + " is not present")
								addLines = True
								x1 = self.domLines[i].PTline.coords[0][0]
								y1 = self.domLines[i].PTline.coords[0][1]
								
								x2= self.domLines[i].PTline.coords[len(self.domLines[i].PTline.coords)-1][0]
								y2 = self.domLines[i].PTline.coords[len(self.domLines[i].PTline.coords)-1][1]
								
								if not (abs(x1-x2) <= EQ_THRESH and abs(y1-y2) <= EQ_THRESH):
									#This is a conditional for the situation where a line has the same first and last point in its sequence
									commonLines.append(self.domLines[i])

									# commonLines.append(self.domLines[i].extrapLine(bothEnds = True))
									print(self.domLines[i].rightSide)
									print(self.domLines[i].PTline)
				else:
					if i < len(self.domLines):
						if checkLeft:
							x1 = self.domLines[i].PTline.coords[0][0]
							y1 = self.domLines[i].PTline.coords[0][1]
							
							x2= self.domLines[i].PTline.coords[len(self.domLines[i].PTline.coords)-1][0]
							y2 = self.domLines[i].PTline.coords[len(self.domLines[i].PTline.coords)-1][1]
							thisField = self.domLines[i].leftSide
							if not (abs(x1-x2) <= EQ_THRESH and abs(y1-y2) <= EQ_THRESH):
								
								commonLines = [self.domLines[i]]
								print(self.domLines[i].leftSide)
								print(self.domLines[i].PTline)
							else:
								commonLines = []
						else:
							thisField = self.domLines[i].rightSide
							commonLines = []
							hasLeft = False
							for j in range(len(self.domLines)):
								if self.domLines[j].leftSide == thisField:
									hasLeft = True

							if hasLeft:
								print(thisField + " is already present")
								addLines = False
							else:
								print(thisField + " is not present")
								addLines = True
								x1 = self.domLines[i].PTline.coords[0][0]
								y1 = self.domLines[i].PTline.coords[0][1]
								
								x2= self.domLines[i].PTline.coords[len(self.domLines[i].PTline.coords)-1][0]
								y2 = self.domLines[i].PTline.coords[len(self.domLines[i].PTline.coords)-1][1]
								
								if not (abs(x1-x2) <= EQ_THRESH and abs(y1-y2) <= EQ_THRESH):
									#This is a conditional for the situation where a line has the same first and last point in its sequence
									commonLines.append(self.domLines[i])

									# commonLines.append(self.domLines[i].extrapLine(bothEnds = True))
									print(self.domLines[i].rightSide)
									print(self.domLines[i].PTline)
			print("Checking right side now")
			checkLeft = False
			addLines = False
			commonLines = []
			self.sortRight(0,len(self.domLines)-1)
			thisField = self.domLines[0].rightSide	

					

					

	def sortLines(self, lineGroup):
		#This method will sort lines based on closeness and shared axis intersections
		#Will return a sorted array of lines
		leftAx = LineString([(self.Tmin,self.Pmin),(self.Tmin,self.Pmax)])
		rightAx = LineString([(self.Tmax,self.Pmin),(self.Tmax,self.Pmax)])
		botAx = LineString([(self.Tmin,self.Pmin),(self.Tmax,self.Pmin)])
		topAx = LineString([(self.Tmin,self.Pmax),(self.Tmax,self.Pmax)])
		axes = [leftAx,rightAx,botAx,topAx]
		multiSortedLines = []
		lineGroup = copy.deepcopy(lineGroup)
		sortedLines = [lineGroup[0]]
		numLines = len(lineGroup)
		lineGroup.pop(0)
		for i in range(numLines-1):
			thisLine = sortedLines[len(sortedLines)-1]
			linkPoint = Point(thisLine.PTline.coords[len(thisLine.PTline.coords)-1])
			#Need to check the end point first to see if it intersects with an axis
			atAxis = False
			# print(thisLine.PTline)
			# print(sortedLines[len(sortedLines)-1].PTline)

			for j in range(len(thisLine.intersectedAx)):
				if thisLine.axInterLoc[j] != 0:
					thisAx = thisLine.intersectedAx[j]
					thisIntersec = thisLine.axIntersec[j]
					atAxis = True
				print(thisLine.axInterLoc[j])
			
			# print("At axis = " + str(atAxis))
			if atAxis:
				#Iterate to find the other intersect, can possibly copy the stuff from above
				#Need to add the extra corner point to the end of thisLine if the corner case
				#Will continue this loop as long as there is another line that has two intersections with the same axis
				#Will run at least once to find the nextLine that intersects the same axis
				#polyPts.extend(list(thisLine.axIntersec))
				numInterAx = 0
				for j in range(len(lineGroup)):
					#Check how many other lines intersect this same axis
					
					nextLine = lineGroup[j]

					for k in range(len(nextLine.intersectedAx)):
						#Check both intersectedAxes and see if they are the same as the axes of the last point added to PolyPts
						if nextLine.intersectedAx[k] == thisAx:
							numInterAx += 1


				if numInterAx == 0:
					#This triggers in the literal corner case
					
					for j in range(len(lineGroup)):
						if j < len(lineGroup):
							nextLine = lineGroup[j]
							#if j != thisIndex:
							for k in range(len(nextLine.intersectedAx)):
								nextAx = nextLine.intersectedAx[k]

								axIntersec = thisAx.intersection(nextAx)

								if isinstance(axIntersec, Point):
									if (nextLine.axInterLoc[k] != 0):
										
										nextLine.PTline.coords = list(nextLine.PTline.coords)[::-1]
										nextLine.axesIntersect(axes)
									#polyPts.append(nextLine.axIntersec[k].coords)
									axCoords = list(axIntersec.coords)
									nextLineCoords = list(nextLine.PTline.coords)
									print(axCoords)
									print(nextLineCoords)
									nextLineCoords.insert(0,axCoords[0])
									print(nextLineCoords)
									nextLine.PTline.coords = nextLineCoords
									
									sortedLines.append(nextLine)

									lineGroup.pop(j)

								#Will this throw an error if >1 axis intersection on other axes? 
								#Not sure I have seen that before



				elif numInterAx == 1:
					#this is the normal case
					for j in range(len(lineGroup)):
					#Checks all commonLines for the same Axis Intersection
						if j < len(lineGroup):
							nextLine = lineGroup[j]

							for k in range(len(nextLine.intersectedAx)):
								
								#Check both intersectedAxes and see if they are the same as the axes of the last point added to PolyPts
								if nextLine.intersectedAx[k] == thisAx:
								#Checks if the intersection is at the beginning or end of nextLine
								#If its at the end then it reverses it
									print("This line:")
									print(thisLine.PTline)
									#print("Axis intersect loc: " + str(thisLine.axInterLoc[lastPt]))
									print("Next line: ")
									print(nextLine.PTline)
									
									print("Axis intersect loc: " + str(nextLine.axInterLoc[k]))
									
									if (nextLine.axInterLoc[k] != 0):
										nextLine.PTline.coords = list(nextLine.PTline.coords)[::-1]
										nextLine.axesIntersect(axes)
									sortedLines.append(nextLine)

									lineGroup.pop(j)

					

				elif numInterAx >1:
					#This is the case when you have multiple polygons poking out different parts of the same axis
					#Im not entirely sure what to do here. I suppose the smartest thing would be to group these as seperate polygons
					#Perhaps what is easier is making an invalid polygon and sort it out later.
					#Question is how to establish the order of operations? 
					#Okay, in this case we will connext it to the line of the next highest pressure that is not already in the array
					

					if thisAx == leftAx or thisAx == rightAx:
						xyIndex = 1
					elif thisAx == botAx or thisAx == topAx:
						xyIndex = 0

					thisInterVal = thisIntersec.coords[0][xyIndex]
					bestLine = None #this is the best next choice for line
					bestInterVal = -1
					bestIndex = -1
					bestPt = -1
					for j in range(len(lineGroup)):
						
						nextLine = lineGroup[j]
						
						for k in range(len(nextLine.intersectedAx)):

							if nextLine.intersectedAx[k] == thisAx:
								nextInterVal = nextLine.axIntersec[k].coords[0][xyIndex]

								if bestLine == None and nextInterVal > thisInterVal:
									bestLine = nextLine
									bestInterVal = nextLine.axIntersec[k].coords[0][xyIndex]
									bestIndex = j
									bestPt = k
								elif nextInterVal > thisInterVal and nextInterVal < bestInterVal:
									bestLine = nextLine
									bestInterVal = nextLine.axIntersec[k].coords[0][xyIndex]
									bestIndex = j
									bestPt = k
					if bestLine == None:
						#Should trigger if there is no nexthighest so will take the lowest value
						for j in range(len(lineGroup)):
							
							for k in range(len(nextLine.intersectedAx)):

								if nextLine.intersectedAx[k] == thisAx:
									nextInterVal = nextLine.axIntersec[k].coords[0][xyIndex]

									if nextInterVal < thisInterVal:
										bestLine = nextLine
										bestInterVal = nextLine.axIntersec[k].coords[0][xyIndex]
										bestIndex = j
										bestPt = k
									elif nextInterVal < thisInterVal and nextInterVal < bestInterVal:
										bestLine = nextLine
										bestInterVal = nextLine.axIntersec[k].coords[0][xyIndex]
										bestIndex = j
										bestPt = k

					if bestLine != None: 
						if bestLine.axInterLoc[bestPt] != 0:
							bestLine.PTline.coords = list(bestLine.PTline.coords)[::-1]
							bestLine.axesIntersect(axes)
					
						
						sortedLines.append(bestLine)

						lineGroup.pop(bestIndex)
				

			else:
				#The normal case where you dont have an axis intersection
				atEnd = False
				minDistance = sys.float_info.max
				nextIndex = 0

				for j in range(len(lineGroup)):
					nextLine = lineGroup[j]
					firstPoint = Point(nextLine.PTline.coords[0])
					lastPoint = Point(nextLine.PTline.coords[len(nextLine.PTline.coords)-1])

					distFirst = firstPoint.distance(linkPoint)
					distLast = lastPoint.distance(linkPoint)

					#iterate the array and find the line with an the end that is closest to the last point in the last array
					if distFirst < minDistance or distLast < minDistance:
						if distLast < distFirst:
							atEnd = True
							minDistance = distLast
						else:
							atEnd = False
							minDistance = distFirst
						nextIndex = j

				if len(lineGroup) > 0:
					bestMatch = lineGroup[nextIndex]
					if atEnd:
						bestMatch.PTline.coords = list(bestMatch.PTline.coords)[::-1]
						bestMatch.axesIntersect(axes)
					if minDistance > DIST_THRESH:
						multiSortedLines.append(sortedLines)
						sortedLines = []
					lineGroup.pop(nextIndex)
					sortedLines.append(bestMatch)

		multiSortedLines.append(sortedLines)
		for group in multiSortedLines:
			firstLine = group[0]
			lastLine = group[len(group)-1]
			print("Checking for corner")

			#For when there is the corner case but the axis intersections are at the beginning and end of the list
			for j in range(len(lastLine.intersectedAx)):
				if lastLine.axInterLoc[j] != 0:
					lastAx = lastLine.intersectedAx[j]
					lastIntersec = lastLine.axIntersec[j]
					print(lastIntersec)
					for k in range(len(firstLine.intersectedAx)):
						firstAx = firstLine.intersectedAx[k]

						axIntersec = lastAx.intersection(firstAx)
						print(axIntersec)
						if isinstance(axIntersec, Point):

						#Should only trigger if they intersect at a point
							if (firstLine.axInterLoc[k] == 0):
								
								
								axCoords = list(axIntersec.coords)
								lastLineCoords = list(lastLine.PTline.coords)
								print(axCoords)
									
								lastLineCoords.append(axCoords[0])
					
								lastLine.PTline.coords = lastLineCoords
					



		
		return multiSortedLines

			
	def printLines(self):
		for i in range(len(self.domLines)):
			print(self.domLines[i].leftSide + ": " + str(list(self.domLines[i].PTline.coords)))