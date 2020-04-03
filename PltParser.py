#Class for parsing out a Plt file from theriak-domino
#Will find each "reaction" in the file and initialize them as DomLines
#Can also connect lines that should be attached

from DomLine import DomLine, EQ_THRESH, MAX_EXTRAP
from DomPoly import DomPoly
import re
from shapely.geometry import Point, LineString, Polygon
from shapely.ops import polygonize
import copy

GRPH_CODE = "gph" #Used to remove nonexistant lines that dont conain graphite in a C saturated rock

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
				if checkLeft and i<len(self.domLines) and self.domLines[i].leftSide == thisField:
					
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

				elif addLines and len(commonLines) > 0:
					#Once a different leftSide is reached, search through the whole array for matching rightSides and add them to the commonLines
					if checkLeft:
						for j in range(len(self.domLines)):
							if self.domLines[j].rightSide == thisField:
								commonLines.append(self.domLines[j])
								# commonLines.append(self.domLines[j].extrapLine(bothEnds = True))
								print(self.domLines[j].rightSide)
								print(self.domLines[j].PTline)

					
					#Start with the smallest line:

					
					lineCount = 0
					thisIndex = 0
					lastIndex = 0
					#print(commonLines)
					thisLine = commonLines[thisIndex]
					# for j in range(len(commonLines)):
					# 	if commonLines[j].PTline.length < thisLine.PTline.length:
					# 		thisIndex = j
					# 		thisLine = commonLines[j]

					indexArray = []
					if checkLeft:
						print(thisLine.leftSide)
					else:
						print(thisLine.rightSide)
					print("Commonline: " + str(thisIndex))
					
					print(thisLine.PTline)
					polyPts = []

					#This is for the edge case that thisLine only intersects the axes
					if len(thisLine.axIntersec) > 1:
						polyPts.extend(list(thisLine.PTline.coords))
						indexArray.append(thisIndex)
						print("Two axes intersect found, adding line:")
						print(thisLine.PTline)
						lineCount += 1
						lastPt = 1
						if thisLine.axInterLoc[lastPt] == 0:
							lastPt = 0
						while len(thisLine.axIntersec) >1 and lineCount < len(commonLines):
							#Will continue this as long as thisLine continues to have two axes intersections
							print("Linecount: " + str(lineCount))
							print("Commonline: " + str(thisIndex))
							for j in range(len(commonLines)):
								#Check all commonLines to see what line intersects with the same axis (Should only be one)
								
								nextLine = commonLines[j]

								for k in range(len(nextLine.intersectedAx)):
									if j != thisIndex and j != lastIndex:
										#Check both intersectedAxes and see if they are the same as the axes of the last point added to PolyPts
										if nextLine.intersectedAx[k] == thisLine.intersectedAx[lastPt]:
											#Checks if the intersection is at the beginning or end of nextLine
											#If its at the end then it reverses it
											print(thisLine.PTline)
											print(thisLine.axInterLoc[lastPt])
											print(nextLine.PTline)
											
											print(nextLine.axInterLoc[k])
											
											#if (nextLine.axInterLoc[k] != 0 and thisLine.axInterLoc[lastPt] != 0) or (nextLine.axInterLoc[k] == 0 and thisLine.axInterLoc[lastPt] == 0) or (nextLine.axInterLoc[k] != 0 and thisLine.axInterLoc[lastPt] == 0:
											if (nextLine.axInterLoc[k] != 0):
												nextLine.PTline.coords = list(nextLine.PTline.coords)[::-1]
												nextLine.axesIntersect(axes)
											
											#polyPts.append(nextLine.axIntersec[k].coords)
											print(nextLine.PTline)
											lastIndex = -1
											#This  is set as such because the last intersection was with an axis
											thisLine = nextLine
											
											if len(thisLine.intersectedAx) > 1:
												if thisLine.axInterLoc[lastPt] == 0:
													if lastPt == 1:
														lastPt = 0
													else:
														lastPt = 1
												polyPts.extend(list(thisLine.PTline.coords))
												indexArray.append(thisIndex)
												print("Other axis intersect found, adding line:")
												print(thisLine.PTline)
												lineCount += 1 
												lastIndex = thisIndex #This is done so that the next time around it doesnt accidentally match with the same line if it can

											thisIndex = j
					#elif len(thisLine.axIntersec) >0:
						#polyPts.extend(thisLine.axIntersec)

					
					print("Total lines = " + str(len(commonLines)))
					failedConnection = 0
					while lineCount < len(commonLines) and failedConnection < 4:
						#Continue until all lines in commonLines are acocunted for in the polygon
						#At this point thisLine should NOT be in polyPts and is either connected to two other domlines or one domline and an axis
						intersectNLine = None
						print("Linecount = " +str(lineCount))
						print(thisLine.PTline)
						nextIndex = -1
						#Loop to check if the extrapolated line on either side matches with one of the other lines
						extrapRatio = 1
						#Start with small extrapolations then go up bigger
						#Making sure it is intersecting with the closest lines
						while extrapRatio <= MAX_EXTRAP and intersectNLine == None:
							print("Extrapolation Ratio = " + str(extrapRatio))
							nextIndex = -1
							while intersectNLine == None and nextIndex < len(commonLines)-1:
								nextIndex += 1
								isInArray= False
								for num in indexArray:
									if num == nextIndex:
										isInArray = True
								# print("Last Index: " + str(lastIndex))
								# print("This Index: " + str(thisIndex))
								# print("Next Index: " + str(nextIndex))
								if nextIndex != thisIndex and ((not isInArray) or (lineCount == len(commonLines)-1) and nextIndex ==indexArray[0]):
									intersectNLine = thisLine.extrapIntersec(commonLines[nextIndex],self.Tmin,self.Tmax,self.Pmin, self.Pmax, extrapRatio) #intersect is thisLine PLUS the intersection coordinate
							
							if extrapRatio < 5:
								extrapRatio += 1	
							elif extrapRatio >= 100:
								extrapRatio += 20
							else:
								extrapRatio += 5
							
						if intersectNLine != None:
							#Here we add thisLine plus the intersection as intersectNLine
							print("Found an intersection! Adding Line:")
							print(intersectNLine)
							polyPts.extend(list(intersectNLine.coords))
							indexArray.append(thisIndex)
							failedConnection = 0
							lineCount += 1
							lastIndex = thisIndex
							thisIndex = nextIndex
							thisLine = commonLines[thisIndex]
							
							
							
						if (nextIndex >= len(commonLines) - 1):
							#This should only trigger if the next endpoint is an intersection with an axis
							
							#So thisLine is still not in polyPts and neither is its intersection with the axis
							#But the intersection with the last line still is in this
							print("No intersection found :(")
							checkTwoIntersec = 2
							lastPt = 0
							
							print("Linecount = " + str(lineCount))
							if len(thisLine.axIntersec)==1 and lineCount > 0:
								#Can just intersect thisLine with the lastLine and attach in reverse order to polyPts (minus the last pt)
								# print("Commonline: " + str(thisIndex))
								# print(thisLine.PTline)
								# print("Commonline: " + str(lastIndex))
								# print(commonLines[lastIndex].PTline)
								intersectNLine = thisLine.extrapIntersec(commonLines[lastIndex],self.Tmin,self.Tmax,self.Pmin, self.Pmax)
								#print(intersectNLine)
								if intersectNLine !=None:
									intersectNLine.coords = list(intersectNLine.coords)[::-1]
									intersectNLine.coords = list(intersectNLine.coords)[1::]
									print("Intersected axis, adding line:")
									print(intersectNLine)
									polyPts.extend(list(intersectNLine.coords))
									indexArray.append(thisIndex)
									failedConnection = 0
									lineCount += 1

									lastIndex = thisIndex
									#Now we are faced with the tricky issue	of finding the other line that intersects with the same axis	
										
									while checkTwoIntersec > 1 and lineCount < len(commonLines):
										#Will continue this loop as long as there is another line that has two intersections with the same axis
										#Will run at least once to find the nextLine that intersects the same axis
										#polyPts.extend(list(thisLine.axIntersec))
										numInterAx = 0
										for j in range(len(commonLines)):
											#Check how many other lines intersect this same axis
											if j != thisIndex:
												nextLine = commonLines[j]

												for k in range(len(nextLine.intersectedAx)):
													#Check both intersectedAxes and see if they are the same as the axes of the last point added to PolyPts
													if nextLine.intersectedAx[k] == thisLine.intersectedAx[lastPt]:
														numInterAx += 1 

										if numInterAx == 0:
											#This triggers in the literal corner case
											thisAx = thisLine.intersectedAx[lastPt]
											for j in range(len(commonLines)):
												nextLine = commonLines[j]
												#if j != thisIndex:
												for k in range(len(nextLine.intersectedAx)):
													nextAx = nextLine.intersectedAx[k]

													axIntersec = thisAx.intersection(nextAx)

													if isinstance(axIntersec, Point):
														if (nextLine.axInterLoc[k] != 0):
															nextLine.PTline.coords = list(nextLine.PTline.coords)[::-1]
															nextLine.axesIntersect(axes)
														#polyPts.append(nextLine.axIntersec[k].coords)
														nextLine.PTline.coords = axIntersec.coords.extend(nextLine.PTline.coords)
														lastIndex = -1
														#This  is set as such because the last intersection was with an axis
														thisLine = nextLine
														checkTwoIntersec = len(thisLine.intersectedAx)
														if len(thisLine.intersectedAx) > 1:
															if k != 0:
															#	polyPts.append(nextLine.axIntersec[0])
																lastPt = 0
															else:
															#	polyPts.append(nextLine.axIntersec[1])
																lastPt = 1
															polyPts.extend(list(thisLine.PTline.coords))
															indexArray.append(thisIndex)
															print("Corner case, adding line:")
															print(thisLine.PTline)
															failedConnection = 0
															lineCount += 1 
															lastIndex = thisIndex #This is done so that the next time around it doesnt accidentally match with the same line if it can

														thisIndex = j



										elif numInterAx == 1:
											#this is the normal case
											for j in range(len(commonLines)):
											#Checks all commonLines for the same Axis Intersection

											#Problem: if at the corner, it could potentially loop forever
											#Other problem: if there are more than one line intersecting the axis, it will attach both of them these issues might be addressed by checking the 
											#number of lines that intersect the same axis in commonLines
												
												nextLine = commonLines[j]

												for k in range(len(nextLine.intersectedAx)):
													if j != thisIndex and j != lastIndex:
														#Check both intersectedAxes and see if they are the same as the axes of the last point added to PolyPts
														if nextLine.intersectedAx[k] == thisLine.intersectedAx[lastPt]:
														#Checks if the intersection is at the beginning or end of nextLine
														#If its at the end then it reverses it
															print("This line:")
															print(thisLine.PTline)
															print("Axis intersect loc: " + str(thisLine.axInterLoc[lastPt]))
															print("Next line: ")
															print(nextLine.PTline)
															
															print("Axis intersect loc: " + str(nextLine.axInterLoc[k]))
															
															if (nextLine.axInterLoc[k] != 0):
																nextLine.PTline.coords = list(nextLine.PTline.coords)[::-1]
														
															#polyPts.append(nextLine.axIntersec[k].coords)
															
															lastIndex = thisIndex
															#This  is set as such because the last intersection was with an axis

															thisLine = nextLine
															
															checkTwoIntersec = len(thisLine.intersectedAx)

															if len(thisLine.intersectedAx) > 1:
																if k != 0:
																#	polyPts.append(nextLine.axIntersec[0])
																	lastPt = 0
																else:
																#	polyPts.append(nextLine.axIntersec[1])
																	lastPt = 1
																polyPts.extend(list(thisLine.PTline.coords))
																indexArray.append(thisIndex)
																print("Other axis intersect found, adding line:")
																print(thisLine.PTline)
																failedConnection = 0
																lineCount += 1 
																lastIndex = thisIndex #This is done so that the next time around it doesnt accidentally match with the same line if it can
															#print(thisLine.PTline)
															thisIndex = j

											if thisIndex != lastIndex:
												lastIndex = thisIndex

										elif numInterAx >1:
											#This is the case when you have multiple polygons poking out different parts of the same axis
											#Im not entirely sure what to do here. I suppose the smartest thing would be to group these as seperate polygons
											#Perhaps what is easier is making an invalid polygon and sort it out later.
											#Question is how to establish the order of operations? 
											#Okay, in this case we will connext it to the line of the next highest pressure that is not already in the array
											thisAx = thisLine.intersectedAx[lastPt]

											if thisAx == leftAx or thisAx == rightAx:
												xyIndex = 1
											elif thisAx == botAx or thisAx == topAx:
												xyIndex = 0

											thisInterVal = thisLine.axIntersec[lastPt].coords[0][xyIndex]
											bestLine = None #this is the best next choice for line
											bestInterVal = -1
											bestIndex = -1
											bestPt = -1
											for j in range(len(commonLines)):
												isInPoly = False
												for val in indexArray:
													if j == val:
														isInPoly = True
												nextLine = commonLines[j]
												
												for k in range(len(nextLine.intersectedAx)):

													if (not isInPoly) and nextLine.intersectedAx[k] == thisAx:
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
												for j in range(len(commonLines)):
													isInPoly = False
													for val in indexArray:
														if j == val:
															isInPoly = True
													nextLine = commonLines[j]
													
													for k in range(len(nextLine.intersectedAx)):

														if (not isInPoly) and nextLine.intersectedAx[k] == thisAx:
															nextInterVal = nextLine.axIntersec[k].coords[0][xyIndex]

															if bestLine == None and nextInterVal < thisInterVal:
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
												if (bestLine.axInterLoc[bestPt] != 0 and thisLine.axInterLoc[lastPt] != 0) or (bestLine.axInterLoc[bestPt]  == 0 and thisLine.axInterLoc[lastPt] == 0):
													bestLine.PTline.coords = list(bestLine.PTline.coords)[::-1]
													bestLine.axesIntersect(axes)
											
												#polyPts.append(nextLine.axIntersec[k].coords)
												
												lastIndex = -1
												#This  is set as such because the last intersection was with an axis
												thisLine = bestLine
												checkTwoIntersec = len(thisLine.intersectedAx)
												if len(thisLine.intersectedAx) > 1:
													if k != 0:
													#	polyPts.append(nextLine.axIntersec[0])
														lastPt = 0
													else:
													#	polyPts.append(nextLine.axIntersec[1])
														lastPt = 1
													polyPts.extend(list(thisLine.PTline.coords))
													indexArray.append(thisIndex)
													print("Other axis intersect found, adding line:")
													print(thisLine.PTline)
													failedConnection = 0
													lineCount += 1 
													lastIndex = thisIndex #This is done so that the next time around it doesnt accidentally match with the same line if it can

												thisIndex = j
										

								else:
									print("Something went wrong making the polygon: trying to add a line without an intersection?")
									failedConnection += 1			



							else:
								print("Something went wrong making the polygon: trying to add a line without an intersection?")
								failedConnection += 1

						
					print("Polygon points in order:")
					for pnt in polyPts:
						print(pnt)
					
					if failedConnection < 4:
						firstLine = commonLines[indexArray[0]]
						lastLine = commonLines[indexArray[len(indexArray)-1]]

						#Need case for when its one line 
						if len(firstLine.axIntersec) == 1 and len(lastLine.axIntersec) ==1:
							print("Begin and end have axis intersection")
							if firstLine.intersectedAx != lastLine.intersectedAx:
								print("Axis intersection dont match")
								lastPt = Point(polyPts[len(polyPts)-1])
								print(lastPt)
								print(lastLine.axIntersec[0])

								x1 = lastPt.coords[0][0]
								y1 = lastPt.coords[0][1]
								
								x2= lastLine.axIntersec[0].coords[0][0]
								y2 = lastLine.axIntersec[0].coords[0][1]
								print(x1)
								print(x2)
								#x3, y3 = self.PTline.xy[len(self.PTline.xy)-1]
								
								
								if abs(x1-x2) <= EQ_THRESH and abs(y1-y2) <= EQ_THRESH:
								
									axIntersec = firstLine.intersectedAx[0].intersection(lastLine.intersectedAx[0])
									print("Found a corner:")
									print(axIntersec)
									if isinstance(axIntersec, Point):
										polyPts.append(axIntersec)
						elif len(commonLines)== 1 and len(firstLine.axIntersec) == 2:
							axIntersec = firstLine.intersectedAx[0].intersection(firstLine.intersectedAx[1])
							print("Found a corner:")
							print(axIntersec)
							if isinstance(axIntersec, Point):
								polyPts.append(axIntersec)


						self.polyList.append(DomPoly(polyPts,thisField))
					else:
						for line in commonLines:
							self.failedPolys.append(line)
							

					if i < len(self.domLines):
						if checkLeft:
							thisField = self.domLines[i].leftSide
							commonLines = [self.domLines[i]]
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
					print("")	
				else:
					if i < len(self.domLines):
						if checkLeft:
							thisField = self.domLines[i].leftSide
							commonLines = [self.domLines[i]]
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




			
	def printLines(self):
		for i in range(len(self.domLines)):
			print(self.domLines[i].leftSide + ": " + str(list(self.domLines[i].PTline.coords)))