#Class for parsing out a Plt file from theriak-domino
#Will find each "reaction" in the file and initialize them as DomLines
#Can also connect lines that should be attached

from DomLine import DomLine
from DomPoly import DomPoly
import re
from shapely.geometry import Point, LineString, Polygon
from shapely.ops import polygonize
import copy

class PltParser:

	def __init__(self, fileName):
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
		thisField = self.domLines[0].leftSide
	
		commonLines = []
		for i in range(len(self.domLines)):

			#Group together all lines that have the same leftside
			if self.domLines[i].leftSide == thisField:
				commonLines.append(self.domLines[i])
				# commonLines.append(self.domLines[i].extrapLine(bothEnds = True))
				print(self.domLines[i].leftSide)
				print(self.domLines[i].PTline)
			
			else:
				#Once a different leftSide is reached, search through the whole array for matching rightSides and add them to the commonLines
				for j in range(len(self.domLines)):
					if self.domLines[j].rightSide == thisField:
						commonLines.append(self.domLines[j])
						# commonLines.append(self.domLines[j].extrapLine(bothEnds = True))
						print(self.domLines[j].rightSide)
						print(self.domLines[j].PTline)
				

				lineCount = 0
				thisIndex = 0
				lastIndex = 0
				thisLine = commonLines[thisIndex]
				print("Commonline: " + str(thisIndex))
				print(thisLine.PTline)
				polyPts = []

				#This is for the edge case that thisLine only intersects the axes
				if len(thisLine.axIntersec) > 1:
					polyPts.extend(list(thisLine.PTline.coords))
					lineCount += 1
					lastPt = 1
					while len(thisLine.axIntersec) >1 and lineCount < len(commonLines):
						#Will continue this as long as thisLine continues to have two axes intersections
						for j in range(len(commonLines)):
							#Check all commonLines to see what line intersects with the same axis (Should only be one)
							if j != thisIndex and j != lastIndex:
								nextLine = commonLines[j]

								for k in range(len(nextLine.intersectedAx)):
									#Check both intersectedAxes and see if they are the same as the axes of the last point added to PolyPts
									if nextLine.intersectedAx[k] == thisLine.intersectedAx[lastPt]:
										#Checks if the intersection is at the beginning or end of nextLine
										#If its at the end then it reverses it
										print(thisLine.axInterLoc[lastPt])
										if (nextLine.axInterLoc[k] != 0 and thisLine.axInterLoc[lastPt] != 0) or (nextLine.axInterLoc[k] == 0 and thisLine.axInterLoc[lastPt] == 0):
											nextLine.PTline.coords = list(nextLine.PTline.coords)[::-1]
										
										#polyPts.append(nextLine.axIntersec[k].coords)
										
										lastIndex = -1
										#This  is set as such because the last intersection was with an axis
										thisLine = nextLine
										
										if len(thisLine.intersectedAx) > 1:
											if k != 0:
											#	polyPts.append(nextLine.axIntersec[0])
												lastPt = 0
											else:
											#	polyPts.append(nextLine.axIntersec[1])
												lastPt = 1
											polyPts.extend(list(thisLine.PTline.coords))
											lineCount += 1 
											lastIndex = thisIndex #This is done so that the next time around it doesnt accidentally match with the same line if it can

										thisIndex = j
				#elif len(thisLine.axIntersec) >0:
					#polyPts.extend(thisLine.axIntersec)

				
				print("Total lines = " + str(len(commonLines)))
				while lineCount < len(commonLines):
					#Continue until all lines in commonLines are acocunted for in the polygon
					#At this point thisLine should NOT be in polyPts and is either connected to two other domlines or one domline and an axis
					intersectNLine = None
					print("Linecount = " +str(lineCount))
					nextIndex = -1
					#Loop to check if the extrapolated line on either side matches with one of the other lines
					while intersectNLine == None and nextIndex < len(commonLines)-1:
						nextIndex += 1
						if nextIndex != thisIndex and nextIndex != lastIndex:
							intersectNLine = thisLine.extrapIntersec(commonLines[nextIndex]) #intersect is thisLine PLUS the intersection coordinate

					if intersectNLine != None:
						#Here we add thisLine plus the intersection as intersectNLine
						print("Found an intersection! Adding Line:")
						print(intersectNLine)
						polyPts.extend(list(intersectNLine.coords))
						lineCount += 1
						lastIndex = thisIndex
						thisIndex = nextIndex
						thisLine = commonLines[thisIndex]
						
						
						
					elif nextIndex >= len(commonLines) - 1:
						#This should only trigger if the next endpoint is an intersection with an axis
						#So thisLine is still not in polyPts and neither is its intersection with the axis
						#But the intersection with the last line still is in this
						print("No intersection found :(")
						checkTwoIntersec = 2
						lastPt = 0

						print("Linecount = " + str(lineCount))
						if len(thisLine.axIntersec)==1:
							#Can just intersect thisLine with the lastLine and attach in reverse order to polyPts (minus the last pt)
							# print("Commonline: " + str(thisIndex))
							# print(thisLine.PTline)
							# print("Commonline: " + str(lastIndex))
							# print(commonLines[lastIndex].PTline)
							intersectNLine = thisLine.extrapIntersec(commonLines[lastIndex])
							print(intersectNLine)
							intersectNLine.coords = list(intersectNLine.coords)[::-1]
							intersectNLine.coords = list(intersectNLine.coords)[1::]
							print("Adding line:")
							print(intersectNLine)
							polyPts.extend(list(intersectNLine.coords))
							lineCount += 1
							lastIndex = thisIndex
						#Now we are faced with the tricky issue	of finding the other line that intersects with the same axis	
						else:
							print("Something went wrong making the polygon: trying to add a line without an intersection?")		
						while checkTwoIntersec > 1 and lineCount < len(commonLines):
							#Will continue this loop as long as there is another line that has two intersections with the same axis
							#Will run at least once to find the nextLine that intersects the same axis
							#polyPts.extend(list(thisLine.axIntersec))
							numInterAx = 0
							for j in range(len(commonLines)):
								#Check how many other lines intersect this same axis
								if j != thisIndex and j != lastIndex:
									nextLine = commonLines[j]

									for k in range(len(nextLine.intersectedAx)):
										#Check both intersectedAxes and see if they are the same as the axes of the last point added to PolyPts
										if nextLine.intersectedAx[k] == thisLine.intersectedAx[lastPt]:
											numInterAx += 1 
							if numInterAx == 0:
								#This triggers in the literal corner case
							elif numInterax == 1:
								#this is the normal case
								for j in range(len(commonLines)):
								#Checks all commonLines for the same Axis Intersection

								#Problem: if at the corner, it could potentially loop forever
								#Other problem: if there are more than one line intersecting the axis, it will attach both of them these issues might be addressed by checking the 
								#number of lines that intersect the same axis in commonLines
									if j != thisIndex and j != lastIndex:
										nextLine = commonLines[j]

										for k in range(len(nextLine.intersectedAx)):
											#Check both intersectedAxes and see if they are the same as the axes of the last point added to PolyPts
											if nextLine.intersectedAx[k] == thisLine.intersectedAx[lastPt]:
											#Checks if the intersection is at the beginning or end of nextLine
											#If its at the end then it reverses it
												if (nextLine.axInterLoc[k] != 0 and thisLine.axInterLoc[lastPt] != 0) or (nextLine.axInterLoc[k] == 0 and thisLine.axInterLoc[lastPt] == 0):
													nextLine.PTline.coords = list(nextLine.PTline.coords)[::-1]
											
												#polyPts.append(nextLine.axIntersec[k].coords)
												matchAxes = True
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
													lineCount += 1 
													lastIndex = thisIndex #This is done so that the next time around it doesnt accidentally match with the same line if it can

												thisIndex = j

							

							elif numInteerAx >1:
								#This is the case when you have multiple polygons poking out different parts of the same axis


				print("Polygon points in order:")
				for pnt in polyPts:
					print(pnt)
			
				self.polyList = DomPoly(polyPts,thisField)
				commonLines =[self.domLines[i]] #Need to reset commonLines to prep next field
				thisField = self.domLines[i].leftSide










			
	def printLines(self):
		for i in range(len(self.domLines)):
			print(self.domLines[i].leftSide + ": " + str(list(self.domLines[i].PTline.coords)))