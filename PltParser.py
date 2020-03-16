#Class for parsing out a Plt file from theriak-domino
#Will find each "reaction" in the file and initialize them as DomLines
#Can also connect lines that should be attached

from DomLine import DomLine
import re
from shapely.geometry import Point, LineString
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



			
	def printLines(self):
		for i in range(len(self.domLines)):
			print(self.domLines[i].leftSide + ": " + str(list(self.domLines[i].PTline.coords)))