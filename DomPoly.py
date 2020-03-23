#Class for handling a Domino Polygon
from shapely.geometry import Point, LineString, MultiLineString, Polygon
from shapely.ops import linemerge, snap
import copy

class DomPoly:

	def __init__(self, polyCoords,name):
		#Coordinates will be defined outside of this class constructor
		self.field = Polygon(polyCoords)
		self.phases = name
