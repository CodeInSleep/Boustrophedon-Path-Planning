from utils import *
from sympy.solvers import solve
from sympy import Symbol
from event import Event
import math

BEGINEVENT = 'BEGINEVENT'
ENDEVENT = 'ENDEVENT'
INEVENT = 'INEVENT'
OUTEVENT = 'OUTEVENT'
MIDEVENTS = [INEVENT, OUTEVENT]

class Point:
	def __init__(self, x, y):
		self.x = x
		self.y = y

	def __eq__(self, other):
		return approxEqual(self.x, other.x) and approxEqual(self.y, other.y)

	def __gt__(self, other):
		if self.x == other.x:
			return self.y >= other.y
		return self.x > other.x

	def __repr__(self):
		return "Point(x={}, y={})".format(self.x, self.y)

class Edge:
	def __init__(self, v1, v2):
		# params:
		# 	v1: start point of the edge
		# 	v2: end point of the edge
		self.v1 = v1
		self.v2 = v2
		self.vertical = False
		if v2.x - v1.x != 0:
			self.slope = (v2.y - v1.y)/(v2.x - v1.x)
			self.intercept = abs(v2.x)*self.slope + v2.y if v2.x < 0 else -v2.x*self.slope + v2.y
			self.eq = "y = {}x + {}".format(self.slope,
				self.intercept)
		else:
			self.eq = "x = {}".format(v1.x)
			self.vertical = True

	def __repr__(self):
		return "Edge({}, {})".format(self.v1, self.v2)

class Polygon:
	def __init__(self, vertices, clockwise=True):
		# params:
		# 	vertices: 2D array of x, y coordinates to build polygon
		# 	clockwise: direction of definition

		print('vertices: ', vertices)
		# TODO: adjust for non-clockwise
		self.vertices = sorted([Point(v[0], v[1]) for v in vertices])
		self.edges = []
		# construct edges
		prev_vertex = self.vertices[-1]
		for v in self.vertices:
			self.edges.append(Edge(prev_vertex, v))
			prev_vertex = v

	def generateEvents(self, sorted=True):
		events = []
		for e in self.edges:
			if e.vertical:
				events.append(Event.edgeEvent(e))
			else:
				events.append(Event.pointEvent(e.v2))
		return events

class SurveyArea(Polygon):
	def __init__(self, vertices, obs_vertices, clockwise=True):
		# params:
		# 	vertices: 2D array of x, y coordinates to build SurveyArea
		super(SurveyArea, self).__init__(vertices, clockwise=clockwise)
		for obs_vertex in obs_vertices:
			self.obstacles = Obstacle(obs_vertex)

	def generateEvents(self):
		events = super().generateEvents(sorted=True)
		# overwrite start and end event names
		for idx, e in enumerate(events):
			if idx == 0:
				eventName = BEGINEVENT
			elif idx == len(events)-1:
				eventName = ENDEVENT
		return events

	def generateBCells(self):
		bcells = []
		events = self.generateEvents()
		for e in events:
			# for each event, find the affected B cells
			if e.eventName == BEGINEVENT:
				bcells.append(BoustrophedonCell(e.x))
			elif e.eventName in MIDEVENTS:
				for cell in bcells:

					cell.update(e)
		return bcells

class Obstacle(Polygon):
	def __init__(self, vertices, clockwise=True):
		super(Obstacle, self).__init__(vertices, clockwise=clockwise)

	def generateEvents(self):
		events = super().generateEvents(sorted=True)
		# overwrite start and end event names
		for idx, e in enumerate(events):
			if idx == 0:
				eventName = INEVENT
			elif idx == len(events)-1:
				eventName = OUTEVENT
		return events

class BoustrophedonCell:
	def __init__(self, startX):
		self.startX = startX
		self.terminated = False
		self.ceilingEdges = []
		self.floorEdges = []

	def toPolygon(self):
		pass

	def terminate(self, endX):
		self.endX = endX

	def update(self, event):
		# return true if updated
		if terminated:
			return False
		# update ceiling and floor edges
		pass
		return True

	def intersects(self, eventLine):
		pass


