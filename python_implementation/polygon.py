from utils import *
from sympy.solvers import solve
from sympy import Symbol
from event import Event
import math
import re

BEGINEVENT = 'BEGIN EVENT'
ENDEVENT = 'END EVENT'
INEVENT = 'IN EVENT'
OUTEVENT = 'OUT EVENT'
CEILINGEVENT = 'CEILING EVENT'
FLOOREVENT = 'FLOOR EVENT'
MIDEVENTS = [INEVENT, OUTEVENT]

def linearEqsToAMF(eqs):
	coeffPattern = '(?:\s?)([-+]?\d+)?([xy])(?:\s?)'
	interceptPattern = '\d+'

	A = []
	b = []
	for eq in eqs:
		lhs, rhs = eq.split('=')
		# eq_pattern = '[\+\-]'.join(pattern for _ in range(len(eq.split('+'))))
		# eq_pattern += '\=(?:\s?)(?P<intercept>\d+)'
		# print('eq_pattern', eq_pattern)
		coeffs = re.findall(coeffPattern, lhs)
		intercept = re.match(interceptPattern, rhs)
		
		coeffs_dict = {var: coeff for coeff, var in coeffs}
		A.append([
			float(coeffs_dict.get('x', 0)),
			float(coeffs_dict.get('y', 0))])
		b.append(intercept.group(0))
	return A, b

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
			if self.slope < 0:
				self.eq = "y + {}x = {}".format(-self.slope,
					self.intercept)
			elif self.slope > 0:
				self.eq = "y - {}x = {}".format(self.slope,
					self.intercept)
			else:
				self.eq = "y = {}".format(self.intercept)
		else:
			self.eq = "x = {}".format(v1.x)
			self.vertical = True

	def __repr__(self):
		return "Edge({}, {}, eq = {})".format(self.v1, self.v2, self.eq)

class Polygon:
	def __init__(self, vertices, clockwise=True):
		# params:
		# 	vertices: 2D array of x, y coordinates to build polygon
		# 	clockwise: direction of definition

		# print('vertices: ', vertices)
		# TODO: adjust for non-clockwise
		self.vertices = [Point(v[0], v[1]) for v in vertices]
		self.edges = []
		# construct edges
		prev_vertex = self.vertices[-1]
		for v in self.vertices:
			self.edges.append(Edge(prev_vertex, v))
			prev_vertex = v

	def generateEvents(self, sort=True):
		events = []
		rightMostVertex = max(self.vertices)

		seenMaxVertex = False
		for e in self.edges:
			eventName = CEILINGEVENT if not seenMaxVertex else FLOOREVENT
			if e.vertical:
				events.append(Event.edgeEvent(e, eventName))
			else:
				events.append(Event.pointEvent(e.v2, eventName))

			if e.v2 == rightMostVertex:
				seenMaxVertex = True

		return sorted(events, key=lambda event: event.x)

class SurveyArea(Polygon):
	def __init__(self, vertices, obs_vertices, clockwise=True):
		# params:
		# 	vertices: 2D array of x, y coordinates to build SurveyArea
		super(SurveyArea, self).__init__(vertices, clockwise=clockwise)
		self.obstacles = []
		for obs_vertex in obs_vertices:
			self.obstacles.append(Obstacle(obs_vertex))

	def generateEvents(self, pp=False):
		# params:
		#   pp: print option
		events = super().generateEvents(sort=True)

		# overwrite start and end event names
		events[0].eventName = BEGINEVENT
		events[-1].eventName = ENDEVENT

		if pp:
			print('survey area events: ')
			for e in events:
				print(e)

		ob_events = []
		print('obstacles events: ')
		for idx, ob in enumerate(self.obstacles):
			_ob_events = ob.generateEvents()
			ob_events.append(_ob_events)

		for idx, ob in enumerate(ob_events):
			ob[0].eventName = INEVENT
			ob[-1].eventName = OUTEVENT
			if pp:
				print('obs %d events: ' % idx)
				for event in ob:
					print(event)

		return events

	def generateBCells(self):
		bcells = []
		events = self.generateEvents()
		for e in events:
			# TODO: process events with same x coordinate simultaneously

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
		events = super().generateEvents(sort=True)
		# overwrite start and end event names
		for idx, e in enumerate(events):
			if idx == 0:
				eventName = INEVENT
			elif idx == len(events)-1:
				eventName = OUTEVENT
		return events

class BoustrophedonCell:
	def __init__(self, startX, startCeiling, startFloor):
		self.startX = startX
		self.terminated = False
		self.ceilingEdges = [startCeiling]
		self.floorEdges = [startFloor]

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
		# params:
		# 	eventLine: equation that represent event line
		if terminated:
			return []

		pass


