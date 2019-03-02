from utils import *
from sympy import *

from sympy.solvers import solve
from event import Event
import math
import re

BEGINEVENT = 'BEGIN EVENT'
ENDEVENT = 'END EVENT'
INEVENT = 'IN EVENT'
OUTEVENT = 'OUT EVENT'
CEILINGEVENT = 'CEILING EVENT'
FLOOREVENT = 'FLOOR EVENT'
MIDEVENTS = [CEILINGEVENT, FLOOREVENT]

def linearEqsToAMF(eqs):
	coeffPattern = '(?:\s?)([-+]?(?:\s?)[\d+\.\d+]+)?([xy])(?:\s?)'
	interceptPattern = '(?:\s?)([-+]?[\d+\.\d+]+)(?:\s?)'

	A = []
	b = []
	for eq in eqs:
		print('eq: ', eq)
		lhs, rhs = eq.split('=')
		
		coeffs = re.findall(coeffPattern, lhs)
		intercept = re.match(interceptPattern, rhs)
		
		coeffs_dict = {var: coeff.replace(' ', '') if coeff else '1' for coeff, var in coeffs}
		# print('coeff_dict: ', type(coeffs_dict.get('x', 0)))
		A.append([
			float(coeffs_dict.get('x', 0)),
			float(coeffs_dict.get('y', 0))])
		b.append(float(intercept.group(0)))
		
	return Matrix(A), Matrix(b)

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
	def __init__(self, v1, v2, edgeType, prevEdge=None, nextEdge=None):
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

		self.prevEdge = prevEdge
		self.nextEdge = nextEdge
		self.edgeType = edgeType

	def __eq__(self, other):
		return (self.v1 == other.v1 and self.v2 == other.v2) or \
			(self.v1 == other.v2 and self.v2 == other.v1)

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
		print('vertices: ', self.vertices)
		self.edges = []
		# construct edges
		prev_vertex = self.vertices[-1]
		rightMostVertex = max(self.vertices)
		edgeType = 'C'
		for v in self.vertices:
			self.edges.append(Edge(prev_vertex, v, edgeType))
			if v == rightMostVertex:
				edgeType = 'F'
			prev_vertex = v

		for e in self.edges:
			print('edge: ', e)
			print('type: ', e.edgeType)

	def generateEvents(self, sort=True):
		events = []
		rightMostVertex = max(self.vertices)

		seenMaxVertex = False
		for idx, e in enumerate(self.edges):
			eventName = CEILINGEVENT if not seenMaxVertex else FLOOREVENT
			nextEdge = self.edges[idx+1] if idx+1 < len(self.edges) else self.edges[0]
			if e.vertical:
				prevEdge = self.edges[idx-1] if idx-1 >= 0 else self.edges[-1]
				events.append(Event.edgeEvent(e, eventName, prevEdge, nextEdge))
			else:
				prevEdge = e
				events.append(Event.pointEvent(e.v2, eventName, prevEdge, nextEdge))

			if e.v2 == rightMostVertex:
				seenMaxVertex = True

		return sorted(events, key=lambda event: event.x)

def intersects(eventLine, edges):
	# params:
	# 	eventLine: equation that represent event line

	def solve_for_intercept(eq1, eq2):
		A, b = linearEqsToAMF([eq1, eq2])
		x, y = symbols('x y')
		
		print('A: ', A)
		print('b: ', b)
		intercept = linsolve((A, b), [x, y])
		return intercept.args[0]

	intercepts = []
	for e in edges:
		intercept = solve_for_intercept(eventLine, e.eq)
		if intercept:
			intercepts.append(Point(*intercept))

	sorted(intercepts)
	return intercepts

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
			for e in events:
				print(e)

		ob_events = []
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
				bcells.append(BoustrophedonCell(e.x, e.prevEdge, e.nextEdge))
			elif e.eventName in MIDEVENTS or e.eventName == ENDEVENT:
				for cell in bcells:
					cell.update(e, self.obstacles)
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

	def terminate(self, endX):
		self.endX = endX

	def intersects(self, eventLine):
		# params:
		# 	eventLine: equation that represent event line
		if self.terminated:
			return []

		_intersects = intersects(eventLine, [self.ceilingEdges[-1], self.floorEdges[-1]])
		print('Update BCell')
		print('ceiling edge: ', self.ceilingEdges[-1])
		print('floor edge: ', self.floorEdges[-1])
		flrIntersect, ceilIntersect = _intersects
		print('ceiling intercept: ', ceilIntersect)
		print('floor intercept: ', flrIntersect)
		return _intersects

	def update(self, event, obstacles):
		# params:
		# 	event: event to update the B. cell with
		# 	obstacles: potential obstacles that may block cell
		#		from updating
		# return true if updated
		if self.terminated:
			return False

		print(event.eventName)

		# last ceiling edge
		lce = self.ceilingEdges[-1]
		# last floor edge
		lfe = self.floorEdges[-1]
		# update ceiling edge if ceilIntercept is same as event x-coordinate
		
		print('event: ', event)
		if event.nextEdge == lce:
			if lce.edgeType == event.prevEdge.edgeType:
				self.ceilingEdges.append(event.prevEdge)

		if event.prevEdge == lfe:
			# print('event->nextEdge.edgeType: ', event.nextEdge.edgeType)
			if lfe.edgeType == event.nextEdge.edgeType:
				self.floorEdges.append(event.nextEdge)

		print('ceilEdge: ', self.ceilingEdges)
		print('floorEdges: ', self.floorEdges)
		return True

	


