from utils import *
from sympy import *

from sympy.solvers import solve
from event import Event
import math
import re

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

	def __ge__(self, other):
		if self.x == other.x:
			return self.y >= other.y
		return self.x > other.x

	def __lt__(self, other):
		return not self.__ge__(other)

	def __add__(self, other):
		return Point(self.x+other.x, self.y+other.y)

	def __mul__(self, other):
		self.x *= other
		self.y *= other
		return self

	def __repr__(self):
		return "Point(x={}, y={})".format(self.x, self.y)

	def __hash__(self):
		return hash((self.x, self.y))

	# Return Euclidean distance between two States
	def dist(self, other):
	    return math.sqrt((self.x - other.x)**2 + (self.y - other.y)**2)

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

	def __gt__(self, other):
		return self.v1 > other.v1

	def __repr__(self):
		return "Edge({}, {}, eq = {})".format(self.v1, self.v2, self.eq)

	def __hash__(self):
		return hash((self.v1, self.v2))

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

def inLine(edge, intercept):
	# if there are infinite solutions (colinear with the eventline)
	if isinstance(intercept.x, Symbol) or isinstance(intercept.y, Symbol):
		return True
	# check if intercept is in the line
	minPoint = min(edge.v1, edge.v2)
	maxPoint = max(edge.v1, edge.v2)

	return intercept >= minPoint and intercept <= maxPoint

def intersects(eventLine, edges):
	# params:
	# 	eventLine: equation that represent event line

	def solve_for_intercept(eq1, eq2, intercepts):
		A, b = linearEqsToAMF([eq1, eq2])
		x, y = symbols('x y')
		
		print('A: ', A)
		print('b: ', b)
		intercept = linsolve((A, b), [x, y])
		print('INTERcept_arg: ', intercept.args)
		
		return set() if not intercept.args else intercept.args[0]

	intercepts = []
	for e in edges:
		_intercept = solve_for_intercept(eventLine, e.eq, intercepts)
		if _intercept:
			intercept_pt = Point(*_intercept)
			if inLine(e, intercept_pt):
				intercepts.append([intercept_pt, e])

	# sorted(intercepts)
	return intercepts

	


