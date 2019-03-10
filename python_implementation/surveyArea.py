from polygon import *
from utils import *

# TODO: make intercept a namedtuple
class Opening:
	def __init__(self, intercept1, intercept2):
		self.intercept1 = intercept1
		self.intercept2 = intercept2

	def __eq__(self, other):
		return self.intercept1[0] == other.intercept1[0] and self.intercept2[0] == other.intercept2[0]

	def __hash__(self):
		return hash((self.intercept1[0], self.intercept2[0]))

	def __repr__(self):
		return 'opening from {} to {}'.format(self.intercept1[0], self.intercept2[0])

def nearestNeighbors(V, target, num=1):
    res = []
    min_distance = None

    # Each point is a State object
    for intercept in V:
        curr_distance = intercept[0].dist(target)
        if min_distance is None:
            res.append(intercept)
            min_distance = curr_distance
        elif curr_distance < min_distance:
            res = [intercept]
            min_distance = curr_distance
        # Can return multiple points if necessary
        elif curr_distance == min_distance:
            res.append(intercept)

    return res[0]

def getOpenings(events, sa):
	# params:
	# 	events: list of events at same x-coordinate
	#	sa: survey area to get openings in
	openings = set()
	_intersects = []
	for e in events:
		_intersects += intersects("x = {}".format(e.x),
			sa.edges+[edge for obs in sa.obstacles for edge in obs.edges])

		# print('event: ', e)
		# print('intersects: ')
		for idx, inter in enumerate(_intersects):
			# replace symbol with average of intersected edge
			if isinstance(inter[0].y, Symbol):
				inter[0] = (inter[1].v1 + inter[1].v2)*(1/2)
			print(inter)
			
		closestAbove = nearestNeighbors(list(filter(lambda pt: float(pt[0].y) > e.y2, _intersects)), Point(e.x, e.y1))
		closestBelow = nearestNeighbors(list(filter(lambda pt: float(pt[0].y) < e.y1, _intersects)), Point(e.x, e.y2))

		# add the opening above the event
		openings.add(Opening(closestBelow, (Point(e.x, e.y1), e.eventObj)))
		# add the opening below the event
		openings.add(Opening((Point(e.x, e.y2), e.eventObj), closestAbove))
	return list(openings)

class SurveyArea(Polygon):
	def __init__(self, vertices, obs_vertices, clockwise=True):
		# params:
		# 	vertices: 2D array of x, y coordinates to build SurveyArea
		super(SurveyArea, self).__init__(vertices, clockwise=clockwise)
		self.obstacles = []
		for obs_vertex in obs_vertices:
			# print('init obstacle: ', obs_vertex)
			self.obstacles.append(Obstacle(obs_vertex))

	def generateEvents(self, pp=False):
		# params:
		#   pp: print option
		events = super().generateEvents()

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
			# event with smallest x-coordinate is the INEVENT
			ob[0].eventName = INEVENT
			# event with largest x-coordinate is the OUTEVENT
			ob[-1].eventName = OUTEVENT
			if pp:
				print('obs %d events: ' % idx)
				for event in ob:
					print(event)

		flatten_ob_events = [e for ob_event in ob_events for e in ob_event]
		return events + flatten_ob_events

	def generateBCells(self):
		bcells = []
		events = self.generateEvents(pp=True)
		eventDict = {}
		for e in events:
			if e.x not in eventDict:
				eventDict[e.x] = [e]
			else:
				eventDict[e.x].append(e)

		for ex in eventDict:
			# TODO: process events with same x coordinate simultaneously
			events = eventDict[ex]
			firstEvent = events[0]
			en = firstEvent.eventName
			# for each event, find the affected B cells
			if en == BEGINEVENT:
				assert len(events) == 1
				bcells.append(BoustrophedonCell(firstEvent.x, firstEvent.prevEdge, firstEvent.nextEdge))
			elif en in INOUTEVENTS:
				openings = getOpenings(events, self)
				# terminate cells interfaced
				closing_rng = []
				# append the events, since they are non-overlapping
				for e in events:
					closing_rng.append([e.y1, e.y2])

				# print('IN/OUT EVENT: ', e)
				# open new cells
			elif en in MIDEVENTS:
				for cell in bcells:
					cell.update(events, self.obstacles)
			elif en == ENDEVENT:
				# print('events: ', events)
				assert len(events) == 1
				for cell in bcells:
					cell.update([firstEvent], self.obstacles)
					cell.terminate(firstEvent)
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

	def terminate(self, event):
		self.endX = event.x
		self.termianted = True

	def intersects(self, eventLine):
		# params:
		# 	eventLine: equation that represent event line
		if self.terminated:
			return []

		_intersects = intersects(eventLine, [self.ceilingEdges[-1], self.floorEdges[-1]])
		# print('Update BCell')
		# print('ceiling edge: ', self.ceilingEdges[-1])
		# print('floor edge: ', self.floorEdges[-1])
		# flrIntersect, ceilIntersect = _intersects
		# print('ceiling intercept: ', ceilIntersect)
		# print('floor intercept: ', flrIntersect)
		return _intersects

	def update(self, events, obstacles):
		# params:
		# 	event: event to update the B. cell with
		# 	obstacles: potential obstacles that may block cell
		#		from updating
		# return true if updated
		if self.terminated:
			return False

		for e in events:
			# last ceiling edge
			lce = self.ceilingEdges[-1]
			# last floor edge
			lfe = self.floorEdges[-1]
			# update ceiling edge if ceilIntercept is same as event x-coordinate
			
			if e.nextEdge == lce:
				if lce.edgeType == e.prevEdge.edgeType:
					self.ceilingEdges.append(e.prevEdge)

			if e.prevEdge == lfe:
				# print('event->nextEdge.edgeType: ', event.nextEdge.edgeType)
				if lfe.edgeType == e.nextEdge.edgeType:
					self.floorEdges.append(e.nextEdge)

		# print('ceilEdge: ', self.ceilingEdges)
		# print('floorEdges: ', self.floorEdges)
		return True