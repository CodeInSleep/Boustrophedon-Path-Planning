class Event:
	def __init__(self, eventObj, eventName, nextEdge, prevEdge):
		# rng of the event
		self.eventObj = eventObj
		self.eventName = eventName

		# nextEdge point to the edge next to the event in counterclkwise direction
		self.nextEdge = nextEdge
		# prevEdge points to the edge previous to the event in counterclkwise direction
		self.prevEdge = prevEdge

	@classmethod
	def pointEvent(cls, *args):
		e = Event(*args)
		e.eventType = 'pointEvent'
		e.x = e.eventObj.x
		e.y1 = e.eventObj.y
		e.y2 = e.eventObj.y
		return e

	@classmethod
	def edgeEvent(cls, *args):
		e = Event(*args)
		e.eventType = 'edgeEvent'
		e.x = e.eventObj.v1.x
		e.y1 = e.eventObj.v1.y
		e.y2 = e.eventObj.v2.y
		if e.y1 > e.y2:
			temp = e.y2
			e.y2 = e.y1
			e.y1 = temp
		return e

	def __repr__(self):
		return 'Event({}, name={})'.format(str(self.eventObj),
			self.eventName, self.nextEdge, self.prevEdge)

	def __ge__(self, other):
		if self.x == other.x:
			return self.y1 >= other.y1
		return self.x > other.x