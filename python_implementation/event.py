class Event:
	def __init__(self, eventObj, eventName, nextEdge, prevEdge):
		# rng of the event
		self.eventObj = eventObj
		self.eventName = eventName
		self.nextEdge = nextEdge
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