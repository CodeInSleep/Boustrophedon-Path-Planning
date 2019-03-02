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
		return e

	@classmethod
	def edgeEvent(cls, *args):
		e = Event(*args)
		e.eventType = 'edgeEvent'
		e.x = e.eventObj.v1.x
		return e

	def __repr__(self):
		return 'Event({}, name={}, nextEdge={}, prevEdge={})'.format(str(self.eventObj),
			self.eventName, self.nextEdge, self.prevEdge)