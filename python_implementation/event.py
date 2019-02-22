class Event:
	def __init__(self, eventObj):
		# rng of the event
		self.eventObj = eventObj

	@classmethod
	def pointEvent(cls, eventObj):
		e = Event(eventObj)
		e.eventType = 'pointEvent'
		e.x = eventObj.x
		return e

	@classmethod
	def edgeEvent(cls, eventObj):
		e = Event(eventObj)
		e.eventType = 'edgeEvent'
		e.x = eventObj.v1.x
		return e

	def __repr__(self):
		return 'Event: {}'.format(str(self.eventObj))