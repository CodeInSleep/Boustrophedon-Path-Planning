#include "Event.hpp"

using namespace std;

IEvent::IEvent(eventTypes evType, string eventName, 
		segment_2d *prevEdge, segment_2d *nextEdge)
	: eventType(evType), eventName(eventName), prevEdge(prevEdge), nextEdge(nextEdge) {}

// eventObj could be a segment_2d or point_2d
PtEvent::PtEvent(enum IEvent::eventTypes eventType, string eventName, const point_2d *eventObj, 
		segment_2d *prevEdge, segment_2d *nextEdge)
	: IEvent(eventType, eventName, prevEdge, nextEdge), eventObj(eventObj) {}

double PtEvent::getxCoord() {
	return get<0>(*eventObj);
}

SegEvent::SegEvent(enum IEvent::eventTypes eventType, string eventName, const segment_2d *eventObj, 
		segment_2d *prevEdge, segment_2d *nextEdge)
	: IEvent(eventType, eventName, prevEdge, nextEdge), eventObj(eventObj) {}

double SegEvent::getxCoord() {
	return get<0>((*eventObj).first);
}
// Event::Event(Event::eventType eventType, string eventName, point_2d *eventPt, segment_2d *prevEdge, segment_2d *nextEdge)
// 	: _eventType(eventType), _eventName(eventName), _eventPt(eventObj),
// 	_prevEdge(prevEdge), _nextEdge(_nextEdge) {}