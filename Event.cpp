#include "Event.hpp"

using namespace std;

IEvent::IEvent(eventTypes evType, string eventName, const point_2d *eventObj, segment_2d *prevEdge, segment_2d *nextEdge)
	: eventType(evType), eventName(eventName), prevEdge(prevEdge), nextEdge(nextEdge) {
		xCoord = get<0>(*eventObj);
		minY = get<1>(*eventObj);
		maxY = get<1>(*eventObj);
	}

IEvent::IEvent(eventTypes evType, string eventName, const segment_2d *eventObj, segment_2d *prevEdge, segment_2d *nextEdge)
	: eventType(evType), eventName(eventName), prevEdge(prevEdge), nextEdge(nextEdge) {
		xCoord = get<0>((*eventObj).first);

		double oneEndY = get<1>((*eventObj).first);
		double otherEndY = get<1>((*eventObj).second);
		int greater = oneEndY > otherEndY ? 1 : 0;
		minY = (greater) ? otherEndY : oneEndY;
		maxY = (greater) ? oneEndY : otherEndY;
	}

// eventObj could be a segment_2d or point_2d
PtEvent::PtEvent(enum IEvent::eventTypes eventType, string eventName, const point_2d *eventObj, 
		segment_2d *prevEdge, segment_2d *nextEdge)
	: IEvent(eventType, eventName, eventObj, prevEdge, nextEdge), eventObj(eventObj) {}

SegEvent::SegEvent(enum IEvent::eventTypes eventType, string eventName, const segment_2d *eventObj, 
		segment_2d *prevEdge, segment_2d *nextEdge)
	: IEvent(eventType, eventName, eventObj, prevEdge, nextEdge), eventObj(eventObj) {}

// Event::Event(Event::eventType eventType, string eventName, point_2d *eventPt, segment_2d *prevEdge, segment_2d *nextEdge)
// 	: _eventType(eventType), _eventName(eventName), _eventPt(eventObj),
// 	_prevEdge(prevEdge), _nextEdge(_nextEdge) {}