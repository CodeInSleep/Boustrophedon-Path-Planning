#ifndef _EVENT_H
#define _EVENT_H

#include "util.hpp"

using namespace std;

string eventTypeNames[] = 
{
	"BEGIN",
	"END",
	"IN", 
	"OUT",
	"CEILING",
	"FLOOR",
};

class IEvent {
public:
	enum eventTypes {BEGIN, END, IN, OUT, CEILING, FLOOR};

	IEvent(eventTypes eventType, string eventName, const point_2d *eventObj, segment_2d *prevEdge, segment_2d *nextEdge);
	IEvent(eventTypes eventType, string eventName, const segment_2d *eventObj, segment_2d *prevEdge, segment_2d *nextEdge);
	const eventTypes eventType;
	const std::string eventName;
	const segment_2d *prevEdge;
	const segment_2d *nextEdge;
	
	inline double getxCoord() { return xCoord; }
	inline double getMinY() { return minY; }
	inline double getMaxY() { return maxY; }
	virtual ~IEvent() {}

protected:
	double xCoord;
	double minY;
	double maxY;
	

};

class PtEvent : public IEvent {
public:
	// Event could be point_2d or segment_2d
	PtEvent(enum IEvent::eventTypes eventType, string eventName, const point_2d *eventObj, 
		segment_2d *prevEdge, segment_2d *nextEdge);

	// virtual point_2d *getEventObj();
	const point_2d *eventObj;	
};

class SegEvent : public IEvent {
public:
	// Event could be point_2d or segment_2d
	SegEvent(enum IEvent::eventTypes eventType, string eventName, const segment_2d *eventObj, 
		segment_2d *prevEdge, segment_2d *nextEdge);

	// virtual segment_2d *getEventObj();
	const segment_2d *eventObj;	
};

#endif
