#ifndef _PATH_PLANNING_H
#define _PATH_PLANNING_H

using namespace boost::geometry;
using namespace std;

typedef model::d2::point_xy<double> point_2d;
typedef model::segment<point_2d> segment_2d;
typedef model::polygon<point_2d, false> polygon_2d;
typedef vector<point_2d> coord_seq_t;
typedef coord_seq_t::iterator vec_iter;

class Event {
public:
	Event(int eventType, std::string eventName, point_2d coordinates);
	Event(int eventType, std::string eventName, segment_2d eventSegment);
	int type();

protected:
	int eventType;
	// generalize event coordinates to segment_2d. If event is a point, then eventSegment.p1 == eventSegment.p2
	segment_2d eventSegment;
	std::string eventName;
};

class InOutEvent : public Event {
public:
	InOutEvent(int eventType, std::string eventName, point_2d coordinates, segment_2d *adjFloorEdgePtr, segment_2d *adjCeilEdgePtr);

private:
	segment_2d *adjFloorEdge;
	segment_2d *adjCeilEdge;
};

class InEvent : public InOutEvent {
public:
	InEvent(std::string eventName, point_2d coordinates, segment_2d *adjFloorEdgePtr, segment_2d *adjCeilEdgePtr);
};

class BeginEvent : public InOutEvent {
public:
	BeginEvent(std::string eventName, point_2d coordinates, segment_2d *adjFloorEdgePtr, segment_2d *adjCeilEdgePtr);
};

class OutEvent : public InOutEvent {
public:
	OutEvent(std::string eventName, point_2d coordinates, segment_2d *adjFloorEdgePtr, segment_2d *adjCeilEdgePtr);
};

class EndEvent : public InOutEvent {
public:
	EndEvent(std::string eventName, point_2d coordinates, segment_2d *adjFloorEdgePtr, segment_2d *adjCeilEdgePtr);
};

class SurveyArea {
public:
	// accepts 2d array of vertices of the survey area
	SurveyArea(double coords[][2], int num_of_pts);
private:
	coord_seq_t vertices;
	polygon_2d poly;
};

class BoustrophedonCell {
public:
	BoustrophedonCell(double leftEndPoint, segment_2d *startCeilSegmentPtr, segment_2d *startFloorSegmentPtr);
	void update(Event event);
	void terminate(double rightEndPoint);

private:
	double leftEndPoint;
	double rightEndPoint;
	std::vector<segment_2d*> ceilSegmentPtrs;
	std::vector<segment_2d*> floorSegmentPtrs;
};

#endif