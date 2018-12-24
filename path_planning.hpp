#ifndef _PATH_PLANNING_H
#define _PATH_PLANNING_H

#include "header.hpp"
#include "Event.hpp"

using namespace std;

class Polygon {
public:
	// accepts 2d array of vertices of the survey area
	Polygon(double coords[][2], int num_of_pts);
	intersect_seq_t line_intersect(double xCoord);
	vector<IEvent *> generateEvents(string startEventName, string endEventName);
	
protected:
	pt_seq_t vertices;
	edge_seq_t edges;

private:
	polygon_2d poly;
};

class SurveyArea : public Polygon {
public:
	SurveyArea(double coords[][2], int num_of_pts);
	vector<IEvent *> generateEvents();
};

class Obstacle : public Polygon {
public:
	Obstacle(double coords[][2], int num_of_pts);
	vector<IEvent *> generateEvents();
};

class BoustrophedonCell {
public:
	BoustrophedonCell(double leftEndPoint, segment_2d *startCeilSegmentPtr, segment_2d *startFloorSegmentPtr);
	void update(IEvent *event);
	void terminate(double rightEndPoint);

private:
	double leftEndPoint;
	double rightEndPoint;
	vector<segment_2d*> ceilSegmentPtrs;
	vector<segment_2d*> floorSegmentPtrs;
};

#endif