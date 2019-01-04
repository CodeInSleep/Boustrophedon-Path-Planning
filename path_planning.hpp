#ifndef _PATH_PLANNING_H
#define _PATH_PLANNING_H

#include "header.hpp"
#include "Event.hpp"

using namespace std;

class BoustrophedonCell {
public:
	BoustrophedonCell(IEvent *startEvent, const segment_2d *startCeilSegPtr, const segment_2d *startFloorSegPtr);
	void update(IEvent *event);
	void terminate(IEvent *endEvent);

private:
	bool terminated;
	IEvent *leftEndEvent;
	IEvent *rightEndEvent;
	vector<const segment_2d*> ceilSegPtrs;
	vector<const segment_2d*> floorSegPtrs;
};

class Polygon {
public:
	enum polyTypes {INNER, OUTER};
	// accepts 2d array of vertices of the survey area
	Polygon(double coords[][2], int num_of_pts, polyTypes polyType);
	intersect_seq_t line_intersect(double xCoord);
	vector<IEvent *> generateEvents(string startEventName, string endEventName);
	
protected:
	pt_seq_t vertices;
	edge_seq_t edges;
	polyTypes polyType;

private:
	polygon_2d poly;
};

class SurveyArea : public Polygon {
public:
	SurveyArea(double coords[][2], int num_of_pts);
	vector<IEvent *> generateEvents();

	void update(IEvent *event);
	segment_2d SurveyArea::eventLine(IEvent &event);
	vector<segment_2d> SurveyArea::newOpenings(IEvent &event);

private:
	vector<BoustrophedonCell> cells;
	vector<Obstacle> obstacles;
};

class Obstacle : public Polygon {
public:
	Obstacle(double coords[][2], int num_of_pts);
	vector<IEvent *> generateEvents();
};
#endif