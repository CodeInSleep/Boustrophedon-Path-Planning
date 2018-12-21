#ifndef _PATH_PLANNING_H
#define _PATH_PLANNING_H

#include "header.hpp"

class Polygon {
public:
	// accepts 2d array of vertices of the survey area
	Polygon(double coords[][2], int num_of_pts);
	intersect_seq_t line_intersect(double xCoord);
private:
	pt_seq_t vertices;
	edge_seq_t edges;
	polygon_2d poly;
};

class SurveyArea : public Polygon {
public:
	SurveyArea(double coords[][2], int num_of_pts);
};

class Obstacle : public Polygon {

};

class BoustrophedonCell {
public:
	BoustrophedonCell(double leftEndPoint, segment_2d *startCeilSegmentPtr, segment_2d *startFloorSegmentPtr);
	void update(Event<point_2d> event);
	void update(Event<segment_2d> event);
	void terminate(double rightEndPoint);

private:
	double leftEndPoint;
	double rightEndPoint;
	std::vector<segment_2d*> ceilSegmentPtrs;
	std::vector<segment_2d*> floorSegmentPtrs;
};

#endif