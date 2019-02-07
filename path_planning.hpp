#ifndef _PATH_PLANNING_H
#define _PATH_PLANNING_H

#include "util.hpp"
#include "Event.hpp"

using namespace std;

class BoustrophedonCell;
class Obstacle;
class Polygon;
class SurveyArea;

segment_2d scanLine(double xCoord, Polygon *poly);
intersect_seq_t line_intersect(Polygon *poly, IEvent *eventPtr);
intersect_seq_t line_intersect(BoustrophedonCell *cell, IEvent *eventPtr);
pair<intersect_t, intersect_t> eventLine(SurveyArea *sa, IEvent *event);


class Polygon {
public:
	enum polyTypes {INNER, OUTER};
	// accepts 2d array of vertices of the survey area
	Polygon(double coords[][2], int num_of_pts, polyTypes polyType);
	vector<IEvent *> generateEvents(string startEventName, string endEventName, bool sorted);
	inline pt_seq_t getVertices() { return vertices; }
	inline edge_seq_t getEdges() { return edges; }
	
protected:
	pt_seq_t vertices;
	edge_seq_t edges;
	polyTypes polyType;

private:
	polygon_2d poly;
};

class SurveyArea : public Polygon {
public:
	SurveyArea(double (*coords)[2], int num_of_pts);
	~SurveyArea();
	vector<IEvent *> generateEvents(bool sorted = true);

	void update(IEvent *event);
	vector<BoustrophedonCell *> openCells(IEvent *event, 
		vector<IEvent *>::iterator botEventIter, vector<IEvent *>::iterator topEventIter);
	void closeCells(IEvent *event);
	void updateCells(IEvent *event);

	vector<BoustrophedonCell *> generateBCells();
	inline vector<Obstacle *> getObstacles() { return obstacles; }

private:
	vector<BoustrophedonCell *> cells;
	vector<Obstacle *> obstacles;
};

class Obstacle : public Polygon {
public:
	Obstacle(double coords[][2], int num_of_pts);
	vector<IEvent *> generateEvents(bool sorted = true);
};

class BoustrophedonCell {
public:
	BoustrophedonCell(double leftStartPoint, const segment_2d *startCeilSegPtr, const segment_2d *startFloorSegPtr, SurveyArea *sa);
	// void update(IEvent *event);
	void terminate(double rightEndPoint);
	bool affected(IEvent *event);
	inline void addCeilSegPtr(const segment_2d* ceilSegPtr) { ceilSegPtrs.push_back(ceilSegPtr); }
	inline void addFloorSegPtr(const segment_2d* floorSegPtr) { floorSegPtrs.push_back(floorSegPtr); }

	inline double fromX() { return start; }
	inline double tillX() { return end; }
	inline vector<const segment_2d*> getCeilSegPtrs() { return ceilSegPtrs; }
	inline vector<const segment_2d*> getFloorSegPtrs() { return floorSegPtrs; }
	inline SurveyArea *within() { return surveyArea; }

private:
	bool terminated;
	double start;
	double end;
	vector<const segment_2d*> ceilSegPtrs;
	vector<const segment_2d*> floorSegPtrs;
	SurveyArea *surveyArea;
};

#endif