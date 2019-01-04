#include <algorithm> // for reverse, unique
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <boost/geometry/algorithms/assign.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
#include "path_planning.hpp"

using namespace std;

std::string boolstr(bool v)
{
    return v ? "true" : "false";
}

// configuration of test drone
struct drone_test_configuration {
    double altitude;
};

BoustrophedonCell::BoustrophedonCell(IEvent *startEvent, 
    const segment_2d *startCeilSegPtr, const segment_2d *startFloorSegPtr)
    : leftEndEvent(startEvent), terminated(false) {
        ceilSegPtrs.push_back(startCeilSegPtr);
        floorSegPtrs.push_back(startFloorSegPtr);
}

void BoustrophedonCell::terminate(IEvent *endEvent) {
    // if the cell has not been terminated, then terminate it. If it
    // has been terminated, we do nothing.
    if (!terminated) {
        terminated = true;
        rightEndEvent = endEvent;
    }
}

void BoustrophedonCell::update(IEvent *event) {
    // TODO: check if the current cell is involved in the event (i.e. in contact with the event line)
}

// function to generate path for a single B. cell
std::vector<point_2d> gen_path(BoustrophedonCell cell) {
    std::vector<point_2d> generatedPath;

    return generatedPath;
}

template <int D, bool increasing>
struct compare {
    bool operator()(const pt_iter_t &a_iter, const pt_iter_t &b_iter) {
        return (increasing) ? get<D>(*a_iter) < get<D>(*b_iter)
                        : get<D>(*a_iter) > get<D>(*b_iter);
    }
};

// function to find minimum and maximum item along dimension D
template <int D>
pair<double, double> minmax(pt_seq_t &seq) {
    double xMin = get<D>(seq[0]);
    double xMax = get<D>(seq[0]);

    // return the iterator(s) that point to the points in seq that have min/max
    // value in the dimension D.
    for (pt_seq_t::iterator iter = seq.begin(); iter < seq.end(); ++iter) {
        double xCoord = get<D>(*iter);
        if (xCoord < xMin)
            xMin = xCoord;

        if (xCoord > xMax)
            xMax = xCoord;
    }

    return pair<double, double>(xMin, xMax);
}

template <int D, int D_prime>
pt_seq_t leftCircularShift(pt_seq_t &pts) {
    // left circular shift so that the leftmost point of
    // the polygon is in the first position, followed by
    // vertices in counterclockwise traversal. Assume poly
    // is in counterclockwise traversal.
    pair<double, double> result = minmax<D>(pts);
    double minVal = result.first;

    vector<pt_iter_t> min_ele_iter_vec;
    for (pt_iter_t iter = pts.begin(); iter != pts.end(); ++iter) {
        if (get<D>(*iter) == minVal)
            min_ele_iter_vec.push_back(iter);
    }

    // sort the iterators in decreasing order according to the D_prime dimension of
    // the points pointed to by the iterators. We do this so that in the return 
    // vector, we know which point is on top.
    sort(min_ele_iter_vec.begin(), min_ele_iter_vec.end(), compare<D_prime, false>());

    // take the point with the smallest x-coordinate with the largest y-coordinate
    pt_const_iter_t startPtr = min_ele_iter_vec[0];
    
    cout << "left shift from: " << get<0>(*startPtr) << ", " << get<1>(*startPtr) << endl;

    // store shifted result
    pt_seq_t shifted_poly;
    
    pt_const_iter_t vtx_iter;
    for (vtx_iter = startPtr; vtx_iter < pts.end(); vtx_iter++) {
        shifted_poly.push_back(*vtx_iter);
    }
    for (vtx_iter = pts.begin(); vtx_iter < startPtr; vtx_iter++) {
        shifted_poly.push_back(*vtx_iter);
    }

    return shifted_poly;
}

template <typename T>
void printSeq(const vector<T> &seq) {
    // loop through coordinates
    for (typename vector<T>::const_iterator iter = seq.begin(); iter != seq.end(); ++iter) {
        cout << typeid(seq[0]).name() << dsv(*iter);
    }

    cout << "\n";
}

template<>
void printSeq<intersect_t>(const vector<intersect_t> &seq) {
    for (vector<intersect_t>::const_iterator iter = seq.begin(); iter != seq.end(); ++iter) {
        cout << typeid(seq[0]).name() << dsv((*iter).first) << ", " << dsv((*iter).second) << endl;
    }
}

Polygon::Polygon(double coords[][2], int num_of_pts, Polygon::polyTypes polyType)
    : polyType(polyType) {
    if (num_of_pts < 3)
        throw invalid_argument("received a polygon with less than 3 vertices..");
    // left circular shift so that the leftmost vertex is in the first position
    for (size_t i = 0; i < num_of_pts; ++i) {
        vertices.push_back(make<point_2d>(coords[i][0], coords[i][1]));
    }

    // the dimension in which to perform left circular shift
    const int D = 0;

    // the dimension in which to compare points that have same value in D
    // to find the top left vertex
    const int D_prime = 1;
    vertices = leftCircularShift<D, D_prime>(vertices);

    // ending point is the starting point
    vertices.push_back(vertices[0]);

    // create edges
    for (pt_const_iter_t iter = vertices.begin(); iter != vertices.end()-1; ++iter) {
        pt_const_iter_t nextPtIter = iter+1;
        edges.push_back(segment_2d(*iter, *nextPtIter));
    }

    cout << "vertices: " << endl;
    printSeq(vertices);
    cout << "edges: " << endl;
    printSeq(edges);

    // make a boost polygon representation from vertices
    assign_points<polygon_2d, vector<point_2d> >(poly, vertices);
}

segment_2d scanLine(double xCoord, double ymin, double ymax) {
    return segment_2d(
        point_2d(xCoord, ymin),
        point_2d(xCoord, ymax)
    );
}

// given x-coordinate of scanLine, reutrn the edges of the polygon intersected
// at xCoord and the corresponding points of intersection as pairs
intersect_seq_t Polygon::line_intersect(double xCoord) {
    pair<double, double> result = minmax<1>(vertices);

    // create vertical line that span the height of the polygon
    segment_2d line = scanLine(xCoord, result.first, result.second);

    intersect_seq_t intersects;
    for (edge_const_iter_t e_iter = edges.begin();e_iter != edges.end(); ++e_iter) {
        vector<point_2d> output;
        intersection(*e_iter, line, output);

        if (output.size() != 0)
            intersects.push_back(intersect_t(output[0], *e_iter));
    }

    return intersects;
}

vector<IEvent *> Polygon::generateEvents(string startEventName, string endEventName) {
    pair<double, double> result = minmax<0>(vertices);
    double minX = result.first;
    double maxX = result.second;

    vector<IEvent *> events;

    for (pt_const_iter_t iter = vertices.begin(); iter != vertices.end()-1; ++iter) {
        const point_2d *point_ptr = &(*iter);
        size_t vtx_idx = iter - vertices.begin();

        string eventName;
        enum IEvent::eventTypes eventType;
        double eventX = get<0>(*iter);
        if (eventX == minX) {
            eventName = startEventName;
            eventType = IEvent::BEGIN;
        }
        else if (eventX == maxX) {
            eventName = endEventName;
            eventType = IEvent::END;
        }
        else {
            eventName = "MIDDLE event";
            eventType = IEvent::MIDDLE;
        }

        size_t next_edge_idx = vtx_idx;
        size_t prev_edge_idx = (vtx_idx == 0) ? vertices.size()-1 : vtx_idx - 1;
        PtEvent *newPtEventPtr = new PtEvent(eventType, eventName, point_ptr, 
                &(edges[prev_edge_idx]), &(edges[next_edge_idx]));
        events.push_back(newPtEventPtr);

        // check if the event object is an edge
        if (iter - vertices.begin() > 0 && (get<0>(*iter) == get<0>(*(iter-1)))) {
            // pop last point event and add a new edge event
            events.pop_back();
            events.pop_back();
            // get the edge referred between current and previous points
            size_t event_edge_idx = prev_edge_idx;
            prev_edge_idx = (event_edge_idx == 0) ? edges.size()-1 : event_edge_idx - 1;
            const segment_2d *edge_ptr = &edges[event_edge_idx];

            SegEvent *newSegEventPtr = new SegEvent(eventType, eventName, edge_ptr,
                &edges[prev_edge_idx], &edges[next_edge_idx]);
            events.push_back(newSegEventPtr);
        }
    }

    return events;
}

// given xCoord and vector of obstacles generate the eventLine, which spans from the
// closest intersect (either with an obstacle or a border) above the event location to
// the close intersect (either with an obstacle or a border) below the event location
segment_2d SurveyArea::eventLine(IEvent &event) {

}

// find new openings to open new B. cells with. "opening lines" = eventLines - eventObjs
vector<segment_2d> SurveyArea::newOpenings(IEvent &event) {

}

SurveyArea::SurveyArea(double coords[][2], int num_of_pts)
    : Polygon(coords, num_of_pts, Polygon::OUTER) {}

void SurveyArea::update(IEvent *eventPtr) {
    enum IEvent::eventTypes eventType = eventPtr->eventType;

    if (eventType == IEvent::BEGIN) {
        // begin a new cell with the begin event
        cells.push_back(BoustrophedonCell(eventPtr, eventPtr->prevEdge, eventPtr->nextEdge));
    } else if (eventType == IEvent::END) {
        // end the cells that are in contact with the eventLine = scanLine - obstacles

        // segment_2d eventLine = eventLine()
        
    } else if (eventType == IEvent::IN) {
        // subtract obstacles that are not events from 
        // close off any opened cells with floor and ceiling edges included
        // in scan line.

    }
}

vector<IEvent *> SurveyArea::generateEvents() {
    return Polygon::generateEvents("BEGIN event", "END event");
}

Obstacle::Obstacle(double coords[][2], int num_of_pts)
    : Polygon(coords, num_of_pts, Polygon::INNER) {}

vector<IEvent *> Obstacle::generateEvents() {
    return Polygon::generateEvents("IN event", "OUT event");
}

